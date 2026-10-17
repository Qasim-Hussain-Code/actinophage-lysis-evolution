#!/usr/bin/env python3
"""Archive frozen cohort-reference search evidence without changing any calls."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'src'))

from Bio import SeqIO
from actinolysis.common import NA, budget_check, load_config, read_tsv, sha256, utcnow, write_json, write_tsv
from actinolysis.identification import annotation_category, evaluate_hit
import actinolysis.identification as identification

HIT_FIELDS = ['protein_id', 'reference_id', 'identity', 'alignment_length',
              'query_coverage', 'target_coverage', 'evalue', 'bits']
QUERY_FIELDS = ['protein_id', 'protein_accession', 'accession', 'cds_index', 'start', 'end',
                'location_parts', 'strand', 'partial', 'product', 'note']
HIT_OUTPUT_FIELDS = QUERY_FIELDS + [
    'reference_id', 'reference_category', 'reference_evidence_type', 'reference_source_accession',
    'reference_protein_accession', 'reference_url', 'identity', 'alignment_length', 'query_coverage',
    'target_coverage', 'evalue', 'bits', 'confidence_before_partial_cap', 'score_confidence',
    'negative_control', 'accepted_under_frozen_rules', 'selected_for_category',
    'selected_evidence_confidence', 'min_identity', 'min_query_coverage', 'min_target_coverage',
    'max_evalue', 'passes_identity', 'passes_query_coverage', 'passes_target_coverage', 'passes_evalue']
STATUS_FIELDS = QUERY_FIELDS + [
    'searched', 'reported_hit_count', 'reported_nonnegative_hit_count', 'negative_control_hit_count',
    'reported_nonnegative_category_count', 'reported_nonnegative_categories',
    'accepted_hit_count', 'accepted_category_count', 'accepted_categories', 'selected_reference_hit_count',
    'selected_category_count', 'selected_categories', 'supported_selected_category_count',
    'supported_selected_categories', 'annotation_category', 'detection_status']
REQUIRED_OUTPUTS = {
    '2': ['data/processed/genomes.tsv', 'data/processed/cds.tsv', 'data/processed/proteins.faa'],
    '3': ['data/processed/reference_hits.tsv', 'data/manifests/reference_manifest.tsv', 'data/processed/references.faa'],
    '4': ['data/processed/lysis_evidence.tsv'],
}


def scientific_fingerprint(cfg):
    fields = ['mode', 'seed', 'sample_sizes', 'sampling', 'identification', 'architecture', 'neighborhood', 'statistics']
    return hashlib.sha256(json.dumps({key: cfg[key] for key in fields}, sort_keys=True).encode()).hexdigest()


def contained_path(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f'Receipt path escapes isolated root: {relative}')
    return path


def verify_receipts(root, cfg):
    """Check completed search, query provenance and selected evidence before export."""
    path = root / f"results/reports/{cfg['mode']}_stage_receipts.json"
    receipts = json.loads(path.read_text(encoding='utf-8'))
    checked = {str(path.relative_to(root)): sha256(path), 'config/analysis.yaml': sha256(root / 'config/analysis.yaml')}
    fingerprint = scientific_fingerprint(cfg)
    for number, required in REQUIRED_OUTPUTS.items():
        receipt = receipts.get(number)
        if not receipt or receipt.get('scientific_config_sha256') != fingerprint:
            raise ValueError(f'Stage {number} is incomplete or its scientific configuration changed')
        outputs = receipt.get('outputs', {})
        if not set(required).issubset(outputs):
            raise ValueError(f'Stage {number} receipt lacks required frozen outputs')
        sources = receipt.get('source_sha256', {})
        required_sources = {'src/actinolysis/common.py', 'scripts/run_all.py',
                            f"src/actinolysis/{'quality' if number == '2' else 'identification'}.py"}
        if number in {'3', '4'}:
            required_sources.add('references/protein_controls.tsv')
        normalized_sources = {name.replace('\\', '/'): digest for name, digest in sources.items()}
        if not required_sources.issubset(normalized_sources):
            raise ValueError(f'Stage {number} receipt lacks required frozen source hashes')
        for name, digest in {**outputs, **sources}.items():
            relative = name.replace('\\', '/')
            target = contained_path(root, relative)
            if not target.is_file() or sha256(target) != digest:
                raise ValueError(f'Stage {number} frozen input changed: {relative}')
            if relative in checked and checked[relative] != digest:
                raise ValueError(f'Stage receipts disagree on input: {relative}')
            checked[relative] = digest
    evaluator = Path(identification.__file__).resolve()
    if sha256(evaluator) != checked['src/actinolysis/identification.py']:
        raise ValueError('Current evaluate_hit source differs from the frozen stage 3/4 identification source')
    return checked, {number: receipts[number] for number in REQUIRED_OUTPUTS}


def unique_rows(path, key):
    rows = read_tsv(path)
    result = {row[key]: row for row in rows}
    if len(result) != len(rows) or any(not value for value in result):
        raise ValueError(f'Duplicate or empty {key}: {path}')
    return result


def read_hits(path, cds, references, searched):
    result = []; seen = set()
    with path.open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, fieldnames=HIT_FIELDS, delimiter='\t'):
            if None in row or any(row[key] is None for key in HIT_FIELDS):
                raise ValueError('Malformed eight-column reference hit row')
            pair = row['protein_id'], row['reference_id']
            if pair in seen or pair[0] not in cds or pair[0] not in searched or pair[1] not in references:
                raise ValueError(f'Duplicate, unsearched or unknown reference hit identifiers: {pair}')
            seen.add(pair)
            for field in ['identity', 'query_coverage', 'target_coverage', 'evalue', 'bits']:
                value = float(row[field])
                if not math.isfinite(value) or value < 0 or (field in {'identity', 'query_coverage', 'target_coverage'} and value > 1):
                    raise ValueError(f'Invalid alignment {field}: {pair}')
            if int(row['alignment_length']) <= 0:
                raise ValueError(f'Invalid alignment length: {pair}')
            result.append(row)
    return sorted(result, key=lambda row: (row['protein_id'], row['reference_id']))


def build_tables(cds, references, hits, evidence, searched, cfg):
    selected = {}; evidence_by_protein = defaultdict(list)
    hit_pairs = {(row['protein_id'], row['reference_id']) for row in hits}
    for row in evidence:
        pid = row['protein_id']
        if pid not in cds:
            raise ValueError(f'Selected evidence contains an unknown protein: {pid}')
        evidence_by_protein[pid].append(row)
        if row['source'] != 'curated_reference_homology':
            continue
        rid = row['reference_id']; key = (pid, row['category'], rid)
        if rid not in references or references[rid]['category'] != row['category'] or (pid, rid) not in hit_pairs or key in selected:
            raise ValueError(f'Selected evidence does not identify a unique reported category hit: {key}')
        selected[key] = row
    limits = cfg['identification']; hit_rows = []; hits_by_protein = defaultdict(list)
    for hit in hits:
        pid = hit['protein_id']; ref = references[hit['reference_id']]; query = cds[pid]
        negative = ref['category'] == 'negative_control'
        confidence = evaluate_hit(hit, ref, cfg)
        final = 'low' if query['partial'] == 'true' and not negative else confidence
        selected_row = selected.get((pid, ref['category'], hit['reference_id']))
        if selected_row and selected_row['confidence'] != final:
            raise ValueError(f'Recomputed confidence disagrees with frozen selected evidence: {pid}/{hit["reference_id"]}')
        row = {key: query.get(key, NA) for key in QUERY_FIELDS}
        row.update(reference_id=hit['reference_id'], reference_category=ref['category'],
                   reference_evidence_type=ref['evidence_type'],
                   reference_source_accession=ref.get('resolved_source_accession', ref.get('source_accession', NA)),
                   reference_protein_accession=ref.get('resolved_protein_accession', ref.get('protein_accession', NA)),
                   reference_url=ref.get('reference_url', NA),
                   **{key: hit[key] for key in HIT_FIELDS if key not in {'protein_id', 'reference_id'}},
                   confidence_before_partial_cap=confidence, score_confidence=final, negative_control=negative,
                   accepted_under_frozen_rules=not negative and final in {'high', 'moderate'},
                   selected_for_category=selected_row is not None,
                   selected_evidence_confidence=selected_row['confidence'] if selected_row else NA,
                   **{key: limits[key] for key in ['min_identity', 'min_query_coverage', 'min_target_coverage', 'max_evalue']},
                   passes_identity=float(hit['identity']) >= limits['min_identity'],
                   passes_query_coverage=float(hit['query_coverage']) >= limits['min_query_coverage'],
                   passes_target_coverage=float(hit['target_coverage']) >= limits['min_target_coverage'],
                   passes_evalue=float(hit['evalue']) <= limits['max_evalue'])
        hit_rows.append(row); hits_by_protein[pid].append(row)
    status_rows = []
    for pid, query in sorted(cds.items()):
        local = hits_by_protein[pid]; ledger = evidence_by_protein[pid]
        accepted = [row for row in local if row['accepted_under_frozen_rules']]
        reported_categories = {row['reference_category'] for row in local if not row['negative_control']}
        accepted_categories = {row['reference_category'] for row in accepted}
        selected_categories = {row['category'] for row in ledger}
        supported = {row['category'] for row in ledger if row['source'] == 'curated_reference_homology' and row['confidence'] in {'high', 'moderate'}}
        if pid not in searched:
            status = 'not_searched_no_eligible_translation'
        elif supported:
            status = 'selected_high_or_moderate_homology_assignment'
        elif accepted:
            status = 'reported_qualifying_hit_without_selected_support'
        elif not local:
            status = 'unresolved_no_reported_reference_hit'
        elif all(row['negative_control'] for row in local):
            status = 'unresolved_only_negative_control_hits'
        else:
            status = 'unresolved_reported_hits_below_criteria_or_evidence_tier'
        status_rows.append({**{key: query.get(key, NA) for key in QUERY_FIELDS}, 'searched': pid in searched,
                           'reported_hit_count': len(local),
                           'reported_nonnegative_hit_count': sum(not row['negative_control'] for row in local),
                           'negative_control_hit_count': sum(row['negative_control'] for row in local),
                           'reported_nonnegative_category_count': len(reported_categories),
                           'reported_nonnegative_categories': ';'.join(sorted(reported_categories)) or NA,
                           'accepted_hit_count': len(accepted),
                           'accepted_category_count': len(accepted_categories),
                           'accepted_categories': ';'.join(sorted(accepted_categories)) or NA,
                           'selected_reference_hit_count': sum(row['selected_for_category'] for row in local),
                           'selected_category_count': len(selected_categories),
                           'selected_categories': ';'.join(sorted(selected_categories)) or NA,
                           'supported_selected_category_count': len(supported),
                           'supported_selected_categories': ';'.join(sorted(supported)) or NA,
                           'annotation_category': annotation_category(query['product'] + ' ' + query.get('note', '')),
                           'detection_status': status})
    return hit_rows, status_rows


def serialize_booleans(rows):
    return [{key: str(value).lower() if isinstance(value, bool) else value for key, value in row.items()} for row in rows]


def archive(root, out, project_root=ROOT):
    root = Path(root).resolve(); out = Path(out).resolve(); project_root = Path(project_root).resolve()
    if not out.is_relative_to(project_root) or out == root or out.is_relative_to(root):
        raise ValueError('Supplement output must be inside the outer project and outside the isolated analysis root')
    cfg = load_config(root / 'config/analysis.yaml')
    inputs, receipts = verify_receipts(root, cfg)
    cds = unique_rows(root / 'data/processed/cds.tsv', 'protein_id')
    references = unique_rows(root / 'data/manifests/reference_manifest.tsv', 'reference_id')
    search_report_path = root / 'results/reports/reference_search.json'
    search_report = json.loads(search_report_path.read_text(encoding='utf-8'))
    if (search_report.get('hit_file_sha256') != inputs['data/processed/reference_hits.tsv']
            or search_report.get('reference_count') != len(references)):
        raise ValueError('Search execution report disagrees with the frozen stage 3 outputs')
    inputs['results/reports/reference_search.json'] = sha256(search_report_path)
    records = list(SeqIO.parse(root / 'data/processed/proteins.faa', 'fasta'))
    searched = {record.id for record in records}
    if len(searched) != len(records) or not searched.issubset(cds):
        raise ValueError('Frozen query FASTA has duplicate or unknown protein identifiers')
    if any(row['partial'] not in {'true', 'false'} for row in cds.values()):
        raise ValueError('Frozen CDS partial flags must be true or false')
    hits = read_hits(root / 'data/processed/reference_hits.tsv', cds, references, searched)
    hit_rows, status_rows = build_tables(cds, references, hits, read_tsv(root / 'data/processed/lysis_evidence.tsv'), searched, cfg)
    code_sha = sha256(Path(__file__))
    frozen = {'isolated_root': str(root), 'inputs_sha256': inputs, 'archive_code_sha256': code_sha,
              'evaluator_sha256': sha256(Path(identification.__file__))}
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    report_path = out / 'archive_report.json'
    if out.exists() and any(out.iterdir()):
        if not report_path.is_file():
            raise ValueError('Nonempty output directory lacks an archive report; choose a dedicated supplement directory')
        old = json.loads(report_path.read_text(encoding='utf-8'))
        if old.get('archive_fingerprint_sha256') != fingerprint:
            raise ValueError('Archive inputs or code changed; preserve this output and choose a new directory')
        for name, digest in old.get('outputs_sha256', {}).items():
            if sha256(contained_path(out, name)) != digest:
                raise ValueError(f'Existing supplement output changed: {name}')
    # Complete all validation, including a second input check, before publishing.
    rechecked = verify_receipts(root, cfg)[0]
    rechecked['results/reports/reference_search.json'] = sha256(search_report_path)
    if rechecked != inputs:
        raise ValueError('Frozen inputs changed while preparing the supplement')
    outer_cfg = load_config(project_root / 'config/analysis.yaml')
    budget_check(outer_cfg, 10_000_000, root=project_root)
    write_tsv(out / 'all_cohort_reference_hits.tsv', serialize_booleans(hit_rows), HIT_OUTPUT_FIELDS)
    write_tsv(out / 'protein_detection_status.tsv', serialize_booleans(status_rows), STATUS_FIELDS)
    report = {**frozen, 'archive_fingerprint_sha256': fingerprint, 'created_utc': utcnow(),
              'verified_stage_receipts': receipts, 'identification_thresholds': cfg['identification'],
              'search_execution_report': search_report,
              'counts': {'cds': len(cds), 'searched_proteins': len(searched), 'reported_reference_hits': len(hit_rows),
                         'accepted_reference_hits': sum(row['accepted_under_frozen_rules'] for row in hit_rows),
                         'selected_reference_hits': sum(row['selected_for_category'] for row in hit_rows),
                         'negative_control_hits': sum(row['negative_control'] for row in hit_rows)},
              'detection_status_counts': dict(sorted(Counter(row['detection_status'] for row in status_rows).items())),
              'reference_evidence_types': dict(sorted(Counter(row['evidence_type'] for row in references.values()).items())),
              'outputs_sha256': {name: sha256(out / name) for name in ['all_cohort_reference_hits.tsv', 'protein_detection_status.tsv']},
              'interpretation': [
                  'Read-only supplement: no primary assignment, function, architecture or confidence was changed.',
                  'Acceptance is the frozen evaluate_hit result with partial proteins capped at low; negative controls are never accepted.',
                  'Selected flags reproduce the frozen category ledger and are separate from each hit scoring above the cutoffs.',
                  'No reported hit means unresolved under this reference panel and search reporting settings, never biological absence.',
                  'The hit table preserves every reported cohort-reference alignment, not unreported alignments below the MMseqs search/reporting gates.',
                  'Input and table hashes are recorded; the report does not hash itself.']}
    write_json(report_path, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='Frozen isolated analysis root')
    parser.add_argument('--out', type=Path, required=True, help='Dedicated supplement directory in the outer project')
    args = parser.parse_args(argv)
    report = archive(args.root, args.out)
    print(json.dumps({'out': str(args.out.resolve()), 'counts': report['counts']}, sort_keys=True))


if __name__ == '__main__':
    main()
