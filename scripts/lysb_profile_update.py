#!/usr/bin/env python3
"""Evaluate fourteen assay references using the original frozen HMM functions."""
from __future__ import annotations
import argparse
import importlib.util
from importlib.metadata import version
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import budget_check, load_config, read_tsv, sha256, write_json, write_tsv

spec = importlib.util.spec_from_file_location('lysb_original_engine', ROOT / 'scripts/lysb_profile_candidates.py')
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
ADDED = {'MissWhite_LysB_2026', 'Bxb1_LysB_2026', 'Constella_LysB_2026', 'Bongo_LysB_2026', 'phiGD57_1_LysB_2026'}
BROAD = engine.BROAD | ADDED
PROTOCOL = ROOT / 'docs/lysb_profile_update_protocol.md'

def normalize_metadata(hmm, path):
    """Omit execution metadata while preserving all numerical model parameters."""
    hmm.creation_time = None
    hmm.command_line = None
    with Path(path).open('wb') as stream: hmm.write(stream)

def validate_training(references):
    eligible = {i for i, row in references.items() if row['category'] == 'lysin_b' and
                row['evidence_type'] in {'experimental', 'experimental_generic_lipolytic', 'experimental_antimycobacterial'}}
    if eligible != BROAD or len(BROAD) != 14:
        raise ValueError('Training must contain exactly the fourteen prespecified assay references')
    if {i for i, row in references.items() if row['category'] == 'lysin_b' and row['evidence_type'] == 'published_profile_prediction'} != engine.PREDICTED:
        raise ValueError('The three Gordonia computational predictions must remain evaluation-only')

def run(root, validation, out):
    import pyfamsa
    import pyhmmer
    root, validation, out = (Path(p).resolve() for p in (root, validation, out))
    if not out.is_relative_to(ROOT) or out == root or root.is_relative_to(out) or out.is_relative_to(validation) or validation.is_relative_to(out):
        raise ValueError('Choose a separate output directory inside the project')
    for protected in ('data', 'references', 'src', 'scripts', 'figures', 'results/tables', 'results/reports'):
        if out.is_relative_to(root / protected): raise ValueError('Output overlaps a protected pipeline directory')
    versions = {p: version(p) for p in ('pyfamsa', 'pyhmmer', 'scoring-matrices', 'psutil')}
    if versions != {'pyfamsa': '0.7.0', 'pyhmmer': '0.12.3', 'scoring-matrices': '0.3.5', 'psutil': '7.2.2'}:
        raise ValueError('Use the pinned profile dependencies')
    paths = {'references': root / 'data/processed/references.faa', 'reference_manifest': root / 'data/manifests/reference_manifest.tsv',
        'cds': root / 'data/processed/cds.tsv', 'queries': validation / 'queries.faa', 'query_manifest': validation / 'queries.tsv',
        'validation_freeze': validation / 'frozen_inputs.json', 'script': Path(__file__), 'protocol': PROTOCOL,
        'engine': ROOT / 'scripts/lysb_profile_candidates.py', 'original_protocol': ROOT / 'docs/lysb_profile_protocol.md'}
    engine.authenticate_receipts({**paths, 'stage_receipts': root / 'results/reports/pilot_stage_receipts.json'})
    input_hashes = {key: sha256(path) for key, path in paths.items()}
    freeze = {'input_sha256': input_hashes, 'builder': engine.BUILDER, 'famsa': engine.FAMSA, 'search': engine.SEARCH,
              'versions': versions, 'broad_reference_ids': sorted(BROAD), 'narrow_reference_ids': sorted(engine.NARROW), 'hmm_workers': 1}
    freeze_path = out / 'frozen_inputs.json'
    if out.exists() and any(out.iterdir()) and not freeze_path.exists():
        raise ValueError('Existing output has no input freeze; choose a new directory')
    if freeze_path.exists() and json.loads(freeze_path.read_text()) != freeze:
        raise ValueError('Changed inputs, code or libraries require a new diagnostic freeze')
    report_path = out / 'profile_report.json'
    if report_path.exists():
        report = json.loads(report_path.read_text())
        if not freeze_path.exists() or any(sha256(out / p) != digest for p, digest in report['artifacts'].items()):
            raise ValueError('Previously frozen model artifacts changed')
        return report
    reference_rows = read_tsv(paths['reference_manifest'])
    references = {row['reference_id']: row for row in reference_rows}
    if len(references) != len(reference_rows): raise ValueError('Duplicate reference identifiers')
    validate_training(references)
    sequences = engine.fasta(paths['references'])
    if set(sequences) != set(references): raise ValueError('Reference FASTA and manifest disagree')
    for identifier, row in references.items():
        if engine.digest(sequences[identifier]) != row['protein_sha256'] or len(sequences[identifier]) != int(row['length_aa']):
            raise ValueError('Reference sequence integrity failure')
    frozen_validation = json.loads(paths['validation_freeze'].read_text())
    for key, field in [('references', 'reference_fasta_sha256'), ('reference_manifest', 'reference_manifest_sha256'),
                       ('queries', 'query_fasta_sha256'), ('query_manifest', 'query_manifest_sha256')]:
        if input_hashes[key] != frozen_validation[field]: raise ValueError('Validation query freeze does not match this panel')
    query_rows = read_tsv(paths['query_manifest'])
    query_sequences = engine.fasta(paths['queries'])
    if {row['query_id'] for row in query_rows} != set(query_sequences) or len(query_rows) != len(query_sequences):
        raise ValueError('Validation query identifiers disagree')
    for row in query_rows:
        sequence = query_sequences[row['query_id']]
        if engine.digest(sequence) != row['protein_sha256'] or len(sequence) != int(row['length_aa']):
            raise ValueError('Validation query sequence integrity failure')
    challenges = [row for row in query_rows if row['query_kind'] == 'specificity_challenge']
    if len(challenges) != 3: raise ValueError('Expected the three fixed cellular-enzyme challenges')
    cds = read_tsv(paths['cds'])
    targets = {row['protein_id']: row['translation'].upper().removesuffix('*') for row in cds}
    if len(targets) != len(cds) or len(cds) != engine.FIXED_Z or any(not seq or set(seq) - set('ACDEFGHIKLMNPQRSTVWYBXZJUO') for seq in targets.values()):
        raise ValueError('Expected all 2460 valid cohort translations')
    budget_check(load_config(ROOT / 'config/analysis.yaml'), additional_bytes=50_000_000)
    groups, edges = engine.reference_groups(BROAD, sequences)
    write_json(freeze_path, freeze)
    write_tsv(out / 'training_references.tsv', [references[i] for i in sorted(BROAD)])
    write_tsv(out / 'reference_group_edges.tsv', edges)
    write_tsv(out / 'reference_groups.tsv', [{'reference_id': i, 'group': groups[i],
        'source': references[i]['resolved_source_accession'], 'sequence_sha256': engine.digest(sequences[i])} for i in sorted(BROAD)])
    metadata = {'REF__' + i: {'reference_id': i, 'expected_category': row['category'], 'query_kind':
        'assay_B_reference' if i in BROAD else 'computational_B_prediction' if i in engine.PREDICTED else 'non_B_curated_reference'}
        for i, row in references.items()}
    validation_sequences = {'REF__' + i: seq for i, seq in sequences.items()}
    for row in challenges:
        validation_sequences[row['query_id']] = query_sequences[row['query_id']]
        metadata[row['query_id']] = {'query_kind': 'specificity_challenge', 'expected_category': row['expected_category']}
    target_metadata = {i: {'query_kind': 'cohort_CDS', 'expected_category': 'unknown_not_assigned'} for i in targets}
    hit_rows, domain_rows, summaries, folds = [], [], [], []
    for panel, ids in [('broad_assay14', BROAD), ('narrow_substrate2', engine.NARROW)]:
        definitions = [('full', sorted(ids), [])]
        for kind, mapping in [('source_out', {i: references[i]['resolved_source_accession'].split('.')[0] for i in ids}),
                              ('family_out', {i: groups[i] for i in ids})]:
            for group in sorted(set(mapping.values())):
                held = sorted(i for i in ids if mapping[i] == group)
                training, _ = engine.holdout_training(ids, held, references, sequences, kind == 'source_out')
                definitions.append((kind + '_' + group, training, held))
        for fold, training, held in definitions:
            label = panel + '__' + fold
            folds.append({'model': label, 'panel': panel, 'fold': fold, 'train_ids': ';'.join(training), 'held_ids': ';'.join(held),
                'excluded_ids': ';'.join(sorted(set(ids) - set(training))), 'n_unique_training': len(engine.unique_training(training, sequences)),
                'status': 'testable' if training else 'untestable_no_training_sequences'})
            hmm = engine.build_model(training, sequences, out / 'models' / label, label, pyfamsa, pyhmmer) if training else None
            if hmm is not None: normalize_metadata(hmm, out / 'models' / label / 'model.hmm')
            scans = [('reference_validation', validation_sequences, metadata), ('cohort', targets, target_metadata)] if not held else [
                ('held_out_reference', {'REF__' + i: sequences[i] for i in held}, metadata)]
            for scan, queries, meta in scans:
                hh, dd = engine.search_model(hmm, queries, label, scan, pyhmmer) if hmm else ([], [])
                hit_rows.extend(hh)
                domain_rows.extend(dd)
                summaries.extend(engine.summarize_queries(label, scan, queries, meta, training, sequences, hh, dd, hmm is None))
    write_tsv(out / 'folds.tsv', folds)
    write_tsv(out / 'full_sequence_hits.tsv', hit_rows, engine.HIT_FIELDS)
    write_tsv(out / 'domain_hits.tsv', domain_rows, engine.DOMAIN_FIELDS)
    write_tsv(out / 'query_summary.tsv', summaries, engine.SUMMARY_FIELDS)
    write_tsv(out / 'cohort_candidates.tsv', [row for row in summaries if row['scan'] == 'cohort' and row['passes_candidate_cuts']], engine.SUMMARY_FIELDS)
    recovery = []
    for fold in folds:
        if not fold['held_ids']: continue
        rows = [row for row in summaries if row['model'] == fold['model'] and row['scan'] == 'held_out_reference']
        recovery.append({'model': fold['model'], 'n_queries': len(rows), 'n_unique_query_sequences': len({row['sequence_sha256'] for row in rows}),
            'n_recovered': sum(row['passes_candidate_cuts'] for row in rows),
            'n_unique_sequences_recovered': len({row['sequence_sha256'] for row in rows if row['passes_candidate_cuts']}), 'status': fold['status']})
    write_tsv(out / 'held_out_recovery.tsv', recovery)
    if input_hashes != {key: sha256(path) for key, path in paths.items()}: raise ValueError('A scientific input changed during scoring')
    report = {'diagnostic_only': True, 'primary_calls_modified': False, 'n_cohort_cds': len(cds), 'fixed_Z': engine.FIXED_Z,
        'fixed_domZ': engine.FIXED_Z, 'n_training_references': 14, 'n_reference_similarity_groups': len(set(groups.values())),
        'n_models': sum(row['status'] == 'testable' for row in folds), 'execution_settings': {'hmm_workers': 1, 'famsa_threads': 4},
        'model_metadata_normalized': ['creation_time', 'command_line'],
        'cohort_candidate_counts': {p + '__full': sum(row['model'] == p + '__full' and row['scan'] == 'cohort' and row['passes_candidate_cuts'] for row in summaries)
                                   for p in ('broad_assay14', 'narrow_substrate2')},
        'cohort_non_training_candidate_counts': {p + '__full': sum(row['model'] == p + '__full' and row['scan'] == 'cohort' and row['passes_candidate_cuts'] and not row['exact_training_alias'] for row in summaries)
                                                for p in ('broad_assay14', 'narrow_substrate2')},
        'cross_family_flags': [row for row in summaries if row['potential_cross_family_flag']],
        'computational_transfer_matches': [row for row in summaries if row['computational_transfer_only'] and row['passes_candidate_cuts']],
        'interpretation': 'Candidate homology only; primary functional confidence, native substrate specificity and evolutionary origins are not inferred.',
        'artifacts': {p.relative_to(out).as_posix(): sha256(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'profile_report.json'}}
    write_json(report_path, report)
    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'validation-dir', 'out'): parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    result = run(args.root, args.validation_dir, args.out)
    print(json.dumps({key: result[key] for key in ('n_models', 'n_training_references', 'cohort_candidate_counts', 'cohort_non_training_candidate_counts')}))
