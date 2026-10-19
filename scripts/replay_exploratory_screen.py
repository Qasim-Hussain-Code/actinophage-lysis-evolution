#!/usr/bin/env python3
"""Replay the fixed exploratory screen from verified caches, then optionally export science."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
from actinolysis.architecture import ARCH_FIELDS, PAIR_FIELDS
from actinolysis.common import NA, read_tsv, sha256, write_json, write_tsv
from scripts import run_exploratory_screen as runner
from scripts.archive_cohort_hits import (HIT_OUTPUT_FIELDS, STATUS_FIELDS, build_tables,
                                         read_hits, serialize_booleans, unique_rows)
from scripts.prepare_exploratory_screen import new_workspace

PUBLIC = Path('results/exploratory50_v1')
RETAINED_VALIDATION_SHA256 = '78374ce14082036a0d035025c3a954217a224c8659c7638e6310e904e795a9a7'
SCIENCE = sorted({'data/manifests/intermediate_inputs.tsv', 'data/manifests/acquisition_ledger.tsv'} | {
    name for paths in runner.OUTPUTS.values() for name in paths
    if name.endswith(('.tsv', '.faa'))} | {
    'results/reports/quality_control.json', 'results/reports/functional_identification.json'})
QUERY = ['protein_id', 'protein_accession', 'accession', 'cds_index', 'start', 'end',
         'location_parts', 'strand', 'partial']
ALIAS_FIELDS = ['reference_source_accession_match', 'source_reference_ids',
                'exact_reference_sequence_alias', 'exact_reference_ids', 'protein_sha256']
GENOME_FIELDS = ['accession', 'phage', 'host_genus', 'cluster', 'subcluster', 'inclusion_status',
    'exclusion_reason', 'length', 'gc_percent', 'cds_count', 'topology', 'sequence_sha256',
    'similarity_group', 'representative']
EVIDENCE_FIELDS = ['protein_id', 'protein_accession', 'accession', 'start', 'end', 'strand',
    'partial', 'category', 'confidence', 'annotation_category', 'source', 'reference_id',
    'identity', 'query_coverage', 'target_coverage', 'evalue', 'reference_evidence_type',
    'evidence_url', 'context_used']


def project(rows, fields):
    return [{name: row.get(name, NA) for name in fields} for row in rows]


def verify_stage_receipts(runtime):
    """Require all five complete stages and authenticate every declared output."""
    inputs_sha = sha256(runtime / 'execution_inputs.json')
    execution = json.loads((runtime / 'results/reports/exploratory_execution.json').read_text(encoding='utf-8'))
    if (execution.get('status') != 'exploratory_complete_awaiting_replay'
            or execution.get('completed_stages') != [1, 2, 3, 4, 5]
            or execution.get('h2_analysis_performed') is not False):
        raise ValueError('A completed five-stage exploratory execution is required')
    receipts = json.loads((runtime / 'results/reports/exploratory_stage_receipts.json').read_text(encoding='utf-8'))
    if set(receipts) != {'1', '2', '3', '4', '5'}:
        raise ValueError('Exactly five exploratory stage receipts are required')
    checked = {}
    for number in range(1, 6):
        receipt = receipts[str(number)]
        if (receipt.get('status') != 'exploratory_stage_complete'
                or receipt.get('execution_inputs_sha256') != inputs_sha
                or set(receipt.get('outputs', {})) != set(runner.stage_outputs(runtime, number))):
            raise ValueError(f'Invalid exploratory stage {number} receipt')
        for name, expected in receipt['outputs'].items():
            if sha256(runner.inside(runtime, name)) != expected:
                raise ValueError(f'Exploratory stage output checksum mismatch: {name}')
            if name in checked and checked[name] != expected:
                raise ValueError('Exploratory stage receipts disagree')
            checked[name] = expected
    return checked


def verified_run(root, workspace):
    prepared, plan, ledger = runner.verify_preparation(root, workspace)
    runtime = runner.inside(prepared, 'runtime')
    inputs, cfg = runner.verify_runtime(root, prepared, runtime, plan)
    verify_stage_receipts(runtime)
    runner.check_acquisition(runtime, ledger)
    runner.check_controls(root, runtime, plan)
    return prepared, runtime, plan, ledger, inputs, cfg


def compare_science(original, replay):
    comparisons = [{'path': name, 'original_sha256': sha256(original / name),
                    'replay_sha256': sha256(replay / name)} for name in SCIENCE]
    if any(row['original_sha256'] != row['replay_sha256'] for row in comparisons):
        raise ValueError('Exploratory scientific replay differs; public export is withheld')
    return comparisons


def cache_records(original):
    manifest = read_tsv(original / 'data/manifests/intermediate_inputs.tsv')
    if len(manifest) != 50:
        raise ValueError('Replay requires all 50 frozen selected source records')
    for row in manifest:
        if sha256(runner.inside(original, row['path'])) != row['sha256']:
            raise ValueError('Exploratory source cache checksum mismatch')
    return manifest


def make_replay(root, workspace, replay_workspace):
    prepared, original, plan, ledger, _, cfg = verified_run(root, workspace)
    manifest = cache_records(original)
    destination = new_workspace(root, replay_workspace)
    destination.mkdir(parents=True, exist_ok=False)
    for name in ('preparation.json', 'selection.tsv'):
        shutil.copyfile(prepared / name, destination / name)
    runtime = runner.prepare_runtime(root, destination, plan, ledger, offline=True)
    runner.root_guard(cfg, root, sum(runner.inside(original, row['path']).stat().st_size for row in manifest))
    for row in manifest:
        target = runner.inside(runtime, row['path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(runner.inside(original, row['path']), target)
    shutil.copyfile(original / 'data/manifests/intermediate_inputs.tsv',
                    runtime / 'data/manifests/intermediate_inputs.tsv')
    shutil.copyfile(original / 'data/manifests/acquisition_ledger.tsv',
                    runtime / 'data/manifests/acquisition_ledger.tsv')
    # No generated proteome, score, search database or temporary search state is reused.
    if (runtime / 'data/processed').exists() or (runtime / '.tmp/reference_search').exists():
        raise ValueError('Replay requires an uncached reference-search workspace')
    write_json(destination / 'replay_origin.json', {
        'method': 'checksum_verified_cache_only', 'original_preparation_sha256': sha256(prepared / 'preparation.json'),
        'original_execution_inputs_sha256': sha256(original / 'execution_inputs.json'),
        'frozen_source_manifest_sha256': sha256(original / 'data/manifests/intermediate_inputs.tsv'),
        'replay_driver_sha256': sha256(__file__), 'reference_search_cache_reused': False})
    subprocess.run([sys.executable, '-B', str(runtime / 'scripts/run_exploratory_screen.py'),
        '--worker', '--source-root', str(root), '--workspace', replay_workspace, '--resume-from', '1'],
        cwd=runtime, check=True)


def alias_flags(cds, references):
    flags = {}
    for pid, query in cds.items():
        digest = hashlib.sha256(query.get('translation', '').encode()).hexdigest() if query.get('translation') else NA
        source = sorted(rid for rid, ref in references.items()
            if ref['resolved_source_accession'].split('.')[0] == query['accession'].split('.')[0])
        exact = sorted(rid for rid, ref in references.items() if digest != NA and ref['protein_sha256'] == digest)
        flags[pid] = {'reference_source_accession_match': bool(source), 'source_reference_ids': ';'.join(source) or NA,
            'exact_reference_sequence_alias': bool(exact), 'exact_reference_ids': ';'.join(exact) or NA,
            'protein_sha256': digest}
    return flags


def group_summary(selected, genomes, cds, evidence, architectures):
    support = defaultdict(set)
    for row in evidence:
        if row['source'] == 'curated_reference_homology' and row['confidence'] in {'high', 'moderate'}:
            support[row['accession']].add(row['category'])
    cds_counts = Counter(row['accession'] for row in cds)
    groups = [('all', 'all', selected)]
    for key in ('host_genus', 'cluster'):
        for value in sorted({row.get(key, NA) for row in selected}):
            groups.append((key, value, [row for row in selected if row.get(key, NA) == value]))
    output = []
    for group_type, group, members in groups:
        accessions = {row['accession'] for row in members}
        included = accessions & {row['accession'] for row in genomes}
        counts = Counter(row['architecture'] for row in architectures if row['accession'] in accessions)
        output.append({'group_type': group_type, 'group': group, 'selected_genomes': len(members),
            'included_genomes': len(included), 'cds': sum(cds_counts[a] for a in accessions),
            'supported_A_genomes': sum('lysin_a' in support[a] for a in included),
            'supported_B_genomes': sum('lysin_b' in support[a] for a in included),
            'supported_A_and_B_genomes': sum({'lysin_a', 'lysin_b'} <= support[a] for a in included),
            **{name: counts[name] for name in ('closely_linked', 'separated', 'complex', 'unresolved')}})
    return output


def public_tables(runtime, ledger, cfg):
    """Build scientific tables with explicit column whitelists, before writing any export."""
    manifest = read_tsv(runtime / 'data/manifests/intermediate_inputs.tsv')
    pilot_by_base = {row['accession'].split('.')[0]: row['pilot'] for row in ledger}
    selected = [{**row, 'rank': index, 'pilot': pilot_by_base[row['accession'].split('.')[0]]}
                for index, row in enumerate(manifest, 1)]
    genomes = read_tsv(runtime / 'data/processed/genomes.tsv')
    evidence = read_tsv(runtime / 'data/processed/lysis_evidence.tsv')
    architectures = read_tsv(runtime / 'results/tables/architectures.tsv')
    cds = unique_rows(runtime / 'data/processed/cds.tsv', 'protein_id')
    references = unique_rows(runtime / 'data/manifests/reference_manifest.tsv', 'reference_id')
    records = list(SeqIO.parse(runtime / 'data/processed/proteins.faa', 'fasta'))
    searched = {record.id for record in records}
    if len(searched) != len(records) or not searched <= cds.keys():
        raise ValueError('Duplicate or unknown searchable protein identifier')
    hits = read_hits(runtime / 'data/processed/reference_hits.tsv', cds, references, searched)
    hit_rows, statuses = build_tables(cds, references, hits, evidence, searched, cfg)
    flags = alias_flags(cds, references)
    for row in hit_rows + statuses:
        row.update(flags[row['protein_id']])
    for row in hit_rows:
        ref = references[row['reference_id']]
        row['same_reference_source_accession'] = ref['resolved_source_accession'].split('.')[0] == row['accession'].split('.')[0]
        row['exact_reference_sequence_match'] = row['protein_sha256'] == ref['protein_sha256']
    tables = {}
    def add(name, rows, fields):
        tables[name] = (project(serialize_booleans(rows), fields), fields)
    add('selected_genomes.tsv', selected,
        ['rank', 'accession', 'source_url', 'sha256', 'phage', 'host_genus', 'cluster', 'subcluster', 'pilot'])
    add('genome_qc.tsv', genomes, GENOME_FIELDS)
    add('exclusions.tsv', read_tsv(runtime / 'results/tables/exclusions.tsv'), GENOME_FIELDS[:7])
    add('near_identical_pairs.tsv', read_tsv(runtime / 'results/tables/near_identical_pairs.tsv'),
        ['accession_a', 'accession_b', 'sketch_jaccard', 'threshold', 'k', 'sketch_size'])
    add('architectures.tsv', architectures, ARCH_FIELDS)
    add('architecture_pairs.tsv', read_tsv(runtime / 'results/tables/architecture_pairs.tsv'), PAIR_FIELDS)
    add('architecture_sensitivity.tsv', read_tsv(runtime / 'results/tables/architecture_sensitivity.tsv'), ARCH_FIELDS)
    add('annotation_coverage.tsv', read_tsv(runtime / 'results/tables/annotation_coverage.tsv'),
        ['accession', 'phage', 'host_genus', 'cluster', 'cds', 'with_specific_annotation',
         'supported_lysis_proteins', 'annotation_only_lysis_proteins', 'hypothetical_proteins'])
    add('lysis_evidence.tsv', evidence, EVIDENCE_FIELDS)
    add('all_reference_hits.tsv', hit_rows,
        [name for name in HIT_OUTPUT_FIELDS if name not in {'product', 'note'}] + ALIAS_FIELDS +
        ['same_reference_source_accession', 'exact_reference_sequence_match'])
    add('protein_detection_status.tsv', statuses,
        [name for name in STATUS_FIELDS if name not in {'product', 'note'}] + ALIAS_FIELDS)
    add('references.tsv', list(references.values()), ['reference_id', 'category', 'evidence_type',
        'resolved_source_accession', 'resolved_protein_accession', 'source_url', 'reference_url',
        'source_sha256', 'protein_sha256', 'length_aa'])
    groups = group_summary(selected, genomes, list(cds.values()), evidence, architectures)
    add('coverage_by_group.tsv', groups, list(groups[0]))
    tiers = []
    for tier, confidence in sorted({(r['reference_evidence_type'], r['confidence']) for r in evidence}):
        local = [r for r in evidence if (r['reference_evidence_type'], r['confidence']) == (tier, confidence)]
        tiers.append({'evidence_type': tier, 'confidence': confidence, 'assignment_rows': len(local),
                      'unique_proteins': len({r['protein_id'] for r in local}),
                      'unique_genomes': len({r['accession'] for r in local})})
    add('coverage_by_evidence_tier.tsv', tiers,
        ['evidence_type', 'confidence', 'assignment_rows', 'unique_proteins', 'unique_genomes'])
    summary = {'analysis': 'exploratory_50_genome_detection_screen', 'scope': 'exploratory_only',
        'counts': {**groups[0], 'searched_proteins': len(searched), 'reported_reference_hits': len(hit_rows),
                  'exact_reference_sequence_alias_proteins': sum(f['exact_reference_sequence_alias'] for f in flags.values()),
                  'reference_source_accession_match_genomes': len({row['accession'] for pid, row in cds.items()
                      if flags[pid]['reference_source_accession_match']})},
        'detection_status_counts': dict(sorted(Counter(r['detection_status'] for r in statuses).items())),
        'identification_thresholds': cfg['identification'], 'h2_analysis_performed': False,
        'primary_cohort_modified': False, 'main_expansion_gate_modified': False,
        'interpretation': [
            'Detection yields are conditional on the selected cohort, frozen references and unchanged evidence tiers.',
            'Exact reference aliases and source-genome matches are identified, not counted as independent functional validation.',
            'Source-accession flags compare resolved accession bases; cross-accession biological aliases can remain unrecognized.',
            'No reported reference hit means unresolved, not gene absence. Annotation-only labels retain low confidence.',
            'Cluster labels are coarse metadata groups. These counts do not establish independent lineage replication.',
            'No H2 estimate, p-value, population prevalence or evolutionary mechanism is inferred; no automatic 300-genome progression.']}
    return tables, summary


def retained_panel_validation(root, runtime, cfg):
    """Authenticate an existing fixed-panel control result; do not claim a rerun."""
    base = root / 'results/lysterase_update/reference_validation'
    summary_path = base / 'validation_summary.json'
    if sha256(summary_path) != RETAINED_VALIDATION_SHA256:
        raise ValueError('Retained reference-panel validation summary changed')
    summary = json.loads(summary_path.read_text(encoding='utf-8'))
    for name, expected in summary['artifacts'].items():
        if sha256(runner.inside(base, name)) != expected:
            raise ValueError('Retained reference-panel validation artifact changed')
    frozen = json.loads((base / 'frozen_inputs.json').read_text(encoding='utf-8'))
    if (summary['status'] != 'success' or summary['reference_count'] != 39
            or summary['identification_thresholds'] != cfg['identification']
            or frozen['reference_manifest_sha256'] != sha256(runtime / 'data/manifests/reference_manifest.tsv')
            or frozen['reference_fasta_sha256'] != sha256(runtime / 'data/processed/references.faa')
            or frozen['source_sha256']['identification'] != sha256(runtime / 'src/actinolysis/identification.py')
            or frozen['identification_thresholds'] != cfg['identification']
            or frozen['search_settings'] != summary['search_settings']
            or frozen['search_settings']['threads'] != cfg['threads']
            or frozen['search_settings']['split_memory_limit'] != cfg['resources']['mmseqs_memory']):
        raise ValueError('Retained validation does not match the frozen screen panel/settings')
    return {'scope': 'retained_unchanged_panel_validation_not_rerun',
        'source_report_sha256': RETAINED_VALIDATION_SHA256,
        'frozen_inputs_sha256': summary['artifacts']['frozen_inputs.json'],
        'reference_count': 39, 'search_settings': summary['search_settings'],
        'specificity_challenges': summary['specificity_challenges'],
        'interpretation': 'Three selected cellular-enzyme challenges; no population specificity estimate or proven native-substrate inactivity.'}


def replay_and_export(workspace, replay_workspace, *, verify_only=False, export=False, root=ROOT):
    root = Path(root).resolve()
    if not verify_only:
        make_replay(root, workspace, replay_workspace)
    prepared, original, plan, ledger, _, cfg = verified_run(root, workspace)
    replay_prepared, replay, replay_plan, _, replay_inputs, replay_cfg = verified_run(root, replay_workspace)
    origin = json.loads((replay_prepared / 'replay_origin.json').read_text(encoding='utf-8'))
    if (origin.get('method') != 'checksum_verified_cache_only'
            or origin['original_preparation_sha256'] != sha256(prepared / 'preparation.json')
            or origin['original_execution_inputs_sha256'] != sha256(original / 'execution_inputs.json')
            or origin['frozen_source_manifest_sha256'] != sha256(original / 'data/manifests/intermediate_inputs.tsv')
            or origin['replay_driver_sha256'] != sha256(__file__)
            or origin.get('reference_search_cache_reused') is not False
            or replay_plan != plan or replay_cfg.get('offline') is not True
            or replay_inputs['planned_configuration_sha256'] != plan['effective_configuration_sha256']):
        raise ValueError('Replay provenance or frozen configuration differs')
    comparisons = compare_science(original, replay)
    controls = retained_panel_validation(root, original, cfg)
    report = {'status': 'pass', 'method': 'checksum_verified_cache_only', 'completed_stages': [1, 2, 3, 4, 5],
        'compared_artifacts': len(comparisons), 'artifacts': comparisons, 'reference_search_cache_reused': False,
        'source_acquisition_scope': 'Reused checksum-verified frozen source records; no fresh retrieval.',
        'software_scope': 'Reused the same pinned software and binary; no independent reinstall.',
        'protocol_sha256': plan['protocol_sha256'], 'reference_panel_sha256': plan['panel_files_sha256'],
        'frozen_provenance': {
            'configuration_file_sha256': plan['configuration_file_sha256'],
            'effective_configuration_sha256': plan['effective_configuration_sha256'],
            'private_selection_ledger_sha256': plan['selection_sha256'],
            'execution_runner_sha256': plan.get('runner_sha256', sha256(root / 'scripts/run_exploratory_screen.py')),
            'preparation_code_sha256': plan['code_sha256'],
            'execution_source_tree_sha256': plan['source_tree_sha256']},
        'replay_driver_sha256': sha256(__file__), 'retained_panel_validation': controls}
    write_json(replay_prepared / 'replay_validation.json', report)
    if export:
        target = runner.inside(root, PUBLIC)
        if target.exists():
            raise ValueError('Preserve the existing public freeze; choose a new analysis version')
        tables, summary = public_tables(original, ledger, cfg)
        runner.root_guard(cfg, root, 20_000_000)
        # Recheck the frozen source and replay receipts before writing the public record.
        verified_run(root, workspace)
        verified_run(root, replay_workspace)
        compare_science(original, replay)
        target.mkdir(parents=True, exist_ok=False)
        for name, (rows, fields) in tables.items():
            write_tsv(target / name, rows, fields)
        write_json(target / 'summary.json', summary)
        write_json(target / 'replay_validation.json', {**report, 'public_cohort_manifest': {
            'artifact': 'selected_genomes.tsv', 'sha256': sha256(target / 'selected_genomes.tsv')}})
        qc = json.loads((original / 'results/reports/quality_control.json').read_text(encoding='utf-8'))
        write_json(target / 'quality_control.json', {k: qc[k] for k in
            ('selected', 'included', 'excluded', 'cds', 'translation_status', 'near_identical_screen')})
        write_json(target / 'artifact_manifest.json', {'analysis': 'exploratory50_v1',
            'scope': 'exploratory_only', 'replay_identical': True,
            'export_driver_sha256': sha256(__file__),
            'hit_ledger_builder_sha256': sha256(root / 'scripts/archive_cohort_hits.py'),
            'artifacts': {p.name: sha256(p) for p in sorted(target.iterdir()) if p.is_file()}})
        runner.root_guard(cfg, root)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', required=True, help='Completed original preparation workspace')
    parser.add_argument('--replay-workspace', required=True, help='Separate private replay workspace')
    parser.add_argument('--verify-only', action='store_true', help='Verify a previously completed cache-only replay')
    parser.add_argument('--export', action='store_true', help='Export scientific tables only after identical replay')
    args = parser.parse_args()
    report = replay_and_export(args.workspace, args.replay_workspace, verify_only=args.verify_only, export=args.export)
    print(json.dumps({'status': report['status'], 'compared_artifacts': report['compared_artifacts'],
                      'public_scientific_export_performed': args.export}))


if __name__ == '__main__':
    main()
