#!/usr/bin/env python3
"""Compute descriptive annotation-conditioned contrasts without primary promotion."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis import annotation_sensitivity as sensitivity
from actinolysis import architecture, neighborhoods
from actinolysis.common import budget_check, load_config, read_tsv, sha256, write_json, write_tsv

def run(root, out):
    root, out = Path(root).resolve(), Path(out).resolve()
    if not out.is_relative_to(ROOT) or out == root or root.is_relative_to(out): raise ValueError('Use a separate output inside the project')
    if any(out.is_relative_to(root / p) for p in ('data','src','scripts','references','figures','results/tables','results/reports')):
        raise ValueError('Output overlaps protected scientific inputs')
    paths = {'genomes': root / 'data/processed/genomes.tsv', 'cds': root / 'data/processed/cds.tsv',
        'evidence': root / 'data/processed/lysis_evidence.tsv', 'memberships': root / 'results/tables/family_memberships.tsv',
        'architectures': root / 'results/tables/architecture_sensitivity.tsv', 'primary_pairs': root / 'results/tables/neighborhood_pairs.tsv',
        'configuration': root / 'config/analysis.yaml', 'protocol': ROOT / 'docs/annotation_neighborhood_protocol.md',
        'module': Path(sensitivity.__file__), 'architecture_code': Path(architecture.__file__),
        'matching_code': Path(neighborhoods.__file__), 'script': Path(__file__)}
    hashes = {key: sha256(path) for key,path in paths.items()}
    genomes, cds = read_tsv(paths['genomes']), read_tsv(paths['cds'])
    if (len(genomes), len(cds)) != (20, 2460): raise ValueError('This sensitivity requires the fixed 20-genome, 2460-CDS cohort')
    published = json.loads((ROOT / 'results/lysterase_update/panel/analysis_manifest.json').read_text())['scientific_sha256']
    if hashes['memberships'] != published['results/tables/family_memberships.tsv']:
        raise ValueError('Protein-family membership differs from the frozen cohort analysis')
    receipts = json.loads((root / 'results/reports/pilot_stage_receipts.json').read_text())
    for stage, key, name in [('2','genomes','data/processed/genomes.tsv'), ('2','cds','data/processed/cds.tsv'),
        ('4','evidence','data/processed/lysis_evidence.tsv'), ('5','architectures','results/tables/architecture_sensitivity.tsv'),
        ('6','primary_pairs','results/tables/neighborhood_pairs.tsv')]:
        if receipts[stage]['outputs'].get(name) != hashes[key]: raise ValueError('Scientific input does not match its pipeline receipt')
    frozen = {'input_sha256': hashes, 'evidence_scope': sensitivity.SCOPE, 'primary_calls_modified': False,
              'formal_inference_performed': False}
    marker = out / 'frozen_inputs.json'
    if out.exists() and any(out.iterdir()) and not marker.exists(): raise ValueError('Existing output lacks a freeze')
    if marker.exists() and json.loads(marker.read_text()) != frozen: raise ValueError('Changed inputs require a separate output freeze')
    report_path = out / 'sensitivity_report.json'
    if report_path.exists():
        report = json.loads(report_path.read_text())
        if any(sha256(out / name) != digest for name,digest in report['artifacts'].items()): raise ValueError('Frozen output changed')
        return report
    cfg = load_config(paths['configuration'])
    budget_check(load_config(ROOT / 'config/analysis.yaml'), additional_bytes=20_000_000)
    write_json(marker, frozen)
    pairs, controls, neighbors = sensitivity.analyze(cfg, genomes, cds, *(read_tsv(paths[key]) for key in
        ('evidence', 'memberships', 'architectures')))
    summary = sensitivity.summarize(pairs)
    write_tsv(out / 'conditional_pair_scores.tsv', pairs, sensitivity.PAIR_FIELDS)
    write_tsv(out / 'conditional_controls.tsv', controls, ['evidence_scope', 'pair_id', 'family_identity', 'window'] +
              [name for name in neighborhoods.CONTROL_FIELDS if name not in {'pair_id','family_identity','window','confidence_set'}])
    write_tsv(out / 'conditional_neighborhoods.tsv', neighbors,
        ['evidence_scope','family_identity','window','accession','anchor','anchor_family','neighbor_proteins','neighbor_families','signed_offsets','relative_strands'])
    write_tsv(out / 'conditional_cluster_summary.tsv', summary)
    def base(row):
        return float(row['family_identity']) == float(cfg['neighborhood']['primary_identity']) and int(row['window']) == int(cfg['neighborhood']['primary_window']) and int(row['architecture_threshold']) == int(cfg['architecture']['separation_threshold']) and row['genome_subset'] == 'all' and row['metric'] == 'jaccard'
    conditional = [row for row in pairs if base(row)]
    primary = [row for row in read_tsv(paths['primary_pairs']) if base(row) and row['confidence_set'] == 'high_moderate']
    comparative = []
    for scope, rows in [('independent_high_moderate', primary), (sensitivity.SCOPE, conditional)]:
        for row in rows:
            comparative.append({'evidence_scope': scope, **{key: row[key] for key in ('pair_id','cluster','pair_architecture','target_score','control_mean','delta','n_controls','status','reason')}})
    write_tsv(out / 'base_scope_comparison.tsv', comparative)
    if hashes != {key: sha256(path) for key,path in paths.items()}: raise ValueError('Scientific inputs changed during sensitivity analysis')
    eligible = [row for row in conditional if row['status'] == 'descriptive_eligible']
    report = {**frozen, 'n_genomes': len(read_tsv(paths['genomes'])), 'n_base_pairs': len(conditional),
        'base_available_pairs': len(eligible), 'base_available_cluster_labels': len({row['cluster'] for row in eligible}),
        'base_separated_available_pairs': sum(row['pair_architecture'] == 'separated' for row in eligible),
        'base_exclusion_counts': dict(Counter(row['reason'] for row in conditional if row['reason'])),
        'interpretation': 'Descriptive contrasts conditional on deposited labels; these labels may use context and do not independently validate lysis function or the primary H2 hypothesis.',
        'artifacts': {p.relative_to(out).as_posix(): sha256(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'sensitivity_report.json'}}
    write_json(report_path, report)
    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    result = run(args.root, args.out)
    print(json.dumps({key: result[key] for key in ('n_base_pairs','base_available_pairs','base_available_cluster_labels','base_separated_available_pairs','base_exclusion_counts')}))
