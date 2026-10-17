#!/usr/bin/env python3
"""Compare fixed-cohort reference panels without changing original pilot outputs."""
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, write_json, write_tsv
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def verify_panel(root, original, baseline_inputs):
    architecture = read_tsv(root / 'results/tables/architectures.tsv')
    if len(architecture) != 20 or len({r['accession'] for r in architecture}) != 20:
        raise ValueError('Expected exactly 20 unique genome architecture rows')
    if set(r['accession'] for r in architecture) != set(original):
        raise ValueError('Panel comparison requires the same frozen genomes')
    inputs = read_tsv(root / 'data/manifests/pilot_inputs.tsv')
    if inputs != baseline_inputs:
        raise ValueError('Amended cohort, sequence hashes or metadata differ from the original freeze')
    if root != ROOT:
        provenance = json.loads((root / 'amendment_provenance.json').read_text())
        if provenance['catalogue_sha256'] != sha256(root / 'references/protein_controls.tsv'):
            raise ValueError('Amended reference catalogue failed provenance check')
        manifest = json.loads((root / 'artifact_manifest.json').read_text())
        for path, digest in manifest['artifacts'].items():
            if not (root / path).is_file() or sha256(root / path) != digest:
                raise ValueError(f'Amended exported artifact changed: {path}')
        receipts = json.loads((root / 'results/reports/pilot_stage_receipts.json').read_text())
        if set(receipts) != {str(i) for i in range(1, 10)}:
            raise ValueError('Amended panel has incomplete stage receipts')
        execution = json.loads((root / 'results/reports/pilot_execution.json').read_text())
        validation = json.loads((root / 'results/reports/validation.json').read_text())
        if execution['status'] != 'success' or validation['checksum_status'] != 'pass':
            raise ValueError('Amended execution or source-integrity audit failed')
    return architecture


def summarize():
    roots = {'original': ROOT, **{name: ROOT / 'results/amendment_v1' / name
                                  for name in ['experimental', 'profile']}}
    summary, changes, proteins = [], [], []
    original = {r['accession']: r for r in read_tsv(ROOT / 'results/tables/architectures.tsv')}
    baseline_inputs = read_tsv(ROOT / 'data/manifests/pilot_inputs.tsv')
    baseline_proteins = {(r['protein_id'], r['category']): r
                        for r in read_tsv(ROOT / 'results/tables/lysis_evidence.tsv')}
    for panel, root in roots.items():
        architecture = verify_panel(root, original, baseline_inputs)
        evidence = read_tsv(root / 'results/tables/lysis_evidence.tsv')
        if len({(r['protein_id'], r['category']) for r in evidence}) != len(evidence):
            raise ValueError('Duplicate protein/category evidence would hide a reconciliation conflict')
        if any(r['confidence'] == 'high' and r['reference_evidence_type'] in
               {'experimental_generic_lipolytic', 'published_profile_prediction'} for r in evidence):
            raise ValueError('Capped reference evidence incorrectly received high confidence')
        supported = [r for r in evidence if r['confidence'] in {'high', 'moderate'}]
        assigned = defaultdict(set)
        for row in supported:
            assigned[row['category']].add(row['accession'])
        categories = Counter(r['architecture'] for r in architecture)
        functional = json.loads((root / 'results/reports/functional_identification.json').read_text())
        evolution = json.loads((root / 'results/reports/evolutionary_tests.json').read_text())
        h2 = evolution['primary_H2']
        summary.append({'panel': panel, 'genomes': len(architecture),
            'supported_lysin_a_genomes': len(assigned['lysin_a']),
            'supported_lysin_b_genomes': len(assigned['lysin_b']),
            'both_a_and_b_detected_genomes': len(assigned['lysin_a'] & assigned['lysin_b']),
            'supported_proteins': len({r['protein_id'] for r in supported}),
            'experimental_reference_supported_proteins': len({r['protein_id'] for r in supported if r['reference_evidence_type'] == 'experimental'}),
            'generic_lipolytic_reference_supported_proteins': len({r['protein_id'] for r in supported if r['reference_evidence_type'] == 'experimental_generic_lipolytic'}),
            'profile_reference_supported_proteins': len({r['protein_id'] for r in supported if r['reference_evidence_type'] == 'published_profile_prediction'}),
            **{f'architecture_{key}': categories[key] for key in ['closely_linked', 'separated', 'complex', 'unresolved']},
            'primary_h2_lineages': h2['n_lineages'], 'primary_h2_pairs': h2['n_pairs'],
            'primary_h2_status': h2['status'], 'effect_mean': h2['effect_mean'],
            'ci_low': h2['ci_low'], 'ci_high': h2['ci_high'], 'p_value': h2['p_value'],
            'evidence_report_sha256': sha256(root / 'results/reports/functional_identification.json')})
        for row in architecture:
            changes.append({'panel': panel, 'accession': row['accession'], 'phage': row['phage'],
                'cluster': row['cluster'], 'original_architecture': original[row['accession']]['architecture'],
                'revised_architecture': row['architecture'], 'lysin_a_count': row['lysin_a_count'],
                'lysin_b_count': row['lysin_b_count'], 'intervening_cds': row['intervening_cds'],
                'gap_bp': row['gap_bp'], 'classification_reason': row['classification_reason']})
        for row in supported:
            previous = baseline_proteins.get((row['protein_id'], row['category']), {})
            proteins.append({'panel': panel, **{key: row[key] for key in ['accession', 'protein_id',
                'protein_accession', 'category', 'confidence', 'reference_id', 'reference_evidence_type',
                'identity', 'query_coverage', 'target_coverage', 'evalue', 'evidence_url']},
                'original_confidence': previous.get('confidence', 'not_detected'),
                'original_reference_id': previous.get('reference_id', ''),
                'original_reference_evidence_type': previous.get('reference_evidence_type', ''),
                'newly_supported': str(previous.get('confidence') not in {'high', 'moderate'}).lower()})
        current = {(r['protein_id'], r['category']): r for r in evidence}
        for key, previous in baseline_proteins.items():
            if previous['confidence'] not in {'high', 'moderate'} or current.get(key, {}).get('confidence') in {'high', 'moderate'}:
                continue
            row = current.get(key, previous)
            proteins.append({'panel': panel, **{field: row[field] for field in ['accession', 'protein_id',
                'protein_accession', 'category', 'reference_id', 'reference_evidence_type', 'identity',
                'query_coverage', 'target_coverage', 'evalue', 'evidence_url']},
                'confidence': current.get(key, {}).get('confidence', 'not_detected'),
                'original_confidence': previous['confidence'],
                'original_reference_id': previous['reference_id'],
                'original_reference_evidence_type': previous['reference_evidence_type'],
                'newly_supported': 'false'})
    output = ROOT / 'results/amendment_v1/comparison'
    write_tsv(output / 'panel_summary.tsv', summary)
    write_tsv(output / 'architecture_changes.tsv', changes)
    write_tsv(output / 'supported_protein_changes.tsv', proteins)
    labels = ['Original', 'Expanded\nassay panel', 'Assay panel +\npublished profiles']
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), layout='constrained')
    positions = list(range(len(summary)))
    for offset, field, label, color in [(-.24, 'supported_lysin_a_genomes', 'Lysin A', '#3b7c99'),
        (0, 'supported_lysin_b_genomes', 'Lysin B', '#d99b45'),
        (.24, 'both_a_and_b_detected_genomes', 'Both detected', '#7a9d63')]:
        axes[0].bar([x + offset for x in positions], [r[field] for r in summary], width=.23, label=label, color=color)
    axes[0].set(ylabel='Genomes with supported assignment', ylim=(0, 21), title='A. Homology-assignment coverage')
    axes[0].set_yticks(range(0, 21, 5)); axes[0].legend(frameon=False, fontsize=9)
    bottom = [0] * len(summary)
    for key, label, color in [('closely_linked', 'Linked', '#3b7c99'), ('separated', 'Separated', '#d99b45'),
        ('complex', 'Complex', '#967eaf'), ('unresolved', 'Unresolved', '#c7c7c7')]:
        values = [r[f'architecture_{key}'] for r in summary]
        axes[1].bar(positions, values, bottom=bottom, label=label, color=color, width=.62)
        bottom = [a + b for a, b in zip(bottom, values)]
    axes[1].set(ylabel='Genomes', ylim=(0, 24), title='B. Unique-pair architecture classification')
    axes[1].set_yticks(range(0, 21, 5)); axes[1].legend(frameon=False, ncols=2, fontsize=8, loc='upper center')
    for ax in axes:
        ax.set_xticks(positions, labels); ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Post-pilot reference expansion: same 20 genomes\nProfile-supported assignments remain computational', fontsize=12)
    for extension in ['png', 'pdf']:
        metadata = {'Software': 'actinolysis'} if extension == 'png' else {'CreationDate': None, 'ModDate': None}
        fig.savefig(output / f'reference_expansion.{extension}', dpi=300, metadata=metadata)
    plt.close(fig)
    files = sorted(p for p in output.iterdir() if p.is_file() and p.name != 'comparison_manifest.json')
    write_json(output / 'comparison_manifest.json', {'amendment': 'exploratory post-pilot analysis',
        'summary_script_sha256': sha256(Path(__file__)),
        'functional_thresholds': 'unchanged identity0.40, both coverages0.70, E<=1e-6',
        'scientific_interpretation': 'Expanded sequence homology to published predictions does not establish experimental activity or independent evolutionary origins.',
        'artifacts': {p.name: sha256(p) for p in files}})
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    summarize()
