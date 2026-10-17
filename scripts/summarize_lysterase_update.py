#!/usr/bin/env python3
"""Compare functional coverage and test eligibility in the fixed genome cohort."""
from __future__ import annotations
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, write_json, write_tsv

OLD = ROOT / 'results/amendment_v1/profile'
NEW = ROOT / 'results/lysterase_update/panel'
OUT = ROOT / 'results/lysterase_update/comparison'
ACCEPTED = {'high', 'moderate'}

def measure(root):
    genomes = read_tsv(root / 'results/tables/dataset_manifest.tsv')
    evidence = read_tsv(root / 'results/tables/lysis_evidence.tsv')
    architectures = read_tsv(root / 'results/tables/architectures.tsv')
    a = {r['accession'] for r in evidence if r['category'] == 'lysin_a' and r['confidence'] in ACCEPTED}
    b = {r['accession'] for r in evidence if r['category'] == 'lysin_b' and r['confidence'] in ACCEPTED}
    primary = [r for r in read_tsv(root / 'results/tables/evolutionary_tests.tsv') if r['primary'] == 'true']
    if len(primary) != 1: raise ValueError('Expected exactly one primary H2 result')
    values = {'genomes': len(genomes), 'cds': sum(int(r['cds_count']) for r in genomes),
        'references': len(read_tsv(root / 'data/manifests/reference_manifest.tsv')),
        'supported_a_genomes': len(a), 'supported_b_genomes': len(b), 'supported_ab_genomes': len(a & b),
        'supported_proteins': len({r['protein_id'] for r in evidence if r['confidence'] in ACCEPTED}),
        'primary_h2_pairs': int(primary[0]['n_pairs']), 'primary_h2_lineages': int(primary[0]['n_lineages'])}
    values.update({name: sum(r['architecture'] == name for r in architectures)
                   for name in ('closely_linked', 'separated', 'complex', 'unresolved')})
    return values, {r['accession']: r for r in genomes}, {(r['protein_id'], r['category']): r for r in evidence}, a, b

def run():
    old, old_genomes, old_evidence, old_a, old_b = measure(OLD)
    new, genomes, evidence, a, b = measure(NEW)
    for filename in ('data/manifests/pilot_inputs.tsv', 'results/tables/family_memberships.tsv'):
        if (OLD / filename).read_bytes() != (NEW / filename).read_bytes():
            raise ValueError('Cohort or protein-family memberships changed')
    if (new['genomes'], new['cds'], old['genomes'], old['cds']) != (20, 2460, 20, 2460):
        raise ValueError('Expected the unchanged 20-genome, 2460-CDS cohort')
    if a != old_a or not old_b <= b: raise ValueError('An earlier supported category was lost')
    rows = [{'panel': label, **values} for label, values in [('combined_34', old), ('assay_updated_39', new)]]
    write_tsv(OUT / 'panel_summary.tsv', rows)
    changed = []
    sequences = {r['protein_id']: r['translation'] for r in read_tsv(ROOT / '.experiments/lysterase_2026/data/processed/cds.tsv')}
    for key in sorted(set(old_evidence) | set(evidence)):
        prior, current = old_evidence.get(key, {}), evidence.get(key, {})
        if prior.get('confidence', 'unresolved') == current.get('confidence', 'unresolved'): continue
        accession = current.get('accession', prior.get('accession'))
        changed.append({'phage': genomes[accession]['phage'], 'cluster': genomes[accession]['cluster'],
            'protein_id': key[0], 'category': key[1], 'previous_confidence': prior.get('confidence', 'unresolved'),
            'updated_confidence': current.get('confidence', 'unresolved'), 'reference_id': current.get('reference_id', ''),
            'reference_evidence_type': current.get('reference_evidence_type', ''),
            'identity': current.get('identity', ''), 'query_coverage': current.get('query_coverage', ''),
            'target_coverage': current.get('target_coverage', ''), 'evalue': current.get('evalue', ''),
            'protein_sha256': hashlib.sha256(sequences[key[0]].encode()).hexdigest()})
    fields = ['phage', 'cluster', 'protein_id', 'category', 'previous_confidence', 'updated_confidence', 'reference_id',
              'reference_evidence_type', 'identity', 'query_coverage', 'target_coverage', 'evalue', 'protein_sha256']
    write_tsv(OUT / 'confidence_changes.tsv', changed, fields)
    coverage = [{'phage': row['phage'], 'accession': accession, 'host_genus': row['host_genus'], 'cluster': row['cluster'],
                 'a_before': int(accession in old_a), 'b_before': int(accession in old_b),
                 'a_after': int(accession in a), 'b_after': int(accession in b)}
                for accession, row in sorted(genomes.items(), key=lambda item: (item[1]['host_genus'], item[1]['cluster'], item[1]['phage']))]
    write_tsv(OUT / 'genome_coverage.tsv', coverage)
    pairs = read_tsv(NEW / 'results/tables/neighborhood_pairs.tsv')
    selected = [r for r in pairs if r['family_identity'] == '0.4' and r['window'] == '5' and
                r['architecture_threshold'] == '4' and r['confidence_set'] == 'high_moderate' and
                r['genome_subset'] == 'all' and r['metric'] == 'jaccard']
    write_tsv(OUT / 'primary_pair_eligibility.tsv', selected, list(pairs[0]))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.8), gridspec_kw={'width_ratios': [1.1, 1]})
    metrics = ['supported_a_genomes', 'supported_b_genomes', 'supported_ab_genomes']
    y = np.arange(3)
    axes[0].barh(y+.17, [old[k] for k in metrics], height=.29, color='#b8c2cc', label='34 references')
    axes[0].barh(y-.17, [new[k] for k in metrics], height=.29, color='#326c96', label='39 references')
    axes[0].set(yticks=y, yticklabels=['Supported A', 'Supported B', 'Supported A and B'],
                xlim=(0, 20), xticks=[0, 5, 10, 15, 20], xlabel='Genomes (fixed cohort, n = 20)',
                title='A  Functional identification')
    for offset, panel in ((.17, old), (-.17, new)):
        for pos, metric in enumerate(metrics): axes[0].text(panel[metric]+.2, pos+offset, str(panel[metric]), va='center')
    axes[0].invert_yaxis()
    axes[0].legend(loc='lower right', frameon=False)
    reason_counts = Counter(r['reason'] if r['status'] != 'eligible' else 'eligible' for r in selected)
    labels = {'requires_one_supported_lysin_b_per_genome': 'B anchor unavailable',
              'architecture_complex_or_unresolved': 'Architecture unresolved'}
    reasons = sorted(reason_counts, key=lambda key: -reason_counts[key])
    axes[1].barh(np.arange(len(reasons)), [reason_counts[key] for key in reasons], color=['#b8c2cc', '#da8b43'][:len(reasons)])
    axes[1].set(yticks=np.arange(len(reasons)), yticklabels=[labels.get(key, key) for key in reasons],
                xlim=(0, 9), xticks=[0, 3, 6, 9], xlabel='Within-cluster genome pairs (n = 9)',
                title='B  Comparison eligibility after update')
    axes[1].invert_yaxis()
    for pos, key in enumerate(reasons): axes[1].text(reason_counts[key]+.08, pos, str(reason_counts[key]), va='center')
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', alpha=.15)
        ax.set_axisbelow(True)
    fig.text(.04, .045, 'Two M1 CDS match one assayed Bongo sequence. Complete A/B pairs remain 1/20; H2 has no eligible lineage contrasts.', fontsize=9)
    fig.tight_layout(rect=[0, .10, 1, 1], w_pad=2.5)
    for suffix in ('png', 'pdf'):
        metadata = {'Software': 'actinolysis'} if suffix == 'png' else {'Creator': 'actinolysis', 'CreationDate': None, 'ModDate': None}
        fig.savefig(OUT / f'reference_update.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    caption = 'Adding five independently assayed references increases supported LysB-family coverage from one to three genomes at unchanged numerical thresholds. The two added M1 CDS encode the same peptide as Bongo LysB. LysA remains unresolved in both M1 genomes. Eight within-cluster pairs lack supported B anchors and the remaining pair lacks a resolved architecture; no primary H2 contrast is available.'
    files = ['panel_summary.tsv', 'confidence_changes.tsv', 'genome_coverage.tsv', 'primary_pair_eligibility.tsv',
             'reference_update.png', 'reference_update.pdf']
    write_json(OUT / 'comparison_manifest.json', {'caption': caption,
        'changed_cds': len(changed), 'changed_unique_sequences': len({r['protein_sha256'] for r in changed}),
        'primary_pair_exclusion_counts': dict(reason_counts), 'artifacts': {name: sha256(OUT / name) for name in files},
        'summary_script_sha256': sha256(__file__), 'baseline': 'v0.3.0-profile-diagnostic',
        'updated_panel_manifest_sha256': sha256(NEW / 'analysis_manifest.json')})
    validation_source = ROOT / '.experiments/lysterase_2026_validation'
    if (validation_source / 'validation.json').exists():
        report = json.loads((validation_source / 'validation.json').read_text())
        if report['status'] != 'success' or report['reference_count'] != 39:
            raise ValueError('Updated reference validation is incomplete')
        target = ROOT / 'results/lysterase_update/reference_validation'
        names = ['queries.faa', 'queries.tsv', 'raw_hits.tsv', 'alignment_evidence.tsv', 'query_summary.tsv',
                 'challenge_inputs.tsv', 'frozen_inputs.json']
        names += [asset['raw_path'] for asset in json.loads((validation_source / 'frozen_inputs.json').read_text())['challenge_assets']]
        target.mkdir(parents=True, exist_ok=True)
        for name in names:
            (target / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(validation_source / name, target / name)
        summary = {key: report[key] for key in ('status', 'mmseqs_version', 'identification_thresholds', 'search_settings',
            'reference_count', 'query_count', 'query_origins', 'positive_recoverability', 'specificity_challenges')}
        summary['exclusion_policy'] = 'Exclude self, resolved nucleotide-accession base and exact protein-sequence aliases before interpretation.'
        summary['artifacts'] = {name: sha256(target / name) for name in names}
        write_json(target / 'validation_summary.json', summary)
    search_source = ROOT / '.experiments/lysterase_2026_search_evidence'
    if (search_source / 'archive_report.json').exists():
        report = json.loads((search_source / 'archive_report.json').read_text())
        target = ROOT / 'results/lysterase_update/search_evidence'
        target.mkdir(parents=True, exist_ok=True)
        for name, digest in report['outputs_sha256'].items():
            if sha256(search_source / name) != digest: raise ValueError('Search supplement integrity failure')
            shutil.copyfile(search_source / name, target / name)
        write_json(target / 'evidence_manifest.json', {key: report[key] for key in
            ('identification_thresholds', 'counts', 'detection_status_counts', 'reference_evidence_types', 'outputs_sha256', 'interpretation')})
    print(json.dumps({'panels': rows, 'changed_cds': len(changed), 'changed_unique_sequences': len({r['protein_sha256'] for r in changed})}))

if __name__ == '__main__': run()
