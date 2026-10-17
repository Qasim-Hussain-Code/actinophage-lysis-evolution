#!/usr/bin/env python3
"""Compare held-out reference recovery without treating training aliases as discoveries."""
from __future__ import annotations
from collections import defaultdict
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, write_json, write_tsv

OLD = ROOT / 'results/amendment_v1/lysb_profile_candidates'
NEW = ROOT / 'results/lysb_profile_update'
OUT = NEW / 'comparison'

def truth(value): return str(value).lower() == 'true'

def heldout(root, prefix):
    result = {}
    for row in read_tsv(root / 'query_summary.tsv'):
        if row['scan'] != 'held_out_reference' or not row['model'].startswith(prefix): continue
        kind = 'source_out' if '__source_out_' in row['model'] else 'family_out'
        key = row['reference_id'], kind
        if key in result: raise ValueError('Reference appears twice in one holdout regime')
        if truth(row['exact_training_alias']): raise ValueError('A held-out query leaked into training')
        result[key] = row
    return result

def run():
    old, new = heldout(OLD, 'broad_assay9'), heldout(NEW, 'broad_assay14')
    references = read_tsv(NEW / 'training_references.tsv')
    rows = []
    for ref in references:
        identifier = ref['reference_id']
        for panel, data in [('assay9', old), ('assay14', new)]:
            for regime in ('source_out', 'family_out'):
                row = data.get((identifier, regime))
                rows.append({'reference_id': identifier, 'reference_evidence_type': ref['evidence_type'],
                    'panel': panel, 'regime': regime, 'tested': row is not None,
                    'passes_candidate_cuts': truth(row['passes_candidate_cuts']) if row else '',
                    'model': row['model'] if row else '', 'sequence_sha256': ref['protein_sha256']})
    write_tsv(OUT / 'paired_reference_recovery.tsv', rows)
    grouped = defaultdict(list)
    for row in rows:
        if row['tested']: grouped[row['panel'], row['regime']].append(row)
    recovery = [{'panel': panel, 'regime': regime, 'n_tested': len(group),
                 'n_recovered': sum(row['passes_candidate_cuts'] for row in group)}
                for (panel, regime), group in sorted(grouped.items())]
    write_tsv(OUT / 'recovery_summary.tsv', recovery)
    domains = read_tsv(NEW / 'domain_hits.tsv')
    failures = []
    for (identifier, regime), row in new.items():
        if truth(row['passes_candidate_cuts']): continue
        matches = [d for d in domains if d['model'] == row['model'] and d['query_id'] == 'REF__' + identifier]
        best = min(matches, key=lambda d: float(d['independent_domain_evalue'])) if matches else None
        failures.append({'reference_id': identifier, 'regime': regime, 'model': row['model'],
            **{key: best[key] if best else '' for key in ('full_evalue', 'independent_domain_evalue', 'model_coverage',
               'target_coverage', 'residue_pairs', 'hmm_length', 'target_length')},
            'domain_selection': 'Smallest independent-domain E-value; descriptive only; no change to candidate rules'})
    fields = ['reference_id', 'regime', 'model', 'full_evalue', 'independent_domain_evalue', 'model_coverage',
              'target_coverage', 'residue_pairs', 'hmm_length', 'target_length', 'domain_selection']
    write_tsv(OUT / 'held_out_failures.tsv', failures, fields)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    import numpy as np
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42})
    identifiers = [row['reference_id'] for row in references]
    columns = [('assay9', 'family_out'), ('assay14', 'source_out'), ('assay14', 'family_out')]
    matrix = np.full((len(identifiers), len(columns)), -1)
    lookup = {(row['reference_id'], row['panel'], row['regime']): row for row in rows}
    for y, identifier in enumerate(identifiers):
        for x, (panel, regime) in enumerate(columns):
            row = lookup[identifier, panel, regime]
            if row['tested']: matrix[y, x] = int(row['passes_candidate_cuts'])
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 7.6), gridspec_kw={'width_ratios': [1.15, 1]})
    axes[0].imshow(matrix, vmin=-1, vmax=1, cmap=ListedColormap(['#e9ecef', '#da8b43', '#326c96']), aspect='auto')
    axes[0].set(yticks=np.arange(len(identifiers)), yticklabels=[i.replace('_LysB_2026', '').replace('_LysB', '') for i in identifiers],
                xticks=np.arange(3), xticklabels=['9 refs\nGroup out', '14 refs\nSource out', '14 refs\nGroup out'],
                title='A  Held-out reference recovery')
    for y in range(len(identifiers)):
        for x in range(3): axes[0].text(x, y, '–' if matrix[y, x] == -1 else str(matrix[y, x]), ha='center', va='center', color='white' if matrix[y, x] == 1 else '#222222')
    axes[0].tick_params(length=0)
    failure = next(row for row in failures if row['reference_id'] == 'phiGD57_1_LysB_2026' and row['regime'] == 'source_out')
    values = [100 * float(failure['model_coverage']), 100 * float(failure['target_coverage'])]
    axes[1].barh([0, 1], values, color=['#326c96', '#da8b43'], height=.42)
    axes[1].axvline(70, color='#333333', linestyle='--', linewidth=1.2, label='70% criterion')
    axes[1].set(yticks=[0, 1], yticklabels=['Complete model', 'Complete target'], xlim=(0, 100),
                ylim=(-.8, 2.9), xlabel='Paired-residue coverage (%)', title='B  phiGD57-1 source-held-out domain')
    axes[1].invert_yaxis()
    for y, value in enumerate(values): axes[1].text(value+.8, y, f'{value:.1f}%', va='center')
    domain_note = (f"Domain E = {float(failure['independent_domain_evalue']):.2e}\n"
                   f"{failure['residue_pairs']} paired residues / {failure['target_length']} target residues\n"
                   'Coverage criterion fails; threshold retained.')
    axes[1].text(0, 2.0, domain_note, fontsize=9, va='center',
                 bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 2})
    axes[1].spines[['top', 'right']].set_visible(False)
    axes[1].legend(frameon=False, loc='lower right')
    fig.legend(handles=[Patch(color='#326c96', label='Pass'), Patch(color='#da8b43', label='Fail'), Patch(color='#e9ecef', label='Not tested')],
               loc='lower center', ncol=3, bbox_to_anchor=(.5, .07), frameon=False)
    fig.text(.035, .025, 'Groups are reference-sequence proxies, not independent clades. Strong profile similarity alone does not satisfy full-target coverage.', fontsize=9)
    fig.tight_layout(rect=[0, .14, 1, 1], w_pad=3)
    for suffix in ('png', 'pdf'):
        metadata = {'Software': 'actinolysis'} if suffix == 'png' else {'Creator': 'actinolysis', 'CreationDate': None, 'ModDate': None}
        fig.savefig(OUT / f'profile_recovery.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    filenames = ['paired_reference_recovery.tsv', 'recovery_summary.tsv', 'held_out_failures.tsv', 'profile_recovery.png', 'profile_recovery.pdf']
    write_json(OUT / 'comparison_manifest.json', {'caption': 'Reference recovery in the original and expanded assay-profile panels. The expanded broad model recovers 13/14 queries in both holdout regimes. phiGD57-1 fails target coverage despite strong profile similarity. Holdout groups are defined separately in each training panel; no population sensitivity or specificity is estimated.',
        'artifacts': {name: sha256(OUT / name) for name in filenames}, 'summary_script_sha256': sha256(__file__),
        'original_profile_report_sha256': sha256(OLD / 'profile_report.json'), 'updated_profile_report_sha256': sha256(NEW / 'profile_report.json')})
    print(json.dumps({'recovery': recovery, 'held_out_failures': len(failures)}))

if __name__ == '__main__': run()
