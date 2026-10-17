#!/usr/bin/env python3
"""Display descriptive contrasts conditional on deposited annotations."""
from __future__ import annotations
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, write_json, write_tsv

SOURCE = ROOT / 'results/annotation_neighborhood_sensitivity'
OUT = SOURCE / 'comparison'

def run():
    rows = read_tsv(SOURCE / 'conditional_pair_scores.tsv')
    base = [row for row in rows if float(row['family_identity']) == .4 and int(row['window']) == 5 and
            int(row['architecture_threshold']) == 4 and row['genome_subset'] == 'all']
    write_tsv(OUT / 'base_conditional_scores.tsv', base, list(rows[0]))
    scope = read_tsv(SOURCE / 'base_scope_comparison.tsv')
    availability = []
    for name in ('independent_high_moderate', 'deposited_annotation_only'):
        selected = [row for row in scope if row['evidence_scope'] == name]
        eligible = [row for row in selected if row['status'] in {'eligible', 'descriptive_eligible'}]
        availability.append({'evidence_scope': name, 'candidate_pairs': len(selected), 'available_pairs': len(eligible),
            'cluster_labels': len({row['cluster'] for row in eligible}),
            'separated_available_pairs': sum(row['pair_architecture'] == 'separated' for row in eligible)})
    write_tsv(OUT / 'scope_availability.tsv', availability)
    separated = [row for row in rows if int(row['architecture_threshold']) == 4 and row['genome_subset'] == 'all' and
                 row['pair_architecture'] == 'separated' and row['metric'] == 'jaccard']
    write_tsv(OUT / 'separated_sensitivities.tsv', separated, list(rows[0]))
    summary = [row for row in read_tsv(SOURCE / 'conditional_cluster_summary.tsv') if float(row['family_identity']) == .4 and
               int(row['window']) == 5 and int(row['architecture_threshold']) == 4 and row['genome_subset'] == 'all']
    write_tsv(OUT / 'base_cluster_summary.tsv', summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={'width_ratios': [1.1, 1.55, 1.0]})
    for index, row in enumerate(availability):
        count, total = row['available_pairs'], row['candidate_pairs']
        axes[0].barh(index, count, color='#326c96', height=.42)
        axes[0].barh(index, total-count, left=count, color='#d6dbe0', height=.42)
        axes[0].text(count/2 if count else .1, index, str(count), color='white' if count else '#222222', ha='center' if count else 'left', va='center')
    axes[0].set(yticks=[0, 1], yticklabels=['Independent\nevidence', 'Deposited\nlabels'], xlim=(0, 9), xticks=[0, 3, 6, 9],
                xlabel='Within-cluster pairs (n = 9)', title='A  Scope and comparison availability')
    axes[0].invert_yaxis()
    axes[0].spines[['top', 'right']].set_visible(False)
    available = sorted([row for row in base if row['metric'] == 'jaccard' and row['status'] == 'descriptive_eligible'],
                       key=lambda row: (row['pair_architecture'] != 'separated', row['cluster']))
    for index, row in enumerate(available):
        target, background = float(row['target_score']), float(row['control_mean'])
        axes[1].plot([background, target], [index, index], color='#999999', linewidth=1.5)
        axes[1].scatter(background, index, marker='s', s=38, facecolor='white', edgecolor='#555555', zorder=3)
        axes[1].scatter(target, index, marker='o', s=30, color='#326c96', zorder=4)
        axes[1].text(1.035, index, f"Δ {float(row['delta']):+.3f}\n{row['n_controls']} anchors", va='center', fontsize=8)
    labels = [row['cluster'] + (' (separated)' if row['pair_architecture'] == 'separated' else ' (linked)') for row in available]
    axes[1].set(yticks=np.arange(len(available)), yticklabels=labels, xlim=(-.02, 1.35), xticks=[0, .5, 1],
                xlabel='Flanking-family Jaccard score', title='B  Conditional target and background')
    axes[1].set_ylim(-.65, max(.5, len(available)-.35))
    axes[1].invert_yaxis()
    axes[1].spines[['top', 'right']].set_visible(False)
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='none', color='#326c96', label='Putative B neighborhood'),
                           Line2D([], [], marker='s', linestyle='none', markerfacecolor='white', color='#555555', label='Matched background mean')],
                   loc='lower center', bbox_to_anchor=(.61, .075), ncol=2, frameon=False, fontsize=8)
    fig.legend(handles=[Patch(color='#326c96', label='Contrast available'), Patch(color='#d6dbe0', label='Unavailable')],
               loc='lower left', bbox_to_anchor=(.035, .075), ncol=2, frameon=False, fontsize=8)
    sep_labels = sorted({row['cluster'] for row in separated})
    matrix = np.full((len(sep_labels), 3), np.nan)
    for row in separated:
        if float(row['family_identity']) == .4 and row['status'] == 'descriptive_eligible':
            matrix[sep_labels.index(row['cluster']), [3, 5, 10].index(int(row['window']))] = float(row['delta'])
    view = axes[2].imshow(matrix, cmap='RdBu', vmin=-.5, vmax=.5, aspect='auto')
    axes[2].set(yticks=np.arange(len(sep_labels)), yticklabels=sep_labels, xticks=np.arange(3), xticklabels=[3, 5, 10],
                xlabel='Flanking CDS per side', title='C  Separated annotations\nWindow sensitivity (40% families)')
    for y in range(len(sep_labels)):
        for x in range(3):
            value = matrix[y, x]
            axes[2].text(x, y, f'{value:+.3f}' if np.isfinite(value) else 'NA', ha='center', va='center',
                         color='white' if np.isfinite(value) and abs(value) > .3 else '#222222', fontsize=8)
    fig.colorbar(view, ax=axes[2], fraction=.06, pad=.04, label='Target minus background')
    fig.text(.03, .018, 'Conditional on deposited labels. Four descriptive comparisons are available; two have separated annotations. Primary H2 remains untestable.', fontsize=9)
    fig.tight_layout(rect=[0, .19, 1, 1], w_pad=2)
    for suffix in ('png','pdf'):
        metadata = {'Software': 'actinolysis'} if suffix == 'png' else {'Creator':'actinolysis','CreationDate':None,'ModDate':None}
        fig.savefig(OUT / f'annotation_contrasts.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    files = ['base_conditional_scores.tsv','scope_availability.tsv','separated_sensitivities.tsv','base_cluster_summary.tsv',
             'annotation_contrasts.png','annotation_contrasts.pdf']
    write_json(OUT / 'figure_manifest.json', {'caption': 'Neighborhood contrasts conditional on deposited lysis labels, using unchanged gene families and matching rules. At the base configuration four pairs have sufficient matched background anchors. The two separated pairs have Jaccard contrasts of 0.000 (AA) and -0.413 (KA). AA becomes slightly positive at the widest window; all family/window/cutoff sensitivities are retained. These observations do not independently validate lysis function or test primary H2.',
        'artifacts': {name: sha256(OUT / name) for name in files}, 'source_report_sha256': sha256(SOURCE / 'sensitivity_report.json'),
        'summary_script_sha256': sha256(__file__)})
    print(json.dumps({'availability': availability, 'conditional_base_comparisons': len(available)}))

if __name__ == '__main__': run()
