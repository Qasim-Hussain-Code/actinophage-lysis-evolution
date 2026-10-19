#!/usr/bin/env python3
"""Plot the five resolved local architectures from the published coordinate ledger."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, write_json, write_tsv

SOURCE = ROOT / 'results/exploratory50_v1'
INPUTS = ('architectures.tsv', 'protein_detection_status.tsv')


def run(output_path=None):
    manifest = json.loads((SOURCE / 'artifact_manifest.json').read_text(encoding='utf-8'))
    hashes = {name: sha256(SOURCE / name) for name in INPUTS}
    if any(manifest['artifacts'].get(name) != digest for name, digest in hashes.items()):
        raise ValueError('Published coordinate or architecture checksum changed')
    pairs = [r for r in read_tsv(SOURCE / 'architectures.tsv') if r['architecture'] == 'closely_linked']
    cds = read_tsv(SOURCE / 'protein_detection_status.tsv')
    by_id = {r['protein_id']: r for r in cds}
    if len(by_id) != len(cds) or len(pairs) != 5 or len({r['accession'] for r in pairs}) != 5:
        raise ValueError('Expected five unique resolved pairs and a unique CDS ledger')
    contexts, rows = [], []
    for pair in pairs:
        if (pair['scope'] != 'high_moderate' or pair['threshold'] != '4'
                or pair['lysin_a_count'] != '1' or pair['lysin_b_count'] != '1'):
            raise ValueError('The displayed pair is not uniquely resolved in the primary scope')
        a, b = by_id[pair['lysin_a']], by_id[pair['lysin_b']]
        if a['accession'] != pair['accession'] or b['accession'] != pair['accession']:
            raise ValueError('Pair and CDS source accession differ')
        lo, hi = sorted((int(a['cds_index']), int(b['cds_index'])))
        if hi - lo - 1 != int(pair['intervening_cds']):
            raise ValueError('CDS order differs from the primary architecture')
        anchor = int(a['start'])
        local = sorted((r for r in cds if r['accession'] == pair['accession']
            and lo - 3 <= int(r['cds_index']) <= hi + 3), key=lambda r: int(r['cds_index']))
        context = []
        for gene in local:
            parts = json.loads(gene['location_parts'])
            if not parts or any(int(end) <= int(start) for start, end in parts):
                raise ValueError('Invalid source CDS interval')
            role = 'A' if gene['protein_id'] == a['protein_id'] else 'B' if gene['protein_id'] == b['protein_id'] else 'flanking_CDS'
            row = {'accession': pair['accession'], 'phage': pair['phage'], 'cluster': pair['cluster'],
                'pair_confidence': pair['confidence'], 'protein_id': gene['protein_id'],
                'cds_index': int(gene['cds_index']), 'start': int(gene['start']), 'end': int(gene['end']),
                'strand': int(gene['strand']), 'location_parts': gene['location_parts'],
                'anchor_start': anchor, 'relative_start': int(gene['start']) - anchor,
                'relative_end': int(gene['end']) - anchor, 'role': role,
                'intervening_cds': int(pair['intervening_cds']), 'gap_bp': int(pair['gap_bp'])}
            rows.append(row)
            context.append(row)
        contexts.append((pair, context))
    output = SOURCE / 'linked_contexts' if output_path is None else (ROOT / output_path).resolve()
    if not output.is_relative_to(ROOT.resolve()) or output in {ROOT.resolve(), SOURCE.resolve()}:
        raise ValueError('Use a named figure directory inside the project')
    if output.exists():
        raise ValueError('Preserve the existing figure freeze before regeneration')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Polygon
    from matplotlib.ticker import MaxNLocator

    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42})
    colors = {'A': '#326c96', 'B': '#da8b43', 'flanking_CDS': '#c4cbd1'}
    fig, axes = plt.subplots(len(contexts), 1, figsize=(11.3, 6.8), sharex=True)
    low = min(r['relative_start'] for r in rows) - 150
    high = max(r['relative_end'] for r in rows) + 150
    for ax, (pair, genes) in zip(axes, contexts):
        ax.axhline(0, color='#d9dfe3', lw=.6, zorder=0)
        for gene in genes:
            parts = json.loads(gene['location_parts'])
            for start, end in parts:
                left, right = int(start) - gene['anchor_start'], int(end) - gene['anchor_start']
                head = min(100, (right - left) * .4)
                if gene['strand'] >= 0:
                    vertices = [(left, -.15), (right - head, -.15), (right, 0), (right - head, .15), (left, .15)]
                else:
                    vertices = [(right, -.15), (left + head, -.15), (left, 0), (left + head, .15), (right, .15)]
                ax.add_patch(Polygon(vertices, facecolor=colors[gene['role']], edgecolor='white', lw=.45))
            if gene['role'] != 'flanking_CDS':
                ax.text((gene['relative_start'] + gene['relative_end']) / 2, .22,
                        gene['role'], ha='center', va='bottom', weight='bold', color=colors[gene['role']])
        ax.set_title(f"{pair['phage']} | {pair['accession']} | cluster {pair['cluster']} | {pair['confidence']}",
                     loc='left', fontsize=9.5, pad=2)
        ax.set(xlim=(low, high), ylim=(-.30, .54), yticks=[])
        ax.spines[['top', 'left', 'right']].set_visible(False)
        ax.spines['bottom'].set_color('#d9dfe3')
        ax.tick_params(axis='x', labelsize=8)
    axes[-1].xaxis.set_major_locator(MaxNLocator(nbins=9, integer=True))
    axes[-1].set_xlabel('Position relative to A-family CDS start (bp)')
    fig.suptitle('Resolved closely linked A/B architectures', y=.985, fontsize=13)
    fig.legend(handles=[Patch(facecolor=colors[k], label=label) for k, label in
        [('A', 'A-family assignment'), ('B', 'B-family assignment'), ('flanking_CDS', 'Flanking CDS')]],
        loc='lower center', bbox_to_anchor=(.5, .025), frameon=False, ncol=3)
    fig.subplots_adjust(left=.055, right=.985, bottom=.14, top=.92, hspace=.66)
    if hashes != {name: sha256(SOURCE / name) for name in INPUTS}:
        raise ValueError('Scientific source tables changed during plotting')
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / 'figure_source.tsv', rows)
    for suffix in ('png', 'pdf'):
        metadata = {'Software': 'actinolysis'} if suffix == 'png' else {
            'Creator': 'actinolysis', 'CreationDate': None, 'ModDate': None}
        fig.savefig(output / f'linked_architectures.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    caption = ('Five uniquely resolved high/moderate A/B pairs in the 50-genome screen. '
        'Three flanking annotated CDS on either side are displayed in original order and strand, '
        'using coordinates relative to the A-family CDS start and a shared scale. '
        'D29 has one intervening CDS; the other pairs overlap by one nucleotide. '
        'Colors denote reference-supported categories, with the pair confidence shown per row; '
        'matching category colors do not imply identical sequences, independent ancestry or demonstrated native activity.')
    names = ('figure_source.tsv', 'linked_architectures.png', 'linked_architectures.pdf')
    write_json(output / 'figure_manifest.json', {'caption': caption, 'png_dpi': 300,
        'pdf_rendering': 'vector', 'source_artifacts_sha256': hashes,
        'plot_script_sha256': sha256(__file__),
        'artifacts': {name: sha256(output / name) for name in names}})
    print(json.dumps({'status': 'pass', 'displayed_pairs': len(pairs), 'displayed_cds': len(rows)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='New figure directory inside the project')
    run(parser.parse_args().output)
