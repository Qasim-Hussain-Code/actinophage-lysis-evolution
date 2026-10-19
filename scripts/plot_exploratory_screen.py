#!/usr/bin/env python3
"""Plot validated counts from the exported exploratory screen, without rerunning analysis."""
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
INPUTS = ('coverage_by_group.tsv', 'architectures.tsv', 'genome_qc.tsv', 'summary.json')
GENERA = ('Mycobacterium', 'Gordonia')
ASSIGNMENTS = ('supported_A_genomes', 'supported_B_genomes', 'supported_A_and_B_genomes')
ARCHITECTURES = ('closely_linked', 'separated', 'complex', 'unresolved')


def figure_data(source):
    groups = read_tsv(source / 'coverage_by_group.tsv')
    architecture = read_tsv(source / 'architectures.tsv')
    genomes = read_tsv(source / 'genome_qc.tsv')
    summary = json.loads((source / 'summary.json').read_text(encoding='utf-8'))
    if (summary.get('analysis') != 'exploratory_50_genome_detection_screen'
            or summary.get('scope') != 'exploratory_only'
            or summary.get('h2_analysis_performed') is not False):
        raise ValueError('Expected completed exploratory-screen scientific tables')
    by_genus = {r['group']: r for r in groups if r['group_type'] == 'host_genus'}
    totals = [r for r in groups if r['group_type'] == 'all' and r['group'] == 'all']
    if (set(by_genus) != set(GENERA) or len(by_genus) != sum(r['group_type'] == 'host_genus' for r in groups)
            or len(totals) != 1 or not genomes):
        raise ValueError('Missing or duplicated genus/overall summaries')
    if (len({r['accession'] for r in genomes}) != len(genomes)
            or any(r['host_genus'] not in GENERA or r['inclusion_status'] != 'included' for r in genomes)
            or len({r['accession'] for r in architecture}) != len(architecture)
            or {r['accession'] for r in architecture} != {r['accession'] for r in genomes}
            or any(r['scope'] != 'high_moderate' or r['threshold'] != '4' for r in architecture)):
        raise ValueError('QC and primary architecture records disagree')
    if any(r['architecture'] not in ARCHITECTURES for r in architecture):
        raise ValueError('Unrecognized primary architecture class')
    rows = []
    for genus in GENERA:
        cohort = [r for r in genomes if r['host_genus'] == genus]
        accessions = {r['accession'] for r in cohort}
        local = [r for r in architecture if r['accession'] in accessions]
        if any(r['host_genus'] != genus for r in local):
            raise ValueError('Architecture host genus differs from the QC record')
        observed = {
            'supported_A_genomes': sum(int(r['lysin_a_count']) > 0 for r in local),
            'supported_B_genomes': sum(int(r['lysin_b_count']) > 0 for r in local),
            'supported_A_and_B_genomes': sum(int(r['lysin_a_count']) > 0 and int(r['lysin_b_count']) > 0 for r in local),
            **{key: sum(r['architecture'] == key for r in local) for key in ARCHITECTURES}}
        if (int(by_genus[genus]['included_genomes']) != len(cohort)
                or int(by_genus[genus]['selected_genomes']) < len(cohort)
                or int(by_genus[genus]['cds']) != sum(int(r['cds_count']) for r in cohort)
                or any(int(by_genus[genus][key]) != count for key, count in observed.items())):
            raise ValueError('Genus coverage summary differs from the source records')
        for row in local:
            if row['architecture'] in {'closely_linked', 'separated'}:
                if int(row['lysin_a_count']) != 1 or int(row['lysin_b_count']) != 1:
                    raise ValueError('A resolved pair requires unique supported A and B')
                if (int(row['intervening_cds']) > 4) != (row['architecture'] == 'separated'):
                    raise ValueError('Resolved architecture contradicts the separation threshold')
        for metric in ASSIGNMENTS + ARCHITECTURES:
            rows.append({'panel': 'A' if metric in ASSIGNMENTS else 'B', 'host_genus': genus,
                'selected_genomes': int(by_genus[genus]['selected_genomes']),
                'included_genomes': len(cohort), 'metric': metric, 'count': observed[metric],
                'denominator': len(cohort), 'scope': 'high_moderate_curated_reference_homology'})
    for key in ('selected_genomes', 'included_genomes', 'cds') + ASSIGNMENTS + ARCHITECTURES:
        total = sum(int(by_genus[g][key]) for g in GENERA)
        if int(totals[0][key]) != total or summary['counts'][key] != total:
            raise ValueError('Overall summary does not reconcile with genus counts')
    if summary['counts']['selected_genomes'] != 50 or summary['counts']['included_genomes'] != len(genomes):
        raise ValueError('Expected 50 selected records and reconciled QC denominators')
    return rows


def run(output_path=None):
    source = SOURCE
    inputs = {name: sha256(source / name) for name in INPUTS}
    rows = figure_data(source)
    output = source / 'figures' if output_path is None else (ROOT / output_path).resolve()
    if not output.is_relative_to(ROOT.resolve()) or output in {ROOT.resolve(), source.resolve()}:
        raise ValueError('Use a named figure directory inside the project')
    if output.exists():
        raise ValueError('Figure output exists; preserve the freeze before regenerating it')
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.cache/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import numpy as np

    plt.rcParams.update({'font.size': 9.5, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    lookup = {(r['host_genus'], r['metric']): r for r in rows}
    included = {g: lookup[g, ASSIGNMENTS[0]]['denominator'] for g in GENERA}
    selected = {g: lookup[g, ASSIGNMENTS[0]]['selected_genomes'] for g in GENERA}
    maximum = max(1, max(included.values()))
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.6))
    specifications = [
        (ASSIGNMENTS, ['A category', 'B category', 'Both categories'],
         ['#326c96', '#da8b43', '#47877c'], 'A  A/B category detection'),
        (ARCHITECTURES, ['Unique linked pair', 'Unique separated pair', 'Complex', 'Unresolved'],
         ['#326c96', '#da8b43', '#8c78a3', '#b8c2cc'], 'B  A/B architecture')]
    for ax, (metrics, labels, colors, title) in zip(axes, specifications):
        width = .72 / len(metrics)
        positions = np.arange(len(GENERA))
        for index, (metric, label, color) in enumerate(zip(metrics, labels, colors)):
            offset = (index - (len(metrics) - 1) / 2) * width
            counts = [lookup[g, metric]['count'] for g in GENERA]
            ax.barh(positions + offset, counts, height=width * .82, color=color, label=label)
            for y, count in zip(positions + offset, counts):
                ax.text(count + maximum * .018, y, str(count), va='center', fontsize=9)
        ax.set(yticks=positions, yticklabels=[f'{g}\n(n = {included[g]} included)' for g in GENERA],
               xlim=(0, maximum * 1.14), ylim=(-.6, len(GENERA) - .4),
               xlabel='Number of genomes', title=title)
        ax.invert_yaxis()
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='x', alpha=.15)
        ax.set_axisbelow(True)
        ax.legend(loc='upper center', bbox_to_anchor=(.5, -.17), frameon=False,
                  ncol=3 if len(metrics) == 3 else 2, fontsize=8.5)
    fig.text(.025, .075, 'Supported = high/moderate curated homology. A and B counts overlap; both categories does not require a unique pair.', fontsize=9)
    fig.text(.025, .025, 'Exploratory cohort counts; unresolved is not absence. Neighborhood conservation was not tested.', fontsize=9)
    fig.subplots_adjust(left=.12, right=.99, bottom=.32, top=.89, wspace=.60)
    if inputs != {name: sha256(source / name) for name in INPUTS}:
        raise ValueError('Scientific inputs changed during figure preparation')
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / 'figure_source.tsv', rows)
    for suffix in ('png', 'pdf'):
        metadata = {'Software': 'actinolysis'} if suffix == 'png' else {
            'Creator': 'actinolysis', 'CreationDate': None, 'ModDate': None}
        fig.savefig(output / f'exploratory_detection_architecture.{suffix}', dpi=300, metadata=metadata)
    plt.close(fig)
    caption = (
        'Reference-based detection and architecture in the exploratory 50-record screen. '
        'Panel A shows high/moderate curated-homology A-family, B-family and both-category genome counts; '
        'A and B overlap, and both-category detection does not require a unique A/B pair. '
        'Panel B classifies unique linked or separated A/B pairs, complex architectures and unresolved genomes. '
        'Separation requires more than four intervening CDS. Included-genome denominators are '
        + '; '.join(f'{g}, n={included[g]} of {selected[g]} selected' for g in GENERA) + '. '
        'Counts retain the frozen evidence-tier caps and do not estimate population prevalence, diagnostic sensitivity or neighborhood conservation.')
    names = ('figure_source.tsv', 'exploratory_detection_architecture.png', 'exploratory_detection_architecture.pdf')
    write_json(output / 'figure_manifest.json', {'caption': caption,
        'source_artifacts_sha256': inputs, 'plot_script_sha256': sha256(__file__),
        'png_dpi': 300, 'pdf_rendering': 'vector',
        'artifacts': {name: sha256(output / name) for name in names}})
    print(json.dumps({'status': 'pass', 'included_denominators': included}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='New figure directory inside the project')
    run(parser.parse_args().output)
