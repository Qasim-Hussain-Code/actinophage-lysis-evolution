#!/usr/bin/env python3
"""Reproduce the published exploratory cohort from accession-pinned public inputs.

This reproduces stages 1–5 conditional on the published cohort. It does not
reconstruct selection from the original eligibility frame. A new workspace is
required; interrupted or changed sources are retained for investigation.

Run in the pinned Linux/WSL environment:
    python scripts/reproduce_exploratory_screen.py
    python scripts/reproduce_exploratory_screen.py --workspace cached_reproduction --source-cache /path/to/genbank_cache
"""
from __future__ import annotations

import argparse
import copy
import importlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.parse

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
from actinolysis.acquisition import download, validate_accession
from actinolysis.common import NA, load_config, read_tsv, sha256, write_json, write_tsv
from scripts import replay_exploratory_screen as publication
from scripts import run_exploratory_screen as runner
from scripts.prepare_exploratory_screen import digest_object, new_workspace

PUBLIC = Path('results/exploratory50_v1')
PANEL = Path('results/lysterase_update/panel')
LIMIT = 3_000_000
SCRIPTS = ('reproduce_exploratory_screen.py', 'replay_exploratory_screen.py',
           'run_exploratory_screen.py', 'archive_cohort_hits.py', 'prepare_exploratory_screen.py')


def verify_artifacts(directory):
    manifest = json.loads((directory / 'artifact_manifest.json').read_text(encoding='utf-8'))
    artifacts = manifest['artifacts']
    for name, expected in artifacts.items():
        if (Path(name).name != name or '\\' in name or name in {'', '.', '..'}
                or not re.fullmatch(r'[0-9a-f]{64}', expected)
                or sha256(runner.inside(directory, name)) != expected):
            raise ValueError('Published scientific artifact checksum or path is invalid')
    return manifest


def verify_cohort(selected, pilot):
    if (len(selected) != 50 or len(pilot) != 20
            or [str(r['rank']) for r in selected] != [str(i) for i in range(1, 51)]
            or [r['pilot'] for r in selected] != ['true'] * 20 + ['false'] * 30
            or len({r['accession'].split('.')[0] for r in selected}) != 50
            or len({r['phage'] for r in selected}) != 50):
        raise ValueError('Expected exactly 50 ordered unique genomes nested in the 20-genome pilot')
    for row in selected:
        validate_accession(row['accession'], version_required=True)
        if not re.fullmatch(r'[0-9a-f]{64}', row['sha256']):
            raise ValueError('Invalid frozen genome source checksum')
    for selected_row, pilot_row in zip(selected[:20], pilot):
        if any(selected_row[k] != pilot_row[k] for k in ('phage', 'accession', 'sha256')):
            raise ValueError('Published cohort does not preserve the frozen pilot versions/checksums')


def verify_public(root):
    public = runner.inside(root, PUBLIC)
    manifest = verify_artifacts(public)
    if (manifest.get('analysis') != 'exploratory50_v1'
            or manifest.get('scope') != 'exploratory_only' or manifest.get('replay_identical') is not True):
        raise ValueError('Expected the published exploratory50 scientific freeze')
    report = json.loads((public / 'replay_validation.json').read_text(encoding='utf-8'))
    cohort = report['public_cohort_manifest']
    if (report.get('status') != 'pass' or report.get('completed_stages') != [1, 2, 3, 4, 5]
            or cohort != {'artifact': 'selected_genomes.tsv', 'sha256': manifest['artifacts']['selected_genomes.tsv']}):
        raise ValueError('Published cohort replay provenance is incomplete')
    selected = read_tsv(public / 'selected_genomes.tsv')
    verify_cohort(selected, read_tsv(root / 'data/manifests/pilot_inputs.tsv'))
    provenance = report['frozen_provenance']
    config_path = root / 'config/analysis.yaml'
    if sha256(config_path) != provenance['configuration_file_sha256']:
        raise ValueError('Published configuration file changed')
    cfg = load_config(config_path)
    cfg['mode'] = 'intermediate'
    if digest_object(cfg) != provenance['effective_configuration_sha256']:
        raise ValueError('Published effective configuration changed')
    if runner.tree_hashes(root, ('src',)) != provenance['execution_source_tree_sha256']:
        raise ValueError('Published execution source tree changed')
    for name, expected in provenance['preparation_code_sha256'].items():
        if sha256(runner.inside(root, name)) != expected:
            raise ValueError('Published preparation source changed')
    code = {'scripts/run_exploratory_screen.py': provenance['execution_runner_sha256'],
            'scripts/replay_exploratory_screen.py': manifest['export_driver_sha256'],
            'scripts/archive_cohort_hits.py': manifest['hit_ledger_builder_sha256'],
            'docs/exploratory_screen_protocol.md': report['protocol_sha256']}
    for name, expected in code.items():
        if sha256(root / name) != expected:
            raise ValueError('Published scientific methods changed')
    panel = root / PANEL
    panel_artifacts = json.loads((panel / 'artifact_manifest.json').read_text(encoding='utf-8'))['artifacts']
    for name, expected in report['reference_panel_sha256'].items():
        if panel_artifacts.get(name) != expected or sha256(runner.inside(panel, name)) != expected:
            raise ValueError('Published reference panel changed')
    references = read_tsv(panel / 'data/manifests/reference_manifest.tsv')
    catalogue = read_tsv(panel / 'references/protein_controls.tsv')
    if (len(references) != 39 or len(catalogue) != 39
            or len({r['reference_id'] for r in references}) != 39
            or {r['reference_id'] for r in catalogue} != {r['reference_id'] for r in references}):
        raise ValueError('Expected the unchanged 39-reference catalogue')
    return selected, references, cfg, manifest


def pinned_url(accession):
    validate_accession(accession, version_required=True)
    return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(
        {'db': 'nuccore', 'id': accession, 'rettype': 'gbwithparts', 'retmode': 'text', 'tool': 'actinolysis'})


def source_specs(selected, references):
    sources = {}
    def add(relative, accession, digest, original_url):
        validate_accession(accession, version_required=True)
        if not re.fullmatch(r'[0-9a-f]{64}', digest):
            raise ValueError('Invalid frozen source checksum')
        row = {'path': relative, 'accession': accession, 'sha256': digest,
               'original_source_url': original_url, 'effective_retrieval_url': pinned_url(accession)}
        if relative in sources and sources[relative] != row:
            raise ValueError('Conflicting frozen source identities')
        sources[relative] = row
    for row in selected:
        add(f"data/raw/genbank/{row['accession']}.gb", row['accession'], row['sha256'], row['source_url'])
    for row in references:
        validate_accession(row['source_accession'])
        add(f"data/raw/references/{row['source_accession']}.gb", row['resolved_source_accession'],
            row['source_sha256'], row['source_url'])
    return list(sources.values())


def acquire_source(spec, runtime, cfg, source_cache=None):
    target = runner.inside(runtime, spec['path'])
    if target.exists():
        raise ValueError('Source destination already exists; use a fresh workspace')
    if source_cache is None:
        download(spec['effective_retrieval_url'], target, cfg, maximum=LIMIT)
    else:
        cache = Path(source_cache).resolve()
        names = (Path(spec['path']), Path(spec['accession'] + '.gb'), Path(spec['path']).name,
                 Path('genbank') / (spec['accession'] + '.gb'),
                 Path('references') / Path(spec['path']).name,
                 Path('lysterase_update') / (spec['accession'] + '.gb'),
                 Path('data/raw/lysterase_update') / (spec['accession'] + '.gb'))
        source = next((runner.inside(cache, name) for name in names if runner.inside(cache, name).is_file()), None)
        if source is None:
            raise FileNotFoundError('Frozen source is unavailable in the supplied cache')
        if source.stat().st_size > LIMIT or sha256(source) != spec['sha256']:
            raise ValueError('Frozen source cache size/checksum mismatch')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    if target.stat().st_size > LIMIT or sha256(target) != spec['sha256']:
        raise ValueError('Retrieved source size/checksum mismatch')
    if SeqIO.read(target, 'genbank').id != spec['accession']:
        raise ValueError('Retrieved source accession-version mismatch')
    return target


def direct_manifest(selected):
    """Use the frozen cohort directly; omitted metadata cannot affect sampling."""
    return [{**{k: row[k] for k in ('phage', 'host_genus', 'cluster', 'subcluster', 'accession', 'sha256', 'source_url')},
             'host': NA, 'metadata_source': 'published_accession_frozen_cohort', 'topology_evidence': NA,
             'path': f"data/raw/genbank/{row['accession']}.gb", 'inclusion_status': 'pending_qc',
             'exclusion_reason': 'not_applicable'} for row in selected]


def compare_public_science(public, output, manifest):
    names = sorted(name for name in manifest['artifacts'] if name.endswith('.tsv')
                   or name in {'summary.json', 'quality_control.json'})
    if not names or {p.name for p in output.iterdir() if p.is_file()} != set(names):
        raise ValueError('Reproduced scientific artifact set differs from the publication')
    rows = [{'path': name, 'expected_sha256': manifest['artifacts'][name],
             'actual_sha256': sha256(output / name)} for name in names]
    if any(row['expected_sha256'] != row['actual_sha256'] for row in rows):
        raise ValueError('Reproduced public scientific artifacts differ')
    return rows


def worker(root, workspace, source_cache):
    runtime = runner.inside(runner.inside(root, '.experiments'), workspace)
    if ROOT.resolve() != runtime:
        raise ValueError('Replay worker must execute its isolated source copy')
    selected, references, planned_cfg, manifest = verify_public(root)
    inputs = json.loads((runtime / 'reproduction_inputs.json').read_text(encoding='utf-8'))
    cfg = load_config(runtime / 'config/analysis.yaml')
    if (runner.tree_hashes(runtime, ('src',)) != runner.tree_hashes(root, ('src',))
            or digest_object(cfg) != inputs['execution_configuration_sha256']
            or sha256(Path(cfg['tools']['mmseqs'])) != inputs['mmseqs_sha256']):
        raise ValueError('Isolated execution inputs changed')
    for name, expected in inputs['snapshot_sha256'].items():
        if sha256(runner.inside(runtime, name)) != expected:
            raise ValueError('Isolated execution input checksum mismatch')
    sources = source_specs(selected, references)
    acquired = []
    report = {'status': 'running', 'method': 'public_accession_frozen_cohort_reproduction',
              'source_acquisition_scope': 'checksum_verified_cache_only' if source_cache else 'fresh_version_pinned_downloads',
              'completed_stages': [], 'selection_frame_reconstructed': False,
              'h2_analysis_performed': False, 'software_scope': 'Active environment and pinned binary; no independent installation claim.'}
    try:
        for spec in sources:
            runner.root_guard(cfg, root, LIMIT)
            acquire_source(spec, runtime, cfg, source_cache)
            acquired.append(spec)
            write_tsv(runtime / 'data/manifests/reproduction_sources.tsv', acquired)
        write_tsv(runtime / 'data/manifests/intermediate_inputs.tsv', direct_manifest(selected))
        report['completed_stages'].append(1)
        for number in range(2, 6):
            verify_public(root)
            runner.root_guard(cfg, root, 600_000_000 if number == 3 else 0)
            module_name, function = runner.STAGES[number]
            getattr(importlib.import_module('actinolysis.' + module_name), function)(cfg)
            runner.root_guard(cfg, root)
            report['completed_stages'].append(number)
        runner.check_controls(root, runtime, {'panel_root': PANEL.as_posix()})
        publication.retained_panel_validation(root, runtime, cfg)
        tables, summary = publication.public_tables(runtime, selected, planned_cfg)
        output = runtime / 'results/reproduced_science'
        runner.root_guard(cfg, root, 20_000_000)
        for name, (rows, fields) in tables.items():
            write_tsv(output / name, rows, fields)
        write_json(output / 'summary.json', summary)
        qc = json.loads((runtime / 'results/reports/quality_control.json').read_text(encoding='utf-8'))
        write_json(output / 'quality_control.json', {k: qc[k] for k in
            ('selected', 'included', 'excluded', 'cds', 'translation_status', 'near_identical_screen')})
        comparisons = compare_public_science(root / PUBLIC, output, manifest)
        report.update(status='pass', compared_artifacts=len(comparisons), artifacts=comparisons,
                      wrapper_sha256=sha256(__file__), mmseqs_sha256=inputs['mmseqs_sha256'])
        runner.root_guard(cfg, root)
    except Exception as error:
        report.update(status='failed', error_type=type(error).__name__, acquired_sources=len(acquired))
        raise
    finally:
        write_json(runtime / 'results/reports/public_accession_reproduction.json', report)
    return report


def prepare(root, workspace, mmseqs=None):
    _, _, cfg, _ = verify_public(root)
    runtime = new_workspace(root, workspace)
    binary = Path(mmseqs).resolve() if mmseqs else runner.inside(root, cfg['tools']['mmseqs'])
    if not binary.is_file():
        raise FileNotFoundError('MMseqs executable is unavailable')
    cfg = copy.deepcopy(cfg)
    cfg['offline'] = True  # All sources are acquired explicitly before unchanged stages 2–5.
    cfg['tools']['mmseqs'] = str(binary)
    cfg['resources']['effective_budget_bytes'] = min(int(cfg['resources']['effective_budget_bytes']),
        shutil.disk_usage(root).free - int(cfg['resources']['reserve_bytes']))
    runner.root_guard(cfg, root, 10_000_000)
    shutil.copytree(root / 'src', runtime / 'src', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in SCRIPTS:
        target = runtime / 'scripts' / name
        target.parent.mkdir(exist_ok=True)
        shutil.copyfile(root / 'scripts' / name, target)
    for relative in ('references/protein_controls.tsv', 'data/manifests/reference_manifest.tsv'):
        target = runtime / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / PANEL / relative, target)
    import yaml
    (runtime / 'config').mkdir()
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    write_json(runtime / 'reproduction_inputs.json', {'execution_configuration_sha256': digest_object(cfg),
        'mmseqs_sha256': sha256(binary), 'wrapper_sha256': sha256(__file__),
        'snapshot_sha256': runner.tree_hashes(runtime, ('src', 'scripts', 'config', 'references'),
                                             ('data/manifests/reference_manifest.tsv',))})
    runner.root_guard(cfg, root)
    return runtime


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', default='public_exploratory50_reproduction', help='New directory beneath .experiments')
    parser.add_argument('--source-cache', type=Path, help='Use verified GenBank cache only; never download missing records')
    parser.add_argument('--mmseqs', type=Path)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--source-root', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.source_root is None:
            parser.error('Worker requires the source repository root')
        report = worker(args.source_root.resolve(), args.workspace, args.source_cache)
    else:
        runtime = prepare(ROOT, args.workspace, args.mmseqs)
        command = [sys.executable, '-B', str(runtime / 'scripts/reproduce_exploratory_screen.py'),
                   '--worker', '--source-root', str(ROOT), '--workspace', args.workspace]
        if args.source_cache:
            command += ['--source-cache', str(args.source_cache.resolve())]
        subprocess.run(command, cwd=runtime, check=True)
        report = json.loads((runtime / 'results/reports/public_accession_reproduction.json').read_text(encoding='utf-8'))
    print(json.dumps({'status': report['status'], 'completed_stages': report['completed_stages'],
                      'source_acquisition_scope': report['source_acquisition_scope']}))


if __name__ == '__main__':
    main()
