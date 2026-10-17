#!/usr/bin/env python3
"""Run the frozen-cohort assay-reference follow-up in an isolated workspace."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import urllib.parse

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.acquisition import download
from actinolysis.common import load_config, read_tsv, sha256, write_json, write_tsv
from replay_cohort import cohort_frame, copy_cached
import yaml

BASE = ROOT / 'results/amendment_v1/profile'
PUBLIC = ROOT / 'results/lysterase_update/panel'
PROTOCOL = ROOT / 'docs/lysterase_update_protocol.md'

def scientific_files(root):
    return {p.relative_to(root).as_posix(): sha256(p) for folder in ('results/tables', 'figures')
            for p in sorted((root / folder).rglob('*')) if p.is_file()}

def catalogue():
    old = read_tsv(BASE / 'references/protein_controls.tsv')
    if len(old) != 34: raise ValueError('Expected the frozen 34-reference combined panel')
    columns = list(old[0])
    new = read_tsv(ROOT / 'references/lysterase_2026_expansion.tsv')
    if len(new) != 5 or any(r['evidence_type'] != 'experimental_antimycobacterial' or r['panel_action'] != 'add' for r in new):
        raise ValueError('Expected exactly five eligible assayed references')
    specifications = list(old)
    for row in new:
        item = {field: row.get(field, '') for field in columns}
        item['description'] = '; '.join(row[key] for key in ('assay', 'evidence_locator', 'caveat'))
        specifications.append(item)
    if len({r['reference_id'] for r in specifications}) != 39:
        raise ValueError('Duplicate or missing updated reference identifiers')
    return specifications

def fetch_frozen(source, destination, cfg, offline):
    accession, digest = source['accession'], source['sha256']
    if destination.exists():
        if sha256(destination) != digest: raise ValueError('Cached source checksum mismatch')
        return
    source_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode({
        'db': 'nuccore', 'id': accession, 'rettype': 'gbwithparts', 'retmode': 'text', 'tool': 'actinolysis'})
    download(source_url, destination, cfg, offline, maximum=3_000_000)
    if sha256(destination) != digest: raise ValueError('Downloaded source differs from the frozen record')

def prepare(runtime, cfg, cache, offline):
    if runtime.exists() and any(runtime.iterdir()):
        raise ValueError('Workspace exists; preserve it and select a new workspace or resume')
    specifications = catalogue()
    for directory in ('src', 'references'):
        shutil.copytree(ROOT / directory, runtime / directory, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    (runtime / 'scripts').mkdir()
    shutil.copyfile(ROOT / 'scripts/run_all.py', runtime / 'scripts/run_all.py')
    write_tsv(runtime / 'references/protein_controls.tsv', specifications)
    cohort = read_tsv(ROOT / 'data/manifests/pilot_inputs.tsv')
    (runtime / 'data/manifests').mkdir(parents=True)
    shutil.copyfile(ROOT / 'data/manifests/pilot_inputs.tsv', runtime / 'data/manifests/pilot_inputs.tsv')
    frame = runtime / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(frame, cohort_frame(cohort, cfg))
    write_json(runtime / 'data/manifests/phagesdb_snapshot.json', {'pages': [], 'rows': 20,
        'selection_snapshot': {'path': 'data/manifests/eligibility_snapshot.tsv', 'rows': 20, 'sha256': sha256(frame),
                               'basis': 'Frozen study cohort; no population resampling'},
        'method': 'accession_frozen_cohort_replay'})
    frozen_sources = {}
    for row in cohort:
        frozen_sources[row['path']] = {'accession': row['accession'], 'sha256': row['sha256']}
    references = read_tsv(BASE / 'data/manifests/reference_manifest.tsv')
    for row in references:
        relative = f"data/raw/references/{row['source_accession']}.gb"
        item = {'accession': row['resolved_source_accession'], 'sha256': row['source_sha256']}
        if relative in frozen_sources and frozen_sources[relative] != item: raise ValueError('Conflicting reference sources')
        frozen_sources[relative] = item
    curated = {row['protein_accession']: row for row in read_tsv(ROOT / 'results/lysterase_update/curation/source_sequences.tsv')}
    for row in specifications[34:]:
        source = curated[row['protein_accession']]
        relative = f"data/raw/references/{row['source_accession']}.gb"
        frozen_sources[relative] = {'accession': source['source_accession'], 'sha256': source['source_sha256']}
    for relative, source in frozen_sources.items():
        target = runtime / relative
        if not copy_cached(cache, relative, source['sha256'], target):
            curated_source = ROOT / 'data/raw/lysterase_update' / f"{source['accession']}.gb"
            if curated_source.is_file() and sha256(curated_source) == source['sha256']:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(curated_source, target)
            else: fetch_frozen(source, target, cfg, offline)
    (runtime / 'config').mkdir()
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    snapshot = {p.relative_to(runtime).as_posix(): sha256(p)
                for directory in ('src', 'references', 'scripts') for p in sorted((runtime / directory).rglob('*')) if p.is_file()}
    write_json(runtime / 'update_inputs.json', {'protocol_sha256': sha256(PROTOCOL), 'snapshot_sha256': snapshot,
        'frozen_sources': frozen_sources, 'cohort_sha256': sha256(ROOT / 'data/manifests/pilot_inputs.tsv'),
        'curation_sha256': sha256(ROOT / 'references/lysterase_2026_expansion.tsv'),
        'driver_sha256': sha256(__file__)})

def export(runtime):
    if PUBLIC.exists(): raise ValueError('Preserve the existing published panel; export requires a new analysis version')
    execution = json.loads((runtime / 'results/reports/pilot_execution.json').read_text())
    receipts = json.loads((runtime / 'results/reports/pilot_stage_receipts.json').read_text())
    if execution['status'] != 'success' or set(receipts) != {str(i) for i in range(1, 10)}:
        raise ValueError('All nine stages must complete before publication')
    for receipt in receipts.values():
        for path, digest in receipt['outputs'].items():
            if sha256(runtime / path) != digest: raise ValueError('Completed stage output changed before export')
    scientific = scientific_files(runtime)
    if len(scientific) != 39: raise ValueError('Expected 39 tables and figure artifacts')
    for directory in ('results/tables', 'figures', 'src', 'references'):
        shutil.copytree(runtime / directory, PUBLIC / directory, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    (PUBLIC / 'scripts').mkdir()
    shutil.copyfile(runtime / 'scripts/run_all.py', PUBLIC / 'scripts/run_all.py')
    (PUBLIC / 'data/manifests').mkdir(parents=True)
    for name in ('pilot_inputs.tsv', 'reference_manifest.tsv'):
        shutil.copyfile(runtime / 'data/manifests' / name, PUBLIC / 'data/manifests' / name)
    for name in ('quality_control.json', 'functional_identification.json', 'evolutionary_tests.json',
                 'neighborhoods.json', 'validation.json', 'figure_manifest.json'):
        source = runtime / 'results/reports' / name
        target = PUBLIC / 'results/reports' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    cfg = load_config(runtime / 'config/analysis.yaml')
    cfg['tools']['mmseqs'] = '.tools/mmseqs/bin/mmseqs'
    (PUBLIC / 'analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    inputs = json.loads((runtime / 'update_inputs.json').read_text())
    write_json(PUBLIC / 'analysis_manifest.json', {'analysis': 'assay_reference_update_2026', 'completed_stages': list(range(1, 10)),
        'n_genomes': 20, 'n_cds': 2460, 'n_references': 39, 'scientific_sha256': scientific,
        'input_provenance': inputs, 'source_acquisition': 'Accession-versioned records verified against frozen source checksums'})
    artifact_hashes = {p.relative_to(PUBLIC).as_posix(): sha256(p) for p in sorted(PUBLIC.rglob('*')) if p.is_file()}
    write_json(PUBLIC / 'artifact_manifest.json', {'artifacts': artifact_hashes})

def run(args):
    runtime = ROOT / '.experiments' / args.workspace
    if not runtime.resolve().is_relative_to((ROOT / '.experiments').resolve()): raise ValueError('Workspace must remain inside .experiments')
    before = {name: scientific_files(root) for name, root in [('original', ROOT), ('experimental', ROOT / 'results/amendment_v1/experimental'), ('profile', BASE)]}
    if not args.resume:
        cfg = load_config(ROOT / 'config/analysis.yaml')
        cfg['tools']['mmseqs'] = str((ROOT / cfg['tools']['mmseqs']).resolve())
        prepare(runtime, cfg, args.cache_root.resolve() if args.cache_root else ROOT, args.offline)
    inputs = json.loads((runtime / 'update_inputs.json').read_text())
    if inputs['driver_sha256'] != sha256(__file__) or inputs['protocol_sha256'] != sha256(PROTOCOL):
        raise ValueError('Driver or protocol changed after preparing this analysis')
    for relative, digest in inputs['snapshot_sha256'].items():
        if sha256(runtime / relative) != digest: raise ValueError('Isolated source snapshot changed')
    command = [sys.executable, '-B', str(runtime / 'scripts/run_all.py'), '--threads', '4', '--offline']
    if args.resume: command += ['--resume-from', str(args.resume)]
    subprocess.run(command, cwd=runtime, check=True)
    if not args.no_export: export(runtime)
    after = {name: scientific_files(root) for name, root in [('original', ROOT), ('experimental', ROOT / 'results/amendment_v1/experimental'), ('profile', BASE)]}
    if before != after: raise ValueError('An earlier scientific panel changed during follow-up')
    if args.no_export and PUBLIC.exists():
        expected = json.loads((PUBLIC / 'analysis_manifest.json').read_text())['scientific_sha256']
        if scientific_files(runtime) != expected: raise ValueError('Replayed scientific results differ from the published panel')
        print('All 39 follow-up scientific artifacts match the published panel.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', default='lysterase_2026')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--cache-root', type=Path)
    parser.add_argument('--resume', type=int, choices=range(1, 10))
    parser.add_argument('--no-export', action='store_true', help='Replay without replacing published results')
    run(parser.parse_args())
