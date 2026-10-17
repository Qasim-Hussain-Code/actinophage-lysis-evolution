#!/usr/bin/env python3
"""Reproduce an accession-frozen cohort without resampling database metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import load_config, read_tsv, sha256, utcnow, write_json, write_tsv


def cohort_frame(manifest, cfg):
    """Preserve the cohort and verify its deterministic genomic-analysis order."""
    from actinolysis.acquisition import select_sample
    if len(manifest) != 20 or len({r['accession'] for r in manifest}) != 20:
        raise ValueError('Expected 20 unique accession-versioned cohort genomes')
    rows = [{**row, 'seq_finished': True, 'length': 'not_available'} for row in manifest]
    expected = [(r['phage'], r['accession']) for r in manifest]
    if [(r['phage'], r['accession']) for r in select_sample(rows, cfg)] != expected:
        raise ValueError('The cohort order differs from the frozen sampling configuration')
    return rows


def copy_cached(cache, relative, digest, destination):
    if cache is None: return False
    for prefix in ('', '.experiments/amendment_experimental_v1', '.experiments/amendment_profile_v1'):
        source = cache / prefix / relative
        if source.is_file():
            if sha256(source) != digest:
                raise ValueError(f'Cached source checksum mismatch: {relative}')
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            return True
    return False


def compare_science(source, runtime):
    paths = sorted(p.relative_to(source) for directory in ('results/tables', 'figures')
                   for p in (source / directory).rglob('*') if p.is_file())
    if len(paths) != 39:
        raise ValueError('Expected 39 scientific tables and figure artifacts')
    rows = [{'path': p.as_posix(), 'expected_sha256': sha256(source / p),
             'actual_sha256': sha256(runtime / p)} for p in paths]
    return rows, [r['path'] for r in rows if r['actual_sha256'] != r['expected_sha256']]


def run(args):
    source = ROOT if args.panel == 'original' else ROOT / 'results/amendment_v1' / args.panel
    runtime = ROOT / '.experiments' / f'public_replay_{args.panel}'
    if runtime.exists() and any(runtime.iterdir()):
        raise ValueError('Replay workspace already exists; preserve it and use a fresh checkout')
    cfg = load_config(source / ('config/analysis.yaml' if args.panel == 'original' else 'analysis.yaml'))
    manifest = read_tsv(source / 'data/manifests/pilot_inputs.tsv')
    frame = cohort_frame(manifest, cfg)
    binary = args.mmseqs.resolve() if args.mmseqs else (ROOT / '.tools/mmseqs/bin/mmseqs').resolve()
    if not binary.is_file(): raise FileNotFoundError('MMseqs executable not found; run the bootstrap or supply --mmseqs')
    cache = args.cache_root.resolve() if args.cache_root else None
    for directory in ('src', 'references'):
        shutil.copytree(source / directory, runtime / directory, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    (runtime / 'scripts').mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source / 'scripts/run_all.py', runtime / 'scripts/run_all.py')
    (runtime / 'data/manifests').mkdir(parents=True, exist_ok=True)
    for name in ('pilot_inputs.tsv', 'reference_manifest.tsv'):
        shutil.copyfile(source / 'data/manifests' / name, runtime / 'data/manifests' / name)
    import yaml
    cfg['tools']['mmseqs'] = str(binary)
    (runtime / 'config').mkdir()
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    selection = runtime / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(selection, frame)
    write_json(runtime / 'data/manifests/phagesdb_snapshot.json', {
        'pages': [], 'rows': 20, 'selection_snapshot': {
            'path': 'data/manifests/eligibility_snapshot.tsv', 'rows': 20,
            'sha256': sha256(selection), 'basis': 'Frozen study cohort; no population resampling'},
        'method': 'accession_frozen_cohort_replay'})
    sources = {row['path']: row['sha256'] for row in manifest}
    for row in read_tsv(runtime / 'data/manifests/reference_manifest.tsv'):
        relative = f"data/raw/references/{row['source_accession']}.gb"
        if relative in sources and sources[relative] != row['source_sha256']:
            raise ValueError('Conflicting frozen source checksums')
        sources[relative] = row['source_sha256']
    cached = sum(copy_cached(cache, relative, digest, runtime / relative)
                 for relative, digest in sources.items())
    command = [sys.executable, '-B', str(runtime / 'scripts/run_all.py'), '--threads', str(args.threads)]
    if args.offline: command.append('--offline')
    subprocess.run(command, cwd=runtime, check=True)
    artifacts, differences = compare_science(source, runtime)
    execution = json.loads((runtime / 'results/reports/pilot_execution.json').read_text())
    if execution['status'] != 'success' or len(execution['stages']) != 9:
        raise ValueError('The replay did not complete all nine stages')
    report = {'status': 'pass' if not differences else 'fail', 'verified_utc': utcnow(),
              'panel': args.panel, 'method': 'Nine-stage accession-frozen cohort replay; no population resampling',
              'source_acquisition': 'Checksum-verified local cache' if args.offline else 'Version-pinned downloads with optional verified cache',
              'software': 'Active Python environment and specified MMseqs executable; no independent installation claim',
              'wrapper_sha256': sha256(__file__), 'mmseqs_sha256': sha256(binary),
              'copied_cache_records': cached, 'completed_stages': 9, 'compared_artifacts': len(artifacts),
              'differences': differences, 'artifacts': artifacts}
    write_json(ROOT / 'results/reports' / f'public_cohort_replay_{args.panel}.json', report)
    print(json.dumps({k:v for k,v in report.items() if k != 'artifacts'}))
    if differences: raise ValueError('Scientific outputs differ from the published analysis')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', choices=('original', 'experimental', 'profile'), default='original')
    parser.add_argument('--threads', type=int, choices=range(1, 5), default=4)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--cache-root', type=Path)
    parser.add_argument('--mmseqs', type=Path)
    run(parser.parse_args())


if __name__ == '__main__': main()
