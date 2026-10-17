#!/usr/bin/env python3
"""Authenticate both cached-source amendment replays against a Git publication."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import load_config, sha256, utcnow, write_json

PANELS = ('experimental', 'profile')
CONFIG_FIELDS = ['mode', 'seed', 'sample_sizes', 'sampling', 'identification', 'architecture', 'neighborhood', 'statistics']
VERSION = re.compile(rb'^VERSION\s+(\S+)', re.MULTILINE)
ACCESSION = re.compile(r'[A-Za-z][A-Za-z0-9_]*\.\d+')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def git(root, *args):
    result = subprocess.run(['git', '-C', str(root), *args], capture_output=True, check=False)
    if result.returncode:
        raise ValueError('A Git publication or clone precondition failed')
    return result.stdout


def inside(root, relative):
    path = (root / relative).resolve()
    if path == root or not path.is_relative_to(root):
        raise ValueError('An evidence path leaves its declared analysis directory')
    return path


def tsv(data):
    return list(csv.DictReader(io.StringIO(data.decode('utf-8')), delimiter='\t'))


def published_snapshot(root, commit):
    resolved = git(root, 'rev-parse', '--verify', commit + '^{commit}').decode().strip()
    archive = git(root, 'archive', '--format=tar', resolved, 'results/amendment_v1', 'results/tables', 'figures')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        for member in stream.getmembers():
            if member.isfile():
                files[member.name] = stream.extractfile(member).read()
    return resolved, files


def clone_identity(clone, commit):
    head = git(clone, 'rev-parse', 'HEAD').decode().strip()
    tracking = git(clone, 'rev-parse', 'refs/remotes/origin/main').decode().strip()
    remote = git(clone, 'remote', 'get-url', 'origin').decode().strip()
    host = urlparse(remote).hostname
    if remote.startswith('git@github.com:'):
        host = 'github.com'
    if head != commit or tracking != commit or host != 'github.com':
        raise ValueError('Replay must use the declared GitHub-cloned publication commit')
    return {'clone_commit': head, 'origin_main_at_clone': tracking, 'origin_host': host,
            'remote_url_policy': 'URL and account identifiers are omitted from this report'}


def scientific_paths(paths):
    return {path for path in paths if path.startswith(('results/tables/', 'figures/'))}


def original_science(published, clone, workspace):
    paths = scientific_paths(published)
    expected = {path: digest(published[path]) for path in paths}
    if len(expected) != 39:
        raise ValueError('Published original pilot does not contain exactly 39 scientific artifacts')
    artifacts = [{'path': path, 'expected_sha256': checksum, 'clone_sha256': sha256(inside(clone, path)),
                  'workspace_sha256': sha256(inside(workspace, path))} for path, checksum in sorted(expected.items())]
    differences = [row['path'] for row in artifacts if row['clone_sha256'] != row['expected_sha256'] or row['workspace_sha256'] != row['expected_sha256']]
    return {'status': 'pass' if not differences else 'fail', 'compared_artifacts': 39,
            'differences': differences, 'artifacts': artifacts}


def manifest_files(base):
    """Authenticate an export, including scientific files, before comparing it."""
    manifest = json.loads((base / 'artifact_manifest.json').read_text(encoding='utf-8'))
    for relative, checksum in manifest['artifacts'].items():
        if sha256(inside(base, relative)) != checksum:
            raise ValueError(f'Export artifact checksum mismatch: {relative}')
    return manifest['artifacts']


def fasta(path):
    result = {}; identifier = None
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            identifier = line[1:].split()[0]
            if identifier in result:
                raise ValueError('Duplicate reference FASTA identifier')
            result[identifier] = ''
        elif line:
            if identifier is None:
                raise ValueError('Malformed reference FASTA')
            result[identifier] += line.strip()
    return result


def cached_source(cache, panel, relative):
    for prefix in ['', f'.experiments/amendment_{panel}_v1/']:
        path = inside(cache, prefix + relative)
        if path.is_file():
            return path
    raise ValueError(f'Verified source cache record is unavailable: {relative}')


def receipt_checks(runtime, cfg):
    receipts = json.loads((runtime / 'results/reports/pilot_stage_receipts.json').read_text())
    execution = json.loads((runtime / 'results/reports/pilot_execution.json').read_text())
    validation = json.loads((runtime / 'results/reports/validation.json').read_text())
    if (set(receipts) != {str(number) for number in range(1, 10)}
            or execution.get('status') != 'success'
            or [row['stage'] for row in execution.get('stages', [])] != list(range(1, 10))
            or execution.get('config', {}).get('offline') is not True
            or validation.get('checksum_status') != 'pass'):
        raise ValueError('Replay requires successful offline stages 1-9 and passing source checksums')
    fingerprint = digest(json.dumps({key: cfg[key] for key in CONFIG_FIELDS}, sort_keys=True).encode())
    count = 0
    for stage in execution['stages']:
        receipt = receipts[str(stage['stage'])]
        if {key: value for key, value in stage.items() if key != 'stage'} != receipt:
            raise ValueError('Execution log and completed stage receipts disagree')
        if receipt.get('scientific_config_sha256') != fingerprint or not receipt.get('outputs') or not receipt.get('source_sha256'):
            raise ValueError('Replay stage receipt lacks scientific configuration or source provenance')
        for relative, checksum in {**receipt['outputs'], **receipt['source_sha256']}.items():
            if sha256(inside(runtime, relative)) != checksum:
                raise ValueError(f'Replay receipt checksum mismatch: {relative}')
            count += 1
    return execution, count


def verify_panel(panel, published, clone, cache):
    prefix = f'results/amendment_v1/{panel}/'
    expected = {path[len(prefix):]: data for path, data in published.items() if path.startswith(prefix)}
    if 'artifact_manifest.json' not in expected:
        raise ValueError(f'Published panel missing from the declared Git commit: {panel}')
    manifest = json.loads(expected['artifact_manifest.json'])['artifacts']
    for path, checksum in manifest.items():
        if path not in expected or digest(expected[path]) != checksum:
            raise ValueError(f'Committed export manifest mismatch: {panel}/{path}')
    science = scientific_paths(manifest)
    if len(science) != 39:
        raise ValueError(f'Published panel does not contain exactly 39 scientific artifacts: {panel}')
    runtime = clone / '.experiments' / f'amendment_{panel}_v1'
    export = clone / 'results/amendment_v1' / panel
    actual_export = manifest_files(export)
    actual_runtime_science = {path.relative_to(runtime).as_posix() for directory in ['results/tables', 'figures']
                              for path in (runtime / directory).rglob('*') if path.is_file()}
    if scientific_paths(actual_export) != science or actual_runtime_science != science:
        raise ValueError(f'Replay scientific file set differs from the publication: {panel}')
    cfg = load_config(runtime / 'config/analysis.yaml')
    import yaml
    expected_cfg = yaml.safe_load(expected['analysis.yaml'])
    def normalized(value):
        value = json.loads(json.dumps(value)); value['tools']['mmseqs'] = '<RUNTIME_MMSEQS_PATH>'
        return value
    if normalized(cfg) != normalized(expected_cfg):
        raise ValueError(f'Replay configuration changed beyond the allowed executable path: {panel}')
    execution, checked = receipt_checks(runtime, cfg)
    if execution['environment']['python'] != json.loads(expected['results/reports/pilot_execution.json'])['environment']['python']:
        raise ValueError('Replay Python version differs from the published pinned environment')
    expected_environment = json.loads(expected['results/reports/pilot_execution.json'])['environment']
    if execution['environment']['packages'] != expected_environment['packages']:
        raise ValueError('Replay core package versions differ from the published pinned environment')
    provenance = json.loads(expected['amendment_provenance.json'])
    if (json.loads((runtime / 'amendment_provenance.json').read_text()) != provenance
            or json.loads((export / 'amendment_provenance.json').read_text()) != provenance):
        raise ValueError(f'Replay provenance differs from the published preparation freeze: {panel}')
    if not provenance.get('snapshot_sha256'):
        raise ValueError('Published preparation provenance lacks its frozen source/input snapshot')
    for relative, checksum in provenance['snapshot_sha256'].items():
        if sha256(inside(runtime, relative)) != checksum:
            raise ValueError(f'Replay preparation snapshot changed: {relative}')
    for relative in ['data/manifests/pilot_inputs.tsv', 'data/manifests/reference_manifest.tsv']:
        if inside(runtime, relative).read_bytes() != expected[relative]:
            raise ValueError(f'Replay frozen manifest changed: {relative}')
    tool = Path(cfg['tools']['mmseqs'])
    tool = tool if tool.is_absolute() else runtime / tool
    if sha256(tool.resolve()) != provenance['tool_binary_sha256']:
        raise ValueError('Replay MMseqs executable differs from the pinned published binary')
    source_specs = {}
    for row in tsv(expected['data/manifests/pilot_inputs.tsv']):
        source_specs[row['path']] = (row['accession'], row['sha256'])
    references = tsv(expected['data/manifests/reference_manifest.tsv'])
    for row in references:
        relative = f"data/raw/references/{row['source_accession']}.gb"
        value = row['resolved_source_accession'], row['source_sha256']
        if relative in source_specs and source_specs[relative] != value:
            raise ValueError('Conflicting frozen reference-source provenance')
        source_specs[relative] = value
    source_checks = []
    for relative, (accession, checksum) in sorted(source_specs.items()):
        actual = inside(runtime, relative); cached = cached_source(cache, panel, relative)
        match = VERSION.search(actual.read_bytes())
        if (not ACCESSION.fullmatch(accession) or not match or match.group(1).decode() != accession
                or sha256(actual) != checksum or sha256(cached) != checksum):
            raise ValueError(f'Cached or replayed version-pinned source mismatch: {relative}')
        source_checks.append({'path': relative, 'accession_version': accession, 'expected_sha256': checksum,
                              'cache_path': cached.relative_to(cache).as_posix(),
                              'cached_sha256': checksum, 'replayed_sha256': checksum, 'pass': True})
    proteins = fasta(runtime / 'data/processed/references.faa')
    if set(proteins) != {row['reference_id'] for row in references}:
        raise ValueError('Replay reference protein set differs from the published freeze')
    for row in references:
        sequence = proteins[row['reference_id']]
        if (not ACCESSION.fullmatch(row.get('resolved_protein_accession', ''))
                or len(sequence) != int(row['length_aa']) or digest(sequence.encode()) != row['protein_sha256']):
            raise ValueError('Replay reference protein sequence differs from the published freeze')
    clustering = json.loads((runtime / 'results/reports/neighborhoods.json').read_text())['family_clustering']
    if not clustering or any(row.get('cache_reused') is not False for row in clustering):
        raise ValueError('Replay must recompute family clustering without a pre-existing family cache')
    artifacts = [{'path': path, 'expected_sha256': manifest[path],
                  'runtime_sha256': sha256(inside(runtime, path)), 'export_sha256': sha256(inside(export, path))}
                 for path in sorted(science)]
    differences = [row['path'] for row in artifacts if row['runtime_sha256'] != row['expected_sha256'] or row['export_sha256'] != row['expected_sha256']]
    return {'status': 'pass' if not differences else 'fail', 'compared_artifacts': 39,
            'differences': differences, 'artifacts': artifacts, 'source_checks': source_checks,
            'reference_proteins_checked': len(references), 'receipt_file_checks': checked,
            'completed_stages': list(range(1, 10)), 'execution_seconds': execution['elapsed_seconds'],
            'peak_child_rss_bytes': execution['peak_child_rss_bytes'],
            'core_environment': execution['environment'], 'mmseqs_binary_sha256': provenance['tool_binary_sha256'],
            'family_clustering_cache_reused': False}


def verify(clone, published_root, cache_root, out, expected_commit=None):
    clone = Path(clone).resolve(); published_root = Path(published_root).resolve(); cache_root = Path(cache_root).resolve()
    out = Path(out).resolve()
    if clone == published_root or not out.is_relative_to(published_root) or out.is_relative_to(clone):
        raise ValueError('Use a separate clone and an output report beneath the outer publication workspace')
    commit, published = published_snapshot(published_root, expected_commit or 'HEAD')
    identity = clone_identity(clone, commit)
    reports = {panel: verify_panel(panel, published, clone, cache_root) for panel in PANELS}
    baseline = original_science(published, clone, published_root)
    report = {'status': 'pass' if baseline['status'] == 'pass' and all(row['status'] == 'pass' for row in reports.values()) else 'fail',
              'verified_utc': utcnow(), **identity, 'published_commit': commit,
              'method': 'GitHub clone replay; checksum-verified raw GenBank cache copied locally; both nine-stage panels rerun offline with uncached family clustering.',
              'software': 'The existing pinned WSL Python environment and exact MMseqs executable were reused. Software was not independently installed.',
              'scope': 'Seventy-eight amendment scientific artifacts and thirty-nine unchanged original-pilot artifacts, plus version-pinned raw sources, reference proteins, stage outputs/source hashes and successful execution receipts. Timing, resource and sanitized publication-metadata bytes need not match.',
              'source_acquisition': 'Copied verified cache, not fresh downloads.',
              'original_pilot': baseline, 'panels': reports, 'verifier_sha256': sha256(Path(__file__))}
    write_json(out, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clone', type=Path, required=True)
    parser.add_argument('--published-root', type=Path, default=ROOT)
    parser.add_argument('--cache-root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/reports/amendment_replay.json')
    parser.add_argument('--expected-commit')
    args = parser.parse_args(argv)
    try:
        report = verify(args.clone, args.published_root, args.cache_root, args.out, args.expected_commit)
        print(json.dumps({'status': report['status'], 'published_commit': report['published_commit'],
                          'amendment_artifacts': sum(row['compared_artifacts'] for row in report['panels'].values()),
                          'original_pilot_artifacts': report['original_pilot']['compared_artifacts']}))
        return 0 if report['status'] == 'pass' else 1
    except Exception as error:
        print(json.dumps({'status': 'blocked', 'error_type': type(error).__name__,
                          'message': str(error) if isinstance(error, ValueError) else 'Replay verification failed; inspect the declared input set'}))
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
