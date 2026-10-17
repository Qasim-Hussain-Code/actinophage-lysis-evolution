#!/usr/bin/env python3
"""Run an isolated amendment while preserving the original pilot and freezes."""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.parse
import uuid
import yaml

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import budget_check, load_config, read_tsv, utcnow

BASELINE_COMMIT = 'v0.1.0-pilot'
BASELINE_OUTPUT_COUNT = 39
UNUSED_MODULES = {'amendment_feasibility.py', 'reference_validation.py'}
TEXT_SUFFIXES = {'.py', '.tsv', '.json', '.yaml', '.yml', '.md', '.svg'}
CORE_TOP_LEVEL = {'results', 'figures', 'data', 'references', 'src', 'scripts',
                  'amendment_provenance.json', 'analysis.yaml', 'artifact_manifest.json'}


def _bytes(path):
    """Match Git LF text checkouts without changing raw GenBank or binaries."""
    path = Path(path)
    data = path.read_bytes()
    return data.replace(b'\r\n', b'\n') if path.suffix in TEXT_SUFFIXES else data


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _hash(path):
    return _digest(_bytes(path))


def _config_signature(cfg):
    stable = json.loads(json.dumps(cfg))
    stable['tools']['mmseqs'] = '<runtime executable path>'
    return _digest(json.dumps(stable, sort_keys=True).encode('utf-8'))


def _json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode('utf-8')


def _tsv_bytes(rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode('utf-8')


def _inside(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    if path == root or not path.is_relative_to(root):
        raise ValueError(f'Path is outside the intended directory: {relative}')
    return path


def verify_baseline():
    """Compare all 39 scientific outputs with one in-memory Git archive."""
    result = subprocess.run(['git', 'archive', '--format=tar', BASELINE_COMMIT,
                             'results/tables', 'figures'], cwd=ROOT, check=True, capture_output=True)
    expected = {}
    with tarfile.open(fileobj=io.BytesIO(result.stdout)) as archive:
        commit = archive.pax_headers.get('comment', BASELINE_COMMIT)
        for member in archive.getmembers():
            if member.isfile():
                data = archive.extractfile(member).read()
                if Path(member.name).suffix in TEXT_SUFFIXES:
                    data = data.replace(b'\r\n', b'\n')
                expected[member.name] = _digest(data)
    if len(expected) != BASELINE_OUTPUT_COUNT:
        raise ValueError('Unexpected scientific-output count in the original Git baseline')
    actual_paths = {p.relative_to(ROOT).as_posix() for directory in ['results/tables', 'figures']
                    for p in (ROOT / directory).rglob('*') if p.is_file()}
    if actual_paths != set(expected):
        raise ValueError('Original pilot scientific-output paths differ from the Git baseline')
    for relative, digest in expected.items():
        if _hash(_inside(ROOT, relative)) != digest:
            raise ValueError(f'Original pilot differs from Git {BASELINE_COMMIT}: {relative}')
    return {'commit': commit, 'scientific_output_count': len(expected), 'scientific_sha256': expected}


def _catalogue(panel):
    base = ROOT / 'references/protein_controls.tsv'
    specifications = read_tsv(base)
    fields = list(specifications[0])
    curation = [base, ROOT / 'references/experimental_lysA_expansion.tsv']
    optional = ROOT / 'references/experimental_lysB_expansion.tsv'
    if optional.exists():
        curation.append(optional)
    if panel == 'profile':
        curation.append(ROOT / 'references/profile_supported_expansion.tsv')
    for path in curation[1:]:
        for row in read_tsv(path):
            if row.get('panel_action', 'add') != 'add':
                continue
            item = {key: row.get(key, '') for key in fields}
            if not item['description']:
                item['description'] = '; '.join(row.get(key, '') for key in
                    ['assay', 'assay_outcome', 'evidence_locator', 'caveat'])
            specifications.append(item)
    if len({r['reference_id'] for r in specifications}) != len(specifications):
        raise ValueError('Duplicate reference IDs in amended panel')
    if any(r['source_accession'].split('.')[0] in {'MZ028627', 'KU963248'} for r in specifications):
        raise ValueError('Held-out VanLee and Yvonnetastic must not enter new reference panels')
    return specifications, curation


def _reference_freeze(path, specifications):
    rows = read_tsv(path)
    if len(rows) != len(specifications) or len({r['reference_id'] for r in rows}) != len(rows):
        raise ValueError('Frozen reference catalogue has missing or duplicate records')
    for spec, row in zip(specifications, rows):
        if any(row.get(key) != value for key, value in spec.items()):
            raise ValueError('Frozen reference catalogue does not match curated specifications')
        for supplied, resolved in [('source_accession', 'resolved_source_accession'),
                                   ('protein_accession', 'resolved_protein_accession')]:
            wanted, found = spec[supplied], row.get(resolved, '')
            if not re.fullmatch(r'[A-Za-z0-9_]+\.\d+', found):
                raise ValueError(f'Frozen reference lacks an accession version: {resolved}')
            if wanted and (wanted != found if '.' in wanted else wanted != found.split('.')[0]):
                raise ValueError(f'Frozen reference accession mismatch: {wanted} vs {found}')
        if any(not re.fullmatch(r'[0-9a-f]{64}', row.get(key, ''))
               for key in ['source_sha256', 'protein_sha256']):
            raise ValueError('Frozen reference lacks valid sequence/source hashes')
        url_id = urllib.parse.parse_qs(urllib.parse.urlparse(row['source_url']).query).get('id')
        if url_id != [row['resolved_source_accession']]:
            raise ValueError('Frozen reference retrieval URL is not version pinned')
    return rows


def _verify_artifacts(root):
    manifest = json.loads((root / 'artifact_manifest.json').read_text(encoding='utf-8'))
    for relative, digest in manifest['artifacts'].items():
        path = _inside(root, relative)
        if not path.is_file() or _hash(path) != digest:
            raise ValueError(f'Amendment export checksum mismatch: {relative}')
    return manifest


def _check_snapshot(root, provenance):
    legacy = 'snapshot_sha256' not in provenance
    ledger = provenance.get('snapshot_sha256', provenance.get('source_sha256', {}))
    for relative, digest in ledger.items():
        path = _inside(root, relative)
        valid = path.is_file() and (_hash(path) == digest or
                                   (legacy and _digest(path.read_bytes()) == digest))
        if not valid:
            raise ValueError(f'Existing amendment snapshot was modified: {relative}')
    if provenance.get('config_scientific_sha256'):
        if _config_signature(load_config(root / 'config/analysis.yaml')) != provenance['config_scientific_sha256']:
            raise ValueError('Existing amendment scientific configuration was modified')


def prepare(panel, offline=False, refresh_prefunctional=False):
    if panel not in {'experimental', 'profile'}:
        raise ValueError('Unknown amendment panel')
    baseline = verify_baseline()
    cfg = load_config(ROOT / 'config/analysis.yaml')
    if cfg['mode'] != 'pilot' or cfg['sample_sizes']['pilot'] != 20:
        raise ValueError('The amendment requires the unchanged 20-genome pilot design')
    budget_check(cfg, 300_000_000)
    destination = ROOT / '.experiments' / f'amendment_{panel}_v1'
    exported = ROOT / 'results/amendment_v1' / panel
    specifications, curation = _catalogue(panel)
    sources = sorted(p for p in (ROOT / 'src/actinolysis').glob('*.py') if p.name not in UNUSED_MODULES)
    sources += [ROOT / 'scripts/run_all.py', ROOT / 'scripts/run_amendment.py',
                ROOT / 'references/architecture_controls.tsv']
    sources += [ROOT / 'data/manifests' / p for p in
                ['pilot_inputs.tsv', 'phagesdb_snapshot.json', 'eligibility_snapshot.tsv']]
    planned = {p.relative_to(ROOT).as_posix(): _bytes(p) for p in sources}
    if any(not data.strip() for relative, data in planned.items() if relative.endswith('.py')):
        raise ValueError('An amendment source module is empty; finish source changes before preparing')
    for source in curation:
        planned[f'references/curation_inputs/{source.name}'] = _bytes(source)
    planned['references/protein_controls.tsv'] = _tsv_bytes(specifications)
    cfg['resources']['max_project_bytes'] = 8_000_000_000
    cfg['resources']['effective_budget_bytes'] = 8_000_000_000
    executable = (ROOT / cfg['tools']['mmseqs']).resolve()
    tool_hash = _digest(executable.read_bytes())
    cfg['tools']['mmseqs'] = Path(os.path.relpath(executable, destination)).as_posix()
    planned['config/analysis.yaml'] = yaml.safe_dump(cfg, sort_keys=False).encode('utf-8')
    all_hashes = {relative: _digest(data) for relative, data in planned.items()}
    wanted = {relative: digest for relative, digest in all_hashes.items() if relative != 'config/analysis.yaml'}
    config_signature = _config_signature(cfg)

    prior_path = destination / 'amendment_provenance.json'
    prior = json.loads(prior_path.read_text(encoding='utf-8')) if prior_path.exists() else None
    if destination.exists() and any(destination.iterdir()) and prior is None:
        raise ValueError('Existing amendment directory lacks provenance; do not overwrite it')
    receipts_path = destination / 'results/reports/pilot_stage_receipts.json'
    receipts = json.loads(receipts_path.read_text()) if receipts_path.exists() else {}
    reference_path = destination / 'data/manifests/reference_manifest.tsv'
    functional_started = reference_path.exists() or any(int(key) >= 3 for key in receipts)
    if prior:
        _check_snapshot(destination, prior)
    changes = [relative for relative, digest in all_hashes.items()
               if not (destination / relative).is_file() or _hash(destination / relative) != digest]
    if prior and prior.get('tool_binary_sha256', tool_hash) != tool_hash:
        changes.append('tools/mmseqs:binary_identity')
    refresh = prior is not None and (bool(changes) or 'snapshot_sha256' not in prior)
    if refresh and (functional_started or not refresh_prefunctional):
        raise ValueError('Amendment snapshot changed: frozen functional work needs a new version; '
                         'before stage3 use --refresh-prefunctional explicitly')

    seed_reference, seed_provenance = None, None
    if exported.exists():
        artifacts = _verify_artifacts(exported)['artifacts']
        for required in ['references/protein_controls.tsv', 'data/manifests/reference_manifest.tsv',
                         'amendment_provenance.json']:
            if required not in artifacts:
                raise ValueError(f'Completed amendment export lacks a frozen input: {required}')
        if read_tsv(exported / 'references/protein_controls.tsv') != specifications:
            raise ValueError('Published amendment catalogue differs; use a new amendment version')
        seed_provenance = json.loads((exported / 'amendment_provenance.json').read_text())
        if (seed_provenance.get('snapshot_sha256') != wanted or
                seed_provenance.get('catalogue_sha256') != wanted['references/protein_controls.tsv'] or
                seed_provenance.get('config_scientific_sha256') != config_signature or
                seed_provenance.get('tool_binary_sha256') != tool_hash):
            raise ValueError('Published amendment source/config/input snapshot differs; use a new version')
        if prior is not None and prior != seed_provenance:
            raise ValueError('Local and published amendment provenance differ')
        seed_reference = exported / 'data/manifests/reference_manifest.tsv'
    frozen_rows = _reference_freeze(reference_path if reference_path.exists() else seed_reference,
                                     specifications) if reference_path.exists() or seed_reference else []
    if reference_path.exists() and seed_reference and _bytes(reference_path) != _bytes(seed_reference):
        raise ValueError('Local reference freeze differs from the published amendment')

    # Check every raw copy before writing any snapshot file.
    cohort = read_tsv(ROOT / 'data/manifests/pilot_inputs.tsv')
    if (len(cohort) != 20 or len({r['accession'] for r in cohort}) != 20 or
            any(not re.fullmatch(r'[A-Z][A-Z0-9_]*\.\d+', r['accession']) for r in cohort)):
        raise ValueError('Amendment requires exactly 20 unique frozen accession versions')
    metadata = json.loads(planned['data/manifests/phagesdb_snapshot.json'])
    if _digest(planned['data/manifests/eligibility_snapshot.tsv']) != metadata['selection_snapshot']['sha256']:
        raise ValueError('Frozen selection metadata checksum mismatch')
    raw_plan = []
    raw_specs = [(r['path'], r['sha256'], r['accession']) for r in cohort]
    pinned = {r['source_accession']: r for r in frozen_rows}
    for source in sorted({r['source_accession'] for r in specifications}):
        row = pinned.get(source, {})
        raw_specs.append((f'data/raw/references/{source}.gb', row.get('source_sha256'),
                          row.get('resolved_source_accession', source)))
    for relative, expected_hash, accession in raw_specs:
        source, target = _inside(ROOT, relative), _inside(destination, relative)
        if source.exists() and target.exists() and source.read_bytes() != target.read_bytes():
            raise ValueError(f'Parent and isolated raw records differ: {relative}')
        available = target if target.exists() else source
        if available.exists():
            data = available.read_bytes()
            if expected_hash and _digest(data) != expected_hash:
                raise ValueError(f'Frozen raw-input checksum mismatch: {relative}')
            version = re.search(rb'^VERSION\s+(\S+)', data, flags=re.MULTILINE)
            actual = version.group(1).decode() if version else ''
            if (not version or actual.split('.')[0] != accession.split('.')[0] or
                    ('.' in accession and actual != accession)):
                raise ValueError(f'Raw source accession/version mismatch: {relative}')
            if not target.exists():
                raw_plan.append((target, data))
        elif offline:
            raise FileNotFoundError(f'Offline amendment raw input missing: {relative}')

    provenance = prior or seed_provenance
    if provenance is None or refresh:
        history = list((prior or {}).get('prefunctional_refreshes', []))
        if refresh:
            history.append({'refreshed_utc': utcnow(), 'changed_paths': changes,
                'previous_source_sha256': prior.get('source_sha256', {}),
                'previous_snapshot_sha256': prior.get('snapshot_sha256', {}),
                'previous_curation_sha256': prior.get('curation_sha256', {}),
                'previous_catalogue_sha256': prior.get('catalogue_sha256'),
                'previous_runtime_config_sha256': _hash(destination / 'config/analysis.yaml'),
                'previous_stage_receipts': receipts})
        provenance = {'prepared_utc': (prior or {}).get('prepared_utc', utcnow()),
            'amendment': 'reference_expansion_v1', 'preparation_protocol': 2, 'panel': panel,
            'original_baseline_commit': baseline['commit'], 'baseline_verification': baseline,
            'cohort': 'unchanged frozen 20-genome pilot',
            'selection': 'published assay positives; optional representatives of published Gordonia candidate families',
            'context_used_for_function': False,
            'context_scope': 'Target assignment does not read genomic position; published seed curation used annotation and context.',
            'profile_tier': 'computational seed homology only, capped moderate; no direct target profile search',
            'generic_lipolytic_tier': 'assayed generic substrates only, capped moderate including exact matches',
            'source_sha256': {p.relative_to(ROOT).as_posix(): _digest(_bytes(p)) for p in sources},
            'curation_sha256': {p.relative_to(ROOT).as_posix(): _digest(_bytes(p)) for p in curation},
            'driver_sha256': wanted['scripts/run_amendment.py'], 'snapshot_sha256': wanted,
            'catalogue_sha256': wanted['references/protein_controls.tsv'],
            'initial_runtime_config_sha256': all_hashes['config/analysis.yaml'],
            'config_scientific_sha256': config_signature,
            'tool_binary_sha256': tool_hash,
            'tool_identity_policy': 'Exact executable hash freezes the build; its CLI version is recorded by stage6. Runtime executable paths may differ in a fresh clone.',
            'offline_at_first_preparation': offline,
            'prefunctional_refreshes': history}

    # All checks above are read-only. Only now materialize the approved plan.
    for relative, data in planned.items():
        target = _inside(destination, relative)
        if not target.exists() or target.read_bytes() != data:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
    for target, data in raw_plan:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    if seed_reference and not reference_path.exists():
        reference_path.parent.mkdir(parents=True, exist_ok=True)
        reference_path.write_bytes(_bytes(seed_reference))
    if not prior_path.exists() or refresh:
        prior_path.write_bytes(_json_bytes(provenance))
    if refresh:
        for name in UNUSED_MODULES:
            path = _inside(destination, f'src/actinolysis/{name}')
            if path.exists():
                path.unlink()
        if receipts_path.exists():
            receipts_path.write_bytes(_json_bytes({}))
        print('Explicit pre-functional refresh recorded; stages1-2 receipts invalidated. Rerun from stage1.', flush=True)
    return destination


def _science(root):
    return {p.relative_to(root).as_posix(): _hash(p) for directory in ['results/tables', 'figures']
            for p in (root / directory).rglob('*') if p.is_file()}


def export(destination, panel):
    verify_baseline()
    destination = Path(destination)
    provenance = json.loads((destination / 'amendment_provenance.json').read_text())
    _check_snapshot(destination, provenance)
    _reference_freeze(destination / 'data/manifests/reference_manifest.tsv',
                      read_tsv(destination / 'references/protein_controls.tsv'))
    receipts = json.loads((destination / 'results/reports/pilot_stage_receipts.json').read_text())
    execution = json.loads((destination / 'results/reports/pilot_execution.json').read_text())
    validation = json.loads((destination / 'results/reports/validation.json').read_text())
    if (set(receipts) != {str(i) for i in range(1, 10)} or execution.get('status') != 'success'
            or validation.get('checksum_status') != 'pass'):
        raise ValueError('Export requires successful stages1-9 and passing input checksums')
    for receipt in receipts.values():
        if not receipt.get('outputs'):
            raise ValueError('Stage receipt has no recorded outputs')
        for relative, digest in receipt['outputs'].items():
            if _hash(_inside(destination, relative)) != digest:
                raise ValueError(f'Completed stage output changed before export: {relative}')
    target = ROOT / 'results/amendment_v1' / panel
    if target.exists():
        _verify_artifacts(target)
        if _science(target) != _science(destination):
            raise ValueError('Amendment scientific tables/figures changed; export a new version')
        if json.loads((target / 'amendment_provenance.json').read_text()) != provenance:
            raise ValueError('Published amendment provenance cannot be overwritten')
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f'.{panel}.staging-', dir=target.parent))
    backup = target.parent / f'.{panel}.previous-{uuid.uuid4().hex}'
    _inside(ROOT, staging.relative_to(ROOT))
    _inside(ROOT, backup.relative_to(ROOT))
    try:
        for directory in ['results/tables', 'results/reports', 'figures', 'data/manifests',
                          'references', 'src', 'scripts']:
            for source in (destination / directory).rglob('*'):
                if source.is_file():
                    copied = staging / source.relative_to(destination)
                    copied.parent.mkdir(parents=True, exist_ok=True)
                    copied.write_bytes(_bytes(source))
        (staging / 'amendment_provenance.json').write_bytes(_bytes(destination / 'amendment_provenance.json'))
        (staging / 'analysis.yaml').write_bytes(_bytes(destination / 'config/analysis.yaml'))
        artifacts = {p.relative_to(staging).as_posix(): _hash(p)
                     for p in sorted(staging.rglob('*')) if p.is_file()}
        (staging / 'artifact_manifest.json').write_bytes(_json_bytes({
            'scope': 'core isolated amendment export only; ancillary top-level analyses excluded',
            'artifacts': artifacts}))
        if target.exists():
            for extra in target.iterdir():
                if extra.name not in CORE_TOP_LEVEL:
                    if extra.is_dir():
                        shutil.copytree(extra, staging / extra.name)
                    else:
                        shutil.copyfile(extra, staging / extra.name)
        _verify_artifacts(staging)
        if target.exists():
            os.replace(target, backup)
        try:
            os.replace(staging, target)
        except BaseException:
            if backup.exists():
                os.replace(backup, target)
            raise
        if backup.exists():
            shutil.rmtree(backup)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', choices=['experimental', 'profile'], required=True)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--refresh-prefunctional', action='store_true')
    parser.add_argument('--stop-after', type=int, choices=range(1, 10), default=9)
    parser.add_argument('--resume-from', type=int, choices=range(1, 10), default=1)
    args = parser.parse_args()
    if args.refresh_prefunctional and args.resume_from != 1:
        parser.error('--refresh-prefunctional requires --resume-from 1')
    destination = prepare(args.panel, args.offline, args.refresh_prefunctional)
    if args.prepare_only:
        print(destination)
        return
    command = [sys.executable, '-B', str(destination / 'scripts/run_all.py'),
               '--threads', '4', '--stop-after', str(args.stop_after), '--resume-from', str(args.resume_from)]
    if args.offline:
        command.append('--offline')
    subprocess.run(command, cwd=destination, check=True)
    if args.stop_after == 9:
        export(destination, args.panel)
    budget_check(load_config(ROOT / 'config/analysis.yaml'))


if __name__ == '__main__':
    main()
