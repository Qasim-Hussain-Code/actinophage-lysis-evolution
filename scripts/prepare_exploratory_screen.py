#!/usr/bin/env python3
"""Freeze a private 50-record exploratory selection without acquiring sequences."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.acquisition import select_sample, validate_accession
from actinolysis.common import NA, load_config, read_tsv, sha256, write_json, write_tsv

SCREEN_SIZE, PILOT_SIZE, SEED = 50, 20, 20261016
GENERA = ['Mycobacterium', 'Gordonia']
CONTROLS = ['D29', 'Yvonnetastic']
PANEL = Path('results/lysterase_update/panel')
PROTOCOL = Path('docs/exploratory_screen_protocol.md')


def digest_object(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def new_workspace(root, workspace):
    root = Path(root).resolve()
    experiments = (root / '.experiments').resolve()
    destination = (experiments / workspace).resolve()
    if (not experiments.is_relative_to(root) or destination == experiments
            or not destination.is_relative_to(experiments)):
        raise ValueError('Workspace must be a new directory beneath .experiments')
    if destination.exists():
        raise ValueError('Workspace already exists; preserve it and choose a new name')
    return destination


def validate_panel(root):
    """Authenticate the current catalogue and reference manifest to its freeze."""
    panel = Path(root) / PANEL
    artifact_path = panel / 'artifact_manifest.json'
    artifacts = json.loads(artifact_path.read_text(encoding='utf-8'))['artifacts']
    names = ('references/protein_controls.tsv', 'data/manifests/reference_manifest.tsv', 'analysis_manifest.json')
    hashes = {}
    for name in names:
        expected = artifacts.get(name, '')
        if not re.fullmatch(r'[0-9a-f]{64}', expected) or sha256(panel / name) != expected:
            raise ValueError(f'Frozen reference-panel checksum mismatch: {name}')
        hashes[name] = expected
    catalogue = read_tsv(panel / names[0])
    references = read_tsv(panel / names[1])
    analysis = json.loads((panel / names[2]).read_text(encoding='utf-8'))
    if (len(catalogue) != 39 or len(references) != 39 or analysis.get('n_references') != 39
            or len({r['reference_id'] for r in catalogue}) != 39
            or len({r['reference_id'] for r in references}) != 39
            or {r['reference_id'] for r in catalogue} != {r['reference_id'] for r in references}):
        raise ValueError('Expected the unchanged 39-reference catalogue and manifest')
    frozen_catalogue = analysis['input_provenance']['snapshot_sha256']['references/protein_controls.tsv']
    if frozen_catalogue != hashes[names[0]]:
        raise ValueError('Reference catalogue differs from its frozen analysis provenance')
    return hashes, sha256(artifact_path)


def prepare(workspace, root=ROOT):
    root = Path(root).resolve()
    destination = new_workspace(root, workspace)
    cfg_path = root / 'config/analysis.yaml'
    cfg = load_config(cfg_path)
    if (cfg['seed'] != SEED or cfg['sample_sizes'].get('pilot') != PILOT_SIZE
            or cfg['sample_sizes'].get('intermediate') != SCREEN_SIZE
            or cfg['sampling']['host_genera'] != GENERA
            or cfg['sampling']['controls'] != CONTROLS):
        raise ValueError('Exploratory selection requires the frozen seed, 50/20 sizes, genera and controls')
    cfg['mode'] = 'intermediate'
    snapshot = root / 'data/manifests/eligibility_snapshot.tsv'
    provenance_path = root / 'data/manifests/phagesdb_snapshot.json'
    provenance = json.loads(provenance_path.read_text(encoding='utf-8'))
    source = provenance['selection_snapshot']
    expected = source.get('sha256', '')
    if (source.get('path') != 'data/manifests/eligibility_snapshot.tsv'
            or not re.fullmatch(r'[0-9a-f]{64}', expected) or sha256(snapshot) != expected):
        raise ValueError('Frozen eligibility metadata checksum mismatch')
    rows = read_tsv(snapshot)
    if 'rows' in source and source['rows'] != len(rows):
        raise ValueError('Frozen eligibility metadata row count mismatch')
    for row in rows:
        row['seq_finished'] = str(row['seq_finished']).lower() == 'true'
    selected = select_sample(rows, cfg)
    if (len(selected) != SCREEN_SIZE or len({r['accession'].split('.')[0] for r in selected}) != SCREEN_SIZE
            or len({r['phage'] for r in selected}) != SCREEN_SIZE):
        raise ValueError('Selection must contain exactly 50 unique records')
    pilot_path = root / 'data/manifests/pilot_inputs.tsv'
    pilot = read_tsv(pilot_path)
    identity = lambda row: (row['phage'], row['accession'].split('.')[0])
    if len(pilot) != PILOT_SIZE or [identity(r) for r in selected[:PILOT_SIZE]] != [identity(r) for r in pilot]:
        raise ValueError('Selected 50-record cohort is not nested in the frozen 20-record pilot')
    panel_hashes, panel_artifact_hash = validate_panel(root)
    frozen = {identity(row): row for row in pilot}
    ledger = []
    for rank, row in enumerate(selected, 1):
        accession = validate_accession(row['accession'])
        known = frozen.get(identity(row))
        cache_status = 'not_pilot'
        if known:
            validate_accession(known['accession'], version_required=True)
            if not re.fullmatch(r'[0-9a-f]{64}', known['sha256']):
                raise ValueError('Invalid frozen pilot checksum')
            cached = (root / known['path']).resolve()
            if not cached.is_relative_to(root):
                raise ValueError('Frozen pilot cache path escapes the project')
            cache_status = 'not_present'
            if cached.exists():
                if not cached.is_file() or sha256(cached) != known['sha256']:
                    raise ValueError('Cached pilot source checksum mismatch')
                cache_status = 'verified'
        ledger.append({'rank': rank, 'accession': accession, 'phage': row['phage'],
            'host_genus': row['host_genus'], 'cluster': row.get('cluster', NA),
            'subcluster': row.get('subcluster', NA), 'pilot': 'true' if known else 'false',
            'source_metadata_frame_sha256': expected,
            'frozen_pilot_accession': known['accession'] if known else NA,
            'frozen_pilot_sha256': known['sha256'] if known else NA,
            'pilot_cache_status': cache_status})
    code_names = ('scripts/prepare_exploratory_screen.py', 'src/actinolysis/acquisition.py', 'src/actinolysis/common.py')
    report = {'status': 'awaiting_sequence_acquisition', 'planned_selected_records': SCREEN_SIZE,
        'planned_frozen_pilot_records': PILOT_SIZE, 'sequence_acquisition_performed': False,
        'completed_analysis_stages': [], 'primary_cohort_modified': False, 'main_expansion_gate_modified': False,
        'protocol_sha256': sha256(root / PROTOCOL), 'configuration_file': 'config/analysis.yaml',
        'configuration_file_sha256': sha256(cfg_path), 'effective_configuration': cfg,
        'effective_configuration_sha256': digest_object(cfg), 'panel_root': PANEL.as_posix(),
        'panel_files_sha256': panel_hashes, 'panel_artifact_manifest_sha256': panel_artifact_hash,
        'eligibility_snapshot': 'data/manifests/eligibility_snapshot.tsv', 'eligibility_snapshot_sha256': expected,
        'phagesdb_snapshot_sha256': sha256(provenance_path), 'pilot_manifest_sha256': sha256(pilot_path),
        'code_sha256': {name: sha256(root / name) for name in code_names},
        'source_tree_sha256': {p.relative_to(root).as_posix(): sha256(p)
            for p in sorted((root / 'src').rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'},
        'interpretation': 'Private preparation only. No new genomes acquired or biological analysis performed.'}
    runner = root / 'scripts/run_exploratory_screen.py'
    if runner.is_file():
        report['runner_sha256'] = sha256(runner)
    destination.mkdir(parents=True, exist_ok=False)
    write_tsv(destination / 'selection.tsv', ledger)
    report['selection_sha256'] = sha256(destination / 'selection.tsv')
    write_json(destination / 'preparation.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', required=True, help='New private directory name beneath .experiments')
    args = parser.parse_args()
    report = prepare(args.workspace)
    print(json.dumps({'status': report['status'], 'planned_selected_records': SCREEN_SIZE,
                      'sequence_acquisition_performed': False}))


if __name__ == '__main__':
    main()
