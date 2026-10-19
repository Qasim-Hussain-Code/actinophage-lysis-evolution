"""Synthetic-only checks for the bounded exploratory selection freeze."""
import hashlib
import json

import pytest
import yaml

from actinolysis.acquisition import select_sample
from actinolysis.common import read_tsv, sha256, write_json, write_tsv
from scripts.prepare_exploratory_screen import PANEL, prepare


@pytest.fixture
def screen_root(tmp_path):
    cfg = {'mode': 'pilot', 'seed': 20261016,
        'sample_sizes': {'pilot': 20, 'intermediate': 50, 'full': 300},
        'sampling': {'host_genera': ['Mycobacterium', 'Gordonia'], 'controls': ['D29', 'Yvonnetastic']}}
    (tmp_path / 'config').mkdir()
    (tmp_path / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg), encoding='utf-8')
    rows = []
    for index in range(64):
        rows.append({'phage': ['D29', 'Yvonnetastic'][index] if index < 2 else f'SyntheticOnly{index}',
            'accession': f'SYN{index:06d}', 'host_genus': ['Mycobacterium', 'Gordonia'][index % 2],
            'cluster': f'SYNTHETIC_ONLY_{index // 4}', 'subcluster': 'not_available', 'seq_finished': True})
    frame = tmp_path / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(frame, rows)
    write_json(tmp_path / 'data/manifests/phagesdb_snapshot.json', {
        'selection_snapshot': {'path': 'data/manifests/eligibility_snapshot.tsv',
                               'sha256': sha256(frame), 'rows': len(rows)}})
    cfg['mode'] = 'intermediate'
    selected = select_sample(rows, cfg)
    payload = b'SYNTHETIC TEST ONLY; this is not a genome record'
    pilot = [{**r, 'accession': r['accession'] + '.1', 'sha256': hashlib.sha256(payload).hexdigest(),
              'path': f"data/raw/genbank/{r['accession']}.1.gb"} for r in selected[:20]]
    write_tsv(tmp_path / 'data/manifests/pilot_inputs.tsv', pilot)
    cached = tmp_path / pilot[0]['path']
    cached.parent.mkdir(parents=True)
    cached.write_bytes(payload)
    catalogue = [{'reference_id': f'SYNTHETIC_ONLY_REF_{i}'} for i in range(39)]
    panel = tmp_path / PANEL
    write_tsv(panel / 'references/protein_controls.tsv', catalogue)
    write_tsv(panel / 'data/manifests/reference_manifest.tsv', catalogue)
    write_json(panel / 'analysis_manifest.json', {'n_references': 39,
        'input_provenance': {'snapshot_sha256': {
            'references/protein_controls.tsv': sha256(panel / 'references/protein_controls.tsv')}}})
    names = ('references/protein_controls.tsv', 'data/manifests/reference_manifest.tsv', 'analysis_manifest.json')
    write_json(panel / 'artifact_manifest.json', {'artifacts': {name: sha256(panel / name) for name in names}})
    for name in ('docs/exploratory_screen_protocol.md', 'scripts/prepare_exploratory_screen.py',
                 'src/actinolysis/acquisition.py', 'src/actinolysis/common.py'):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('SYNTHETIC TEST ONLY provenance fixture\n', encoding='utf-8')
    return tmp_path


def test_preparation_freezes_exactly_50_without_acquiring_or_changing_sources(screen_root):
    before = {p.relative_to(screen_root): sha256(p) for p in screen_root.rglob('*') if p.is_file()}
    report = prepare('synthetic_screen', root=screen_root)
    workspace = screen_root / '.experiments/synthetic_screen'
    ledger = read_tsv(workspace / 'selection.tsv')
    assert len(ledger) == 50
    assert [r['rank'] for r in ledger] == [str(i) for i in range(1, 51)]
    assert [r['pilot'] for r in ledger] == ['true'] * 20 + ['false'] * 30
    assert ledger[0]['pilot_cache_status'] == 'verified'
    assert ledger[1]['pilot_cache_status'] == 'not_present'
    assert report['status'] == 'awaiting_sequence_acquisition'
    assert report['effective_configuration']['mode'] == 'intermediate'
    assert report['completed_analysis_stages'] == []
    assert report['sequence_acquisition_performed'] is False
    assert set(p.name for p in workspace.iterdir()) == {'selection.tsv', 'preparation.json'}
    assert report['selection_sha256'] == sha256(workspace / 'selection.tsv')
    assert json.loads((workspace / 'preparation.json').read_text()) == report
    assert all(sha256(screen_root / name) == digest for name, digest in before.items())


def test_preparation_rejects_changed_eligibility_frame_before_output(screen_root):
    frame = screen_root / 'data/manifests/eligibility_snapshot.tsv'
    frame.write_bytes(frame.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='eligibility metadata checksum'):
        prepare('changed_frame', root=screen_root)
    assert not (screen_root / '.experiments/changed_frame').exists()


def test_preparation_rejects_nonnested_frozen_pilot(screen_root):
    path = screen_root / 'data/manifests/pilot_inputs.tsv'
    pilot = read_tsv(path)
    pilot[2], pilot[3] = pilot[3], pilot[2]
    write_tsv(path, pilot)
    with pytest.raises(ValueError, match='not nested'):
        prepare('changed_pilot', root=screen_root)
    assert not (screen_root / '.experiments/changed_pilot').exists()


@pytest.mark.parametrize('workspace', ['../escape', '../.experiments', '.'])
def test_preparation_rejects_workspace_escape_or_root(screen_root, workspace):
    with pytest.raises(ValueError, match='beneath .experiments'):
        prepare(workspace, root=screen_root)


def test_preparation_preserves_even_an_empty_existing_destination(screen_root):
    existing = screen_root / '.experiments/existing'
    existing.mkdir(parents=True)
    with pytest.raises(ValueError, match='already exists'):
        prepare('existing', root=screen_root)
    assert not list(existing.iterdir())


def test_preparation_rejects_an_absolute_destination_outside_experiments(screen_root):
    with pytest.raises(ValueError, match='beneath .experiments'):
        prepare(str(screen_root / 'outside'), root=screen_root)


def test_preparation_rejects_changed_reference_panel(screen_root):
    path = screen_root / PANEL / 'references/protein_controls.tsv'
    path.write_bytes(path.read_bytes().replace(b'REF_0', b'REF_CHANGED'))
    with pytest.raises(ValueError, match='reference-panel checksum'):
        prepare('changed_panel', root=screen_root)


def test_preparation_requires_the_full_planned_selection(screen_root):
    frame = screen_root / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(frame, read_tsv(frame)[:49])
    provenance_path = screen_root / 'data/manifests/phagesdb_snapshot.json'
    provenance = json.loads(provenance_path.read_text())
    provenance['selection_snapshot'].update(sha256=sha256(frame), rows=49)
    write_json(provenance_path, provenance)
    with pytest.raises(ValueError, match='exactly 50|Requested 50 genomes'):
        prepare('too_small', root=screen_root)


def test_preparation_rejects_a_changed_cached_pilot_source(screen_root):
    pilot = read_tsv(screen_root / 'data/manifests/pilot_inputs.tsv')
    (screen_root / pilot[0]['path']).write_bytes(b'SYNTHETIC CHANGED SOURCE')
    with pytest.raises(ValueError, match='Cached pilot source checksum'):
        prepare('changed_cache', root=screen_root)


@pytest.mark.parametrize('field,value', [('seed', 1), ('intermediate', 49), ('pilot', 19)])
def test_preparation_rejects_changed_sampling_settings(screen_root, field, value):
    path = screen_root / 'config/analysis.yaml'
    cfg = yaml.safe_load(path.read_text())
    if field == 'seed':
        cfg[field] = value
    else:
        cfg['sample_sizes'][field] = value
    path.write_text(yaml.safe_dump(cfg), encoding='utf-8')
    with pytest.raises(ValueError, match='frozen seed, 50/20 sizes'):
        prepare('changed_settings', root=screen_root)
