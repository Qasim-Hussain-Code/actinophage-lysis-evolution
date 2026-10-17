"""Amendment integrity tests use synthetic records and never execute a pipeline."""
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tarfile
from types import SimpleNamespace

import pytest
import yaml

from actinolysis.common import write_json, write_tsv


@pytest.fixture
def amendment(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('synthetic_amendment_driver', repo / 'scripts/run_amendment.py')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    monkeypatch.setattr(driver, 'ROOT', tmp_path)
    monkeypatch.setattr(driver, 'budget_check', lambda *a, **kw: {})
    monkeypatch.setattr(driver, 'verify_baseline', lambda: {
        'commit': '0' * 40, 'scientific_output_count': 39, 'scientific_sha256': {}})
    cfg = yaml.safe_load((repo / 'config/analysis.yaml').read_text())
    (tmp_path / 'config').mkdir()
    (tmp_path / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg))
    tool = tmp_path / cfg['tools']['mmseqs']
    tool.parent.mkdir(parents=True)
    tool.write_bytes(b'synthetic tool; never executed')
    for relative in ['scripts/run_amendment.py', 'scripts/run_all.py',
                     'src/actinolysis/common.py', 'src/actinolysis/statistics.py',
                     'src/actinolysis/reference_validation.py', 'src/actinolysis/amendment_feasibility.py']:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'# synthetic source only\r\n')
    reference = {'reference_id': 'Synthetic_LysA', 'category': 'lysin_a',
        'protein_accession': 'ABC00001.1', 'source_accession': 'NC_900001', 'gene': '1',
        'reference_url': 'https://example.invalid/synthetic-test-only', 'evidence_type': 'experimental',
        'expected_category': 'lysin_a', 'description': 'Synthetic test only; never biological evidence'}
    write_tsv(tmp_path / 'references/protein_controls.tsv', [reference])
    write_tsv(tmp_path / 'references/experimental_lysA_expansion.tsv', [], ['panel_action'])
    write_tsv(tmp_path / 'references/architecture_controls.tsv', [{'phage': 'Synthetic test only'}])
    raw = tmp_path / 'data/raw/references/NC_900001.gb'
    raw.parent.mkdir(parents=True)
    raw.write_bytes(b'LOCUS       synthetic\nVERSION     NC_900001.1\n//\n')
    cohort = []
    for number in range(20):
        accession = f'SYN{number:04}.1'
        path = tmp_path / f'data/raw/genbank/{accession}.gb'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f'LOCUS       synthetic\nVERSION     {accession}\n//\n'.encode())
        cohort.append({'accession': accession, 'path': path.relative_to(tmp_path).as_posix(),
                       'sha256': driver._digest(path.read_bytes())})
    write_tsv(tmp_path / 'data/manifests/pilot_inputs.tsv', cohort)
    selection = tmp_path / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(selection, [{'phage': 'Synthetic test only'}])
    write_json(tmp_path / 'data/manifests/phagesdb_snapshot.json', {
        'selection_snapshot': {'sha256': driver._hash(selection)}, 'pages': []})
    return driver


def freeze_reference(driver, destination):
    specs, _ = driver._catalogue('experimental')
    row = {**specs[0], 'resolved_source_accession': 'NC_900001.1',
           'resolved_protein_accession': 'ABC00001.1',
           'source_sha256': driver._digest((destination / 'data/raw/references/NC_900001.gb').read_bytes()),
           'protein_sha256': driver._digest(b'MGP'),
           'source_url': 'https://example.invalid/efetch?id=NC_900001.1'}
    write_tsv(destination / 'data/manifests/reference_manifest.tsv', [row])


def complete_synthetic_run(driver, destination, value='synthetic observation'):
    path = destination / 'results/tables/synthetic_test_only.tsv'
    write_tsv(path, [{'value': value}])
    receipts = {str(i): {'outputs': {'results/tables/synthetic_test_only.tsv': driver._hash(path)}}
                for i in range(1, 10)}
    write_json(destination / 'results/reports/pilot_stage_receipts.json', receipts)
    write_json(destination / 'results/reports/pilot_execution.json', {'status': 'success'})
    write_json(destination / 'results/reports/validation.json', {'checksum_status': 'pass'})


def snapshot(path):
    return {p.relative_to(path).as_posix(): (p.read_bytes(), p.stat().st_mtime_ns)
            for p in path.rglob('*') if p.is_file()}


def test_frozen_source_change_rejected_before_any_existing_file_mutates(amendment):
    driver = amendment
    destination = driver.prepare('experimental', offline=True)
    freeze_reference(driver, destination)
    before = snapshot(destination)
    (driver.ROOT / 'src/actinolysis/statistics.py').write_text('# changed source\n')
    with pytest.raises(ValueError, match='frozen functional work'):
        driver.prepare('experimental', offline=True, refresh_prefunctional=True)
    assert snapshot(destination) == before


def test_prefunctional_refresh_is_explicit_preserves_time_and_old_ledger(amendment):
    driver = amendment
    destination = driver.prepare('experimental', offline=True)
    provenance_path = destination / 'amendment_provenance.json'
    old = json.loads(provenance_path.read_text())
    write_json(destination / 'results/reports/pilot_stage_receipts.json', {'1': {}, '2': {}})
    (driver.ROOT / 'src/actinolysis/statistics.py').write_bytes(b'# changed source\r\n')
    before = snapshot(destination)
    with pytest.raises(ValueError, match='refresh-prefunctional'):
        driver.prepare('experimental', offline=True)
    assert snapshot(destination) == before
    driver.prepare('experimental', offline=True, refresh_prefunctional=True)
    new = json.loads(provenance_path.read_text())
    assert new['prepared_utc'] == old['prepared_utc']
    assert new['prefunctional_refreshes'][-1]['previous_source_sha256'] == old['source_sha256']
    assert json.loads((destination / 'results/reports/pilot_stage_receipts.json').read_text()) == {}
    assert (destination / 'src/actinolysis/statistics.py').read_bytes() == b'# changed source\n'
    assert not (destination / 'src/actinolysis/reference_validation.py').exists()
    assert not (destination / 'src/actinolysis/amendment_feasibility.py').exists()


def test_fresh_workspace_rehydrates_pinned_manifest_and_original_provenance(amendment):
    driver = amendment
    destination = driver.prepare('experimental', offline=True)
    freeze_reference(driver, destination)
    complete_synthetic_run(driver, destination)
    driver.export(destination, 'experimental')
    expected = (destination / 'amendment_provenance.json').read_bytes()
    assert destination.resolve().is_relative_to(driver.ROOT.resolve())
    shutil.rmtree(destination)
    # A clone may share an identical executable through a different path.
    shared = driver.ROOT / 'shared_tools/mmseqs'
    shared.parent.mkdir()
    shared.write_bytes((driver.ROOT / '.tools/mmseqs/bin/mmseqs').read_bytes())
    cfg_path = driver.ROOT / 'config/analysis.yaml'
    cfg = yaml.safe_load(cfg_path.read_text())
    cfg['tools']['mmseqs'] = str(shared)
    cfg_path.write_text(yaml.safe_dump(cfg))
    destination = driver.prepare('experimental', offline=True)
    rows = driver._reference_freeze(destination / 'data/manifests/reference_manifest.tsv',
                                     driver._catalogue('experimental')[0])
    assert rows[0]['resolved_source_accession'] == 'NC_900001.1'
    assert (destination / 'amendment_provenance.json').read_bytes() == expected


def test_published_catalogue_change_rejects_before_creating_destination(amendment):
    driver = amendment
    destination = driver.prepare('experimental', offline=True)
    freeze_reference(driver, destination)
    complete_synthetic_run(driver, destination)
    driver.export(destination, 'experimental')
    assert destination.resolve().is_relative_to(driver.ROOT.resolve())
    shutil.rmtree(destination)
    path = driver.ROOT / 'references/protein_controls.tsv'
    rows = driver._catalogue('experimental')[0]
    rows[0]['description'] = 'Changed scientific evidence'
    write_tsv(path, rows)
    with pytest.raises(ValueError, match='catalogue differs'):
        driver.prepare('experimental', offline=True)
    assert not destination.exists()


def test_export_preserves_ancillary_analysis_and_rejects_scientific_overwrite(amendment):
    driver = amendment
    destination = driver.prepare('experimental', offline=True)
    freeze_reference(driver, destination)
    complete_synthetic_run(driver, destination)
    driver.export(destination, 'experimental')
    target = driver.ROOT / 'results/amendment_v1/experimental'
    extra = target / 'reference_validation/synthetic_test_only.txt'
    extra.parent.mkdir()
    extra.write_text('separate analysis')
    driver.export(destination, 'experimental')
    assert extra.read_text() == 'separate analysis'
    manifest = driver._verify_artifacts(target)
    assert not any(p.startswith('reference_validation/') for p in manifest['artifacts'])
    before = snapshot(target)
    complete_synthetic_run(driver, destination, 'different scientific observation')
    with pytest.raises(ValueError, match='scientific tables/figures changed'):
        driver.export(destination, 'experimental')
    assert snapshot(target) == before


def test_raw_version_failure_precedes_any_preparation_write(amendment):
    driver = amendment
    raw = driver.ROOT / 'data/raw/genbank/SYN0000.1.gb'
    raw.write_bytes(b'VERSION     SYN0000.2\n')
    with pytest.raises(ValueError, match='raw-input checksum'):
        driver.prepare('experimental', offline=True)
    assert not (driver.ROOT / '.experiments').exists()


def test_git_baseline_lf_normalization_and_binary_integrity(tmp_path, monkeypatch):
    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('synthetic_baseline_check', repo / 'scripts/run_amendment.py')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    monkeypatch.setattr(driver, 'ROOT', tmp_path)
    archive = io.BytesIO()
    with tarfile.open(fileobj=archive, mode='w') as tar:
        for number in range(39):
            relative = f'results/tables/test{number}.tsv' if number else 'figures/test.png'
            data = b'column\nvalue\n' if number else b'binary\r\nimage'
            info = tarfile.TarInfo(relative)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
            path = tmp_path / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data.replace(b'\n', b'\r\n') if number else data)
    calls = []
    def git_archive(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(stdout=archive.getvalue())
    monkeypatch.setattr(driver.subprocess, 'run', git_archive)
    assert driver.verify_baseline()['scientific_output_count'] == 39
    assert len(calls) == 1
    (tmp_path / 'figures/test.png').write_bytes(b'binary\nimage')
    with pytest.raises(ValueError, match='Original pilot differs'):
        driver.verify_baseline()
