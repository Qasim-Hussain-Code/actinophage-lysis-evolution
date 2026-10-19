"""Synthetic-only integrity tests for the isolated five-stage driver."""
from types import SimpleNamespace
import urllib.error

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from actinolysis.common import sha256, write_json, write_tsv
from scripts import run_exploratory_screen as screen


def cohort():
    ledger, acquired = [], []
    for index in range(50):
        base = f'SYN{index:06d}'
        ledger.append({'phage': f'SyntheticTestOnly{index}', 'accession': base,
            'frozen_pilot_accession': base + '.1', 'frozen_pilot_sha256': 'a' * 64})
        acquired.append({'phage': f'SyntheticTestOnly{index}', 'accession': base + '.1',
                         'sha256': 'a' * 64, 'path': f'data/raw/genbank/{base}.1.gb'})
    return ledger, acquired


def test_acquisition_requires_all_frozen_records_and_preserves_pilot_versions(tmp_path):
    ledger, acquired = cohort()
    manifest = tmp_path / 'data/manifests/intermediate_inputs.tsv'
    write_tsv(manifest, acquired)
    screen.check_acquisition(tmp_path, ledger)
    acquired[0]['accession'] = acquired[0]['accession'].replace('.1', '.2')
    write_tsv(manifest, acquired)
    with pytest.raises(ValueError, match='original pilot accession version'):
        screen.check_acquisition(tmp_path, ledger)


def test_acquisition_rejects_replacements_or_underfilled_selection(tmp_path):
    ledger, acquired = cohort()
    write_tsv(tmp_path / 'data/manifests/intermediate_inputs.tsv', acquired[:-1])
    with pytest.raises(ValueError, match='frozen 50-record selection'):
        screen.check_acquisition(tmp_path, ledger)


def test_outer_budget_check_accounts_for_original_root_and_download_bound(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(screen, 'budget_check', lambda cfg, proposed, root: calls.append((cfg, proposed, root)))
    cfg = {'resources': {'synthetic_test_only': True}}
    screen.root_guard(cfg, tmp_path, 90_000_000)
    assert calls == [(cfg, 90_000_000, tmp_path)]
    write_json(tmp_path / 'results/reports/software_storage.json', {'immutable_directories': {'.experiments': 100}})
    with pytest.raises(ValueError, match='cannot exclude the exploratory workspace'):
        screen.root_guard(cfg, tmp_path)


def test_worker_cannot_execute_from_primary_source_tree(tmp_path, monkeypatch):
    prepared = tmp_path / '.experiments/synthetic_test_only'
    monkeypatch.setattr(screen, 'verify_preparation', lambda root, workspace: (prepared, {}, []))
    with pytest.raises(ValueError, match='isolated runtime copy'):
        screen.worker(tmp_path, 'synthetic_test_only', 1)


def test_runtime_detects_changed_input_before_any_stage(tmp_path):
    prepared = tmp_path / '.experiments/synthetic_test_only'
    runtime = prepared / 'runtime'
    code = runtime / 'src/actinolysis/quality.py'
    code.parent.mkdir(parents=True)
    code.write_text('SYNTHETIC TEST ONLY frozen source\n')
    write_json(prepared / 'preparation.json', {'synthetic_test_only': True})
    write_json(runtime / 'execution_inputs.json', {
        'preparation_sha256': sha256(prepared / 'preparation.json'),
        'selection_sha256': 'a' * 64, 'driver_sha256': sha256(screen.__file__),
        'snapshot_sha256': {'src/actinolysis/quality.py': sha256(code)}})
    code.write_text('SYNTHETIC TEST ONLY changed source\n')
    with pytest.raises(ValueError, match='Isolated execution input changed'):
        screen.verify_runtime(tmp_path, prepared, runtime, {'selection_sha256': 'a' * 64})


def test_worker_scope_has_no_neighborhood_or_inference_stages():
    assert set(screen.STAGES) == {1, 2, 3, 4, 5}
    assert screen.STAGES[4] == ('identification', 'evidence')
    assert all(module not in {'neighborhoods', 'statistics', 'validation', 'figures'}
               for module, _ in screen.STAGES.values())


def control_fixture(tmp_path, complex_control=False):
    runtime = tmp_path / '.experiments/synthetic_test_only/runtime'
    panel = tmp_path / 'synthetic_panel'
    artifacts = {}
    for relative, phage, scope, architecture in [
            ('results/tables/architectures.tsv', 'D29', 'high_moderate', 'closely_linked'),
            ('results/tables/architecture_sensitivity.tsv', 'Yvonnetastic', 'annotation', 'separated')]:
        row = {'phage': phage, 'scope': scope, 'threshold': '4', 'accession': 'SYN000001.1',
               'architecture': architecture, 'confidence': 'high', 'lysin_a': 'SYNTHETIC_A',
               'lysin_b': 'SYNTHETIC_B', 'intervening_cds': '1', 'gap_bp': '100',
               'lysin_a_count': '1', 'lysin_b_count': '1'}
        if phage == 'Yvonnetastic' and complex_control:
            row.update(architecture='complex', lysin_a='not_available', lysin_b='not_available',
                       intervening_cds='not_available', gap_bp='not_available', lysin_a_count='2')
        write_tsv(panel / relative, [row])
        write_tsv(runtime / relative, [row])
        artifacts[relative] = sha256(panel / relative)
    pairs = [{'phage': phage, 'scope': scope, 'lysin_a': 'SYNTHETIC_A', 'lysin_b': 'SYNTHETIC_B',
              'intervening_cds': distance, 'gap_bp': '100', 'a_start': '100', 'b_start': '500'}
             for phage, scope, distance in [('D29', 'high_moderate', '1'), ('Yvonnetastic', 'annotation', '113')]]
    if complex_control:
        pairs.append({**pairs[-1], 'lysin_a': 'SYNTHETIC_A2', 'intervening_cds': '112', 'a_start': '200'})
    relative = 'results/tables/architecture_pairs.tsv'
    write_tsv(panel / relative, pairs)
    write_tsv(runtime / relative, pairs)
    artifacts[relative] = sha256(panel / relative)
    write_json(panel / 'artifact_manifest.json', {'artifacts': artifacts})
    return runtime


def test_control_contradiction_stops_exploratory_completion(tmp_path):
    runtime = control_fixture(tmp_path)
    screen.check_controls(tmp_path, runtime, {'panel_root': 'synthetic_panel'})
    path = runtime / 'results/tables/architectures.tsv'
    path.write_text(path.read_text().replace('closely_linked', 'unresolved'))
    with pytest.raises(ValueError, match='D29 control architecture is contradicted'):
        screen.check_controls(tmp_path, runtime, {'panel_root': 'synthetic_panel'})


def test_complex_control_cannot_pass_after_losing_one_component_pair(tmp_path):
    from actinolysis.common import read_tsv
    runtime = control_fixture(tmp_path, complex_control=True)
    screen.check_controls(tmp_path, runtime, {'panel_root': 'synthetic_panel'})
    path = runtime / 'results/tables/architecture_pairs.tsv'
    pairs = read_tsv(path)
    write_tsv(path, pairs[:-1])
    with pytest.raises(ValueError, match='Yvonnetastic control component pairs are contradicted'):
        screen.check_controls(tmp_path, runtime, {'panel_root': 'synthetic_panel'})


def write_source(path, accession):
    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(Seq('ACGT' * 25), id=accession,
                       description='Synthetic test-only phage, complete genome')
    record.annotations['molecule_type'] = 'DNA'
    SeqIO.write(record, path, 'genbank')
    return path


def acquisition_fixture(tmp_path):
    selection = []
    for index in range(50):
        base = f'SYN{index:06d}'
        row = {'phage': f'SyntheticTestOnly{index}', 'accession': base,
               'pilot': 'true' if index < 20 else 'false',
               'frozen_pilot_accession': base + '.1', 'frozen_pilot_sha256': ''}
        if index < 20:
            path = write_source(tmp_path / 'data/raw/genbank' / (base + '.gb'), base + '.1')
            row['frozen_pilot_sha256'] = sha256(path)
        selection.append(row)
    return selection


def source_request(module, runtime, selected):
    accession = selected['accession']
    return module.download('https://example.invalid/efetch?id=' + accession,
                           runtime / 'data/raw/genbank' / (accession + '.gb'), {}, maximum=3_000_000)


def test_source_is_frozen_immediately_before_later_acquisition_failure(tmp_path):
    from actinolysis.common import read_tsv
    selection = acquisition_fixture(tmp_path)
    def download(url, path, cfg, offline, maximum):
        if 'SYN000021' in url:
            raise urllib.error.HTTPError(url, 404, 'Synthetic test-only unavailable source', None, None)
        return write_source(path, path.stem + '.1')
    module = SimpleNamespace(download=download, SeqIO=SeqIO)
    with screen.AcquisitionLedger(module, tmp_path, selection):
        first = source_request(module, tmp_path, selection[20])
        assert first.name == 'SYN000020.1.gb'
        with pytest.raises(RuntimeError, match='private acquisition ledger'):
            source_request(module, tmp_path, selection[21])
    ledger = read_tsv(tmp_path / 'data/manifests/acquisition_ledger.tsv')
    assert len(ledger) == 50
    assert ledger[20]['status'] == 'acquired'
    assert ledger[20]['resolved_accession'] == 'SYN000020.1'
    assert ledger[20]['source_sha256'] == sha256(first)
    assert ledger[21]['status'] == 'failed'
    assert ledger[21]['failure_reason'] == 'http_404'
    assert all(r['status'] == 'not_attempted' for r in ledger[22:])


def test_selected_accession_version_must_resolve_exactly(tmp_path):
    from actinolysis.common import read_tsv
    selection = acquisition_fixture(tmp_path)
    selection[20]['accession'] += '.1'
    module = SimpleNamespace(SeqIO=SeqIO,
        download=lambda url, path, cfg, offline, maximum: write_source(path, 'SYN000020.2'))
    with screen.AcquisitionLedger(module, tmp_path, selection):
        with pytest.raises(RuntimeError, match='private acquisition ledger'):
            source_request(module, tmp_path, selection[20])
    row = read_tsv(tmp_path / 'data/manifests/acquisition_ledger.tsv')[20]
    assert row['status'] == 'failed'
    assert row['failure_reason'] == 'resolved_accession_version_mismatch'
    assert not row['resolved_accession']


def test_resume_rejects_changed_or_unregistered_new_source_cache(tmp_path):
    selection = acquisition_fixture(tmp_path)
    module = SimpleNamespace(SeqIO=SeqIO,
        download=lambda url, path, cfg, offline, maximum: write_source(path, path.stem + '.1'))
    with screen.AcquisitionLedger(module, tmp_path, selection):
        frozen = source_request(module, tmp_path, selection[20])
    frozen.write_bytes(frozen.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='source cache checksum mismatch'):
        screen.AcquisitionLedger(module, tmp_path, selection)


def test_unregistered_source_cannot_be_accepted_as_a_cache(tmp_path):
    selection = acquisition_fixture(tmp_path)
    write_source(tmp_path / 'data/raw/genbank/SYN000020.gb', 'SYN000020.1')
    module = SimpleNamespace(SeqIO=SeqIO, download=lambda *args, **kwargs: pytest.fail('No download expected'))
    with pytest.raises(ValueError, match='Unregistered acquisition source cache'):
        screen.AcquisitionLedger(module, tmp_path, selection)


def test_resume_reuses_registered_sources_and_finishes_remaining_records(tmp_path):
    from actinolysis.common import read_tsv
    selection = acquisition_fixture(tmp_path)
    def first_download(url, path, cfg, offline, maximum):
        if 'SYN000021' in url:
            raise urllib.error.HTTPError(url, 503, 'Synthetic test-only temporary failure', None, None)
        return write_source(path, path.stem + '.1')
    module = SimpleNamespace(SeqIO=SeqIO, download=first_download)
    with screen.AcquisitionLedger(module, tmp_path, selection):
        frozen = source_request(module, tmp_path, selection[20])
        before = sha256(frozen)
        with pytest.raises(RuntimeError):
            source_request(module, tmp_path, selection[21])
    downloads = []
    def resumed_download(url, path, cfg, offline, maximum):
        downloads.append(path.stem)
        return write_source(path, path.stem + '.1')
    module.download = resumed_download
    with screen.AcquisitionLedger(module, tmp_path, selection):
        for selected in selection:
            source_request(module, tmp_path, selected)
    rows = read_tsv(tmp_path / 'data/manifests/acquisition_ledger.tsv')
    assert all(r['status'] == 'acquired' for r in rows)
    assert len(downloads) == 29
    assert 'SYN000020' not in downloads
    assert sha256(frozen) == before
