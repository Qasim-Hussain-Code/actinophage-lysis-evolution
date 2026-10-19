"""Synthetic-only public replay integrity tests; no network access or biological claims."""
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest
import yaml

from actinolysis.common import sha256, write_json, write_tsv
from scripts import reproduce_exploratory_screen as replay


def cohort():
    selected = [{'rank': str(i + 1), 'accession': f'SYN{i:06d}.1', 'sha256': 'a' * 64,
                 'phage': f'SyntheticOnly{i}', 'host_genus': 'Mycobacterium', 'cluster': 'SYNTHETIC_ONLY',
                 'subcluster': 'not_available', 'pilot': 'true' if i < 20 else 'false',
                 'source_url': 'https://example.invalid/synthetic-only'} for i in range(50)]
    return selected, [dict(row) for row in selected[:20]]


@pytest.fixture
def public_root(tmp_path):
    selected, pilot = cohort()
    write_tsv(tmp_path / 'data/manifests/pilot_inputs.tsv', pilot)
    cfg = {'mode': 'pilot', 'seed': 20261016, 'sample_sizes': {'pilot': 20, 'intermediate': 50},
           'identification': {'min_identity': .4}}
    config = tmp_path / 'config/analysis.yaml'
    config.parent.mkdir()
    config.write_text(yaml.safe_dump(cfg))
    code_names = ['src/synthetic_only.txt', 'scripts/prepare_exploratory_screen.py',
                  'scripts/run_exploratory_screen.py', 'scripts/replay_exploratory_screen.py',
                  'scripts/archive_cohort_hits.py', 'docs/exploratory_screen_protocol.md']
    for name in code_names:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('SYNTHETIC TEST ONLY\n')
    references = [{'reference_id': f'SYNTHETIC_REF_{i}'} for i in range(39)]
    panel = tmp_path / replay.PANEL
    write_tsv(panel / 'references/protein_controls.tsv', references)
    write_tsv(panel / 'data/manifests/reference_manifest.tsv', references)
    panel_hashes = {name: sha256(panel / name) for name in
                    ('references/protein_controls.tsv', 'data/manifests/reference_manifest.tsv')}
    write_json(panel / 'artifact_manifest.json', {'artifacts': panel_hashes})
    public = tmp_path / replay.PUBLIC
    write_tsv(public / 'selected_genomes.tsv', selected)
    cfg['mode'] = 'intermediate'
    report = {'status': 'pass', 'completed_stages': [1, 2, 3, 4, 5],
              'public_cohort_manifest': {'artifact': 'selected_genomes.tsv', 'sha256': sha256(public / 'selected_genomes.tsv')},
              'protocol_sha256': sha256(tmp_path / code_names[-1]), 'reference_panel_sha256': panel_hashes,
              'frozen_provenance': {'configuration_file_sha256': sha256(config),
                  'effective_configuration_sha256': replay.digest_object(cfg),
                  'execution_source_tree_sha256': replay.runner.tree_hashes(tmp_path, ('src',)),
                  'preparation_code_sha256': {code_names[1]: sha256(tmp_path / code_names[1])},
                  'execution_runner_sha256': sha256(tmp_path / code_names[2])}}
    write_json(public / 'replay_validation.json', report)
    write_json(public / 'artifact_manifest.json', {'analysis': 'exploratory50_v1', 'scope': 'exploratory_only',
        'replay_identical': True, 'export_driver_sha256': sha256(tmp_path / code_names[3]),
        'hit_ledger_builder_sha256': sha256(tmp_path / code_names[4]),
        'artifacts': {name: sha256(public / name) for name in ('selected_genomes.tsv', 'replay_validation.json')}})
    return tmp_path


def test_public_inputs_need_no_private_selection_frame(public_root):
    selected, references, cfg, _ = replay.verify_public(public_root)
    assert len(selected) == 50 and len(references) == 39 and cfg['mode'] == 'intermediate'
    assert not (public_root / 'data/manifests/eligibility_snapshot.tsv').exists()
    assert not (public_root / '.experiments').exists()


@pytest.mark.parametrize('mutation', ['rank', 'version', 'duplicate', 'pilot_hash'])
def test_cohort_fails_closed(mutation):
    selected, pilot = cohort()
    if mutation == 'rank':
        selected[21]['rank'] = '1'
    elif mutation == 'version':
        selected[21]['accession'] = selected[21]['accession'].split('.')[0]
    elif mutation == 'duplicate':
        selected[21]['accession'] = selected[22]['accession']
    else:
        selected[0]['sha256'] = 'b' * 64
    with pytest.raises(ValueError):
        replay.verify_cohort(selected, pilot)


def test_configuration_drift_is_rejected(public_root):
    path = public_root / 'config/analysis.yaml'
    path.write_bytes(path.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='configuration file'):
        replay.verify_public(public_root)


def test_public_artifact_paths_cannot_escape(tmp_path):
    write_json(tmp_path / 'artifact_manifest.json', {'artifacts': {'../outside.tsv': 'a' * 64}})
    with pytest.raises(ValueError, match='path is invalid'):
        replay.verify_artifacts(tmp_path)


def write_record(path, accession):
    record = SeqRecord(Seq('ATG' * 10), id=accession, description='Synthetic test only')
    record.annotations['molecule_type'] = 'DNA'
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(record, path, 'genbank')


def source_spec(path):
    return {'path': 'data/raw/genbank/SYN000001.1.gb', 'accession': 'SYN000001.1',
            'sha256': sha256(path), 'effective_retrieval_url': replay.pinned_url('SYN000001.1')}


def test_cache_only_never_downloads_and_checks_checksum(tmp_path, monkeypatch):
    cache = tmp_path / 'cache'
    source = cache / 'SYN000001.1.gb'
    write_record(source, 'SYN000001.1')
    spec = source_spec(source)
    monkeypatch.setattr(replay, 'download', lambda *a, **k: pytest.fail('Cache-only replay tried to download'))
    target = replay.acquire_source(spec, tmp_path / 'runtime', {}, cache)
    assert sha256(target) == spec['sha256']
    source.write_bytes(source.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='checksum mismatch'):
        replay.acquire_source(spec, tmp_path / 'other_runtime', {}, cache)


def test_cache_exact_version_is_verified_even_when_checksum_matches(tmp_path):
    source = tmp_path / 'cache/SYN000001.1.gb'
    write_record(source, 'SYN000001.2')
    with pytest.raises(ValueError, match='accession-version mismatch'):
        replay.acquire_source(source_spec(source), tmp_path / 'runtime', {}, source.parent)


def test_fresh_retrieval_uses_frozen_url_and_bounded_download(tmp_path, monkeypatch):
    source = tmp_path / 'fixture.gb'
    write_record(source, 'SYN000001.1')
    calls = []
    def download(url, target, cfg, maximum):
        calls.append((url, maximum))
        target.parent.mkdir(parents=True)
        target.write_bytes(source.read_bytes())
    monkeypatch.setattr(replay, 'download', download)
    spec = source_spec(source)
    replay.acquire_source(spec, tmp_path / 'runtime', {})
    assert calls == [(replay.pinned_url('SYN000001.1'), 3_000_000)]


def test_direct_manifest_never_resamples_and_retains_source_provenance(monkeypatch):
    import actinolysis.acquisition
    monkeypatch.setattr(actinolysis.acquisition, 'select_sample', lambda *a: pytest.fail('Unexpected resampling'))
    selected, _ = cohort()
    manifest = replay.direct_manifest(selected)
    assert [r['accession'] for r in manifest] == [r['accession'] for r in selected]
    assert all(r['source_url'] == selected[0]['source_url'] for r in manifest)
    assert all(r['host'] == 'not_available' for r in manifest)


def test_source_plan_pins_retrieval_versions_separately_from_original_urls():
    selected, _ = cohort()
    refs = [{'source_accession': 'SYN123456', 'resolved_source_accession': 'SYN123456.2',
             'source_sha256': 'b' * 64, 'source_url': 'https://example.invalid/original-resolution'}]
    specs = replay.source_specs(selected, refs)
    assert specs[-1]['path'] == 'data/raw/references/SYN123456.gb'
    assert specs[-1]['original_source_url'] == refs[0]['source_url']
    assert specs[-1]['effective_retrieval_url'] == replay.pinned_url('SYN123456.2')
    assert specs[0]['effective_retrieval_url'] == replay.pinned_url(selected[0]['accession'])


def test_comparison_uses_scientific_tables_and_excludes_replay_provenance(tmp_path):
    public, output = tmp_path / 'public', tmp_path / 'output'
    public.mkdir(); output.mkdir()
    for name in ('selected_genomes.tsv', 'summary.json', 'quality_control.json'):
        (public / name).write_text('SYNTHETIC TEST ONLY\n')
        (output / name).write_bytes((public / name).read_bytes())
    manifest = {'artifacts': {name: sha256(public / name) for name in
                             ('selected_genomes.tsv', 'summary.json', 'quality_control.json')}}
    manifest['artifacts']['replay_validation.json'] = 'a' * 64
    manifest['artifacts']['figure_manifest.json'] = 'b' * 64
    assert len(replay.compare_public_science(public, output, manifest)) == 3
    (output / 'summary.json').write_text('SYNTHETIC ALTERED RESULT\n')
    with pytest.raises(ValueError, match='artifacts differ'):
        replay.compare_public_science(public, output, manifest)
