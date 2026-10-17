"""Replay authentication uses synthetic files and never runs a search or clone."""
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from actinolysis.common import sha256, write_json, write_tsv

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/verify_amendment_replay.py'
SPEC = importlib.util.spec_from_file_location('verify_amendment_replay', SCRIPT)
replay = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(replay)


@pytest.fixture
def panel(tmp_path):
    clone = tmp_path / 'synthetic_clone'; cache = tmp_path / 'synthetic_cache'
    runtime = clone / '.experiments/amendment_experimental_v1'
    export = clone / 'results/amendment_v1/experimental'
    tool = cache / 'synthetic_mmseqs_never_executed'; tool.parent.mkdir(parents=True)
    tool.write_bytes(b'SYNTHETIC TEST ONLY; NEVER EXECUTED')
    cfg = dict(mode='pilot', seed=1, sample_sizes={}, sampling={}, identification={}, architecture={},
               neighborhood={}, statistics={}, tools={'mmseqs': str(tool)})
    (runtime / 'config').mkdir(parents=True)
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg))
    environment = {'python': 'SYNTHETIC_TEST_ONLY', 'packages': {'synthetic': '1'}}
    manifests = {}
    for relative, accession in [('data/raw/genbank/SYNGEN.1.gb', 'SYNGEN.1'),
                                ('data/raw/references/SYNREF.1.gb', 'SYNREF.1')]:
        path = runtime / relative; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f'LOCUS synthetic\nVERSION {accession}\n//\n'.encode())
        cached = cache / relative; cached.parent.mkdir(parents=True, exist_ok=True); cached.write_bytes(path.read_bytes())
        manifests[relative] = sha256(path)
    write_tsv(runtime / 'data/manifests/pilot_inputs.tsv', [dict(path='data/raw/genbank/SYNGEN.1.gb', accession='SYNGEN.1', sha256=manifests['data/raw/genbank/SYNGEN.1.gb'])])
    write_tsv(runtime / 'data/manifests/reference_manifest.tsv', [dict(reference_id='SYNTHETIC_TEST_ONLY', source_accession='SYNREF.1',
              resolved_source_accession='SYNREF.1', resolved_protein_accession='SYNPROTEIN.1',
              source_sha256=manifests['data/raw/references/SYNREF.1.gb'], length_aa=3, protein_sha256=replay.digest(b'MAA'))])
    (runtime / 'data/processed').mkdir(parents=True)
    (runtime / 'data/processed/references.faa').write_text('>SYNTHETIC_TEST_ONLY\nMAA\n')
    source = runtime / 'src/actinolysis/common.py'; source.parent.mkdir(parents=True)
    source.write_text('# SYNTHETIC TEST ONLY\n')
    science = {}
    for index in range(39):
        relative = f'results/tables/synthetic_test_only_{index}.tsv'
        path = runtime / relative; path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f'SYNTHETIC_TEST_ONLY\n{index}\n'); science[relative] = path.read_bytes()
    fingerprint = replay.digest(json.dumps({key: cfg[key] for key in replay.CONFIG_FIELDS}, sort_keys=True).encode())
    receipts = {str(stage): dict(outputs={f'results/tables/synthetic_test_only_{stage}.tsv': replay.digest(science[f'results/tables/synthetic_test_only_{stage}.tsv'])},
                 source_sha256={'src/actinolysis/common.py': sha256(source)}, scientific_config_sha256=fingerprint, seconds=1)
                for stage in range(1, 10)}
    execution = dict(status='success', stages=[dict(stage=stage, **receipts[str(stage)]) for stage in range(1, 10)],
                     config={**cfg, 'offline': True}, environment=environment, elapsed_seconds=9, peak_child_rss_bytes=1)
    write_json(runtime / 'results/reports/pilot_stage_receipts.json', receipts)
    write_json(runtime / 'results/reports/pilot_execution.json', execution)
    write_json(runtime / 'results/reports/validation.json', {'checksum_status': 'pass'})
    write_json(runtime / 'results/reports/neighborhoods.json', {'family_clustering': [{'cache_reused': False}]})
    provenance = dict(tool_binary_sha256=sha256(tool), snapshot_sha256={'src/actinolysis/common.py': sha256(source)})
    write_json(runtime / 'amendment_provenance.json', provenance)
    expected = {**science, 'analysis.yaml': yaml.safe_dump({**cfg, 'tools': {'mmseqs': 'runtime/path/may/differ'}}).encode(),
                'amendment_provenance.json': (runtime / 'amendment_provenance.json').read_bytes(),
                'results/reports/pilot_execution.json': (runtime / 'results/reports/pilot_execution.json').read_bytes(),
                **{relative: (runtime / relative).read_bytes() for relative in ['data/manifests/pilot_inputs.tsv', 'data/manifests/reference_manifest.tsv']}}
    manifest = {relative: replay.digest(data) for relative, data in expected.items()}
    expected['artifact_manifest.json'] = json.dumps({'artifacts': manifest}).encode()
    for relative, data in expected.items():
        path = export / relative; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
    published = {'results/amendment_v1/experimental/' + relative: data for relative, data in expected.items()}
    return dict(clone=clone, cache=cache, runtime=runtime, export=export, published=published)


def test_every_scientific_artifact_and_versioned_cached_source_authenticated(panel):
    result = replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])
    assert result['status'] == 'pass' and result['compared_artifacts'] == 39
    assert len(result['artifacts']) == 39 and len(result['source_checks']) == 2
    assert result['reference_proteins_checked'] == 1 and result['completed_stages'] == list(range(1, 10))


def test_cloned_committed_exports_are_insufficient_without_actual_runtime(panel):
    path = panel['runtime'] / 'results/reports/pilot_stage_receipts.json'
    path.unlink()
    with pytest.raises(FileNotFoundError):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


@pytest.mark.parametrize('relative', ['data/raw/genbank/SYNGEN.1.gb', 'data/raw/references/SYNREF.1.gb'])
def test_changed_cached_source_is_rejected_even_when_runtime_and_results_match(panel, relative):
    (panel['cache'] / relative).write_bytes(b'TAMPERED SYNTHETIC TEST ONLY\n')
    with pytest.raises(ValueError, match='version-pinned source mismatch'):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


def test_wrong_accession_version_is_rejected_even_if_checksums_are_reforged(panel):
    path = panel['runtime'] / 'data/raw/references/SYNREF.1.gb'
    path.write_bytes(b'LOCUS synthetic\nVERSION SYNREF.2\n//\n')
    (panel['cache'] / 'data/raw/references/SYNREF.1.gb').write_bytes(path.read_bytes())
    relative = 'data/manifests/reference_manifest.tsv'; prefix = 'results/amendment_v1/experimental/'
    rows = replay.tsv(panel['published'][prefix + relative]); rows[0]['source_sha256'] = sha256(path)
    write_tsv(panel['runtime'] / relative, rows)
    data = (panel['runtime'] / relative).read_bytes()
    (panel['export'] / relative).write_bytes(data); panel['published'][prefix + relative] = data
    manifest = json.loads(panel['published'][prefix + 'artifact_manifest.json'])
    manifest['artifacts'][relative] = replay.digest(data)
    manifest_data = json.dumps(manifest).encode()
    (panel['export'] / 'artifact_manifest.json').write_bytes(manifest_data)
    panel['published'][prefix + 'artifact_manifest.json'] = manifest_data
    with pytest.raises(ValueError, match='version-pinned source mismatch'):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


def test_runtime_source_change_is_not_hidden_by_updated_receipt_hashes(panel):
    source = panel['runtime'] / 'src/actinolysis/common.py'; source.write_text('# CHANGED SYNTHETIC TEST ONLY\n')
    receipts_path = panel['runtime'] / 'results/reports/pilot_stage_receipts.json'; receipts = json.loads(receipts_path.read_text())
    execution_path = panel['runtime'] / 'results/reports/pilot_execution.json'; execution = json.loads(execution_path.read_text())
    for stage in range(1, 10):
        receipts[str(stage)]['source_sha256']['src/actinolysis/common.py'] = sha256(source)
        execution['stages'][stage - 1] = {'stage': stage, **receipts[str(stage)]}
    write_json(receipts_path, receipts); write_json(execution_path, execution)
    with pytest.raises(ValueError, match='preparation snapshot changed'):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


def test_resumed_stage9_only_log_does_not_claim_full_replay(panel):
    path = panel['runtime'] / 'results/reports/pilot_execution.json'; execution = json.loads(path.read_text())
    execution['stages'] = execution['stages'][-1:]; write_json(path, execution)
    with pytest.raises(ValueError, match='successful offline stages 1-9'):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


def test_validly_receipted_science_change_is_reported_as_reproduction_failure(panel):
    relative = 'results/tables/synthetic_test_only_0.tsv'
    for base in [panel['runtime'], panel['export']]:
        (base / relative).write_text('CHANGED SYNTHETIC TEST ONLY\n')
    manifest_path = panel['export'] / 'artifact_manifest.json'; manifest = json.loads(manifest_path.read_text())
    manifest['artifacts'][relative] = sha256(panel['export'] / relative); write_json(manifest_path, manifest)
    result = replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])
    assert result['status'] == 'fail' and result['differences'] == [relative]


def test_reused_family_cache_is_not_misreported_as_uncached(panel):
    path = panel['runtime'] / 'results/reports/neighborhoods.json'
    write_json(path, {'family_clustering': [{'cache_reused': True}]})
    with pytest.raises(ValueError, match='without a pre-existing family cache'):
        replay.verify_panel('experimental', panel['published'], panel['clone'], panel['cache'])


def test_clone_must_be_the_declared_github_publication(monkeypatch, tmp_path):
    commit = '1' * 40
    def fake_git(root, *args):
        return (b'https://github.com/synthetic-test-only/repository.git\n' if args[0] == 'remote' else (commit + '\n').encode())
    monkeypatch.setattr(replay, 'git', fake_git)
    assert replay.clone_identity(tmp_path, commit)['origin_host'] == 'github.com'
    with pytest.raises(ValueError, match='declared GitHub-cloned publication'):
        replay.clone_identity(tmp_path, '2' * 40)


def test_original_pilot_is_also_checked_against_the_committed_scientific_bytes(tmp_path):
    clone = tmp_path / 'clone'; workspace = tmp_path / 'workspace'
    published = {f'results/tables/synthetic_test_only_{index}.tsv': f'SYNTHETIC TEST ONLY\n{index}\n'.encode() for index in range(39)}
    for base in [clone, workspace]:
        for relative, data in published.items():
            path = base / relative; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
    assert replay.original_science(published, clone, workspace)['status'] == 'pass'
    relative = next(iter(published)); (workspace / relative).write_text('CHANGED SYNTHETIC TEST ONLY\n')
    result = replay.original_science(published, clone, workspace)
    assert result['status'] == 'fail' and result['differences'] == [relative]
