"""Checks for accession-frozen replay and checksum-verified source reuse."""
import importlib.util
from pathlib import Path

import pytest

from actinolysis.common import load_config, read_tsv, sha256

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('cohort_replay', ROOT / 'scripts/replay_cohort.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def test_frozen_cohort_order_is_preserved():
    manifest = read_tsv(ROOT / 'data/manifests/pilot_inputs.tsv')
    frame = replay.cohort_frame(manifest, load_config(ROOT / 'config/analysis.yaml'))
    assert [r['accession'] for r in frame] == [r['accession'] for r in manifest]


def test_changed_cohort_order_is_rejected():
    manifest = read_tsv(ROOT / 'data/manifests/pilot_inputs.tsv')
    manifest[2], manifest[3] = manifest[3], manifest[2]
    with pytest.raises(ValueError, match='cohort order'):
        replay.cohort_frame(manifest, load_config(ROOT / 'config/analysis.yaml'))


def test_cache_checksum_mismatch_is_rejected(tmp_path):
    (tmp_path / 'source').write_bytes(b'test-only data')
    with pytest.raises(ValueError, match='checksum mismatch'):
        replay.copy_cached(tmp_path, 'source', '0' * 64, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_verified_cache_record_is_copied(tmp_path):
    (tmp_path / 'source').write_bytes(b'test-only data')
    assert replay.copy_cached(tmp_path, 'source', sha256(tmp_path / 'source'), tmp_path / 'out')
    assert (tmp_path / 'out').read_bytes() == b'test-only data'
