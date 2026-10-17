"""Runner integrity tests using only explicit synthetic stage outputs."""
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from actinolysis.common import write_json, write_tsv


def runner_fixture(tmp_path, monkeypatch):
    # resource is Linux-only. This stub supplies accounting fields for the
    # isolated runner tests on native Windows; no real stage is executed.
    try:
        import resource
    except ImportError:
        monkeypatch.setitem(sys.modules, "resource", SimpleNamespace(
            RUSAGE_SELF=0, RUSAGE_CHILDREN=1,
            getrusage=lambda _: SimpleNamespace(ru_maxrss=0)))
    repo = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("synthetic_runner_test", repo/"scripts/run_all.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    monkeypatch.setattr(runner, "ROOT", tmp_path)
    monkeypatch.setattr(runner, "budget_check", lambda *a, **kw: {})
    monkeypatch.setattr(runner, "record_environment", lambda: {"synthetic_test_only": True})
    monkeypatch.setattr(runner, "project_bytes", lambda: 0)
    monkeypatch.chdir(tmp_path)
    cfg = yaml.safe_load((repo/"config/analysis.yaml").read_text())
    (tmp_path/"config").mkdir()
    (tmp_path/"config/analysis.yaml").write_text(yaml.safe_dump(cfg))
    for relative in ("scripts/run_all.py", "src/actinolysis/common.py",
                     "src/actinolysis/acquisition.py", "src/actinolysis/quality.py"):
        path = tmp_path/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Synthetic test-only source fingerprint\n")
    (tmp_path/"references").mkdir()
    (tmp_path/"references/protein_controls.tsv").write_text("reference_id\nsynthetic_test_only\n")
    (tmp_path/"references/architecture_controls.tsv").write_text("phage\nsynthetic_test_only\n")

    def fake_acquisition(cfg):
        raw = tmp_path/"data/raw/synthetic_test_only.txt"
        raw.parent.mkdir(parents=True, exist_ok=True)
        raw.write_text("synthetic input, never biological data")
        write_tsv(tmp_path/"data/manifests/pilot_inputs.tsv", [{"path": str(raw.relative_to(tmp_path))}])
        write_json(tmp_path/"data/manifests/phagesdb_snapshot.json", {"pages": []})
        write_json(tmp_path/"data/raw/phagesdb/eligible_metadata.json", {"synthetic_test_only": True})
        write_tsv(tmp_path/"data/manifests/eligibility_snapshot.tsv", [{"phage":"synthetic_test_only"}])

    def fake_quality(cfg):
        for relative in runner.OUTPUTS[2]:
            path = tmp_path/relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("synthetic stage output, never biological data\n")

    modules = {"actinolysis.acquisition": SimpleNamespace(run=fake_acquisition),
               "actinolysis.quality": SimpleNamespace(run=fake_quality)}
    monkeypatch.setattr(runner.importlib, "import_module", lambda name: modules[name])
    return runner, cfg


def test_resume_rejects_changed_scientific_config(tmp_path, monkeypatch):
    runner, cfg = runner_fixture(tmp_path, monkeypatch)
    runner.main(["--stop-after", "2"])
    cfg["identification"]["min_identity"] = 0.5
    (tmp_path/"config/analysis.yaml").write_text(yaml.safe_dump(cfg))
    with pytest.raises(RuntimeError, match="config or source changed"):
        runner.main(["--resume-from", "3", "--stop-after", "3"])


def test_resume_rejects_changed_source_and_changed_stage_output(tmp_path, monkeypatch):
    runner, _ = runner_fixture(tmp_path, monkeypatch)
    runner.main(["--stop-after", "2"])
    source = tmp_path/"src/actinolysis/quality.py"
    original = source.read_text()
    source.write_text(original+"# altered synthetic source\n")
    with pytest.raises(RuntimeError, match="config or source changed"):
        runner.main(["--resume-from", "3", "--stop-after", "3"])
    source.write_text(original)
    (tmp_path/"data/processed/proteins.faa").write_text("altered synthetic output")
    with pytest.raises(ValueError, match="Resume input changed"):
        runner.main(["--resume-from", "3", "--stop-after", "3"])


def test_rerunning_earlier_stage_invalidates_later_receipts(tmp_path, monkeypatch):
    runner, _ = runner_fixture(tmp_path, monkeypatch)
    runner.main(["--stop-after", "2"])
    receipts = tmp_path/"results/reports/pilot_stage_receipts.json"
    previous = json.loads(receipts.read_text())
    previous["9"] = {"outputs": {"synthetic_later_output": "obsolete"}}
    write_json(receipts, previous)
    runner.main(["--resume-from", "2", "--stop-after", "2"])
    refreshed = json.loads(receipts.read_text())
    assert set(refreshed) == {"1", "2"}
