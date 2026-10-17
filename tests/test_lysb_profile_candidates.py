"""Leakage, coverage and optional real-API tests for the candidate-only script."""
import importlib.util
import json
from pathlib import Path
import random

import pytest

spec = importlib.util.spec_from_file_location("lysb_profile_candidate_test_module", Path(__file__).resolve().parents[1] / "scripts/lysb_profile_candidates.py")
profile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(profile)


def test_source_holdout_removes_source_siblings_and_cross_source_exact_aliases():
    refs = {"a": {"resolved_source_accession": "AB123456.1"}, "b": {"resolved_source_accession": "AB123456.2"},
            "alias": {"resolved_source_accession": "CD123456.1"}, "unrelated": {"resolved_source_accession": "EF123456.1"}}
    seqs = {"a": "MALAMALA", "b": "MGGSMGGS", "alias": "MALAMALA", "unrelated": "TWYRTWYR"}
    train, excluded = profile.holdout_training(set(refs), ["a"], refs, seqs, True)
    assert train == ["unrelated"] and excluded == ["a", "alias", "b"]
    train, excluded = profile.holdout_training(set(refs), ["a"], refs, seqs, False)
    assert train == ["b", "unrelated"] and excluded == ["a", "alias"]
    assert not {profile.digest(seqs[i]) for i in train} & {profile.digest(seqs["a"])}


def test_duplicate_training_aliases_do_not_increase_alignment_weights():
    seqs = {"b": "MALA", "a": "MALA", "c": "MGGS"}
    assert sorted(profile.unique_training(seqs, seqs).values()) == [["a", "b"], ["c"]]


def test_global_coverage_counts_residue_pairs_instead_of_all_consumed_residues():
    metrics = profile.global_metrics("MALAMALAMALA", "MALAMALA")
    assert metrics["coverage_a"] == pytest.approx(8 / 12)
    assert metrics["coverage_b"] == 1
    assert metrics["identity"] == pytest.approx(8 / 12)


def test_reference_groups_are_transitive_proxy_components(monkeypatch):
    def metrics(a, b):
        return dict(identity=.5 if {a, b} != {"A", "C"} else .2, coverage_a=1., coverage_b=1., paired_residues=1,
                    alignment_columns=1, score=1.)
    monkeypatch.setattr(profile, "global_metrics", metrics)
    groups, edges = profile.reference_groups({"a", "b", "c"}, {"a": "A", "b": "B", "c": "C"})
    assert len(set(groups.values())) == 1 and len(edges) == 3
    assert sum(r["group_edge"] == "true" for r in edges) == 2


def test_domain_span_cannot_inflate_actual_coverage_across_deletions_or_insertions():
    # HMM has 10 match positions and 4 insert columns; target has 5 deletions.
    cov = profile.coverage("AAAAAAAAAA....", "AAAAA-----GGGG", 10, 9, 1, 10, 1, 9)
    assert cov["model_span_coverage"] == cov["target_span_coverage"] == 1
    assert cov["residue_pairs"] == 5
    assert cov["model_coverage"] == .5 and cov["target_coverage"] == pytest.approx(5 / 9)
    assert not profile.candidate(1e-20, 1e-20, cov["model_coverage"], cov["target_coverage"])
    with pytest.raises(ValueError):
        profile.coverage("AAA", "AA", 3, 3, 1, 3, 1, 3)


def test_candidate_cuts_are_inclusive_and_do_not_substitute_sequence_identity():
    assert profile.candidate(1e-6, 1e-6, .7, .7)
    assert not profile.candidate(1.01e-6, 1e-20, 1., 1.)
    assert not profile.candidate(1e-20, 1.01e-6, 1., 1.)
    assert not profile.candidate(1e-20, 1e-20, .69999, 1.)
    assert profile.SEARCH["Z"] == profile.SEARCH["domZ"] == 2460
    assert profile.SEARCH["bias_filter"] is profile.SEARCH["null2"] is True
    assert profile.BUILDER["seed"] == 20261016
    assert not profile.BROAD & profile.PREDICTED and profile.NARROW < profile.BROAD


def test_empty_fold_and_no_hit_are_explicit_and_alias_is_not_independent():
    seqs = {"q": "MALAMALA"}; meta = {"q": {"query_kind": "assay_B_reference"}}
    row = profile.summarize_queries("m", "held_out_reference", seqs, meta, [], {}, [], [], True)[0]
    assert row["status"] == "untestable_no_training_sequences" and not row["passes_candidate_cuts"]
    row = profile.summarize_queries("m", "cohort", seqs, meta, ["training"], {"training": "MALAMALA"}, [], [
        {"query_id": "q", "passes_candidate_cuts": True, "domain_reported": True}])[0]
    assert row["exact_training_alias"] and row["status"] == "exact_training_sequence_not_independent"


def test_frozen_stage_receipts_authenticate_first_cds_and_reference_inputs(tmp_path):
    paths = {}
    for name in ("cds", "references", "reference_manifest"):
        path = tmp_path / name
        path.write_text(f"SYNTHETIC TEST INPUT {name}\n")
        paths[name] = path
    receipt = {"2": {"outputs": {"data/processed/cds.tsv": profile.sha256(paths["cds"])}}, "3": {"outputs": {
        "data/processed/references.faa": profile.sha256(paths["references"]),
        "data/manifests/reference_manifest.tsv": profile.sha256(paths["reference_manifest"])}}}
    paths["stage_receipts"] = tmp_path / "receipts.json"
    paths["stage_receipts"].write_text(json.dumps(receipt))
    profile.authenticate_receipts(paths)
    paths["cds"].write_text("TAMPERED TEST INPUT\n")
    with pytest.raises(ValueError, match="stage 2"):
        profile.authenticate_receipts(paths)


def test_source_hashes_are_lf_portable_but_scientific_inputs_are_byte_exact(tmp_path):
    script, protocol, data = (tmp_path / n for n in ("code.py", "protocol.md", "data.tsv"))
    for path in (script, protocol, data):
        path.write_bytes(b"first\r\nsecond\r\n")
    paths = {"script": script, "protocol": protocol, "cds": data}
    original = profile.input_hashes(paths)
    for path in paths.values():
        path.write_bytes(b"first\nsecond\n")
    normalized = profile.input_hashes(paths)
    assert original["script"] == normalized["script"] and original["protocol"] == normalized["protocol"]
    assert original["cds"] != normalized["cds"]


def test_real_api_training_only_msa_fixed_search_space_and_evalue_scaling(tmp_path):
    pyfamsa = pytest.importorskip("pyfamsa", reason="Optional real HMM API test requires requirements-profile.txt in WSL")
    pyhmmer = pytest.importorskip("pyhmmer", reason="Optional real HMM API test requires requirements-profile.txt in WSL")
    rng = random.Random(301)
    sequence = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(140))
    sequences = {"SYNTH_TRAIN_A": sequence, "SYNTH_TRAIN_B": sequence[:50] + "A" + sequence[51:],
                 "SYNTH_HELD_OUT": sequence[:70] + "G" + sequence[71:]}
    hmm = profile.build_model(["SYNTH_TRAIN_A", "SYNTH_TRAIN_B"], sequences, tmp_path, "SYNTH_MODEL", pyfamsa, pyhmmer)
    training = json.loads((tmp_path / "training.json").read_text())
    assert training["train_ids"] == ["SYNTH_TRAIN_A", "SYNTH_TRAIN_B"]
    assert "SYNTH_HELD_OUT" not in (tmp_path / "alignment.sto").read_text()
    assert training["actual_builder_seed"] == 20261016
    assert "DATE" not in (tmp_path / "model.hmm").read_text()
    profile.build_model(["SYNTH_TRAIN_A", "SYNTH_TRAIN_B"], sequences, tmp_path / "repeat", "SYNTH_MODEL", pyfamsa, pyhmmer)
    assert profile.sha256(tmp_path / "model.hmm") == profile.sha256(tmp_path / "repeat/model.hmm")
    target = {"SYNTH_HELD_OUT": sequences["SYNTH_HELD_OUT"]}
    full_hits, full_domains = profile.search_model(hmm, target, "SYNTH_MODEL", "unit_test", pyhmmer)
    small = {**profile.SEARCH, "Z": 1, "domZ": 1}
    small_hits, small_domains = profile.search_model(hmm, target, "SYNTH_MODEL", "unit_test", pyhmmer, small)
    assert full_hits and small_hits and full_domains and small_domains
    assert full_hits[0]["Z"] == full_hits[0]["domZ"] == 2460
    assert full_hits[0]["full_bits"] == pytest.approx(small_hits[0]["full_bits"])
    assert full_hits[0]["full_evalue"] > 0
    assert full_hits[0]["full_evalue"] / small_hits[0]["full_evalue"] == pytest.approx(2460, rel=1e-5)
    assert full_domains[0]["independent_domain_evalue"] / small_domains[0]["independent_domain_evalue"] == pytest.approx(2460, rel=1e-5)
    assert full_domains[0]["model_coverage"] <= full_domains[0]["model_span_coverage"]


def test_real_api_single_remaining_reference_fold_is_buildable(tmp_path):
    pyfamsa = pytest.importorskip("pyfamsa", reason="Optional real HMM API test requires requirements-profile.txt in WSL")
    pyhmmer = pytest.importorskip("pyhmmer", reason="Optional real HMM API test requires requirements-profile.txt in WSL")
    seq = "MALAGGEITKELQKALVEAGADVVVVSGSGNSGALAALAAAAGADVVVVATSGNGTADKAALEAGADVVVVSGHGSTALAALAAAGADVVVVATSGNGTA"
    hmm = profile.build_model(["SYNTH_SINGLE"], {"SYNTH_SINGLE": seq}, tmp_path, "SYNTH_SINGLE_FOLD", pyfamsa, pyhmmer)
    assert hmm.nseq == 1 and hmm.M > 0
    saved = json.loads((tmp_path / "training.json").read_text())
    assert saved["actual_builder_seed"] == 20261016 and saved["model_sequence_count"] == 1
