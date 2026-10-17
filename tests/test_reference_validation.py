"""Query-only validation tests; every local sequence and score is synthetic."""
import hashlib
import json
from pathlib import Path

import pytest

from actinolysis.common import NA, read_tsv, write_tsv
from actinolysis import reference_validation as validation


def config():
    return {"identification": {"min_identity": .4, "min_query_coverage": .7,
                               "min_target_coverage": .7, "max_evalue": 1e-6}}


def reference(identifier, source, category="lysin_b", evidence="experimental", sequence="MAA"):
    return dict(reference_id=identifier, resolved_source_accession=source, resolved_protein_accession="SYNTHETIC_TEST_ONLY_" + identifier,
                category=category, evidence_type=evidence, length_aa=len(sequence),
                protein_sha256=hashlib.sha256(sequence.encode()).hexdigest())


def panel():
    rows = [reference("self", "NC_001900.1"), reference("same_source", "NC_001900.2", sequence="MCC"),
            reference("other", "AF123456.1", sequence="MGG"),
            reference("profile", "AB123456.1", evidence="published_profile_prediction", sequence="MLL"),
            reference("generic", "AC123456.1", evidence="experimental_generic_lipolytic", sequence="MKK"),
            reference("negative", "NC_001900.1", category="negative_control", sequence="MPP")]
    sequences = {row["reference_id"]: sequence for row, sequence in zip(rows, ["MAA", "MCC", "MGG", "MLL", "MKK", "MPP"])}
    return {row["reference_id"]: row for row in rows}, sequences


def challenge_spec():
    return dict(challenge_id="SYNTHETIC_TEST_ONLY", source_id="9SYN_1", pdb_id="9SYN", polymer_entity_id="1",
                author_chain_id="A", equivalent_author_chains="A;B", expected_length_aa="3",
                sequence_url="https://www.rcsb.org/fasta/entry/9SYN/display",
                analysis_use="query_only_exclude_from_functional_reference_panel")


def fake_download(url, path, cfg, offline=False, maximum=None):
    assert url == challenge_spec()["sequence_url"] and maximum == 64_000
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(">9SYN_1|Chains A, B|SYNTHETIC TEST ONLY\nMWW\n", encoding="ascii")
    return path


def prepared_queries(tmp_path):
    assets, challenges = validation.prepare_challenges([challenge_spec()], tmp_path, {}, fetch=fake_download)
    references, sequences = panel()
    queries = validation.write_queries(tmp_path, references, sequences, challenges, assets)
    return queries, references


def hit(query_id, reference_id, **changes):
    row = dict(protein_id=query_id, reference_id=reference_id, identity=".8", alignment_length="3",
               query_coverage=".9", target_coverage=".9", evalue="1e-20", bits="100")
    row.update({key: str(value) for key, value in changes.items()})
    return row


def test_pdb_entity_query_preparation_freezes_identity_chain_length_and_sequence(tmp_path):
    assets, sequences = validation.prepare_challenges([challenge_spec()], tmp_path, {}, fetch=fake_download)
    assert sequences == {"SYNTHETIC_TEST_ONLY": "MWW"}
    assert assets[0]["source_id"] == "9SYN_1" and assets[0]["length_aa"] == 3
    assert assets[0]["protein_sha256"] == hashlib.sha256(b"MWW").hexdigest()
    freeze = tmp_path / "frozen_inputs.json"
    validation.freeze_inputs(freeze, {"challenge_assets": assets})
    before = freeze.read_bytes()
    raw = tmp_path / assets[0]["raw_path"]
    raw.write_text(">9SYN_1|Chains A, B|SYNTHETIC TEST ONLY\nMWA\n")
    changed, _ = validation.prepare_challenges([challenge_spec()], tmp_path, {}, fetch=fake_download)
    with pytest.raises(ValueError, match="provenance changed"):
        validation.freeze_inputs(freeze, {"challenge_assets": changed})
    assert freeze.read_bytes() == before


@pytest.mark.parametrize("payload, message", [
    (">9SYN_2|Chains A, B|SYNTHETIC TEST ONLY\nMWW\n", "exactly one"),
    (">9SYN_1|Chain C|SYNTHETIC TEST ONLY\nMWW\n", "chain labels"),
    (">9SYN_1|Chains A, B|SYNTHETIC TEST ONLY\nMWWW\n", "length mismatch"),
    (">9SYN_1|Chains A, B|SYNTHETIC TEST ONLY\nM*W\n", "Invalid protein"),
])
def test_challenge_entity_mismatch_fails_before_search(tmp_path, payload, message):
    path = tmp_path / "synthetic_test_only.faa"
    path.write_text(payload)
    with pytest.raises(ValueError, match=message):
        validation.parse_challenge_entity(path, challenge_spec())


def test_queries_include_every_nonnegative_reference_and_keep_origin_tiers(tmp_path):
    queries, references = prepared_queries(tmp_path)
    positives = [row for row in queries if row["query_kind"] == "positive_reference"]
    assert len(positives) == len(references) - 1
    assert {row["origin_evidence_type"] for row in positives} == {"experimental", "published_profile_prediction", "experimental_generic_lipolytic"}
    assert "negative" not in {row["source_reference_id"] for row in positives}
    assert len(read_tsv(tmp_path / "queries.tsv")) == len(queries)


def test_source_holdout_excludes_all_same_accession_base_before_recovery(tmp_path):
    queries, references = prepared_queries(tmp_path)
    qid = "POS__self"
    hits = [hit(qid, "self", identity=1, query_coverage=1, target_coverage=1),
            hit(qid, "same_source", identity=1, query_coverage=1, target_coverage=1),
            hit(qid, "other")]
    alignments, summaries = validation.score_queries(queries, references, hits, config())
    by_ref = {row["reference_id"]: row for row in alignments}
    assert by_ref["self"]["score_confidence"] == "high" and not by_ref["self"]["accepted_call"]
    assert by_ref["self"]["exclusion_reason"] == "self_reference"
    assert by_ref["same_source"]["exclusion_reason"] == "same_resolved_source_accession_base"
    assert by_ref["other"]["fold_eligible"] and by_ref["other"]["accepted_expected_category"]
    summary = next(row for row in summaries if row["query_id"] == qid)
    assert summary["expected_category_recovered"] and summary["best_reference_id"] == "other"
    assert summary["n_excluded_hits"] == 2
    assert summary["n_excluded_panel_references"] == 3  # Self, same-source LysB, same-source structural control.
    assert validation.source_base("NC_001900.1") == validation.source_base("NC_001900.2")
    assert validation.source_base("NC_001900.1") != validation.source_base("AF123456.1")


def test_self_matches_cannot_make_unrecoverable_fold_pass_and_no_hit_queries_are_retained(tmp_path):
    queries, references = prepared_queries(tmp_path)
    alignments, summaries = validation.score_queries(queries, references,
        [hit("POS__self", "self", identity=1, query_coverage=1, target_coverage=1),
         hit("POS__self", "same_source", identity=1, query_coverage=1, target_coverage=1)], config())
    assert len(summaries) == len(queries)
    summary = next(row for row in summaries if row["query_id"] == "POS__self")
    assert not summary["expected_category_recovered"]
    assert summary["n_eligible_hits"] == 0 and summary["best_reference_id"] == NA
    assert summary["status"] == "expected_category_not_recovered"
    assert all(not row["accepted_call"] for row in alignments)
    profile = next(row for row in summaries if row["query_id"] == "POS__profile")
    assert profile["origin_evidence_type"] == "published_profile_prediction"
    assert profile["n_reported_hits"] == 0 and not profile["expected_category_recovered"]


def test_exact_protein_alias_from_another_source_cannot_validate_itself(tmp_path):
    queries, references = prepared_queries(tmp_path)
    references['other']['protein_sha256'] = references['self']['protein_sha256']
    alignments, summaries = validation.score_queries(queries, references,
        [hit('POS__self', 'other', identity=1, query_coverage=1, target_coverage=1)], config())
    assert alignments[0]['exclusion_reason'] == 'exact_sequence_alias'
    assert not alignments[0]['accepted_call']
    summary = next(row for row in summaries if row['query_id'] == 'POS__self')
    assert not summary['expected_category_recovered'] and summary['n_excluded_panel_references'] == 4


def test_challenge_call_is_flagged_with_capped_generic_evidence_and_no_mAGP_inactivity_claim(tmp_path):
    queries, references = prepared_queries(tmp_path)
    qid = "CHAL__SYNTHETIC_TEST_ONLY"
    alignments, summaries = validation.score_queries(queries, references,
        [hit(qid, "generic", identity=1, query_coverage=1, target_coverage=1)], config())
    assert alignments[0]["score_confidence"] == "moderate"
    assert alignments[0]["reference_evidence_type"] == "experimental_generic_lipolytic"
    assert alignments[0]["potential_cross_family_lysis_call"]
    summary = next(row for row in summaries if row["query_id"] == qid)
    assert summary["status"] == "potential_cross_family_lysis_label"
    assert "mAGP inactivity is not established" in summary["interpretation"]
    assert summary["minimum_identity"] == .4 and summary["minimum_target_coverage"] == .7


def test_reference_panel_length_and_sequence_hash_are_enforced(tmp_path):
    references, sequences = panel()
    write_tsv(tmp_path / "data/manifests/reference_manifest.tsv", references.values())
    fasta = tmp_path / "data/processed/references.faa"
    fasta.parent.mkdir(parents=True)
    fasta.write_text("".join(f">{key}\n{sequence}\n" for key, sequence in sequences.items()))
    assert validation.read_reference_panel(tmp_path)[1] == sequences
    fasta.write_text(fasta.read_text().replace("MAA", "MGA"))
    with pytest.raises(ValueError, match="checksum mismatch"):
        validation.read_reference_panel(tmp_path)


def test_search_arguments_match_declared_pipeline_limits_without_executing_tool():
    command = validation.search_command("mmseqs", "queries.faa", "references.faa", "hits.tsv", "temporary")
    assert command[1:6] == ["easy-search", "queries.faa", "references.faa", "hits.tsv", "temporary"]
    for flag, expected in [("--threads", "4"), ("-s", "7.5"), ("-k", "5"), ("--alignment-mode", "3"),
                           ("-e", "0.001"), ("--min-seq-id", "0.25"), ("-c", "0.5"),
                           ("--max-seqs", "100"), ("--split-memory-limit", "1G")]:
        assert command[command.index(flag) + 1] == expected


def test_raw_hit_parser_rejects_duplicate_and_unknown_identifiers(tmp_path):
    queries, references = prepared_queries(tmp_path)
    row = hit("POS__self", "other")
    line = "\t".join(row[key] for key in validation.RAW_FIELDS)
    path = tmp_path / "raw_hits.tsv"
    path.write_text(line + "\n")
    assert len(validation.parse_hits(path, queries, references)) == 1
    path.write_text(line + "\n" + line + "\n")
    with pytest.raises(ValueError, match="duplicate"):
        validation.parse_hits(path, queries, references)
    path.write_text(line.replace("POS__self", "POS__unknown") + "\n")
    with pytest.raises(ValueError, match="Unknown"):
        validation.parse_hits(path, queries, references)


def test_isolated_validation_io_cache_and_query_freeze_without_running_mmseqs(tmp_path, monkeypatch):
    """Exercise real files and reports with a synthetic subprocess replacement."""
    from types import SimpleNamespace
    import yaml
    project = tmp_path / "synthetic_test_only_project"
    isolated = project / ".experiments/synthetic_test_only_panel"
    output = project / "results/synthetic_test_only_validation"
    cfg = {"mode": "pilot", "tools": {"mmseqs": "synthetic_test_only_mmseqs"}, **config()}
    for root in (project, isolated):
        (root / "config").mkdir(parents=True)
        (root / "config/analysis.yaml").write_text(yaml.safe_dump(cfg))
    references, sequences = panel()
    write_tsv(isolated / "data/manifests/reference_manifest.tsv", references.values())
    fasta = isolated / "data/processed/references.faa"
    fasta.parent.mkdir(parents=True)
    fasta.write_text("".join(f">{key}\n{sequence}\n" for key, sequence in sequences.items()))
    write_tsv(project / "references/specificity_challenges.tsv", [challenge_spec()])
    originals = {path: path.read_bytes() for path in isolated.rglob("*") if path.is_file()}
    monkeypatch.setattr(validation, "download", fake_download)
    monkeypatch.setattr(validation, "budget_check", lambda cfg, additional_bytes=0, root=None: {
        "project_bytes": 1000, "free_bytes": 9999999999, "cap_bytes": 10000000000, "proposed_bytes": additional_bytes})
    searches = []
    def fake_subprocess(command, **kwargs):
        if command[1] == "version":
            return SimpleNamespace(stdout="synthetic-test-only-version", returncode=0)
        assert command[1] == "easy-search"
        assert kwargs["env"]["OMP_NUM_THREADS"] == "4"
        searches.append(command)
        rows = [hit("POS__self", "self", identity=1, query_coverage=1, target_coverage=1),
                hit("POS__self", "other"),
                hit("CHAL__SYNTHETIC_TEST_ONLY", "generic", identity=1, query_coverage=1, target_coverage=1)]
        Path(command[4]).write_text("".join("\t".join(row[field] for field in validation.RAW_FIELDS) + "\n" for row in rows))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(validation.subprocess, "run", fake_subprocess)
    first = validation.run(isolated, output, project_root=project)
    assert first["status"] == "success" and not first["search_reused"]
    assert first["query_count"] == 6 and first["reference_count"] == 6
    assert first["specificity_challenges"]["n_potential_cross_family_lysis_calls"] == 1
    assert first["specificity_challenges"]["biochemical_mAGP_inactivity_established"] is False
    assert {row["origin_evidence_type"] for row in first["positive_recoverability"]} == {
        "experimental", "published_profile_prediction", "experimental_generic_lipolytic"}
    assert len(read_tsv(output / "query_summary.tsv")) == 6
    assert len(read_tsv(output / "alignment_evidence.tsv")) == 3
    assert all(path.read_bytes() == content for path, content in originals.items())
    frozen = (output / "frozen_inputs.json").read_bytes()
    second = validation.run(isolated, output, project_root=project, offline=True)
    assert second["search_reused"] and len(searches) == 1
    assert (output / "frozen_inputs.json").read_bytes() == frozen
    (output / "queries.faa").write_text((output / "queries.faa").read_text() + "\n")
    with pytest.raises(ValueError, match="Frozen validation query input changed"):
        validation.run(isolated, output, project_root=project, offline=True)
    assert len(searches) == 1
