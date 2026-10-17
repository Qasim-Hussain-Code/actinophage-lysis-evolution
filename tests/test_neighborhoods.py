"""Analytical and explicitly synthetic tests; these are not biological results."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from actinolysis.neighborhoods import (cluster_proteins, compare_neighborhoods, completeness, jaccard, matched_controls,
                                      neighborhood, parse_clusters, read_fasta, signed_order)


def genes(accession="A", count=100):
    return [{"protein_id": f"{accession}_{i}", "accession": accession, "cds_index": str(i),
             "start": str(i * 100), "end": str(i * 100 + 90), "strand": "1", "partial": "false",
             "translation": "M" * 100, "product": "hypothetical protein", "family": f"F{i}"}
            for i in range(count)]


def config():
    return {"neighborhood": {"windows": [3], "primary_window": 3, "primary_identity": 0.4},
            "architecture": {"separation_threshold": 4, "thresholds": [4]},
            "statistics": {"min_controls": 5, "matched_controls": 20}}


def comparison_fixture():
    ga, gb = genes("A"), genes("B")
    genomes = [{"accession": acc, "length": "10000", "topology": "linear", "cluster": "TEST",
                "representative": "true", "similarity_group": acc} for acc in ("A", "B")]
    evidence = [{"protein_id": f"{acc}_50", "category": "lysin_b", "confidence": "high"} for acc in ("A", "B")]
    architectures = [{"accession": acc, "architecture": "separated", "lysin_b": f"{acc}_50", "confidence": "high", "intervening_cds": 40}
                     for acc in ("A", "B")]
    assignments = {g["protein_id"]: g["family"] for g in ga + gb}
    members = [{"family_identity": 0.4, "family_id": f"F{i}", "prevalence": 1.0} for i in range(100)]
    return genomes, ga + gb, evidence, architectures, {0.4: assignments}, members


def test_linear_boundary_is_clipped_and_anchor_excluded():
    data = genes(count=8)
    rows = neighborhood(data, "A_0", 3)
    assert [r["protein_id"] for r in rows] == ["A_1", "A_2", "A_3"]
    assert completeness(rows) == (0, 3)
    assert "A_7" not in {r["protein_id"] for r in rows}
    assert len(neighborhood(data, "A_0", 20)) == 7


def test_explicit_circular_boundary_wraps_once():
    data = genes(count=8)
    rows = neighborhood(data, "A_0", 3, circular=True)
    assert [r["protein_id"] for r in rows] == ["A_5", "A_6", "A_7", "A_1", "A_2", "A_3"]
    assert completeness(rows) == (3, 3)
    assert len(neighborhood(data, "A_0", 20, circular=True)) == 7


def test_whole_genome_reverse_complement_preserves_signed_order():
    original = genes(count=9)
    original[2]["strand"] = "-1"
    reversed_genes = [{**r, "start": str(900 - int(r["end"])), "end": str(900 - int(r["start"])),
                       "strand": str(-int(r["strand"]))} for r in original]
    left = neighborhood(original, "A_4", 3)
    right = neighborhood(reversed_genes, "A_4", 3)
    assert left == right
    assert signed_order(left, right) == 1


def test_short_circular_antipodal_tie_is_reverse_complement_invariant():
    original = genes(count=8)
    reverse = [{**r, "start": str(800 - int(r["end"])), "end": str(800 - int(r["start"])), "strand": "-1"}
               for r in original]
    left = neighborhood(original, "A_0", 20, circular=True)
    right = neighborhood(reverse, "A_0", 20, circular=True)
    assert left == right
    assert signed_order(left, right) == 1


def test_order_orientation_and_side_matter_separately_from_content():
    original = neighborhood(genes(count=9), "A_4", 3)
    inverted = deepcopy(original)
    inverted[0]["relative_strand"] = -1
    reordered = deepcopy(original)
    reordered[0]["family"], reordered[1]["family"] = reordered[1]["family"], reordered[0]["family"]
    assert jaccard([r["family"] for r in original], [r["family"] for r in inverted]) == 1
    assert signed_order(original, inverted) == pytest.approx(5 / 6)
    assert signed_order(original, reordered) < 1
    assert jaccard([], []) is None
    assert signed_order([], []) is None
    assert jaccard(["a", "a", "b"], ["a", "c"]) == pytest.approx(1 / 3)


def test_family_ids_depend_on_sequence_content_not_member_order_or_representative(tmp_path):
    sequences = {"a": "MALA", "b": "MGLA", "c": "MPP"}
    first, second = tmp_path / "one.tsv", tmp_path / "two.tsv"
    first.write_text("a\ta\na\tb\nc\tc\n")
    second.write_text("c\tc\nb\tb\nb\ta\n")
    assignments_a, _ = parse_clusters(first, sequences)
    assignments_b, _ = parse_clusters(second, sequences)
    assert assignments_a == assignments_b
    sequences["b"] = "MGLL"
    assignments_c, _ = parse_clusters(first, sequences)
    assert assignments_c["a"] != assignments_a["a"]


@pytest.mark.parametrize("contents", ["a\ta\na\ta\nb\tb\n", "a\ta\n", "a\ta\nb\tunknown\n", "b\ta\na\tb\n"])
def test_cluster_integrity_rejects_duplicate_missing_unknown_and_nonmember_representative(tmp_path, contents):
    path = tmp_path / "families.tsv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        parse_clusters(path, {"a": "MALA", "b": "MGLA"})


def test_fasta_duplicate_identifier_is_not_silently_overwritten(tmp_path):
    path = tmp_path / "proteins.faa"
    path.write_text(">a\nMALA\n>a\nMGLA\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_fasta(path)


def test_identical_short_peptides_may_have_distinct_significance_filtered_families(tmp_path):
    path = tmp_path / "families.tsv"
    path.write_text("a\ta\nb\tb\n")
    assignments, families = parse_clusters(path, {"a": "MAP", "b": "MAP"})
    assert assignments["a"] != assignments["b"]
    assert len(families) == 2


def test_mmseqs_invocation_limits_resources_and_validates_cached_assignments(tmp_path, monkeypatch):
    import actinolysis.neighborhoods as module
    proteome = tmp_path / "data/processed/proteins.faa"
    proteome.parent.mkdir(parents=True)
    proteome.write_text(">a\nMALA\n>b\nMGLA\n")
    sequences = {"a": "MALA", "b": "MGLA"}
    cds = [{"protein_id": p, "accession": "A", "translation": s, "product": "hypothetical protein"}
           for p, s in sequences.items()]
    # This partial CDS is deliberately excluded from the input FASTA.
    cds.append({"protein_id": "partial", "accession": "A", "translation": "M", "partial": "true"})
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        if command[1] == "version":
            return SimpleNamespace(stdout="TEST_ONLY_MMSEQS")
        prefix = Path(command[3])
        Path(str(prefix) + "_cluster.tsv").write_text("a\ta\na\tb\n")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.setattr(module, "budget_check", lambda *args, **kwargs: None)
    configuration = {"neighborhood": {"family_identities": [0.3, 0.4, 0.5], "primary_identity": 0.4,
                        "coverage": 0.8, "evalue": 0.001}, "threads": 9,
                     "resources": {"mmseqs_memory": "500M"}, "tools": {"mmseqs": "mmseqs"}}
    assignments, memberships, audit = cluster_proteins(configuration, sequences, cds, [{"accession": "A"}], tmp_path)
    assert len(assignments) == 3 and len(memberships) == 6
    assert len(commands) == 4
    for command in commands[1:]:
        for flag, expected in [("--threads", "4"), ("--split-memory-limit", "500M"), ("--cov-mode", "0"),
                               ("--seq-id-mode", "0"), ("--alignment-mode", "3"), ("--single-step-clustering", "1"),
                               ("-k", "5"), ("-c", "0.8"), ("-e", "0.001")]:
            assert command[command.index(flag) + 1] == expected
    assert all(not row["cache_reused"] for row in audit)
    _, _, audit = cluster_proteins(configuration, sequences, cds, [{"accession": "A"}], tmp_path)
    assert all(row["cache_reused"] for row in audit)
    assert len(commands) == 5  # Version checked; the three workflows are reused.
    Path(str(tmp_path / "data/processed/mmseqs_work/040/families") + "_cluster.tsv").write_text("a\ta\n")
    with pytest.raises(ValueError, match="Unassigned"):
        cluster_proteins(configuration, sequences, cds, [{"accession": "A"}], tmp_path)


def test_matched_null_has_zero_contrast_and_uses_homologous_nonlysis_anchors():
    fixture = comparison_fixture()
    pairs, controls, _ = compare_neighborhoods(config(), *fixture)
    eligible = [r for r in pairs if r["status"] == "eligible"]
    assert eligible and controls
    assert all(r["delta"] == 0 for r in eligible)
    assert all(r["n_controls"] >= 5 for r in eligible)
    assert all(r["anchor_a"].split("_")[1] == r["anchor_b"].split("_")[1] for r in controls)
    assert all(r["control_family"] != "F50" for r in controls)
    repeated = compare_neighborhoods(config(), *fixture)
    assert (pairs, controls) == repeated[:2]


def test_control_selection_does_not_use_conservation_scores():
    left, right = genes("A"), genes("B")
    windows = {r["protein_id"]: neighborhood(rows, r["protein_id"], 3)
               for rows in (left, right) for r in rows}
    prevalence = {f"F{i}": 1.0 for i in range(100)}
    genome = {"length": "10000"}
    first, _ = matched_controls(left, right, left[50], right[50], windows, prevalence, {"F50"}, genome, genome, 5)
    changed = deepcopy(windows)
    for key, rows in changed.items():
        if key.startswith("B_"):
            for index, row in enumerate(rows):
                row["family"] = f"new_{key}_{index}"
    second, _ = matched_controls(left, right, left[50], right[50], changed, prevalence, {"F50"}, genome, genome, 5)
    assert [r["control_family"] for r in first] == [r["control_family"] for r in second]
    assert all(r["jaccard"] == 1 for r in first)
    assert all(r["jaccard"] == 0 for r in second)


def test_insufficient_controls_are_retained_as_untestable():
    fixture = list(comparison_fixture())
    # Only the Lysin B family has prevalence near the target, precluding controls.
    fixture[-1] = [{**r, "prevalence": 1.0 if r["family_id"] == "F50" else 0.1} for r in fixture[-1]]
    pairs, _, _ = compare_neighborhoods(config(), *fixture)
    assert pairs and all(r["status"] == "untestable" for r in pairs)
    assert all(r["delta"] == "" for r in pairs)
    assert all(r["reason"] == "fewer_than_minimum_matched_homologous_controls" for r in pairs)


def test_missing_family_cds_are_distinct_and_not_controls():
    fixture = list(comparison_fixture())
    fixture[4] = {0.4: {p: f for p, f in fixture[4][0.4].items() if p not in {"A_49", "B_49"}}}
    _, controls, neighborhoods = compare_neighborhoods(config(), *fixture)
    target = [r for r in neighborhoods if r["anchor"] == "A_50"][0]
    assert target["n_unclustered_neighbors"] == 1
    assert "UNCLUSTERED::A_49" in target["neighbor_families"]
    assert all(not r["control_family"].startswith("UNCLUSTERED") for r in controls)


def test_near_identical_pair_is_not_independent_evidence():
    fixture = list(comparison_fixture())
    fixture[0][0]["similarity_group"] = fixture[0][1]["similarity_group"] = "same"
    pairs, _, _ = compare_neighborhoods(config(), *fixture)
    assert all(r["status"] == "untestable" for r in pairs)
    assert all(r["reason"] == "same_near_identical_similarity_group" for r in pairs)


def test_high_only_rejects_moderate_architecture():
    fixture = list(comparison_fixture())
    fixture[3][0]["confidence"] = "moderate"
    pairs, _, _ = compare_neighborhoods(config(), *fixture)
    assert all(r["status"] == "untestable" for r in pairs if r["confidence_set"] == "high_only")
    assert all(r["status"] == "eligible" for r in pairs if r["confidence_set"] == "high_moderate")


def test_threshold_sensitivities_do_not_change_matched_controls():
    fixture = comparison_fixture()
    settings = config()
    settings["architecture"]["thresholds"] = [2, 4, 6, 10]
    fixture[3][0]["intervening_cds"] = fixture[3][1]["intervening_cds"] = 5
    pairs, controls, _ = compare_neighborhoods(settings, *fixture)
    assert {r["architecture_threshold"] for r in pairs} == {2, 4, 6, 10}
    assert all(r["pair_architecture"] == "separated" for r in pairs if r["architecture_threshold"] <= 4)
    assert all(r["pair_architecture"] == "closely_linked" for r in pairs if r["architecture_threshold"] >= 6)
    assert len(controls) < len(pairs) * 20  # Selection is recorded once per confidence/window, not repeated per cutoff.


def test_high_only_uses_independently_reconstructed_architecture():
    fixture = list(comparison_fixture())
    fixture[3][0]["architecture"] = fixture[3][1]["architecture"] = "complex"
    fixture[3][0]["confidence"] = fixture[3][1]["confidence"] = "unresolved"
    high = [{"accession": acc, "scope": "high", "threshold": 4, "architecture": "separated",
             "lysin_b": f"{acc}_50", "intervening_cds": 20, "confidence": "high"} for acc in ("A", "B")]
    pairs, _, _ = compare_neighborhoods(config(), *fixture, architecture_sensitivities=high)
    assert all(r["status"] == "untestable" for r in pairs if r["confidence_set"] == "high_moderate")
    assert all(r["status"] == "eligible" for r in pairs if r["confidence_set"] == "high_only")
