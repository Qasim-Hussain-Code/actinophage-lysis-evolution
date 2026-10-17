"""Analytical coordinate examples; all records here are synthetic test data."""
import json

import pytest

from actinolysis.architecture import NA, _validate, classify, pair_metrics


def gene(pid, start, end, index=0, strand=1, partial=False, parts=None, product="hypothetical protein"):
    row = dict(protein_id=pid, accession="TEST.1", start=start, end=end,
               cds_index=index, strand=strand, partial=partial, product=product)
    if parts is not None:
        row["location_parts"] = json.dumps(parts)
    return row


def genome(genes, topology="linear", length=1000):
    return dict(accession="TEST.1", phage="Synthetic", host_genus="test_only", cluster="test",
                subcluster=NA, topology=topology, length=length, cds_count=len(genes))


def ev(pid, category, confidence="high", source="reference_homology", **extra):
    return dict(protein_id=pid, category=category, confidence=confidence, source=source, **extra)


def test_linear_forward_reverse_and_intervening():
    a, other, b = gene("a", 100, 200, 0), gene("x", 250, 300, 1), gene("b", 400, 500, 2, -1)
    measured = pair_metrics(a, b, [a, other, b], 1000, "linear")
    assert measured["gap_bp"] == 200
    assert measured["absolute_midpoint_distance_bp"] == 300
    assert measured["midpoint_distance_bp"] == 300
    assert measured["intervening_cds"] == 1
    assert measured["relative_strands"] == "opposite_strands"
    assert measured["orientation"] == "convergent"
    assert measured["a_normalized_position"] == .15
    assert pair_metrics(b, a, [a, other, b], 1000, "linear")["orientation"] == "convergent"


def test_overlap_touching_and_nested_genes():
    a, b = gene("a", 100, 250), gene("b", 200, 300, strand=-1)
    measured = pair_metrics(a, b, [a, b], 1000, "linear")
    assert (measured["gap_bp"], measured["overlap_bp"], measured["intervening_cds"]) == (0, 50, 0)
    assert measured["orientation"] == "overlapping_opposite_strands"
    b = gene("b", 250, 300)
    assert pair_metrics(a, b, [a, b], 1000, "linear")["gap_bp"] == 0
    b = gene("b", 150, 175)
    assert pair_metrics(a, b, [a, b], 1000, "linear")["overlap_bp"] == 25


@pytest.mark.parametrize("topology", ["linear", "unknown", "not_available", "circularly_permuted"])
def test_linear_termini_are_not_adjacent(topology):
    a, x, b = gene("a", 10, 50, 0), gene("x", 400, 450, 1), gene("b", 950, 990, 2)
    measured = pair_metrics(a, b, [a, x, b], 1000, topology)
    assert measured["gap_bp"] == 900
    assert measured["midpoint_distance_bp"] == 940
    assert measured["intervening_cds"] == 1


def test_circular_origin_shortest_path_and_contained_crossing_gene():
    a, x, b = gene("a", 10, 50, 0), gene("x", 400, 450, 1), gene("b", 950, 990, 2)
    measured = pair_metrics(a, b, [a, x, b], 1000, "circular")
    assert measured["gap_bp"] == 20
    assert measured["absolute_midpoint_distance_bp"] == 940
    assert measured["midpoint_distance_bp"] == 60
    assert measured["intervening_cds"] == 0
    a, b = gene("a", 50, 100, 1), gene("b", 900, 950, 2)
    crossing = gene("x", 0, 1000, 0, parts=[[0, 20], [980, 1000]])
    assert pair_metrics(a, b, [crossing, a, b], 1000, "circular")["intervening_cds"] == 1


def test_joined_origin_cds_uses_parts_and_original_envelope():
    a = gene("a", 0, 1000, 0, parts=[[0, 40], [950, 1000]])
    x, b = gene("x", 60, 70, 1), gene("b", 100, 200, 2)
    measured = pair_metrics(a, b, [a, x, b], 1000, "circular")
    assert measured["gap_bp"] == 60
    assert measured["overlap_bp"] == 0
    assert measured["midpoint_distance_bp"] == 155
    assert measured["intervening_cds"] == 1
    assert measured["a_normalized_position"] == .995
    row, _ = classify(genome([a, x, b], "circular"), [a, x, b], [ev("a", "lysin_a"), ev("b", "lysin_b")])
    assert row["a_start"] == 0 and row["a_end"] == 1000
    assert json.loads(row["a_location_parts"]) == [[0, 40], [950, 1000]]


def test_circular_rotation_and_reverse_complement_preserve_distances():
    original = [gene("a", 10, 70, 0), gene("x", 150, 180, 1), gene("b", 300, 400, 2, -1)]
    expected = pair_metrics(original[0], original[2], original, 1000, "circular")
    rotated = [gene("a", 0, 1000, 0, parts=[[0, 20], [960, 1000]]),
               gene("x", 100, 130, 1), gene("b", 250, 350, 2, -1)]
    reverse = [gene("b", 600, 700, 0, 1), gene("x", 820, 850, 1, -1), gene("a", 930, 990, 2, -1)]
    for measured in [pair_metrics(rotated[0], rotated[2], rotated, 1000, "circular"),
                     pair_metrics(reverse[2], reverse[0], reverse, 1000, "circular")]:
        for field in ["gap_bp", "overlap_bp", "midpoint_distance_bp", "intervening_cds", "orientation"]:
            assert measured[field] == expected[field]


def test_linear_discontinuous_join_does_not_count_covered_gene_as_gap():
    a = gene("a", 100, 300, parts=[[100, 150], [250, 300]])
    x, b = gene("x", 310, 330), gene("b", 350, 400)
    measured = pair_metrics(a, b, [a, x, b], 1000, "linear")
    assert (measured["gap_bp"], measured["intervening_cds"]) == (50, 1)


def test_overlap_boundary_genes_count_in_cds_order_and_partial_cds_still_count():
    a, b = gene("a", 100, 200, 0), gene("b", 400, 500, 3)
    overlap = gene("overlap", 190, 250, 1)
    partial = gene("partial", 300, 350, 2, partial=True)
    measured = pair_metrics(a, b, [a, overlap, partial, b], 1000, "linear")
    assert measured["intervening_cds"] == 2
    assert measured["complete_gap_cds"] == 1
    assert measured["gene_order_basis"] == "linear_cds_index_order"


def test_d29_analytical_overlapping_holin_is_intervening_in_gene_order():
    a, holin, b = gene("a", 4696, 6178, 0), gene("holin", 6174, 6600, 1), gene("b", 6596, 7361, 2)
    measured = pair_metrics(a, b, [a, holin, b], 10000, "linear")
    assert measured["intervening_cds"] == 1
    assert measured["complete_gap_cds"] == 0
    assert measured["gap_bp"] == 418
    reverse = [gene("b", 2639, 3404, 0, -1), gene("holin", 3400, 3826, 1, -1), gene("a", 3822, 5304, 2, -1)]
    reversed_metrics = pair_metrics(reverse[2], reverse[0], reverse, 10000, "linear")
    assert reversed_metrics["intervening_cds"] == 1 and reversed_metrics["complete_gap_cds"] == 0
    assert reversed_metrics["gap_bp"] == 418


def test_circular_overlapping_boundary_gene_order_is_rotation_and_reversal_invariant():
    original = [gene("a", 100, 300, 0), gene("holin", 296, 354, 1), gene("b", 350, 500, 2)]
    rotated = [gene("a", 0, 1000, 0, parts=[[0, 100], [900, 1000]]), gene("holin", 96, 154, 1), gene("b", 150, 300, 2)]
    reverse = [gene("b", 500, 650, 0, -1), gene("holin", 646, 704, 1, -1), gene("a", 700, 900, 2, -1)]
    for genes, a_index, b_index in [(original, 0, 2), (rotated, 0, 2), (reverse, 2, 0)]:
        measured = pair_metrics(genes[a_index], genes[b_index], genes, 1000, "circular")
        assert measured["intervening_cds"] == 1
        assert measured["complete_gap_cds"] == 0
        assert measured["gap_bp"] == 50
        assert measured["midpoint_distance_bp"] == 225


def test_nested_linear_anchor_order_is_explicit_and_circular_overlap_does_not_wrap_long_way():
    a, x, b = gene("a", 100, 500, 0), gene("x", 150, 180, 1), gene("b", 200, 300, 2)
    measured = pair_metrics(a, b, [a, x, b], 1000, "linear")
    assert measured["intervening_cds"] == 1 and measured["complete_gap_cds"] == 0
    assert measured["overlap_bp"] == 100
    x, a, b = gene("x", 200, 250, 0), gene("a", 900, 980, 1), gene("b", 940, 995, 2)
    assert pair_metrics(a, b, [x, a, b], 1000, "circular")["intervening_cds"] == 0
    assert pair_metrics(b, a, [x, a, b], 1000, "circular")["intervening_cds"] == 0


def test_circular_tied_paths_preserve_ambiguous_gene_count():
    a, x, b = gene("a", 0, 100, 0), gene("x", 200, 250, 1), gene("b", 500, 600, 2)
    genes = [a, x, b]
    measured = pair_metrics(a, b, genes, 1000, "circular")
    assert measured["gap_bp"] == 400
    assert measured["shortest_path_count"] == 2
    assert measured["intervening_cds"] == NA
    assert (measured["intervening_cds_min"], measured["intervening_cds_max"]) == (0, 1)
    row, _ = classify(genome(genes, "circular"), genes, [ev("a", "lysin_a"), ev("b", "lysin_b")])
    assert row["architecture"] == "unresolved"


def test_duplicate_evidence_does_not_duplicate_candidates():
    a, b = gene("a", 100, 200, 0), gene("b", 400, 500, 1)
    evidence = [ev("a", "lysin_a", "moderate"), ev("a", "lysin_a"), ev("a", "lysin_a"), ev("b", "lysin_b")]
    row, pairs = classify(genome([a, b]), [a, b], evidence)
    assert row["architecture"] == "closely_linked"
    assert row["lysin_a_count"] == 1 and len(pairs) == 1
    assert row["confidence"] == "high"


def test_multiple_copies_and_conflicting_assignments_remain_unselected():
    a, a2, b = gene("a", 100, 200, 0), gene("a2", 250, 300, 1), gene("b", 400, 500, 2)
    evidence = [ev("a", "lysin_a"), ev("a2", "lysin_a"), ev("b", "lysin_b")]
    row, pairs = classify(genome([a, a2, b]), [a, a2, b], evidence)
    assert row["architecture"] == "complex" and row["lysin_a"] == NA
    assert len(pairs) == 2
    evidence.append(ev("a2", "lysin_b"))
    row, _ = classify(genome([a, a2, b]), [a, a2, b], evidence)
    assert row["architecture"] == "unresolved"
    assert row["conflicting_protein_ids"] == "a2"


def test_partial_and_confidence_scope_annotation_sensitivity():
    a, b = gene("a", 100, 200, 0), gene("b", 400, 500, 1)
    evidence = [ev("a", "lysin_a", annotation_category="lysin_a"),
                ev("b", "lysin_b", "moderate", annotation_category="lysin_b")]
    assert classify(genome([a, b]), [a, b], evidence, scope="high")[0]["architecture"] == "unresolved"
    assert classify(genome([a, b]), [a, b], evidence, scope="high_moderate")[0]["confidence"] == "moderate"
    assert classify(genome([a, b]), [a, b], evidence, scope="annotation")[0]["confidence"] == "annotation_only"
    b["partial"] = True
    row, pairs = classify(genome([a, b]), [a, b], evidence)
    assert row["architecture"] == "unresolved" and row["partial_lysis_count"] == 1 and not pairs


def test_threshold_boundary_holin_dedup_and_structural_qualification():
    genes = [gene("a", 0, 50, 0)] + [gene(f"x{i}", 100 * i, 100 * i + 50, i) for i in range(1, 6)]
    genes.append(gene("b", 600, 650, 6))
    genes[3]["product"] = "major capsid protein"
    evidence = [ev("a", "lysin_a"), ev("b", "lysin_b"), ev("x1", "holin_like"), ev("x1", "holin_like")]
    row, _ = classify(genome(genes), genes, evidence, threshold=4)
    assert row["architecture"] == "separated" and row["intervening_cds"] == 5
    assert row["nearby_holin_count"] == 1 and row["other_lysis_count"] == 1
    assert row["structural_distance_cds"] == 3
    assert "not_validated_module" in row["structural_basis"]
    assert classify(genome(genes), genes, evidence, threshold=5)[0]["architecture"] == "closely_linked"


def test_coordinate_input_validation():
    a, b = gene("a", 100, 200, 0), gene("b", 400, 500, 1)
    assert _validate([genome([a, b])], [a, b], [])
    with pytest.raises(ValueError, match="Duplicate protein"):
        _validate([genome([a, b])], [a, a], [])
    bad = gene("b", 400, 500, 2)
    with pytest.raises(ValueError, match="indices"):
        _validate([genome([a, bad])], [a, bad], [])
    with pytest.raises(ValueError, match="unknown protein"):
        _validate([genome([a, b])], [a, b], [ev("missing", "lysin_a")])
    bad = gene("b", 400, 1100, 1)
    with pytest.raises(ValueError, match="coordinates"):
        _validate([genome([a, bad])], [a, bad], [])
    bad = gene("b", 400, 500, 1, parts=[[410, 500]])
    with pytest.raises(ValueError, match="envelope"):
        _validate([genome([a, bad])], [a, bad], [])


def test_stage_writes_deterministic_tables_and_empty_pair_schema(tmp_path, monkeypatch):
    import actinolysis.architecture as architecture
    from actinolysis.common import read_tsv, write_tsv
    a, b = gene("a", 100, 200, 0), gene("b", 400, 500, 1)
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    write_tsv(tmp_path / "data/processed/genomes.tsv", [genome([a, b])])
    write_tsv(tmp_path / "data/processed/cds.tsv", [a, b])
    write_tsv(tmp_path / "data/processed/lysis_evidence.tsv", [ev("a", "lysin_a"), ev("b", "lysin_b", "moderate")])
    cfg = {"architecture": {"separation_threshold": 4, "thresholds": [2, 4, 6, 10]}}
    result = architecture.run(cfg)
    assert result == {"genomes": 1, "eligible_pairs": 1, "sensitivity_rows": 12}
    table = tmp_path / "results/tables/architecture_sensitivity.tsv"
    before = table.read_bytes()
    architecture.run(cfg)
    assert table.read_bytes() == before
    assert len(read_tsv(table)) == 12
    write_tsv(tmp_path / "data/processed/lysis_evidence.tsv", [ev("a", "lysin_a")])
    assert architecture.run(cfg)["eligible_pairs"] == 0
    assert read_tsv(tmp_path / "data/processed/architecture_pairs.tsv") == []


def test_stage_retains_all_scopes_and_all_annotation_pairs_for_complex_genome(tmp_path, monkeypatch):
    import actinolysis.architecture as architecture
    from actinolysis.common import read_tsv, write_tsv
    a, a2, b = gene("a", 100, 200, 0), gene("a2", 250, 300, 1), gene("b", 800, 900, 2)
    monkeypatch.setattr(architecture, "ROOT", tmp_path)
    write_tsv(tmp_path / "data/processed/genomes.tsv", [genome([a, a2, b])])
    write_tsv(tmp_path / "data/processed/cds.tsv", [a, a2, b])
    evidence = [ev("a", "lysin_a", annotation_category="lysin_a"),
                ev("a2", "lysin_a", "low", annotation_category="lysin_a"),
                ev("b", "lysin_b", annotation_category="lysin_b")]
    write_tsv(tmp_path / "data/processed/lysis_evidence.tsv", evidence)
    architecture.run({"architecture": {"separation_threshold": 4, "thresholds": [2, 4]}})
    pairs = read_tsv(tmp_path / "data/processed/architecture_pairs.tsv")
    assert [row["scope"] for row in pairs].count("high") == 1
    assert [row["scope"] for row in pairs].count("high_moderate") == 1
    assert [row["scope"] for row in pairs].count("annotation") == 2
    annotations = [row for row in read_tsv(tmp_path / "results/tables/architecture_sensitivity.tsv")
                   if row["scope"] == "annotation" and row["threshold"] == "4"]
    assert annotations[0]["architecture"] == "complex" and annotations[0]["lysin_a"] == NA
