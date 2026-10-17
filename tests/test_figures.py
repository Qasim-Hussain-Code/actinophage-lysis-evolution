"""Figure-data consistency using explicitly synthetic test-only tables."""
from collections import defaultdict
import json

import matplotlib.pyplot as plt
import pytest

from actinolysis import figures
from actinolysis.common import NA, read_tsv, sha256, write_tsv


def synthetic_inputs():
    genomes, cds, evidence, architectures, sensitivities, families = [], [], [], [], [], []
    cases = [("SYN000001.1", "D29", "A", "closely_linked", "high", 3),
             ("SYN000002.1", "Yvonnetastic", "B", "separated", "moderate", 6),
             ("SYN000003.1", "Synthetic test-only 3", "A", "closely_linked", "high", 2),
             ("SYN000004.1", "Synthetic test-only 4", "A", "unresolved", "unresolved", None)]
    for number, (accession, name, cluster, architecture, confidence, b_index) in enumerate(cases):
        genus = "Mycobacterium" if number % 2 == 0 else "Gordonia"
        genomes.append(dict(accession=accession, phage=name, host_genus=genus, cluster=cluster, length=1000 + 100 * number))
        for index in range(8):
            pid = f"{accession}__cds{index:04d}"
            start, end = index * 100 + 20, index * 100 + 80
            cds.append(dict(protein_id=pid, accession=accession, start=start, end=end, strand=1 if index % 2 else -1, cds_index=index, partial="false", location_parts=json.dumps([[start, end]])))
            families.append(dict(protein_id=pid, accession=accession, family_id=f"FAM_SYNTHETIC_TEST_ONLY_{index}"))
            if index == 0 or index == b_index:
                evidence.append(dict(protein_id=pid, category="lysin_a" if index == 0 else "lysin_b", confidence=confidence if b_index is not None else "low", source="synthetic_test_only"))
        a = f"{accession}__cds0000"
        b = f"{accession}__cds{b_index:04d}" if b_index is not None else NA
        row = dict(accession=accession, phage=name, cluster=cluster, host_genus=genus, architecture=architecture, confidence=confidence,
                   lysin_a=a, lysin_b=b, lysin_a_count=1, lysin_b_count=1 if b_index is not None else 0,
                   intervening_cds=b_index - 1 if b_index is not None else NA, gap_bp=b_index * 100 - 60 if b_index is not None else NA)
        architectures.append(row)
        for scope in ("high", "high_moderate", "annotation"):
            for threshold in (2, 4, 6, 10):
                category = "unresolved" if b_index is None or (scope == "high" and confidence != "high") else "separated" if b_index - 1 > threshold else "closely_linked"
                sensitivities.append({**row, "scope": scope, "threshold": threshold, "architecture": category})
    return dict(genomes=genomes, cds=cds, evidence=evidence, architectures=architectures, sensitivities=sensitivities, families=families, pairs=[])


@pytest.fixture
def cfg():
    return {"architecture": {"separation_threshold": 4}, "neighborhood": {"primary_identity": .4, "primary_window": 5},
            "statistics": {"min_controls": 5, "min_lineages": 5}}


def test_lineage_weighting_preserves_unresolved_and_does_not_duplicate_large_cluster():
    inputs = synthetic_inputs()
    fig, rows, _ = figures.lineage_figure(inputs["genomes"], inputs["architectures"])
    try:
        weighted = {r["architecture"]: r["mean_fraction"] for r in rows if r["panel"] == "equal_weight_clusters"}
        raw = {r["architecture"]: r["n_genomes"] for r in rows if r["panel"] == "raw_counts"}
        assert raw == {"closely_linked": 2, "separated": 1, "complex": 0, "unresolved": 1}
        assert weighted["closely_linked"] == pytest.approx(1 / 3)
        assert weighted["separated"] == .5
        assert weighted["unresolved"] == pytest.approx(1 / 6)
        assert sum(weighted.values()) == pytest.approx(1)
    finally:
        plt.close(fig)


def test_confidence_counts_deduplicate_evidence_and_keep_strongest():
    inputs = synthetic_inputs()
    evidence = inputs["evidence"] + [dict(inputs["evidence"][0]), {**inputs["evidence"][0], "confidence": "moderate"}]
    fig, rows, _ = figures.robustness_figure(inputs["sensitivities"], evidence)
    try:
        counts = {(row["category"], row["confidence"]): row["n_proteins"] for row in rows if row["panel"] == "functional_confidence"}
        assert counts[("lysin_a", "high")] == 2
        assert counts[("lysin_a", "moderate")] == 1
        assert counts[("lysin_a", "low")] == 1
        assert counts[("lysin_b", "high")] == 2
        assert counts[("lysin_b", "moderate")] == 1
    finally:
        plt.close(fig)


def test_zero_neighborhood_pairs_produces_observed_diagnostic(cfg):
    inputs = synthetic_inputs()
    fig, rows, caption = figures.neighborhood_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], inputs["families"], [], cfg)
    try:
        pair_counts = {row["criterion"]: row["n_pairs"] for row in rows if row["panel"] == "pair_eligibility"}
        assert pair_counts == {"Within-cluster genome pairs": 3, "Unique supported B in both genomes": 1, "Shared B homologous family": 0, "Eligible matched comparisons": 0}
        assert "No primary genome-pair comparison" in caption
        assert all(row["panel"] != "observed_vs_background" for row in rows)
        assert "H2 unavailable" in fig._suptitle.get_text()
    finally:
        plt.close(fig)


def test_neighborhood_score_figure_uses_actual_values_and_rejects_inconsistent_delta(cfg):
    inputs = synthetic_inputs()
    pair = dict(pair_id="SYN000001.1|SYN000003.1", cluster="A", family_identity=.4, window=5, architecture_threshold=4,
                confidence_set="high_moderate", genome_subset="all", metric="jaccard", status="eligible", reason="",
                pair_architecture="separated", target_score=.8, control_mean=.3, delta=.5, n_controls=5)
    fig, rows, _ = figures.neighborhood_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], inputs["families"], [pair], cfg)
    try:
        observed = next(row for row in rows if row["panel"] == "observed_vs_background")
        assert observed["target_score"] == .8 and observed["control_mean"] == .3 and observed["delta"] == .5
        assert next(row for row in rows if row["panel"] == "within_cluster_contrast")["mean_delta"] == .5
    finally:
        plt.close(fig)
    with pytest.raises(ValueError, match="delta disagrees"):
        figures.neighborhood_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], inputs["families"], [{**pair, "delta": .6}], cfg)
    plt.close("all")


def test_genome_maps_preserve_original_parts_and_explicit_low_confidence(cfg):
    inputs = synthetic_inputs()
    original = inputs["cds"][0]
    original.update(start=0, end=1000, strand=0, location_parts="[[0,20],[980,1000]]")
    inputs["evidence"][0]["confidence"] = "low"
    fig, rows, _ = figures.genome_map_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], figures.primary_annotations(inputs["sensitivities"], cfg), inputs["families"])
    try:
        parts = [row for row in rows if row["panel"] == "cds_parts" and row["protein_id"] == original["protein_id"]]
        assert [(row["start"], row["end"]) for row in parts] == [(0, 20), (980, 1000)]
        assert all(row["confidence"] == "low" and row["category"] == "lysin_a" for row in parts)
        assert all(row["arrow_type"] == "unoriented" for row in parts)
        label = next(row for row in rows if row.get("label") == "A (low)")
        assert label["plot_x"] == 10
    finally:
        plt.close(fig)


def test_adjacent_same_category_map_labels_use_separate_heights(cfg):
    inputs = synthetic_inputs()
    inputs["cds"][1].update(start=60, end=100, location_parts="[[60,100]]")
    inputs["evidence"].append(dict(protein_id=inputs["cds"][1]["protein_id"], category="lysin_a", confidence="low", source="synthetic_test_only"))
    fig, rows, _ = figures.genome_map_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], figures.primary_annotations(inputs["sensitivities"], cfg), inputs["families"])
    try:
        labels = [row for row in rows if row["panel"] == "lysis_labels" and row["accession"] == "SYN000001.1" and row["category"] == "lysin_a"]
        assert len(labels) == 2 and labels[0]["plot_y"] != labels[1]["plot_y"]
        assert len({row["protein_id"] for row in labels}) == 2
    finally:
        plt.close(fig)


def _write_inputs(root, inputs):
    paths = dict(genomes="data/processed/genomes.tsv", cds="data/processed/cds.tsv", architectures="data/processed/architectures.tsv",
                 evidence="data/processed/lysis_evidence.tsv", sensitivities="results/tables/architecture_sensitivity.tsv",
                 families="data/processed/families.tsv", pairs="results/tables/neighborhood_pairs.tsv")
    for key, path in paths.items():
        fields = list(inputs[key][0]) if inputs[key] else ["pair_id", "family_identity", "window", "architecture_threshold", "confidence_set", "genome_subset", "metric", "status", "reason"]
        write_tsv(root / path, inputs[key], fields)


def test_six_figure_artifacts_and_supporting_values(tmp_path, monkeypatch, cfg):
    from PIL import Image
    inputs = synthetic_inputs()
    _write_inputs(tmp_path, inputs)
    monkeypatch.setattr(figures, "ROOT", tmp_path)
    report = figures.run(cfg)
    assert report["created_count"] == 6
    assert len(report["source_inputs"]) == 7
    assert not plt.get_fignums()
    for figure in report["figures"]:
        assert figure["status"] == "created"
        assert figure["data_rows"] > 0
        for artifact in figure["files"]:
            path = tmp_path / artifact["path"]
            assert path.is_file() and sha256(path) == artifact["sha256"]
            if path.suffix == ".tsv":
                assert len(read_tsv(path)) == figure["data_rows"]
            elif path.suffix == ".png":
                with Image.open(path) as image:
                    assert image.info["dpi"][0] == pytest.approx(300, abs=.01)
                    assert image.width > 2000
            else:
                content = path.read_bytes()
                assert content.startswith(b"%PDF")
                assert b"D:20000101000000" in content
    table = read_tsv(tmp_path / "results/tables/figure_01_dataset_data.tsv")
    assert sum(int(row["n_genomes"]) for row in table if row["panel"] == "host_representation") == 4
    architecture_data = read_tsv(tmp_path / "results/tables/figure_02_architecture_data.tsv")
    assert sum(int(row["n_genomes"]) for row in architecture_data if row["panel"] == "architecture_status" and row["scope"] == "high_moderate") == 4
    manifest = json.loads((tmp_path / "results/reports/figure_manifest.json").read_text())
    assert manifest == report


def test_unsupportable_map_is_skipped_with_explicit_status(tmp_path):
    fig, rows, reason = figures.genome_map_figure([dict(accession="SYN000001.1", phage="Synthetic test only", length=1000)], [], [], [dict(accession="SYN000001.1", architecture="unresolved")], [], [])
    result = figures._save(tmp_path, "test_only_unsupported_map", "Synthetic test only", (fig, rows, reason))
    assert result["status"] == "skipped" and result["files"] == []
    assert "No control genome" in result["reason"]
    fig, rows, reason = figures.architecture_figure(synthetic_inputs()["architectures"], [])
    assert fig is None and not rows and "Annotation-only architecture rows are unavailable" in reason


def test_fixed_metadata_makes_repeated_artifacts_identical(tmp_path):
    checks = []
    for _ in range(2):
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.plot([0, 1], [2, 3])
        result = figures._save(tmp_path, "synthetic_test_only_repeatability", "Synthetic test only", (fig, [dict(x=0, y=2), dict(x=1, y=3)], "Synthetic test only"))
        checks.append([row["sha256"] for row in result["files"]])
    assert checks[0] == checks[1]
