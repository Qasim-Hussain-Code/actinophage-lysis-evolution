"""Synthetic component configurations and analytical control-availability checks."""
import json

import pytest
import yaml

from actinolysis import amendment_feasibility as amendment
from actinolysis import neighborhoods
from actinolysis.common import sha256, write_tsv


def fixture():
    genomes, cds, evidence, primary, families = [], [], [], [], []
    for accession in ("SYNTH_A", "SYNTH_B"):
        genomes.append({"accession": accession, "phage": accession, "cluster": "SYNTH_CLUSTER", "length": "10000",
                        "cds_count": "100", "topology": "linear", "similarity_group": accession})
        primary.append({"accession": accession, "architecture": "complex"})
        for i in range(100):
            protein = f"{accession}__cds{i:04d}"
            cds.append({"protein_id": protein, "accession": accession, "cds_index": str(i), "start": str(i * 100),
                        "end": str(i * 100 + 90), "strand": "1", "partial": "false", "product": "hypothetical protein",
                        "translation": "M" * 100})
            families.append({"family_identity": 0.4, "protein_id": protein, "family_id": f"SYNTH_F{i}", "prevalence": 1.0})
        for i, category in ((5, "lysin_a"), (10, "lysin_a"), (50, "lysin_b")):
            evidence.append({"protein_id": f"{accession}__cds{i:04d}", "category": category, "confidence": "high",
                             "source": "curated_reference_homology", "reference_id": "SYNTH_TEST_REFERENCE",
                             "evidence_url": "https://example.invalid/synthetic-test-only", "context_used": "false"})
    return genomes, cds, evidence, primary, families


def cfg():
    return {"mode": "pilot", "architecture": {"thresholds": [4], "separation_threshold": 4},
            "neighborhood": {"family_identities": [0.4], "windows": [5]},
            "statistics": {"min_controls": 5, "matched_controls": 20, "min_lineages": 5}}


def one_genome(data):
    genomes, cds, evidence, primary, _ = data
    accession = genomes[0]["accession"]
    return genomes[0], [g for g in cds if g["accession"] == accession], [r for r in evidence if r["protein_id"].startswith(accession)], primary[0]


def test_all_components_distal_preserves_complex_primary_and_every_pair():
    result, calls, pairs = amendment.component_classification(*one_genome(fixture()), 4)
    assert result["primary_architecture"] == result["primary_scope_architecture"] == "complex"
    assert result["secondary_classification"] == "all_detected_components_distal"
    assert result["min_intervening_cds"] == 39 and result["max_intervening_cds"] == 44
    assert len(pairs) == result["n_component_pairs"] == 2 and len(calls) == 3
    assert {p["lysin_a"] for p in pairs} == set(result["detected_a"].split(";"))
    assert all(json.loads(p["measurements_json"])["gene_order_basis"] == "linear_cds_index_order" for p in pairs)


def test_nearest_component_determines_mixed_not_a_favorable_distal_pair():
    data = fixture()
    data[2][0]["protein_id"] = "SYNTH_A__cds0049"
    result, _, pairs = amendment.component_classification(*one_genome(data), 4)
    assert result["secondary_classification"] == "mixed_detected_components"
    assert result["min_intervening_cds"] == 0 and result["max_intervening_cds"] == 39 and len(pairs) == 2


def test_missing_a_is_not_vacuously_distal_and_multiple_b_retains_pairs():
    data = fixture()
    args = list(one_genome(data))
    args[2] = [r for r in args[2] if r["category"] == "lysin_b"]
    result, _, _ = amendment.component_classification(*args, 4)
    assert result["status"] == "unresolved" and "no_complete_independently_supported_a" in result["reason"]
    args = list(one_genome(data))
    args[2].append({**args[2][-1], "protein_id": "SYNTH_A__cds0060"})
    result, _, pairs = amendment.component_classification(*args, 4)
    assert result["status"] == "unresolved" and result["n_detected_b"] == 2 and len(pairs) == 4


def test_partial_and_conflicting_calls_cannot_establish_all_components_distal():
    data = fixture()
    data[1][5]["partial"] = "true"
    result, calls, _ = amendment.component_classification(*one_genome(data), 4)
    assert result["status"] == "unresolved" and "partial_independently_supported_component" in result["reason"]
    assert any(c["complete"] == "false" for c in calls)
    data = fixture()
    data[2].append({**data[2][0], "protein_id": "SYNTH_A__cds0050"})
    result, _, pairs = amendment.component_classification(*one_genome(data), 4)
    assert result["status"] == "unresolved" and "conflicting_a_b_assignments" in result["reason"]
    assert any(p["distance_reason"] == "conflicting_same_protein" for p in pairs)


@pytest.mark.parametrize("change,reason", [({"source": "genbank_annotation"}, "source_not_independent_functional_evidence"),
    ({"context_used": "true"}, "contextual_evidence"), ({"reference_id": "not_available"}, "missing_reference_identifier"),
    ({"evidence_url": "not_available"}, "missing_traceable_evidence_url")])
def test_annotation_or_untraceable_evidence_is_not_promoted(change, reason):
    row = {**fixture()[2][0], **change}
    assert amendment.independent_evidence(row, "high_moderate") == (False, reason)


def test_high_only_recomputes_sets_and_counts_only_removed_assignments():
    data = fixture()
    data[2][0]["confidence"] = "moderate"
    data[2].append({**data[2][2], "confidence": "moderate"})
    result, _, pairs = amendment.component_classification(*one_genome(data), 4, "high_only")
    assert result["n_detected_a"] == 1 and len(pairs) == 1 and result["n_omitted_moderate_calls"] == 1


def test_funnel_uses_covariates_without_reading_conservation_and_is_deterministic(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Conservation outcomes must not be calculated in this audit")
    monkeypatch.setattr(neighborhoods, "jaccard", forbidden)
    monkeypatch.setattr(neighborhoods, "signed_order", forbidden)
    data = fixture()
    tables = amendment.analyze(cfg(), *data)
    assert tables == amendment.analyze(cfg(), *data)
    controls = tables["control_availability"]
    assert len(controls) == 2
    assert all(r["n_single_copy_homologs"] == 100 and r["n_nonlysis"] == 97 for r in controls)
    assert all(r["n_position"] == r["n_length"] == r["n_prevalence"] == r["n_selected"] == 19 for r in controls)
    for row in tables["pair_feasibility"]:
        assert row["first_failed_gate"] == ("a_resolved" if row["eligibility_rule"] == "primary_unique_pair" else "eligible")
    final = [r for r in tables["availability_funnel"] if r["stage"] == "minimum_matched_controls" and r["eligibility_rule"] == "all_detected_components"]
    assert all(r["n_pairs"] == r["n_cluster_labels"] == 1 for r in final)


def test_boundary_position_length_prevalence_funnel_matches_existing_selector():
    genomes, cds, evidence, _, memberships = fixture()
    genes, windows = [], {}
    for genome in genomes:
        acc = genome["accession"]
        rows = [{**g, "family": f"SYNTH_F{g['cds_index']}"} for g in cds if g["accession"] == acc]
        genes.append(rows)
        windows.update({g["protein_id"]: neighborhoods.neighborhood(rows, g["protein_id"], 5) for g in rows})
    genes[0][40]["translation"] = "M" * 10
    prevalence = {r["family_id"]: float(r["prevalence"]) for r in memberships}
    prevalence["SYNTH_F41"] = 0.1
    excluded = {f"SYNTH_F{i}" for i in (5, 10, 50)}
    counts = amendment.control_availability(*genes, genes[0][50], genes[1][50], windows, prevalence, excluded, *genomes)
    selected, n_eligible = neighborhoods.matched_controls(*genes, genes[0][50], genes[1][50], windows, prevalence, excluded, *genomes)
    assert counts["n_position"] == 19 and counts["n_length"] == 18 and counts["n_prevalence"] == 17
    assert counts["n_prevalence"] == n_eligible and counts["n_selected"] == len(selected)
    assert counts["eligible_control_families"].split(";") == [r["control_family"] for r in selected]


def test_near_identical_pair_does_not_increase_replication():
    data = fixture()
    data[0][1]["similarity_group"] = data[0][0]["similarity_group"]
    tables = amendment.analyze(cfg(), *data)
    assert all(r["first_failed_gate"] == "not_near_identical" for r in tables["pair_feasibility"] if r["eligibility_rule"] == "all_detected_components")
    assert all(r["n_pairs"] == 0 for r in tables["availability_funnel"] if r["stage"] == "minimum_matched_controls")


def test_partial_component_does_not_rewrite_original_primary_rule():
    data = fixture()
    data[1][5]["partial"] = data[1][105]["partial"] = "true"
    tables = amendment.analyze(cfg(), *data)
    assert all(r["first_failed_gate"] == "eligible" for r in tables["pair_feasibility"] if r["eligibility_rule"] == "primary_unique_pair")
    assert all(r["first_failed_gate"] == "a_resolved" for r in tables["pair_feasibility"] if r["eligibility_rule"] == "all_detected_components")


def test_isolated_io_preserves_inputs_and_rejects_changed_freeze(tmp_path):
    root, out = tmp_path / "input", tmp_path / "amendment"
    root.mkdir()
    paths = ["data/processed/genomes.tsv", "data/processed/cds.tsv", "data/processed/lysis_evidence.tsv",
             "data/processed/architectures.tsv", "results/tables/family_memberships.tsv"]
    for path, rows in zip(paths, fixture()):
        write_tsv(root / path, rows)
    configuration = root / "config/analysis.yaml"
    configuration.parent.mkdir()
    configuration.write_text(yaml.safe_dump(cfg()))
    before = {str(p): sha256(p) for p in root.rglob("*") if p.is_file()}
    report = amendment.run(root, out)
    assert report["outcome_scores_read_or_computed"] is False
    assert report["two_sided_exact_minimum_p_with_five_nonzero_clusters"] == 0.0625
    assert before == {str(p): sha256(p) for p in root.rglob("*") if p.is_file()}
    assert amendment.run(root, out)["outputs"] == report["outputs"]
    with pytest.raises(ValueError, match="protected"):
        amendment.run(root, root / "results/tables")
    configuration.write_text(yaml.safe_dump({**cfg(), "seed": 1}))
    with pytest.raises(ValueError, match="different inputs"):
        amendment.run(root, out)
