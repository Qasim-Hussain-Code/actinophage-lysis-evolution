"""Analytical checks and synthetic null calibration, not study observations."""
from itertools import product

import numpy as np
import pytest

from actinolysis.statistics import (annotation_summaries, architecture_summaries, bh_adjust,
                                   lineage_bootstrap, sign_flip_test, summarize_comparisons)


def cfg():
    return {"seed": 43, "neighborhood": {"primary_identity": 0.4, "primary_window": 5},
            "architecture": {"thresholds": [2, 4, 6, 10], "separation_threshold": 4},
            "statistics": {"min_controls": 5, "min_lineages": 5, "permutations": 999, "bootstrap": 500}}


def pair(cluster, delta, pair_id=None):
    return {"family_identity": 0.4, "window": 5, "architecture_threshold": 4, "confidence_set": "high_moderate", "genome_subset": "all",
            "pair_architecture": "separated", "metric": "jaccard", "pair_id": pair_id or f"{cluster}_A|{cluster}_B", "cluster": cluster,
            "n_controls": 10, "target_score": 0.5 + delta, "control_mean": 0.5, "delta": delta,
            "status": "eligible", "reason": ""}


def test_bh_matches_independent_hand_calculation_and_preserves_missing():
    # Ordered p: .01,.02,.04,.20 -> raw .04,.04,.053333,.20.
    assert bh_adjust([0.04, None, 0.01, 0.2, 0.02]) == pytest.approx([0.04 * 4 / 3, None, 0.04, 0.2, 0.04], nan_ok=True)
    assert bh_adjust([0, 1, 1]) == [0, 1, 1]
    assert bh_adjust([None, ""]) == [None, None]
    with pytest.raises(ValueError):
        bh_adjust([1.01])


def test_exact_sign_flip_probability_is_analytically_known():
    pvalue, method = sign_flip_test([1] * 6)
    assert pvalue == 2 / 64
    assert method == "exact_cluster_sign_flip"
    assert sign_flip_test([0] * 10)[0] == 1
    assert sign_flip_test([1, -1])[0] == 1


def test_exact_null_calibration_over_all_sign_configurations():
    magnitudes = np.asarray([0.11, 0.23, 0.31, 0.47, 0.59, 0.83])
    pvalues = [sign_flip_test(magnitudes * signs)[0] for signs in product([-1, 1], repeat=6)]
    # Exhausting the null sample space gives a deterministic type-I error check.
    for alpha in [0.01, 0.05, 0.10, 0.20]:
        assert np.mean(np.asarray(pvalues) <= alpha) <= alpha + 1e-12


def test_monte_carlo_and_bootstrap_are_deterministic():
    values = np.linspace(-0.5, 0.8, 24)
    assert sign_flip_test(values, 999, 19) == sign_flip_test(values, 999, 19)
    assert 0 < sign_flip_test(values, 999, 19)[0] <= 1
    assert lineage_bootstrap(values, 501, 19) == lineage_bootstrap(values, 501, 19)
    assert lineage_bootstrap([0.3], 20, 1) == (None, None)
    assert lineage_bootstrap([0.2] * 5, 20, 1) == pytest.approx((0.2, 0.2))


def test_lineages_receive_equal_weight_despite_unequal_pair_counts():
    rows = [pair("ONE", 0.4, f"A{i}|B{i}") for i in range(100)] + [pair("TWO", -0.2)]
    tests, lineages, _ = summarize_comparisons(rows, cfg())
    primary = next(r for r in tests if r["primary"] == "true")
    assert primary["effect_mean"] == pytest.approx(0.1)
    assert primary["n_pairs"] == 101
    assert primary["n_lineages"] == 2
    assert primary["p_value"] == ""
    assert len(lineages) == 2


def test_five_lineage_gate_and_leave_one_out_are_conservative():
    rows = [pair(f"C{i}", 0.2) for i in range(5)]
    tests, _, loo = summarize_comparisons(rows, cfg())
    primary = next(r for r in tests if r["primary"] == "true")
    assert primary["status"] == "conditional_inference"
    assert primary["p_value"] == 2 / 32
    assert primary["q_value"] == primary["p_value"]
    assert len(loo) == 5
    assert all(r["status"] == "descriptive" and r["p_value"] == "" for r in loo)


def test_constant_matched_null_is_not_significant():
    tests, _, _ = summarize_comparisons([pair(f"C{i}", 0) for i in range(6)], cfg())
    primary = next(r for r in tests if r["primary"] == "true")
    assert primary["effect_mean"] == 0
    assert primary["p_value"] == 1
    assert primary["ci_low"] == primary["ci_high"] == 0


def test_empty_comparisons_still_report_untestable_primary():
    tests, lineages, loo = summarize_comparisons([], cfg())
    assert len(tests) == 1 and tests[0]["primary"] == "true"
    assert tests[0]["status"] == "untestable"
    assert tests[0]["p_value"] == tests[0]["effect_mean"] == ""
    assert lineages == loo == []


def test_invalid_pair_integrity_fails_instead_of_changing_sample_size():
    bad = pair("C", 0.1)
    bad["n_controls"] = 4
    with pytest.raises(ValueError, match="too few"):
        summarize_comparisons([bad], cfg())
    bad = pair("C", 0.1)
    bad["delta"] = 0.2
    with pytest.raises(ValueError, match="disagrees"):
        summarize_comparisons([bad], cfg())
    with pytest.raises(ValueError, match="Duplicate"):
        summarize_comparisons([pair("C", 0.1), pair("C", 0.1)], cfg())


def test_architecture_thresholds_and_high_only_are_explicit():
    genomes = [{"accession": "A", "cluster": "C"}, {"accession": "B", "cluster": "C"}]
    architectures = [{"accession": "A", "architecture": "separated", "intervening_cds": 5, "confidence": "high"},
                     {"accession": "B", "architecture": "closely_linked", "intervening_cds": 3, "confidence": "moderate"}]
    results = architecture_summaries(genomes, architectures, cfg())
    primary = next(r for r in results if r["threshold"] == 4 and r["confidence_set"] == "high_moderate")
    assert primary["fraction_separated"] == 0.5
    assert primary["within_cluster_replicated"] == "true"
    high = next(r for r in results if r["threshold"] == 4 and r["confidence_set"] == "high_only")
    assert high["n_resolved"] == high["n_unresolved"] == 1
    assert high["fraction_separated"] == 1
    relaxed = next(r for r in results if r["threshold"] == 6 and r["confidence_set"] == "high_moderate")
    assert relaxed["n_separated"] == 0


def test_annotation_counts_cannot_be_called_gene_absence():
    genomes = [{"accession": "A"}, {"accession": "B"}]
    cds = [{"protein_id": "p1", "accession": "A"}, {"protein_id": "p2", "accession": "B"}]
    evidence = [{"protein_id": "p1", "category": "lysin_b", "confidence": "high"},
                {"protein_id": "p2", "category": "lysin_b", "confidence": "low"}]
    results = annotation_summaries(genomes, cds, evidence)
    assert all(r["genomes_without_high_or_moderate_assignment"] == 1 for r in results)
    assert all("does_not_establish_gene_absence" in r["interpretation"] for r in results)
