"""Conservative evolutionary summaries and prespecified lineage-block tests.

Genome pairs and matched controls are not treated as independent replicates.
Each observed cluster contributes one mean contrast. Cluster labels are only
coarse ancestry proxies, so inferential results retain an explicit independence
and sign-symmetry qualification; they cannot establish recurrent origins.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import math

import numpy as np

from .common import ROOT, read_tsv, write_json, write_tsv

MISSING = {"", "not_available", "unknown", "unresolved", "none", "na", "singleton", "unclustered"}
BLOCK = ["family_identity", "window", "architecture_threshold", "confidence_set", "genome_subset", "pair_architecture", "metric"]
TEST_FIELDS = ["hypothesis", "analysis", *BLOCK, "primary", "n_lineages", "n_pairs", "effect_mean",
               "ci_low", "ci_high", "p_value", "q_value", "correction_family", "test_method", "status", "reason"]
LINEAGE_FIELDS = [*BLOCK, "cluster", "n_pairs", "mean_delta", "mean_target", "mean_control", "min_controls"]
LOO_FIELDS = ["analysis", *BLOCK, "omitted_cluster", "n_lineages", "effect_mean", "ci_low", "ci_high",
              "p_value", "test_method", "status"]
ARCH_FIELDS = ["threshold", "confidence_set", "cluster", "n_genomes", "n_resolved", "n_separated", "n_closely_linked",
               "n_complex", "n_unresolved", "fraction_separated", "median_intervening_cds", "within_cluster_replicated"]
ANNOTATION_FIELDS = ["category", "confidence", "n_proteins", "n_genomes", "genomes_without_high_or_moderate_assignment",
                     "interpretation"]


def bh_adjust(pvalues):
    """Benjamini-Hochberg q values in original order, preserving missing values."""
    output = [None] * len(pvalues)
    indexed = []
    for index, pvalue in enumerate(pvalues):
        if pvalue is None or pvalue == "":
            continue
        value = float(pvalue)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("BH input p values must be finite and in [0,1]")
        indexed.append((value, index))
    indexed.sort()
    bound, n = 1.0, len(indexed)
    for rank in range(n, 0, -1):
        value, original_index = indexed[rank - 1]
        bound = min(bound, value * n / rank)
        output[original_index] = min(1.0, bound)
    return output


def sign_flip_test(values, permutations=9999, seed=0):
    """Two-sided test of an equal-weight mean with cluster-wide sign flips.

    Exact enumeration is used through 20 nonzero cluster contrasts. Larger
    analyses use a Monte Carlo plus-one p value, so zero p values cannot occur.
    Validity requires independent clusters and sign-symmetric null contrasts;
    this is not a consequence of matching alone.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Sign-flip values must be a nonempty finite vector")
    if int(permutations) < 1:
        raise ValueError("At least one Monte Carlo permutation is required")
    # Zero contrasts do not change any permutation statistic or its probability.
    values = values[values != 0]
    if not len(values):
        return 1.0, "exact_cluster_sign_flip"
    observed = abs(float(values.sum()))
    tolerance = max(1e-12, observed * 1e-12)
    n, exceedances = len(values), 0
    if n <= 20:
        total = 1 << n
        shifts = np.arange(n, dtype=np.uint64)
        for start in range(0, total, 2048):
            integers = np.arange(start, min(total, start + 2048), dtype=np.uint64)
            signs = 2 * ((integers[:, None] >> shifts) & 1).astype(float) - 1
            exceedances += int(np.count_nonzero(np.abs(np.sum(signs * values, axis=1)) >= observed - tolerance))
        return exceedances / total, "exact_cluster_sign_flip"
    rng = np.random.default_rng(seed)
    total = int(permutations)
    for start in range(0, total, 2048):
        signs = rng.choice([-1.0, 1.0], size=(min(2048, total - start), n))
        exceedances += int(np.count_nonzero(np.abs(np.sum(signs * values, axis=1)) >= observed - tolerance))
    return (exceedances + 1) / (total + 1), "monte_carlo_cluster_sign_flip"


def lineage_bootstrap(values, iterations=5000, seed=0):
    """Percentile interval by resampling whole lineage effects with equal weight."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError("Bootstrap values must be a nonempty finite vector")
    if iterations < 1:
        raise ValueError("At least one bootstrap iteration is required")
    if len(values) < 2:
        return None, None
    rng, bootstrap = np.random.default_rng(seed), []
    for start in range(0, int(iterations), 2048):
        indices = rng.integers(0, len(values), size=(min(2048, int(iterations) - start), len(values)))
        bootstrap.extend(values[indices].mean(axis=1).tolist())
    low, high = np.quantile(bootstrap, [0.025, 0.975])
    return float(low), float(high)


def derived_seed(seed, key):
    return (int(seed) + int(sha256("|".join(map(str, key)).encode()).hexdigest()[:8], 16)) % (2 ** 32)


def _block_key(row):
    return (float(row["family_identity"]), int(row["window"]), int(row["architecture_threshold"]), row["confidence_set"], row["genome_subset"],
            row["pair_architecture"], row["metric"])


def _primary(key, cfg):
    return key == (float(cfg["neighborhood"]["primary_identity"]), int(cfg["neighborhood"]["primary_window"]),
                   int(cfg["architecture"]["separation_threshold"]),
                   "high_moderate", "all", "separated", "jaccard")


def summarize_comparisons(pair_rows, cfg):
    """Average within clusters before estimating uncertainty or testing H2."""
    blocks, reasons, seen = defaultdict(lambda: defaultdict(list)), defaultdict(Counter), set()
    for row in pair_rows:
        key = _block_key(row)
        identifier = (key, row["pair_id"])
        if identifier in seen:
            raise ValueError("Duplicate genome pair in the same analysis block")
        seen.add(identifier)
        # Keep blocks even when none of their rows are testable.
        blocks[key]
        if row["status"] != "eligible":
            reasons[key][row["reason"]] += 1
            continue
        if row["cluster"].lower() in MISSING:
            raise ValueError("An eligible pair cannot lack a usable lineage label")
        if int(row["n_controls"]) < max(5, int(cfg["statistics"]["min_controls"])):
            raise ValueError("Eligible neighborhood pair has too few matched controls")
        expected = float(row["target_score"]) - float(row["control_mean"])
        if not math.isclose(float(row["delta"]), expected, abs_tol=1e-10):
            raise ValueError("Neighborhood pair delta disagrees with observed-minus-control score")
        if not -1 <= expected <= 1:
            raise ValueError("Neighborhood conservation contrast outside [-1,1]")
        blocks[key][row["cluster"]].append(row)
    # The primary test remains explicit when no separated genome pair exists.
    primary_key = (float(cfg["neighborhood"]["primary_identity"]), int(cfg["neighborhood"]["primary_window"]),
                   int(cfg["architecture"]["separation_threshold"]),
                   "high_moderate", "all", "separated", "jaccard")
    blocks[primary_key]
    tests, lineage_rows, loo_rows = [], [], []
    minimum = max(5, int(cfg["statistics"]["min_lineages"]))
    for key, by_cluster in sorted(blocks.items()):
        block = dict(zip(BLOCK, key))
        effects = []
        for cluster, rows in sorted(by_cluster.items()):
            effect = float(np.mean([float(row["delta"]) for row in rows]))
            effects.append((cluster, effect))
            lineage_rows.append({**block, "cluster": cluster, "n_pairs": len(rows), "mean_delta": effect,
                "mean_target": float(np.mean([float(r["target_score"]) for r in rows])),
                "mean_control": float(np.mean([float(r["control_mean"]) for r in rows])),
                "min_controls": min(int(r["n_controls"]) for r in rows)})
        values = [value for _, value in effects]
        seed = derived_seed(cfg.get("seed", 0), key)
        low, high = lineage_bootstrap(values, cfg["statistics"]["bootstrap"], seed) if values else (None, None)
        testable = len(values) >= minimum and key[5] == "separated"
        pvalue, method = sign_flip_test(values, cfg["statistics"]["permutations"], seed) if testable else (None, "none")
        if not values:
            reason = ";".join(f"{reason}:{count}" for reason, count in sorted(reasons[key].items())) or "no_within_cluster_separated_pairs"
            status = "untestable"
        elif key[5] != "separated":
            status, reason = "descriptive", "closely_linked_or_mixed_architecture_comparison_is_exploratory"
        elif not testable:
            status, reason = "descriptive", f"fewer_than_{minimum}_independent_replicated_clusters"
        else:
            status, reason = "conditional_inference", "requires_independent_clusters_and_sign_symmetric_null_contrasts"
        is_primary = _primary(key, cfg)
        tests.append({"hypothesis": "H2", "analysis": "lysin_b_minus_matched_homologous_background", **block,
                      "primary": str(is_primary).lower(), "n_lineages": len(values),
                      "n_pairs": sum(len(rows) for rows in by_cluster.values()),
                      "effect_mean": float(np.mean(values)) if values else "", "ci_low": low if low is not None else "",
                      "ci_high": high if high is not None else "", "p_value": pvalue if pvalue is not None else "",
                      "q_value": "", "correction_family": "single_prespecified_primary" if is_primary else "exploratory_H2_sensitivities",
                      "test_method": method, "status": status, "reason": reason})
        if is_primary:
            for omitted, _ in effects:
                reduced = [value for cluster, value in effects if cluster != omitted]
                reduced_seed = derived_seed(seed, (omitted,))
                lo, hi = lineage_bootstrap(reduced, cfg["statistics"]["bootstrap"], reduced_seed) if reduced else (None, None)
                p, test_method = sign_flip_test(reduced, cfg["statistics"]["permutations"], reduced_seed) if len(reduced) >= minimum else (None, "none")
                loo_rows.append({"analysis": "leave_one_cluster_out", **block, "omitted_cluster": omitted,
                    "n_lineages": len(reduced), "effect_mean": float(np.mean(reduced)) if reduced else "",
                    "ci_low": lo if lo is not None else "", "ci_high": hi if hi is not None else "",
                    "p_value": p if p is not None else "", "test_method": test_method,
                    "status": "conditional_inference" if len(reduced) >= minimum else "descriptive"})
    for correction_family in ("single_prespecified_primary", "exploratory_H2_sensitivities"):
        indices = [i for i, row in enumerate(tests) if row["correction_family"] == correction_family]
        adjusted = bh_adjust([tests[i]["p_value"] for i in indices])
        for index, qvalue in zip(indices, adjusted):
            tests[index]["q_value"] = qvalue if qvalue is not None else ""
    return tests, lineage_rows, loo_rows


def architecture_summaries(genomes, architectures, cfg, sensitivities=None):
    """Describe lineage-balanced architecture frequencies over distance cutoffs."""
    by_accession = {row["accession"]: row for row in architectures}
    sensitivity_map = {(row["accession"], int(row["threshold"]), row["scope"]): row for row in (sensitivities or [])}
    grouped = defaultdict(list)
    for genome in genomes:
        if genome["cluster"].lower() not in MISSING:
            grouped[genome["cluster"]].append(by_accession.get(genome["accession"],
                                                           {"accession": genome["accession"], "architecture": "unresolved"}))
    rows = []
    for threshold in sorted(set(int(t) for t in cfg["architecture"]["thresholds"])):
        for confidence_set in ("high_moderate", "high_only"):
            for cluster, records in sorted(grouped.items()):
                counts, distances = Counter(), []
                for record in records:
                    scope = "high" if confidence_set == "high_only" else "high_moderate"
                    record = sensitivity_map.get((record["accession"], threshold, scope), record)
                    category = record["architecture"]
                    if confidence_set == "high_only" and record.get("confidence") != "high":
                        category = "unresolved"
                    if category in {"separated", "closely_linked"}:
                        try:
                            distance = int(record["intervening_cds"])
                        except (KeyError, ValueError, TypeError):
                            category = "unresolved"
                        else:
                            if distance < 0:
                                raise ValueError("Intervening CDS cannot be negative")
                            distances.append(distance)
                            category = "separated" if distance > threshold else "closely_linked"
                    if category not in {"separated", "closely_linked", "complex", "unresolved"}:
                        raise ValueError(f"Unexpected architecture category: {category}")
                    counts[category] += 1
                resolved = counts["separated"] + counts["closely_linked"]
                rows.append({"threshold": threshold, "confidence_set": confidence_set, "cluster": cluster,
                    "n_genomes": len(records), "n_resolved": resolved, "n_separated": counts["separated"],
                    "n_closely_linked": counts["closely_linked"], "n_complex": counts["complex"],
                    "n_unresolved": counts["unresolved"],
                    "fraction_separated": counts["separated"] / resolved if resolved else "",
                    "median_intervening_cds": float(np.median(distances)) if distances else "",
                    "within_cluster_replicated": str(resolved >= 2).lower()})
    return rows


def annotation_summaries(genomes, cds, evidence):
    protein_genome = {row["protein_id"]: row["accession"] for row in cds}
    grouped, supported = defaultdict(set), defaultdict(set)
    for row in evidence:
        if row["protein_id"] not in protein_genome:
            raise ValueError("Annotation evidence references absent protein")
        grouped[(row["category"], row["confidence"])].add(row["protein_id"])
        if row["confidence"] in {"high", "moderate"}:
            supported[row["category"]].add(protein_genome[row["protein_id"]])
    rows = []
    for (category, confidence), proteins in sorted(grouped.items()):
        rows.append({"category": category, "confidence": confidence, "n_proteins": len(proteins),
            "n_genomes": len({protein_genome[p] for p in proteins}),
            "genomes_without_high_or_moderate_assignment": len(genomes) - len(supported[category]),
            "interpretation": "no_supported_assignment_does_not_establish_gene_absence_or_loss"})
    return rows


def _descriptive(hypothesis, analysis, status, reason, n_lineages=0, n_pairs=0):
    return {"hypothesis": hypothesis, "analysis": analysis, **{field: "" for field in BLOCK},
            "primary": "false", "n_lineages": n_lineages, "n_pairs": n_pairs, "effect_mean": "", "ci_low": "",
            "ci_high": "", "p_value": "", "q_value": "", "correction_family": "none", "test_method": "none",
            "status": status, "reason": reason}


def run(cfg):
    root = Path(ROOT)
    pairs = read_tsv(root / "results/tables/neighborhood_pairs.tsv")
    genomes = read_tsv(root / "data/processed/genomes.tsv")
    cds = read_tsv(root / "data/processed/cds.tsv")
    architectures = read_tsv(root / "data/processed/architectures.tsv")
    sensitivities = read_tsv(root / "results/tables/architecture_sensitivity.tsv")
    evidence = read_tsv(root / "data/processed/lysis_evidence.tsv")
    if not genomes:
        raise ValueError("Evolutionary analysis requires a nonempty validated genome cohort")
    tests, lineage_rows, loo_rows = summarize_comparisons(pairs, cfg)
    architecture_rows = architecture_summaries(genomes, architectures, cfg, sensitivities)
    primary_architecture = [row for row in architecture_rows if row["threshold"] == int(cfg["architecture"]["separation_threshold"])
                            and row["confidence_set"] == "high_moderate"]
    replicated = sum(row["within_cluster_replicated"] == "true" for row in primary_architecture)
    tests.append(_descriptive("H1", "architecture_variation_by_lineage", "descriptive",
        "lineage_balanced_frequencies_only;no_validated_phylogeny_or_exchangeable_architecture_null;"
        + ("insufficient_within_cluster_replication" if replicated < 5 else "replication_available_but_clusters_are_coarse_ancestry_proxies"),
        n_lineages=len(primary_architecture)))
    tests.append(_descriptive("H3", "candidate_recurrent_reorganization", "untestable",
        "similar_architecture_in_distinct_clusters_does_not_identify_independent_origins;requires_gene_and_species_history",
        n_lineages=sum(row["n_separated"] > 0 for row in primary_architecture)))
    tests.append(_descriptive("H4", "annotation_uncertainty", "descriptive",
        "confidence_stratification_cannot_distinguish_undetected_genes_from_true_gene_loss"))
    annotation_rows = annotation_summaries(genomes, cds, evidence)
    write_tsv(root / "results/tables/evolutionary_tests.tsv", tests, TEST_FIELDS)
    write_tsv(root / "results/tables/lineage_effects.tsv", lineage_rows, LINEAGE_FIELDS)
    write_tsv(root / "results/tables/leave_one_cluster_out.tsv", loo_rows, LOO_FIELDS)
    write_tsv(root / "results/tables/architecture_by_lineage.tsv", architecture_rows, ARCH_FIELDS)
    write_tsv(root / "results/tables/annotation_uncertainty.tsv", annotation_rows, ANNOTATION_FIELDS)
    primary = next(row for row in tests if row["primary"] == "true")
    report = {"primary_H2": primary, "n_lineages_with_architecture": len(primary_architecture),
        "n_within_cluster_replicated_architecture_lineages": replicated,
        "equal_weight_lineage_fraction_separated": float(np.mean([row["fraction_separated"] for row in primary_architecture
                if row["fraction_separated"] != ""])) if any(row["fraction_separated"] != "" for row in primary_architecture) else None,
        "n_conditional_H2_tests": sum(r["hypothesis"] == "H2" and r["status"] == "conditional_inference" for r in tests),
        "minimum_independent_lineages": max(5, int(cfg["statistics"]["min_lineages"])),
        "inference": "Pair deltas are averaged within each cluster, then clusters receive equal weight."
            " The primary test is two-sided. Bootstrap intervals resample whole cluster effects."
            " The single prespecified primary test is separate from BH correction of all exploratory H2 sensitivities.",
        "limitations": ["Cluster labels are coarse ancestry proxies, not a reconstructed phylogeny or proof of independent origins.",
            "Sign-flip p values require independent clusters and sign-symmetric contrasts under the null; neither is guaranteed by matching.",
            "Controls share genome pairs and may overlap; they determine a background mean, not additional biological sample size.",
            "Protein family, boundary, position, length and prevalence matching cannot remove unmeasured functional or regional confounding.",
            "Percentile bootstrap intervals may be unstable with few lineages and are descriptive when the inference gate fails.",
            "Sensitivity results and leave-one-cluster-out summaries must be assessed jointly; no favorable configuration is selected."]}
    write_json(root / "results/reports/evolutionary_tests.json", report)
    return report
