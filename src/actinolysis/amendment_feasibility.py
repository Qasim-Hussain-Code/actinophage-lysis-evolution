"""Isolated component classifications and covariate-only availability audits.

No neighborhood-conservation score or inferential result is read or computed.
The original primary classifications and scientific output files are untouched.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
import math
from pathlib import Path
from urllib.parse import urlparse

from . import architecture, neighborhoods
from .common import NA, load_config, read_tsv, sha256, utcnow, write_json, write_tsv

SCOPES = ("high_moderate", "high_only")
RULES = ("primary_unique_pair", "all_detected_components")
CLASS_FIELDS = ["accession", "phage", "cluster", "confidence_set", "architecture_threshold",
    "primary_architecture", "primary_scope_architecture", "secondary_classification", "status", "reason",
    "n_detected_a", "n_detected_b", "detected_a", "detected_b", "partial_supported_ids", "conflicting_ids",
    "min_intervening_cds", "max_intervening_cds", "n_component_pairs", "n_omitted_moderate_calls"]
CALL_FIELDS = ["accession", "phage", "confidence_set", "protein_id", "category", "confidence", "complete",
    "start", "end", "strand", "reference_ids", "evidence_urls", "evidence_json"]
PAIR_DETAIL_FIELDS = ["accession", "phage", "confidence_set", "lysin_a", "lysin_b", "a_confidence", "b_confidence",
    "a_start", "a_end", "a_strand", "a_location_parts", "b_start", "b_end", "b_strand", "b_location_parts",
    "intervening_cds", "gap_bp", "distance_status", "distance_reason", "measurements_json"]
CONTROL_FIELDS = ["pair_id", "cluster", "family_identity", "window", "confidence_set", "lysin_b_a", "lysin_b_b",
    "complete_flanks", "left_a", "right_a", "left_b", "right_b", "n_single_copy_homologs", "n_nonlysis",
    "n_oriented", "n_boundary", "n_position", "n_length", "n_prevalence", "n_selected", "eligible_control_families"]
GATES = ["all_within_cluster_pairs", "supported_unique_b", "homologous_single_copy_b", "a_resolved",
    "separated", "not_near_identical", "background_single_copy", "background_nonlysis", "background_oriented",
    "boundary_controls", "position_controls", "length_controls", "prevalence_controls", "minimum_matched_controls"]
BLOCK = ["family_identity", "window", "confidence_set", "architecture_threshold", "eligibility_rule"]
FEASIBILITY_FIELDS = [*BLOCK, "pair_id", "cluster", "accession_a", "accession_b", "phage_a", "phage_b", *GATES,
    "complete_flanks", "n_controls", "first_failed_gate", "secondary_a", "secondary_b"]
FUNNEL_FIELDS = [*BLOCK, "stage", "diagnostic_branch", "n_pairs", "n_cluster_labels", "n_cluster_labels_at_least_two_pairs"]


def truth(value):
    return str(value).strip().lower() in {"true", "1", "yes"}


def independent_evidence(row, confidence_set):
    """Return a transparent acceptance flag/reason for non-contextual calls."""
    accepted = {"high"} if confidence_set == "high_only" else {"high", "moderate"}
    if row.get("confidence") not in accepted:
        return False, "confidence_outside_scope"
    reference = str(row.get("reference_id", ""))
    if reference.lower() in neighborhoods.MISSING:
        return False, "missing_reference_identifier"
    source = str(row.get("source", "")).lower()
    if truth(row.get("context_used", "false")) or any(x in source for x in ("context", "synteny", "neighborhood")):
        return False, "contextual_evidence"
    if not any(x in source for x in ("experimental", "homology", "alignment", "hmm")):
        return False, "source_not_independent_functional_evidence"
    url = row.get("evidence_url") or row.get("reference_url") or ""
    parsed = urlparse(str(url))
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        return False, "missing_traceable_evidence_url"
    return True, "accepted"


def component_classification(genome, genes, evidence, primary, threshold, confidence_set="high_moderate"):
    """Retain all supported A/B proteins, classify their entire distance set.

    Partial supported calls prevent an all-components claim. All distinct
    complete A/B pairs are exported even when multiple B calls preclude a unique
    secondary anchor. A conflicting same-CDS A/B assignment is never measured as
    a valid pair. No empty A set is interpreted as vacuously distal.
    """
    if confidence_set not in SCOPES or int(threshold) < 0:
        raise ValueError("Invalid component confidence set or threshold")
    by_id = {gene["protein_id"]: gene for gene in genes}
    selected, confidence, omitted = defaultdict(list), {}, set()
    for row in evidence:
        if row.get("category") not in {"lysin_a", "lysin_b"}:
            continue
        if row["protein_id"] not in by_id:
            raise ValueError("Component evidence references a protein outside its genome")
        if confidence_set == "high_only" and row.get("confidence") == "moderate" and independent_evidence(row, "high_moderate")[0]:
            omitted.add((row["protein_id"], row["category"]))
        if independent_evidence(row, confidence_set)[0]:
            key = (row["protein_id"], row["category"])
            selected[key].append(row)
            if row["confidence"] == "high" or key not in confidence:
                confidence[key] = row["confidence"]
    omitted -= set(selected)
    partial = sorted({pid for pid, _ in selected if truth(by_id[pid].get("partial", "false"))})
    cats = defaultdict(set)
    for pid, category in selected:
        cats[pid].add(category)
    conflicts = sorted(pid for pid, categories in cats.items() if len(categories) > 1)
    components, calls = {}, []
    for category in ("lysin_a", "lysin_b"):
        components[category] = sorted([by_id[pid] for pid, cat in selected if cat == category and pid not in partial],
                                      key=lambda g: (int(g["cds_index"]), g["protein_id"]))
    for (pid, category), rows in sorted(selected.items()):
        gene = by_id[pid]
        calls.append({"accession": genome["accession"], "phage": genome.get("phage", NA), "confidence_set": confidence_set,
            "protein_id": pid, "category": category, "confidence": confidence[(pid, category)],
            "complete": str(pid not in partial).lower(), "start": gene["start"], "end": gene["end"], "strand": gene["strand"],
            "reference_ids": ";".join(sorted({r["reference_id"] for r in rows})),
            "evidence_urls": ";".join(sorted({r.get("evidence_url") or r.get("reference_url") for r in rows})),
            "evidence_json": json.dumps(rows, sort_keys=True, separators=(",", ":"))})
    details = []
    for a in components["lysin_a"]:
        for b in components["lysin_b"]:
            pid_a, pid_b = a["protein_id"], b["protein_id"]
            metrics = {} if pid_a == pid_b else architecture.pair_metrics(a, b, genes, int(genome["length"]), genome.get("topology", "linear"))
            valid = bool(metrics) and str(metrics.get("intervening_cds", NA)) != NA
            row = {"accession": genome["accession"], "phage": genome.get("phage", NA), "confidence_set": confidence_set,
                "lysin_a": pid_a, "lysin_b": pid_b, "a_confidence": confidence[(pid_a, "lysin_a")],
                "b_confidence": confidence[(pid_b, "lysin_b")], "intervening_cds": metrics.get("intervening_cds", ""),
                "gap_bp": metrics.get("gap_bp", ""), "distance_status": "resolved" if valid else "unresolved",
                "distance_reason": "all_component_pairs_retained" if valid else "conflicting_same_protein" if pid_a == pid_b else "ambiguous_gene_order_distance",
                "measurements_json": json.dumps(metrics, sort_keys=True, separators=(",", ":"))}
            for prefix, gene in (("a", a), ("b", b)):
                row.update({f"{prefix}_{key}": gene.get(key, "") for key in ("start", "end", "strand", "location_parts")})
            details.append(row)
    reasons = []
    if partial:
        reasons.append("partial_independently_supported_component")
    if conflicts:
        reasons.append("conflicting_a_b_assignments")
    if not components["lysin_a"]:
        reasons.append("no_complete_independently_supported_a")
    if len(components["lysin_b"]) != 1:
        reasons.append("requires_one_complete_independently_supported_b")
    if any(row["distance_status"] != "resolved" for row in details):
        reasons.append("component_distance_unresolved")
    distances = [int(row["intervening_cds"]) for row in details if row["distance_status"] == "resolved"]
    category = "unresolved"
    if not reasons:
        category = ("all_detected_components_distal" if min(distances) > threshold else
                    "all_detected_components_linked" if max(distances) <= threshold else "mixed_detected_components")
    original_scope = "high" if confidence_set == "high_only" else "high_moderate"
    scope_primary, _ = architecture.classify(genome, genes, evidence, threshold, original_scope)
    classification = {"accession": genome["accession"], "phage": genome.get("phage", NA), "cluster": genome.get("cluster", NA),
        "confidence_set": confidence_set, "architecture_threshold": threshold, "primary_architecture": primary.get("architecture", NA),
        "primary_scope_architecture": scope_primary["architecture"], "secondary_classification": category,
        "status": "resolved" if not reasons else "unresolved", "reason": ";".join(reasons) or "all_detected_components_retained",
        "n_detected_a": len(components["lysin_a"]), "n_detected_b": len(components["lysin_b"]),
        "detected_a": ";".join(g["protein_id"] for g in components["lysin_a"]),
        "detected_b": ";".join(g["protein_id"] for g in components["lysin_b"]),
        "partial_supported_ids": ";".join(partial), "conflicting_ids": ";".join(conflicts),
        "min_intervening_cds": min(distances) if distances else "", "max_intervening_cds": max(distances) if distances else "",
        "n_component_pairs": len(details), "n_omitted_moderate_calls": len(omitted)}
    return classification, calls, details


def control_availability(genes_a, genes_b, a, b, windows, prevalence, excluded, genome_a, genome_b, maximum=20):
    """Count the unchanged matching filters without computing outcome scores."""
    grouped = []
    for genes in (genes_a, genes_b):
        families = defaultdict(list)
        for gene in genes:
            families[gene["family"]].append(gene)
        grouped.append(families)
    fa, fb = grouped
    families = sorted(f for f in set(fa) & set(fb) if len(fa[f]) == len(fb[f]) == 1 and f in prevalence)
    counts = {"n_single_copy_homologs": len(families)}
    families = [f for f in families if f != a["family"] and f not in excluded]
    counts["n_nonlysis"] = len(families)
    families = [f for f in families if neighborhoods.strand_value(fa[f][0]["strand"]) and neighborhoods.strand_value(fb[f][0]["strand"])]
    counts["n_oriented"] = len(families)
    targets = windows[a["protein_id"]], windows[b["protein_id"]]
    families = [f for f in families if all(neighborhoods.completeness(windows[g[0]["protein_id"]]) == neighborhoods.completeness(t)
                                                 for g, t in zip((fa[f], fb[f]), targets))]
    counts["n_boundary"] = len(families)
    circular_a, circular_b = (str(g.get("topology", "")).lower() == "circular" for g in (genome_a, genome_b))
    bins = neighborhoods.position_bin(a, genome_a["length"], circular_a), neighborhoods.position_bin(b, genome_b["length"], circular_b)
    families = [f for f in families if (neighborhoods.position_bin(fa[f][0], genome_a["length"], circular_a),
                    neighborhoods.position_bin(fb[f][0], genome_b["length"], circular_b)) == bins]
    counts["n_position"] = len(families)
    lengths = len(a["translation"]), len(b["translation"])
    families = [f for f in families if all(t and c and .5 <= c / t <= 2 for c, t in
                         zip((len(fa[f][0]["translation"]), len(fb[f][0]["translation"])), lengths))]
    counts["n_length"] = len(families)
    families = [f for f in families if abs(prevalence[f] - prevalence[a["family"]]) <= .20 + 1e-12]
    counts["n_prevalence"] = len(families)
    rank = lambda f: (sum(abs(math.log2(c / t)) for c, t in
                           zip((len(fa[f][0]["translation"]), len(fb[f][0]["translation"])), lengths))
                      + abs(prevalence[f] - prevalence[a["family"]]), f, fa[f][0]["protein_id"], fb[f][0]["protein_id"])
    families.sort(key=rank)
    counts["n_selected"] = min(maximum, len(families))
    counts["eligible_control_families"] = ";".join(families)
    return counts


def analyze(cfg, genomes, cds, evidence, primary_architectures, memberships):
    """Return isolated tables; this function does not read or write any files."""
    gm = {g["accession"]: g for g in genomes}
    cm = {g["protein_id"]: g for g in cds}
    pm = {g["accession"]: g for g in primary_architectures}
    if len(gm) != len(genomes) or len(cm) != len(cds) or len(pm) != len(primary_architectures) or set(gm) != set(pm):
        raise ValueError("Unique genomes/proteins and exactly one primary architecture per genome are required")
    genes, ev, clusters = defaultdict(list), defaultdict(list), defaultdict(list)
    for g in cds:
        if g["accession"] not in gm:
            raise ValueError("CDS references an unknown genome")
        genes[g["accession"]].append(g)
    for row in evidence:
        if row["protein_id"] not in cm:
            raise ValueError("Evidence references an unknown protein")
        ev[cm[row["protein_id"]]["accession"]].append(row)
    for g in genomes:
        if g["cluster"].lower() not in neighborhoods.MISSING | {"unassigned"}:
            clusters[g["cluster"]].append(g["accession"])
    thresholds = sorted(set(int(x) for x in cfg["architecture"]["thresholds"]))
    if not thresholds:
        raise ValueError("At least one architecture threshold is required")
    minimum = max(5, int(cfg["statistics"]["min_controls"]))
    maximum = int(cfg["statistics"]["matched_controls"])
    if maximum < minimum:
        raise ValueError("Maximum matched controls is below the unchanged minimum")
    classes, calls, details, class_map = [], [], [], {}
    for accession in sorted(gm):
        for scope in SCOPES:
            for threshold in thresholds:
                result, cc, dd = component_classification(gm[accession], genes[accession], ev[accession], pm[accession], threshold, scope)
                classes.append(result)
                class_map[(accession, scope, threshold)] = result
                if threshold == thresholds[0]:
                    calls.extend(cc)
                    details.extend(dd)
    exclusions = {r["protein_id"] for r in evidence if r["category"] in neighborhoods.LYSIS}
    controls, feasibility = [], []
    for identity in sorted(set(float(x) for x in cfg["neighborhood"]["family_identities"])):
        selected = [r for r in memberships if float(r["family_identity"]) == identity]
        families, prevalence = {}, {}
        for row in selected:
            pid, family = row["protein_id"], row["family_id"]
            value = float(row["prevalence"])
            if pid not in cm or pid in families or not 0 <= value <= 1 or (family in prevalence and prevalence[family] != value):
                raise ValueError("Invalid, duplicate, or inconsistent family membership")
            families[pid], prevalence[family] = family, value
        if not families:
            raise ValueError(f"Missing family memberships for identity {identity}")
        family_genes = {acc: [{**g, "family": families.get(g["protein_id"], "UNCLUSTERED::" + g["protein_id"])} for g in rows]
                        for acc, rows in genes.items()}
        fm = {g["protein_id"]: g for rows in family_genes.values() for g in rows}
        excluded_families = {families[p] for p in exclusions if p in families}
        for window in sorted(set(int(x) for x in cfg["neighborhood"]["windows"])):
            windows = {g["protein_id"]: neighborhoods.neighborhood(rows, g["protein_id"], window,
                       gm[acc].get("topology", "").lower() == "circular") for acc, rows in family_genes.items() for g in rows}
            for cluster, accessions in sorted(clusters.items()):
                for acc_a, acc_b in combinations(sorted(accessions), 2):
                    pair_id = f"{acc_a}|{acc_b}"
                    for scope in SCOPES:
                        c_a, c_b = (class_map[(acc, scope, thresholds[0])] for acc in (acc_a, acc_b))
                        b_ids = [c["detected_b"].split(";") if c["detected_b"] else [] for c in (c_a, c_b)]
                        unique_b = all(len(ids) == 1 for ids in b_ids)
                        anchors = [fm[ids[0]] for ids in b_ids] if unique_b else []
                        homologous = unique_b and all(a["protein_id"] in families for a in anchors) and anchors[0]["family"] == anchors[1]["family"]
                        if homologous:
                            homologous = all(sum(g["family"] == a["family"] for g in family_genes[acc]) == 1
                                             and neighborhoods.strand_value(a["strand"]) != 0 for a, acc in zip(anchors, (acc_a, acc_b)))
                        same_group = str(gm[acc_a].get("similarity_group", "")).lower() not in neighborhoods.MISSING and gm[acc_a].get("similarity_group") == gm[acc_b].get("similarity_group")
                        counts, complete = None, False
                        if homologous:
                            counts = control_availability(family_genes[acc_a], family_genes[acc_b], *anchors, windows, prevalence,
                                                          excluded_families, gm[acc_a], gm[acc_b], maximum)
                            left_a, right_a = neighborhoods.completeness(windows[anchors[0]["protein_id"]])
                            left_b, right_b = neighborhoods.completeness(windows[anchors[1]["protein_id"]])
                            complete = left_a == right_a == left_b == right_b == window
                            controls.append({"pair_id": pair_id, "cluster": cluster, "family_identity": identity, "window": window,
                                "confidence_set": scope, "lysin_b_a": anchors[0]["protein_id"], "lysin_b_b": anchors[1]["protein_id"],
                                "complete_flanks": str(complete).lower(), "left_a": left_a, "right_a": right_a,
                                "left_b": left_b, "right_b": right_b, **counts})
                        for threshold in thresholds:
                            aa, bb = (class_map[(acc, scope, threshold)] for acc in (acc_a, acc_b))
                            for rule in RULES:
                                resolved = all(c["status"] == "resolved" for c in (aa, bb))
                                if rule == "primary_unique_pair":
                                    # The original rule excludes partial calls from its
                                    # component set; the secondary rule conservatively
                                    # leaves their all-components interpretation unresolved.
                                    resolved = all(c["n_detected_a"] == c["n_detected_b"] == 1 and
                                        c["primary_scope_architecture"] in {"separated", "closely_linked"} and
                                        not (set(c["conflicting_ids"].split(";")) & set((c["detected_a"] + ";" + c["detected_b"]).split(";")))
                                        for c in (aa, bb))
                                if rule == "primary_unique_pair":
                                    separated = resolved and all(c["primary_scope_architecture"] == "separated" for c in (aa, bb))
                                else:
                                    separated = resolved and all(c["secondary_classification"] == "all_detected_components_distal" for c in (aa, bb))
                                flags = [True, unique_b, bool(homologous), resolved, separated, not same_group]
                                control_counts = [counts[k] if counts else 0 for k in ("n_single_copy_homologs", "n_nonlysis", "n_oriented", "n_boundary", "n_position", "n_length", "n_prevalence", "n_selected")]
                                flags.extend(n > 0 for n in control_counts[:-1])
                                flags.append(control_counts[-1] >= minimum)
                                cumulative, passed = {}, True
                                for gate, flag in zip(GATES, flags):
                                    passed = passed and flag
                                    cumulative[gate] = passed
                                feasibility.append({**dict(zip(BLOCK, (identity, window, scope, threshold, rule))),
                                    "pair_id": pair_id, "cluster": cluster, "accession_a": acc_a, "accession_b": acc_b,
                                    "phage_a": gm[acc_a].get("phage", NA), "phage_b": gm[acc_b].get("phage", NA),
                                    **{k: str(v).lower() for k, v in cumulative.items()}, "complete_flanks": str(complete).lower(),
                                    "n_controls": counts["n_selected"] if counts else "",
                                    "first_failed_gate": next((g for g in GATES if not cumulative[g]), "eligible"),
                                    "secondary_a": aa["secondary_classification"], "secondary_b": bb["secondary_classification"]})
    grouped = defaultdict(list)
    for row in feasibility:
        grouped[tuple(row[k] for k in BLOCK)].append(row)
    funnels = []
    for key, rows in sorted(grouped.items()):
        for stage in [*GATES, "complete_flanks_after_near_identical_filter"]:
            diagnostic = stage == "complete_flanks_after_near_identical_filter"
            retained = [r for r in rows if truth(r["not_near_identical"]) and truth(r["complete_flanks"])] if diagnostic else [r for r in rows if truth(r[stage])]
            per_cluster = Counter(r["cluster"] for r in retained)
            funnels.append({**dict(zip(BLOCK, key)), "stage": stage, "diagnostic_branch": str(diagnostic).lower(),
                "n_pairs": len(retained), "n_cluster_labels": len(per_cluster),
                "n_cluster_labels_at_least_two_pairs": sum(n >= 2 for n in per_cluster.values())})
    return {"component_classifications": classes, "component_calls": calls, "component_pairs": details,
            "control_availability": controls, "pair_feasibility": feasibility, "availability_funnel": funnels}


def run(root, out, config=None):
    root, out = Path(root).resolve(), Path(out).resolve()
    config = Path(config).resolve() if config else root / "config/analysis.yaml"
    protected = [root / p for p in ("data", "references", "src", "scripts", "figures", "results/tables", "results/reports")]
    if out == root or root.is_relative_to(out) or any(out.is_relative_to(p) for p in protected):
        raise ValueError("Use a separate output directory; baseline scientific/source paths are protected")
    inputs = {"genomes": root / "data/processed/genomes.tsv", "cds": root / "data/processed/cds.tsv",
              "evidence": root / "data/processed/lysis_evidence.tsv", "architectures": root / "data/processed/architectures.tsv",
              "families": root / "results/tables/family_memberships.tsv", "configuration": config}
    checksums = {key: sha256(path) for key, path in inputs.items()}
    prior = out / "feasibility_report.json"
    if prior.exists() and json.loads(prior.read_text(encoding="utf-8"))["input_sha256"] != checksums:
        raise ValueError("This output freeze has different inputs; choose a new output directory")
    cfg = load_config(config)
    data = {key: read_tsv(path) for key, path in inputs.items() if key != "configuration"}
    tables = analyze(cfg, data["genomes"], data["cds"], data["evidence"], data["architectures"], data["families"])
    if {key: sha256(path) for key, path in inputs.items()} != checksums:
        raise ValueError("An input changed during the feasibility calculation; no secondary outputs were written")
    fields = [CLASS_FIELDS, CALL_FIELDS, PAIR_DETAIL_FIELDS, CONTROL_FIELDS, FEASIBILITY_FIELDS, FUNNEL_FIELDS]
    for (name, rows), schema in zip(tables.items(), fields):
        write_tsv(out / f"{name}.tsv", rows, schema)
    report = {"analysis": "secondary_component_and_covariate_availability_amendment", "created_utc": utcnow(),
        "input_root": str(root), "input_paths": {k: str(p) for k, p in inputs.items()}, "input_sha256": checksums,
        "code_sha256": {Path(p).name: sha256(p) for p in (__file__, architecture.__file__, neighborhoods.__file__)},
        "outcome_scores_read_or_computed": False, "n_genomes": len(data["genomes"]),
        "minimum_controls": max(5, int(cfg["statistics"]["min_controls"])),
        "primary_minimum_cluster_gate_unchanged": max(5, int(cfg["statistics"]["min_lineages"])),
        "two_sided_exact_minimum_p_with_five_nonzero_clusters": 2 / 2 ** 5,
        "minimum_nonzero_clusters_to_permit_p_below_0_05": 6,
        "limitations": ["Availability in this cohort is not a population power calculation.",
            "Named clusters are ancestry proxies, not proven independent evolutionary replicates.",
            "All detected components does not mean all biological components or experimentally demonstrated multipart enzymes.",
            "Complete flanks are a diagnostic branch, not an additional primary eligibility requirement.",
            "Exact clipped-boundary matching may intrinsically leave fewer than five distinct control anchors.",
            "High-only filtering can remove a proximal moderate component; all removed calls and component sets must be considered.",
            "Feasibility is calculated using current global family prevalence and sampling, which can change with cohort expansion."],
        "outputs": {name: {"path": str(out / f"{name}.tsv"), "rows": len(rows), "sha256": sha256(out / f"{name}.tsv")}
                    for name, rows in tables.items()}}
    write_json(prior, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Existing isolated pipeline root")
    parser.add_argument("--out", type=Path, required=True, help="Separate directory for secondary tables")
    parser.add_argument("--config", type=Path, help="Defaults to ROOT/config/analysis.yaml")
    args = parser.parse_args(argv)
    report = run(args.root, args.out, args.config)
    print(json.dumps({"output_directory": str(args.out.resolve()), "n_genomes": report["n_genomes"],
                      "outcome_scores_read_or_computed": False}, sort_keys=True))


if __name__ == "__main__":
    main()
