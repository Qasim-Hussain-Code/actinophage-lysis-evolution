#!/usr/bin/env python3
"""Descriptive post-score annotation; profile candidates receive no function call."""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))
from actinolysis import architecture
from actinolysis.amendment_feasibility import independent_evidence, truth
from actinolysis.common import NA, read_tsv, sha256, write_json, write_tsv

FIELDS = ["phage", "accession", "host_genus", "cluster", "candidate_protein_id", "candidate_protein_accession", "candidate_gene",
    "candidate_cds_index", "candidate_start", "candidate_end", "candidate_strand", "candidate_partial", "candidate_product",
    "candidate_note", "passing_models", "existing_b_confidence", "n_supported_complete_a", "supported_partial_a_ids",
    "a_protein_id", "a_protein_accession", "a_gene", "a_product", "a_confidence", "a_reference_ids", "a_reference_evidence_types", "a_evidence_urls",
    "a_start", "a_end", "a_strand", "intervening_cds", "gap_bp", "relative_strands", "distance_status", "measurements_json", "interpretation"]


def annotate(candidates, cds, genomes, evidence):
    by_id = {g["protein_id"]: g for g in cds}
    gm = {g["accession"]: g for g in genomes}
    genes, supported_a, partial_a, b_confidence, models = defaultdict(list), defaultdict(dict), defaultdict(set), {}, defaultdict(set)
    for gene in cds:
        genes[gene["accession"]].append(gene)
    for row in candidates:
        if truth(row["passes_candidate_cuts"]) and not truth(row["exact_training_alias"]):
            if row["query_id"] not in by_id:
                raise ValueError("Candidate not present in frozen CDS table")
            models[row["query_id"]].add(row["model"])
    rank = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0}
    for row in evidence:
        if row["protein_id"] not in by_id:
            raise ValueError("Existing evidence references absent CDS")
        gene = by_id[row["protein_id"]]
        if row["category"] == "lysin_b" and rank[row["confidence"]] > rank.get(b_confidence.get(row["protein_id"]), -1):
            b_confidence[row["protein_id"]] = row["confidence"]
        if row["category"] == "lysin_a" and independent_evidence(row, "high_moderate")[0]:
            if truth(gene.get("partial", "false")):
                partial_a[gene["accession"]].add(gene["protein_id"])
            else:
                supported_a[gene["accession"]].setdefault(gene["protein_id"], []).append(row)
    rows = []
    for pid in sorted(models):
        candidate = by_id[pid]; accession = candidate["accession"]; genome = gm[accession]
        a_ids = sorted(supported_a[accession], key=lambda i: (int(by_id[i]["cds_index"]), i))
        for a_id in a_ids or [None]:
            row = {key: "" for key in FIELDS}
            row.update({k: genome.get(k, NA) for k in ("phage", "accession", "host_genus", "cluster")})
            row.update({"candidate_" + key: candidate.get(key, "") for key in
                ("protein_id", "protein_accession", "gene", "cds_index", "start", "end", "strand", "partial", "product", "note")})
            row.update(passing_models=";".join(sorted(models[pid])), existing_b_confidence=b_confidence.get(pid, "not_detected"),
                n_supported_complete_a=len(a_ids), supported_partial_a_ids=";".join(sorted(partial_a[accession])),
                interpretation="Unvalidated profile candidate; spacing does not establish a functional lysis architecture or gene absence.")
            if a_id is None:
                row.update(a_protein_id=NA, distance_status="no_complete_independently_supported_a_detected")
            else:
                a = by_id[a_id]; source = supported_a[accession][a_id]
                row.update({"a_" + key: a.get(key, "") for key in ("protein_id", "protein_accession", "gene", "product", "start", "end", "strand")})
                row.update(a_confidence=max((r["confidence"] for r in source), key=rank.get),
                    a_reference_ids=";".join(sorted({r["reference_id"] for r in source})),
                    a_reference_evidence_types=";".join(sorted({r.get("reference_evidence_type", NA) for r in source})),
                    a_evidence_urls=";".join(sorted({r.get("evidence_url") or r.get("reference_url") for r in source})))
                if a_id == pid:
                    row["distance_status"] = "candidate_already_supported_as_a_conflicting_function_requires_review"
                else:
                    metrics = architecture.pair_metrics(a, candidate, genes[accession], int(genome["length"]), genome["topology"])
                    row.update({key: metrics[key] for key in ("intervening_cds", "gap_bp", "relative_strands")})
                    row.update(distance_status="observed_coordinate_distance" if metrics["intervening_cds"] != NA else "ambiguous_coordinate_distance",
                               measurements_json=json.dumps(metrics, sort_keys=True, separators=(",", ":")))
            rows.append(row)
    return rows


def run(root, diagnostic, out):
    root, diagnostic, out = (Path(p).resolve() for p in (root, diagnostic, out))
    protected = [root / p for p in ("data", "references", "src", "scripts", "figures", "results/tables", "results/reports")]
    if out == root or root.is_relative_to(out) or any(out.is_relative_to(p) for p in protected) or out.is_relative_to(diagnostic) or diagnostic.is_relative_to(out):
        raise ValueError("Use a separate review output directory")
    report = json.loads((diagnostic / "profile_report.json").read_text())
    for relative, expected in report["artifacts"].items():
        if sha256(diagnostic / relative) != expected:
            raise ValueError("Completed HMM diagnostic artifact changed")
    paths = {"cds": root / "data/processed/cds.tsv", "genomes": root / "data/processed/genomes.tsv",
        "evidence": root / "data/processed/lysis_evidence.tsv", "input_manifest": root / "data/manifests/pilot_inputs.tsv",
        "stage_receipts": root / "results/reports/pilot_stage_receipts.json", "candidates": diagnostic / "cohort_candidates.tsv",
        "profile_report": diagnostic / "profile_report.json"}
    checksums = {k: sha256(p) for k, p in paths.items()}
    receipt = json.loads(paths["stage_receipts"].read_text())
    for stage, key, relative in [("1", "input_manifest", "data/manifests/pilot_inputs.tsv"), ("2", "cds", "data/processed/cds.tsv"),
                                ("2", "genomes", "data/processed/genomes.tsv"), ("4", "evidence", "data/processed/lysis_evidence.tsv")]:
        if receipt[stage]["outputs"].get(relative) != checksums[key]:
            raise ValueError("Frozen cohort metadata or existing evidence changed")
    frozen = json.loads((diagnostic / "frozen_inputs.json").read_text())
    if any(checksums[k] != frozen["input_sha256"][k] for k in ("cds", "stage_receipts")):
        raise ValueError("Review inputs differ from the completed HMM freeze")
    rows = annotate(read_tsv(paths["candidates"]), read_tsv(paths["cds"]), read_tsv(paths["genomes"]), read_tsv(paths["evidence"]))
    if checksums != {k: sha256(p) for k, p in paths.items()}:
        raise ValueError("An input changed during descriptive review")
    write_tsv(out / "candidate_review.tsv", rows, FIELDS)
    write_json(out / "review_provenance.json", {"descriptive_only": True, "primary_calls_modified": False,
        "input_sha256": checksums, "review_sha256": sha256(out / "candidate_review.tsv"),
        "script_sha256": sha256(__file__), "geometry_code_sha256": sha256(architecture.__file__),
        "distinct_non_training_candidate_cds": len({r["candidate_protein_id"] for r in rows}), "component_rows": len(rows)})
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "diagnostic", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps({"review_rows": len(run(args.root, args.diagnostic, args.out)), "descriptive_only": True}))
