"""Protein families and lineage-aware comparisons of local gene order.

Family membership is sequence based. Functional evidence is used only to choose
anchors and exclude lysis-associated controls, never to define neighboring gene
families or their conservation scores.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import math
import shutil
import subprocess

from .common import ROOT, budget_check, read_tsv, write_json, write_tsv

MISSING = {"", "not_available", "unknown", "unresolved", "none", "na", "singleton", "unclustered"}
LYSIS = {"lysin_a", "lysin_b", "holin", "holin_like"}
METRICS = ("jaccard", "signed_order")
PAIR_FIELDS = ["pair_id", "accession_a", "accession_b", "cluster", "family_identity", "window",
               "architecture_threshold", "confidence_set", "genome_subset", "pair_architecture", "metric", "lysin_b_a", "lysin_b_b",
               "anchor_family", "target_score", "control_mean", "delta", "n_controls", "n_eligible_controls",
               "status", "reason"]
CONTROL_FIELDS = ["pair_id", "family_identity", "window", "confidence_set", "control_family", "anchor_a", "anchor_b",
                  "match_distance", "jaccard", "signed_order", "target_length_a", "target_length_b",
                  "control_length_a", "control_length_b", "target_prevalence", "control_prevalence"]
MEMBER_FIELDS = ["family_identity", "protein_id", "family_id", "representative", "accession", "sequence_count",
                 "genome_count", "prevalence", "representative_function", "all_hypothetical"]
FAMILY_FIELDS = ["family_identity", "family_id", "representative", "sequence_count", "genome_count", "prevalence",
                 "unique_sequence_count", "unique_sequence_fraction", "min_length_aa", "median_length_aa", "max_length_aa",
                 "representative_function", "function_annotations", "all_hypothetical", "diversity_basis"]
NEIGHBOR_FIELDS = ["family_identity", "window", "accession", "anchor", "anchor_family", "anchor_confidence",
                   "n_neighbors", "n_unclustered_neighbors", "left_count", "right_count", "complete", "neighbor_proteins",
                   "neighbor_families", "signed_offsets", "relative_strands"]


def read_fasta(path):
    """Read the small study-specific proteome, rejecting duplicated identifiers."""
    sequences, current, chunks = {}, None, []
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks)
                current, chunks = line[1:].split()[0], []
                if current in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
            else:
                if current is None:
                    raise ValueError("FASTA sequence precedes its identifier")
                chunks.append(line.upper())
    if current is not None:
        if current in sequences:
            raise ValueError(f"Duplicate FASTA identifier: {current}")
        sequences[current] = "".join(chunks)
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError("Proteome FASTA is empty or contains empty proteins")
    return sequences


def parse_clusters(path, sequences):
    """Validate exactly-once assignment and name families from sequence content.

    The family hash uses sorted (stable protein identifier, sequence SHA-256)
    pairs. Its identity is independent of MMseqs representative choice and input
    ordering. Adding a member changes the identifier. Including identifiers also
    distinguishes short identical peptides left in separate nonsignificant
    clusters; the naming scheme does not merge them after clustering.
    """
    grouped, seen = defaultdict(list), set()
    with Path(path).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Invalid MMseqs cluster row {line_number}: expected two fields")
            representative, protein = fields
            if representative not in sequences or protein not in sequences:
                raise ValueError(f"Unknown protein in MMseqs output: {fields}")
            if protein in seen:
                raise ValueError(f"Protein assigned more than once: {protein}")
            seen.add(protein)
            grouped[representative].append(protein)
    if seen != set(sequences):
        raise ValueError(f"Unassigned proteins: {sorted(set(sequences) - seen)[:10]}")
    assignments, families = {}, {}
    for representative, members in sorted(grouped.items()):
        if representative not in members:
            raise ValueError(f"MMseqs representative absent from its own family: {representative}")
        digest = "\n".join(p + "\t" + sha256(sequences[p].encode()).hexdigest() for p in sorted(members))
        family = "FAM_" + sha256(digest.encode()).hexdigest()[:20]
        if family in families:
            raise ValueError("Unexpected family content-hash collision")
        families[family] = {"representative": representative, "members": sorted(members)}
        assignments.update({p: family for p in members})
    return assignments, families


def cluster_proteins(cfg, sequences, cds, genomes, root=ROOT):
    """Run three small MMseqs2 clusters with explicit resource and search limits."""
    root = Path(root)
    params = cfg["neighborhood"]
    identities = sorted(set(float(x) for x in params["family_identities"]))
    if float(params["primary_identity"]) not in identities:
        raise ValueError("primary_identity must be included in family_identities")
    tool = cfg["tools"]["mmseqs"]
    command_prefix = [str(x) for x in tool] if isinstance(tool, list) else [str(tool)]
    log_dir = root / "results/logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    version = subprocess.run(command_prefix + ["version"], check=True, capture_output=True, text=True).stdout.strip()
    input_fasta = root / "data/processed/proteins.faa"
    input_hash = sha256(input_fasta.read_bytes()).hexdigest()
    protein_rows = {row["protein_id"]: row for row in cds}
    if len(protein_rows) != len(cds) or not set(sequences) <= set(protein_rows):
        raise ValueError("Proteome FASTA must contain unique identifiers present in the CDS table")
    if any(sequence != protein_rows[protein]["translation"].upper().removesuffix("*") for protein, sequence in sequences.items()):
        raise ValueError("Proteome FASTA and CDS translation content disagree")
    member_rows, summary_rows, assignments_by_identity, audit = [], [], {}, []
    for identity in identities:
        if not 0 < identity <= 1:
            raise ValueError("Family identity thresholds must lie in (0,1]")
        budget_check(cfg, additional_bytes=max(50_000_000, input_fasta.stat().st_size * 20), root=root)
        tag = f"{round(identity * 100):03d}"
        work = root / "data/processed/mmseqs_work" / tag
        work.mkdir(parents=True, exist_ok=True)
        prefix, temp = work / "families", work / "tmp"
        output = Path(str(prefix) + "_cluster.tsv")
        metadata_path = work / "provenance.json"
        metadata = {"input_sha256": input_hash, "mmseqs_version": version, "identity": identity,
                    "coverage": float(params["coverage"]), "evalue": float(params["evalue"]),
                    "coverage_mode": 0, "sequence_identity_mode": 0, "cluster_mode": 0,
                    "alignment_mode": 3, "kmer_length": 5, "single_step_clustering": True,
                    "threads": min(4, int(cfg.get("threads", 4))),
                    "split_memory_limit": str(cfg["resources"].get("mmseqs_memory", "500M"))}
        reusable = output.exists() and metadata_path.exists() and json.loads(metadata_path.read_text()) == metadata
        command = command_prefix + ["easy-cluster", str(input_fasta), str(prefix), str(temp),
                    "--min-seq-id", str(identity), "-c", str(params["coverage"]), "--cov-mode", "0",
                    "--seq-id-mode", "0", "--alignment-mode", "3", "--cluster-mode", "0",
                    "--single-step-clustering", "1", "-k", "5", "-e", str(params["evalue"]),
                    "--threads", str(metadata["threads"]), "--split-memory-limit", metadata["split_memory_limit"],
                    "--remove-tmp-files", "1"]
        if not reusable:
            # A previous failed workflow must not cause MMseqs to reuse stale DBs.
            if any(work.iterdir()):
                resolved = work.resolve()
                if root.resolve() not in resolved.parents or resolved == root.resolve():
                    raise RuntimeError("Refusing to clean an MMseqs directory outside the project")
                shutil.rmtree(work)
                work.mkdir(parents=True)
            with (log_dir / f"mmseqs_families_{tag}.log").open("w", encoding="utf-8") as log:
                log.write("Command arguments: " + json.dumps(command) + "\n")
                log.flush()
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, cwd=root)
        assignments, families = parse_clusters(output, sequences)
        write_json(metadata_path, metadata)
        assignments_by_identity[identity] = assignments
        for family, info in sorted(families.items()):
            members = info["members"]
            genome_count = len({protein_rows[p]["accession"] for p in members})
            products = [protein_rows[p].get("product", "not_available") for p in members]
            informative = [p for p in products if p.lower() not in MISSING and "hypothetical" not in p.lower()]
            function = protein_rows[info["representative"]].get("product", "not_available")
            lengths = sorted(len(sequences[p]) for p in members)
            middle = len(lengths) // 2
            median_length = lengths[middle] if len(lengths) % 2 else (lengths[middle - 1] + lengths[middle]) / 2
            unique_sequence_count = len({sequences[p] for p in members})
            summary_rows.append({"family_identity": identity, "family_id": family, "representative": info["representative"],
                "sequence_count": len(members), "genome_count": genome_count, "prevalence": genome_count / len(genomes),
                "unique_sequence_count": unique_sequence_count, "unique_sequence_fraction": unique_sequence_count / len(members),
                "min_length_aa": min(lengths), "median_length_aa": median_length, "max_length_aa": max(lengths),
                "representative_function": function, "function_annotations": "; ".join(sorted(set(informative))) or "not_available",
                "all_hypothetical": str(not informative).lower(),
                "diversity_basis": "exact_sequence_and_length_variation;not_alignment_identity"})
            for protein in members:
                member_rows.append({"family_identity": identity, "protein_id": protein, "family_id": family,
                                    "representative": info["representative"], "accession": protein_rows[protein]["accession"],
                                    "sequence_count": len(members), "genome_count": genome_count,
                                    "prevalence": genome_count / len(genomes), "representative_function": function,
                                    "all_hypothetical": str(not informative).lower()})
        audit.append({**metadata, "n_proteins": len(assignments), "n_families": len(families),
                      "cache_reused": reusable, "log": f"results/logs/mmseqs_families_{tag}.log"})
        budget_check(cfg, root=root)
    write_tsv(root / "results/tables/family_memberships.tsv", member_rows, MEMBER_FIELDS)
    write_tsv(root / "results/tables/family_summary.tsv", summary_rows, FAMILY_FIELDS)
    primary = float(params["primary_identity"])
    write_tsv(root / "data/processed/families.tsv",
              [row for row in member_rows if row["family_identity"] == primary], MEMBER_FIELDS)
    return assignments_by_identity, member_rows, audit


def strand_value(value):
    value = str(value)
    if value in {"1", "+", "+1"}:
        return 1
    if value in {"-1", "-"}:
        return -1
    if value in {"0", "", "not_available", "None"}:
        return 0
    raise ValueError(f"Unknown CDS strand: {value}")


def neighborhood(genes, anchor_id, window, circular=False):
    """Return distinct flanking CDS in anchor-oriented order, excluding anchor.

    Linear boundaries are clipped. Circular windows wrap only with explicit
    circular topology. In very short genomes, each other CDS occurs once, at its
    nearest signed offset (negative in anchor orientation breaks a distance tie).
    """
    if window < 1:
        raise ValueError("Neighborhood window must be positive")
    genes = sorted(genes, key=lambda row: (int(row["start"]), int(row["end"]), row["protein_id"]))
    anchor_index = next((i for i, gene in enumerate(genes) if gene["protein_id"] == anchor_id), None)
    if anchor_index is None:
        raise ValueError(f"Anchor not in its genome CDS: {anchor_id}")
    anchor_strand, selected = strand_value(genes[anchor_index]["strand"]) or 1, {}
    for offset in sorted(range(-window, window + 1), key=lambda x: (abs(x), x * anchor_strand)):
        if offset == 0:
            continue
        index = anchor_index + offset
        if circular:
            index %= len(genes)
        if not 0 <= index < len(genes) or index == anchor_index or index in selected:
            continue
        gene = genes[index]
        selected[index] = {"protein_id": gene["protein_id"], "family": gene.get("family", ""),
                           "offset": offset * anchor_strand,
                           "relative_strand": strand_value(gene["strand"]) * anchor_strand}
    return sorted(selected.values(), key=lambda row: row["offset"])


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a or b else None


def signed_order(a, b):
    """LCS overlap of family/anchor-side/relative-orientation gene tokens.

    This bounded score penalizes losses, insertions, reorderings and inversions;
    it does not call gene orientation a biochemical mechanism.
    """
    if not a and not b:
        return None
    encode = lambda rows: [(r["family"], -1 if r["offset"] < 0 else 1,
                            r["relative_strand"] if r["relative_strand"] else "unknown::" + r["protein_id"])
                           for r in rows]
    x, y = encode(a), encode(b)
    previous = [0] * (len(y) + 1)
    for left in x:
        current = [0]
        for index, right in enumerate(y, 1):
            current.append(previous[index - 1] + 1 if left == right else max(previous[index], current[-1]))
        previous = current
    return 2 * previous[-1] / (len(x) + len(y))


def scores(a, b):
    return {"jaccard": jaccard([g["family"] for g in a], [g["family"] for g in b]),
            "signed_order": signed_order(a, b)}


def completeness(rows):
    return (sum(row["offset"] < 0 for row in rows), sum(row["offset"] > 0 for row in rows))


def position_bin(gene, genome_length, circular=False):
    # Distance from the nearest terminus is invariant to whole-genome reversal.
    # An explicitly circular record has no privileged origin or terminus.
    if circular:
        return 0
    position = (int(gene["start"]) + int(gene["end"])) / (2 * int(genome_length))
    return min(4, int(max(0, min(position, 1 - position)) * 10))


def matched_controls(genes_a, genes_b, anchor_a, anchor_b, windows, prevalence,
                     excluded_families, genome_a, genome_b, max_controls=20):
    """Choose homologous anchors in the same genome pair using covariates only.

    Eligible controls are single-copy in each genome, nonlysis, matched exactly
    on oriented boundary completeness and nearest-terminus decile, within a
    factor of two on both protein lengths, and within 0.20 family prevalence.
    Sorting uses these covariates and identifiers, never conservation scores.
    """
    by_family_a, by_family_b = defaultdict(list), defaultdict(list)
    for gene in genes_a:
        by_family_a[gene["family"]].append(gene)
    for gene in genes_b:
        by_family_b[gene["family"]].append(gene)
    a, b = anchor_a, anchor_b
    target_family = a["family"]
    target_windows = (windows[a["protein_id"]], windows[b["protein_id"]])
    target_lengths = (len(a["translation"]), len(b["translation"]))
    circular_a = str(genome_a.get("topology", "")).lower() == "circular"
    circular_b = str(genome_b.get("topology", "")).lower() == "circular"
    target_bins = (position_bin(a, genome_a["length"], circular_a), position_bin(b, genome_b["length"], circular_b))
    if not all(target_lengths):
        return [], 0
    eligible = []
    for family in sorted(set(by_family_a) & set(by_family_b)):
        if family == target_family or family in excluded_families or family not in prevalence:
            continue
        if len(by_family_a[family]) != 1 or len(by_family_b[family]) != 1:
            continue
        ca, cb = by_family_a[family][0], by_family_b[family][0]
        if not strand_value(ca["strand"]) or not strand_value(cb["strand"]):
            continue
        control_windows = (windows[ca["protein_id"]], windows[cb["protein_id"]])
        if any(completeness(c) != completeness(t) for c, t in zip(control_windows, target_windows)):
            continue
        if (position_bin(ca, genome_a["length"], circular_a), position_bin(cb, genome_b["length"], circular_b)) != target_bins:
            continue
        control_lengths = (len(ca["translation"]), len(cb["translation"]))
        if not all(control_lengths) or any(not 0.5 <= c / t <= 2 for c, t in zip(control_lengths, target_lengths)):
            continue
        prevalence_distance = abs(prevalence[family] - prevalence[target_family])
        if prevalence_distance > 0.20 + 1e-12:
            continue
        distance = sum(abs(math.log2(c / t)) for c, t in zip(control_lengths, target_lengths)) + prevalence_distance
        eligible.append({"control_family": family, "anchor_a": ca["protein_id"], "anchor_b": cb["protein_id"],
                         "match_distance": distance, **scores(*control_windows),
                         "target_length_a": target_lengths[0], "target_length_b": target_lengths[1],
                         "control_length_a": control_lengths[0], "control_length_b": control_lengths[1],
                         "target_prevalence": prevalence[target_family], "control_prevalence": prevalence[family]})
    eligible.sort(key=lambda row: (row["match_distance"], row["control_family"], row["anchor_a"], row["anchor_b"]))
    return eligible[:max_controls], len(eligible)


def _truth(value):
    return str(value).lower() in {"true", "1", "yes"}


def compare_neighborhoods(cfg, genomes, cds, evidence, architectures, assignments_by_identity, memberships,
                          architecture_sensitivities=None):
    """Construct all prespecified comparison blocks, retaining untestable pairs."""
    genome_map = {r["accession"]: r for r in genomes}
    architecture_map = {r["accession"]: r for r in architectures}
    high_architecture_map = {r["accession"]: r for r in (architecture_sensitivities or [])
                             if r.get("scope") == "high" and int(r["threshold"]) == int(cfg["architecture"]["separation_threshold"])}
    thresholds = sorted(set(int(t) for t in cfg.get("architecture", {}).get("thresholds", [4])))
    cds_map = {r["protein_id"]: r for r in cds}
    confidence_rank = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0}
    b_confidence, all_lysis = {}, set()
    for row in evidence:
        protein, category = row["protein_id"], row["category"]
        if protein not in cds_map:
            raise ValueError(f"Functional evidence references absent CDS: {protein}")
        if category in LYSIS:
            all_lysis.add(protein)
        if category == "lysin_b":
            confidence = row["confidence"]
            if confidence not in confidence_rank:
                raise ValueError(f"Unknown functional confidence: {confidence}")
            if confidence_rank[confidence] > confidence_rank.get(b_confidence.get(protein), -1):
                b_confidence[protein] = confidence
    if len(genome_map) != len(genomes) or len(cds_map) != len(cds):
        raise ValueError("Duplicate genome or protein identifier")
    per_genome = defaultdict(list)
    for gene in cds:
        if gene["accession"] not in genome_map:
            raise ValueError("CDS references an absent genome")
        per_genome[gene["accession"]].append(gene)
    cluster_genomes = defaultdict(list)
    for genome in genomes:
        if genome["cluster"].lower() not in MISSING:
            cluster_genomes[genome["cluster"]].append(genome["accession"])
    pairs, controls, neighbor_rows = [], [], []
    min_controls = max(5, int(cfg["statistics"]["min_controls"]))
    max_controls = int(cfg["statistics"]["matched_controls"])
    if max_controls < min_controls:
        raise ValueError("matched_controls cannot be smaller than min_controls (at least five)")
    has_representatives = all("representative" in g for g in genomes)
    for identity, assignments in sorted(assignments_by_identity.items()):
        if not set(assignments) <= set(cds_map):
            raise ValueError("Family threshold references a protein absent from the CDS table")
        genes = {acc: [{**row, "family": assignments.get(row["protein_id"], "UNCLUSTERED::" + row["protein_id"])} for row in rows]
                 for acc, rows in per_genome.items()}
        prevalence = {row["family_id"]: float(row["prevalence"]) for row in memberships
                      if float(row["family_identity"]) == identity}
        excluded_families = {assignments[p] for p in all_lysis if p in assignments}
        for window in sorted(set(int(w) for w in cfg["neighborhood"]["windows"])):
            windows = {g["protein_id"]: neighborhood(rows, g["protein_id"], window,
                       genome_map[accession].get("topology", "").lower() == "circular")
                       for accession, rows in genes.items() for g in rows}
            for protein in sorted(b_confidence):
                rows = windows[protein]
                left, right = completeness(rows)
                neighbor_rows.append({"family_identity": identity, "window": window,
                    "accession": cds_map[protein]["accession"], "anchor": protein,
                    "anchor_family": assignments.get(protein, "UNCLUSTERED::" + protein), "anchor_confidence": b_confidence[protein],
                    "n_neighbors": len(rows), "n_unclustered_neighbors": sum(r["family"].startswith("UNCLUSTERED::") for r in rows),
                    "left_count": left, "right_count": right,
                    "complete": str(left == right == window).lower(),
                    "neighbor_proteins": ";".join(r["protein_id"] for r in rows),
                    "neighbor_families": ";".join(r["family"] for r in rows),
                    "signed_offsets": ";".join(str(r["offset"]) for r in rows),
                    "relative_strands": ";".join(str(r["relative_strand"]) for r in rows)})
            for cluster, accessions in sorted(cluster_genomes.items()):
                for acc_a, acc_b in combinations(sorted(accessions), 2):
                    ga, gb = genome_map[acc_a], genome_map[acc_b]
                    rep_eligible = has_representatives and _truth(ga["representative"]) and _truth(gb["representative"])
                    for confidence_set, accepted in (("high_moderate", {"high", "moderate"}), ("high_only", {"high"})):
                        candidates_a = [g for g in genes[acc_a] if b_confidence.get(g["protein_id"]) in accepted]
                        candidates_b = [g for g in genes[acc_b] if b_confidence.get(g["protein_id"]) in accepted]
                        pair_id = f"{acc_a}|{acc_b}"
                        reason, anchors = "", None
                        if len(candidates_a) != 1 or len(candidates_b) != 1:
                            reason = "requires_one_supported_lysin_b_per_genome"
                        else:
                            a, b = candidates_a[0], candidates_b[0]
                            anchors = a, b
                            if a["protein_id"] not in assignments or b["protein_id"] not in assignments:
                                reason = "lysin_b_anchor_excluded_from_family_clustering"
                            elif a["family"] != b["family"]:
                                reason = "lysin_b_anchors_not_in_same_homologous_family"
                            elif any(sum(g["family"] == anchor["family"] for g in genes[acc]) != 1
                                     for acc, anchor in ((acc_a, a), (acc_b, b))):
                                reason = "lysin_b_family_has_multiple_copies"
                        chosen_architectures = high_architecture_map if confidence_set == "high_only" and high_architecture_map else architecture_map
                        aa, ab = chosen_architectures.get(acc_a, {}), chosen_architectures.get(acc_b, {})
                        architecture_a, architecture_b = aa.get("architecture", "unresolved"), ab.get("architecture", "unresolved")
                        pair_architecture = architecture_a if architecture_a == architecture_b else "mixed"
                        if not reason and any(x not in {"separated", "closely_linked"} for x in (architecture_a, architecture_b)):
                            reason = "architecture_complex_or_unresolved"
                        if not reason and anchors and any(str(arch.get("lysin_b", "")) != anchor["protein_id"]
                                                              for arch, anchor in zip((aa, ab), anchors)):
                            reason = "architecture_and_supported_anchor_disagree"
                        if not reason and confidence_set == "high_only" and any(arch.get("confidence") != "high" for arch in (aa, ab)):
                            reason = "architecture_not_high_confidence"
                        if not reason and anchors and any(not strand_value(anchor["strand"]) for anchor in anchors):
                            reason = "lysin_b_anchor_strand_unresolved"
                        if not reason and str(ga.get("similarity_group", "")).lower() not in MISSING and ga.get("similarity_group") == gb.get("similarity_group"):
                            reason = "same_near_identical_similarity_group"
                        selected, n_eligible, target = [], 0, {m: None for m in METRICS}
                        if not reason:
                            a, b = anchors
                            target = scores(windows[a["protein_id"]], windows[b["protein_id"]])
                            if any(target[m] is None for m in METRICS):
                                reason = "empty_anchor_neighborhoods"
                            else:
                                selected, n_eligible = matched_controls(genes[acc_a], genes[acc_b], a, b, windows,
                                    prevalence, excluded_families, ga, gb, max_controls)
                                if len(selected) < min_controls:
                                    reason = "fewer_than_minimum_matched_homologous_controls"
                                # Controls do not depend on confidence or representative filtering.
                                controls.extend({"pair_id": pair_id, "family_identity": identity, "window": window,
                                                 "confidence_set": confidence_set, **row} for row in selected)
                        for subset in ("all", "representatives"):
                            subset_reason = reason
                            if subset == "representatives" and not has_representatives:
                                subset_reason = "qc_representative_metadata_unavailable"
                            elif subset == "representatives" and not rep_eligible:
                                subset_reason = "not_both_qc_representatives"
                            for threshold in thresholds:
                                classified, threshold_reason = [], subset_reason
                                for arch in (aa, ab):
                                    category = arch.get("architecture", "unresolved")
                                    if category in {"separated", "closely_linked"}:
                                        try:
                                            intervening = int(arch["intervening_cds"])
                                        except (KeyError, TypeError, ValueError):
                                            category = "unresolved"
                                            threshold_reason = threshold_reason or "missing_intervening_cds_for_threshold"
                                        else:
                                            category = "separated" if intervening > threshold else "closely_linked"
                                    classified.append(category)
                                threshold_architecture = classified[0] if classified[0] == classified[1] else "mixed"
                                for metric in METRICS:
                                    control_values = [row[metric] for row in selected if row[metric] is not None]
                                    mean = sum(control_values) / len(control_values) if control_values else None
                                    eligible = not threshold_reason
                                    pairs.append({"pair_id": pair_id, "accession_a": acc_a, "accession_b": acc_b,
                                        "cluster": cluster, "family_identity": identity, "window": window,
                                        "architecture_threshold": threshold, "confidence_set": confidence_set, "genome_subset": subset,
                                        "pair_architecture": threshold_architecture, "metric": metric,
                                        "lysin_b_a": anchors[0]["protein_id"] if anchors else "",
                                        "lysin_b_b": anchors[1]["protein_id"] if anchors else "",
                                        "anchor_family": anchors[0]["family"] if anchors else "",
                                        "target_score": target[metric] if target[metric] is not None else "",
                                        "control_mean": mean if mean is not None else "",
                                        "delta": target[metric] - mean if eligible else "",
                                        "n_controls": len(selected), "n_eligible_controls": n_eligible,
                                        "status": "eligible" if eligible else "untestable", "reason": threshold_reason})
    return pairs, controls, neighbor_rows


def run(cfg):
    root = Path(ROOT)
    genomes = read_tsv(root / "data/processed/genomes.tsv")
    cds = read_tsv(root / "data/processed/cds.tsv")
    if not genomes or not cds:
        raise ValueError("Neighborhood analysis requires nonempty validated genomes and CDS")
    evidence = read_tsv(root / "data/processed/lysis_evidence.tsv")
    architectures = read_tsv(root / "data/processed/architectures.tsv")
    architecture_sensitivities = read_tsv(root / "results/tables/architecture_sensitivity.tsv")
    sequences = read_fasta(root / "data/processed/proteins.faa")
    assignments, members, clustering_audit = cluster_proteins(cfg, sequences, cds, genomes, root)
    pair_rows, control_rows, neighbor_rows = compare_neighborhoods(cfg, genomes, cds, evidence, architectures, assignments, members,
                                                                architecture_sensitivities)
    write_tsv(root / "results/tables/neighborhood_pairs.tsv", pair_rows, PAIR_FIELDS)
    write_tsv(root / "results/tables/neighborhood_controls.tsv", control_rows, CONTROL_FIELDS)
    write_tsv(root / "results/tables/neighborhoods.tsv", neighbor_rows, NEIGHBOR_FIELDS)
    report = {"family_clustering": clustering_audit, "n_cds": len(cds), "n_clustered_proteins": len(sequences),
              "n_excluded_cds": len(cds) - len(sequences), "n_comparison_rows": len(pair_rows),
              "n_eligible_rows": sum(r["status"] == "eligible" for r in pair_rows),
              "untestable_reasons": dict(Counter(r["reason"] for r in pair_rows if r["reason"])),
              "control_matching": {"same_genome_pair": True, "homologous_single_copy_anchors": True,
                  "exact_boundary_completeness": True, "exact_nearest_terminus_position_bin": True,
                  "protein_length_ratio": [0.5, 2.0], "max_prevalence_difference": 0.20,
                  "minimum_controls": max(5, cfg["statistics"]["min_controls"]),
                  "selection_uses_conservation_score": False},
              "limitations": ["Homologous families do not establish identical biochemical function.",
                  "Matched background comparisons are observational; lysis and control anchors need not be exchangeable.",
                  "Genomic position, annotation, unmeasured family age and regional constraints may leave residual confounding.",
                  "Position matching uses nearest-terminus deciles only for linear records; circular origins are arbitrary.",
                  "Single-step MMseqs clustering uses exact aligned residue identity; a k=5 prefilter reduces its index memory."
                  " Candidate detection is not an exhaustive pairwise homology search, and thresholds apply to representative-member edges."
                  " Family membership does not guarantee the same minimum identity between every pair of members.",
                  "The MMseqs split-memory limit bounds prefilter structures, not total process resident memory."]}
    write_json(root / "results/reports/neighborhoods.json", report)
    return report
