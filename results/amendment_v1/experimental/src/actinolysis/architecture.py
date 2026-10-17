"""Measure lysis architecture without inferring gene absence or genome adjacency.

Distances use original zero-based, half-open CDS coordinates. The principal
intervening count measures CDS order, including genes overlapping the physical
gap boundaries: absolute CDS-index difference minus one on linear records.
The nonnegative edge gap, overlapping bases and the supplemental number of CDS
entirely inside that gap are separate measurements. Only explicitly circular
records may wrap across the origin. Their order is defined by CDS midpoint loci
along the selected shortest edge path, avoiding envelope-start artifacts for
origin-spanning joins. A circular midpoint is the midpoint of the shortest arc
containing the CDS parts, with a deterministic tie rule. Compound locations use
``location_parts`` as a JSON list of [start, end] intervals; an envelope alone
must not represent an origin-spanning CDS.
"""
from __future__ import annotations

from collections import defaultdict
import json
import re

from .common import NA, ROOT, read_tsv, write_tsv

CATEGORIES = {"lysin_a", "lysin_b", "holin", "holin_like"}
CONFIDENCE = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0}
SCOPES = ("high", "high_moderate", "annotation")
PAIR_FIELDS = [
    "accession", "phage", "host_genus", "cluster", "scope", "lysin_a", "lysin_b",
    "confidence", "a_start", "a_end", "a_strand", "a_location_parts",
    "b_start", "b_end", "b_strand", "b_location_parts", "a_normalized_position",
    "b_normalized_position", "gap_bp", "overlap_bp", "absolute_midpoint_distance_bp",
    "midpoint_distance_bp", "intervening_cds", "intervening_cds_min",
    "intervening_cds_max", "complete_gap_cds", "complete_gap_cds_min", "complete_gap_cds_max",
    "gene_order_basis", "relative_strands", "orientation", "distance_basis",
    "shortest_path_count", "path_ambiguous", "nearby_holin_count",
    "nearby_holin_ids", "other_lysis_count", "structural_distance_cds",
    "nearest_structural_protein", "structural_basis",
]
ARCH_FIELDS = [
    "accession", "phage", "host_genus", "cluster", "subcluster", "scope",
    "architecture", "classification_reason", "threshold", "confidence",
    "lysin_a_count", "lysin_b_count", "holin_count", "holin_like_count",
    "lysis_protein_count", "partial_lysis_count", "conflicting_protein_ids",
] + [f for f in PAIR_FIELDS if f not in {"accession", "phage", "host_genus", "cluster", "scope", "confidence"}]


def _integer(value, label):
    try:
        return int(value)
    except (ValueError, TypeError) as error:
        raise ValueError(f"Invalid integer {label}: {value!r}") from error


def _truth(value):
    return str(value).strip().lower() in {"1", "true", "yes"}


def _parts(gene, length):
    """Validate and merge a CDS footprint; retain original coordinates elsewhere."""
    start = _integer(gene["start"], "CDS start")
    end = _integer(gene["end"], "CDS end")
    if not 0 <= start < end <= length:
        raise ValueError(f"Invalid CDS coordinates for {gene['protein_id']}: {start}:{end}/{length}")
    raw = gene.get("location_parts", "")
    try:
        intervals = json.loads(raw) if isinstance(raw, str) and raw not in {"", NA} else raw
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid location_parts for {gene['protein_id']}") from error
    intervals = intervals or [[start, end]]
    if not isinstance(intervals, list) or not all(isinstance(x, (list, tuple)) and len(x) == 2 for x in intervals):
        raise ValueError(f"Invalid location_parts for {gene['protein_id']}")
    intervals = sorted((_integer(s, "part start"), _integer(e, "part end")) for s, e in intervals)
    if any(not 0 <= s < e <= length for s, e in intervals):
        raise ValueError(f"Invalid compound CDS coordinates for {gene['protein_id']}")
    if intervals[0][0] != start or max(e for _, e in intervals) != end:
        raise ValueError(f"CDS envelope does not match location_parts for {gene['protein_id']}")
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(e, merged[-1][1]))
        else:
            merged.append((s, e))
    return merged


def _circular(topology):
    # Circular permutation alone does not establish physical end adjacency.
    return str(topology).strip().lower() == "circular"


def _midpoint(parts, length, circular):
    if not circular or len(parts) == 1:
        return (parts[0][0] + parts[-1][1]) / 2
    gaps = [(parts[(i + 1) % len(parts)][0] - end) % length for i, (_, end) in enumerate(parts)]
    largest = max(gaps)
    # The smallest start after a tied largest gap fixes an otherwise ambiguous
    # midpoint convention; distances still preserve path ties independently.
    start = min(parts[(i + 1) % len(parts)][0] for i, gap in enumerate(gaps) if gap == largest)
    return (start + (length - largest) / 2) % length


def _inside(parts, start, gap, length, circular):
    """Whether every covered base lies in the selected directed gap."""
    if gap <= 0:
        return False
    if circular:
        return all((s - start) % length + (e - s) <= gap for s, e in parts)
    return all(start <= s and e <= start + gap for s, e in parts)


def _orientation(a, b, left_id):
    a_strand = _integer(a["strand"], "strand")
    b_strand = _integer(b["strand"], "strand")
    if a_strand not in {-1, 1} or b_strand not in {-1, 1}:
        return "unresolved", "unresolved"
    if a_strand == b_strand:
        return "same_strand", "tandem_forward" if a_strand == 1 else "tandem_reverse"
    left, right = (a_strand, b_strand) if left_id == a["protein_id"] else (b_strand, a_strand)
    return "opposite_strands", "convergent" if (left, right) == (1, -1) else "divergent"


def pair_metrics(a, b, genes, length, topology):
    """Return CDS-order separation and physical distances for distinct anchors.

    Circular edge-path ties are retained. If their intervening counts differ,
    the single intervening count is unavailable and architecture is unresolved.
    For overlapping circular anchors, the shorter midpoint arc defines order.
    Linear order follows the validated, start-sorted CDS indices, including
    overlaps and nested intervals. Supplemental complete-gap counts require
    every covered base of the other CDS to lie inside the physical edge gap.
    """
    if a["protein_id"] == b["protein_id"]:
        raise ValueError("Lysin A and B cannot be the same CDS")
    circular = _circular(topology)
    parts = {g["protein_id"]: _parts(g, length) for g in genes}
    ap, bp = parts[a["protein_id"]], parts[b["protein_id"]]
    overlap = sum(max(0, min(ae, be) - max(ast, bst)) for ast, ae in ap for bst, be in bp)
    am, bm = _midpoint(ap, length, circular), _midpoint(bp, length, circular)
    absolute = abs(am - bm)
    if overlap:
        if circular:
            directions = [((bm - am) % length, a["protein_id"]), ((am - bm) % length, b["protein_id"])]
            shortest = min(distance for distance, _ in directions)
            paths = [(0, 0, pid) for distance, pid in directions if distance == shortest]
        else:
            paths = [(0, 0, a["protein_id"])]
    else:
        paths = []
        for ast, ae in ap:
            for bst, be in bp:
                if circular:
                    paths.extend([((bst - ae) % length, ae % length, a["protein_id"]),
                                  ((ast - be) % length, be % length, b["protein_id"])])
                elif ae <= bst:
                    paths.append((bst - ae, ae, a["protein_id"]))
                elif be <= ast:
                    paths.append((ast - be, be, b["protein_id"]))
        shortest = min(p[0] for p in paths)
        paths = sorted(set(p for p in paths if p[0] == shortest))
    complete_counts = [sum(_inside(parts[g["protein_id"]], start, gap, length, circular)
                           for g in genes if g["protein_id"] not in {a["protein_id"], b["protein_id"]})
                       for gap, start, _ in paths]
    if circular:
        midpoints = {pid: _midpoint(intervals, length, True) for pid, intervals in parts.items()}
        counts = []
        for _, _, left_id in paths:
            left, right = (am, bm) if left_id == a["protein_id"] else (bm, am)
            span = (right - left) % length
            counts.append(sum(0 < (position - left) % length < span
                              for pid, position in midpoints.items() if pid not in {a["protein_id"], b["protein_id"]}))
    else:
        ordered = sorted(genes, key=lambda g: (_integer(g["cds_index"], "CDS index"), int(g["start"]), int(g["end"]), g["protein_id"]))
        positions = {g["protein_id"]: index for index, g in enumerate(ordered)}
        counts = [abs(positions[a["protein_id"]] - positions[b["protein_id"]]) - 1] * len(paths)
    relative, orientation = _orientation(a, b, paths[0][2])
    if overlap and relative == "opposite_strands":
        orientation = "overlapping_opposite_strands"
    if not overlap and len({p[2] for p in paths}) > 1 and relative == "opposite_strands":
        orientation = "ambiguous_circular_path"
    return {
        "gap_bp": paths[0][0], "overlap_bp": overlap,
        "absolute_midpoint_distance_bp": absolute,
        "midpoint_distance_bp": min(absolute, length - absolute) if circular else absolute,
        "intervening_cds": counts[0] if len(set(counts)) == 1 else NA,
        "intervening_cds_min": min(counts), "intervening_cds_max": max(counts),
        "complete_gap_cds": complete_counts[0] if len(set(complete_counts)) == 1 else NA,
        "complete_gap_cds_min": min(complete_counts), "complete_gap_cds_max": max(complete_counts),
        "gene_order_basis": "circular_midpoint_locus_order_along_selected_edge_path" if circular else "linear_cds_index_order",
        "relative_strands": relative, "orientation": orientation,
        "distance_basis": "circular_shortest_edge_arc" if circular else "linear_record_bounds",
        "shortest_path_count": len(paths), "path_ambiguous": len(paths) > 1,
        "a_normalized_position": am / length, "b_normalized_position": bm / length,
    }


def _assignments(genes, evidence, scope):
    labels = defaultdict(dict)
    annotations = defaultdict(set)
    for gene in genes:
        category = gene.get("annotation_category", "")
        if category in CATEGORIES:
            annotations[gene["protein_id"]].add(category)
    for row in evidence:
        category = row.get("annotation_category", "")
        if category in CATEGORIES:
            annotations[row["protein_id"]].add(category)
        elif str(row.get("source", "")).lower() in {"annotation", "genbank", "genbank_annotation"} and row["category"] in CATEGORIES:
            annotations[row["protein_id"]].add(row["category"])
        category = row["category"]
        confidence = row["confidence"]
        if category not in CATEGORIES or confidence not in CONFIDENCE:
            continue
        if CONFIDENCE[confidence] > CONFIDENCE.get(labels[row["protein_id"]].get(category, "unresolved"), 0):
            labels[row["protein_id"]][category] = confidence
    selected = defaultdict(dict)
    partial = set()
    for gene in genes:
        pid = gene["protein_id"]
        categories = {c: "annotation_only" for c in annotations[pid]} if scope == "annotation" else {
            c: confidence for c, confidence in labels[pid].items()
            if CONFIDENCE[confidence] >= (3 if scope == "high" else 2)}
        if _truth(gene.get("partial", "false")):
            if categories:
                partial.add(pid)
            continue
        selected[pid].update(categories)
    conflicts = {pid for pid, categories in selected.items() if {"lysin_a", "lysin_b"} <= categories.keys()}
    return selected, partial, conflicts


def _index_distance(a, b, count, circular):
    distance = abs(_integer(a["cds_index"], "CDS index") - _integer(b["cds_index"], "CDS index"))
    return min(distance, count - distance) if circular else distance


def _pair_row(genome, a, b, genes, selected, scope):
    row = {key: genome.get(key, NA) for key in ("accession", "phage", "host_genus", "cluster")}
    row.update(scope=scope, lysin_a=a["protein_id"], lysin_b=b["protein_id"])
    row["confidence"] = "annotation_only" if scope == "annotation" else min(
        [selected[a["protein_id"]]["lysin_a"], selected[b["protein_id"]]["lysin_b"]], key=CONFIDENCE.get)
    for prefix, gene in (("a", a), ("b", b)):
        for name in ("start", "end", "strand"):
            row[f"{prefix}_{name}"] = gene[name]
        row[f"{prefix}_location_parts"] = json.dumps(_parts(gene, int(genome["length"])), separators=(",", ":"))
    row.update(pair_metrics(a, b, genes, int(genome["length"]), genome["topology"]))
    circular = _circular(genome["topology"])
    holins = [g for g in genes if {"holin", "holin_like"}.intersection(selected[g["protein_id"]])
              and min(_index_distance(g, a, len(genes), circular), _index_distance(g, b, len(genes), circular)) <= 5]
    row["nearby_holin_count"] = len(holins)
    row["nearby_holin_ids"] = ";".join(sorted(g["protein_id"] for g in holins)) or NA
    row["other_lysis_count"] = sum(bool(c) for pid, c in selected.items() if pid not in {a["protein_id"], b["protein_id"]})
    structural = [g for g in genes if re.search(r"\b(?:capsid|tail)\b", g.get("product", ""), re.IGNORECASE)
                  and g["protein_id"] not in {a["protein_id"], b["protein_id"]}]
    if structural:
        distance, pid = min((min(_index_distance(g, a, len(genes), circular), _index_distance(g, b, len(genes), circular)),
                             g["protein_id"]) for g in structural)
        row.update(structural_distance_cds=distance, nearest_structural_protein=pid,
                   structural_basis="explicit_capsid_or_tail_annotation_proxy_not_validated_module")
    else:
        row.update(structural_distance_cds=NA, nearest_structural_protein=NA, structural_basis="no_explicit_structural_annotation")
    return row


def classify(genome, genes, evidence, threshold=4, scope="high_moderate"):
    """Return one genome classification and all distinct eligible A/B pairs."""
    if scope not in SCOPES or int(threshold) < 0:
        raise ValueError("Invalid architecture scope or separation threshold")
    selected, partial, conflicts = _assignments(genes, evidence, scope)
    a_genes = [g for g in genes if "lysin_a" in selected[g["protein_id"]] and g["protein_id"] not in conflicts]
    b_genes = [g for g in genes if "lysin_b" in selected[g["protein_id"]] and g["protein_id"] not in conflicts]
    pairs = [_pair_row(genome, a, b, genes, selected, scope) for a in a_genes for b in b_genes]
    row = {key: NA for key in ARCH_FIELDS}
    row.update({key: genome.get(key, NA) for key in ("accession", "phage", "host_genus", "cluster", "subcluster")})
    row.update(scope=scope, threshold=int(threshold), confidence="unresolved",
               lysin_a_count=len(a_genes), lysin_b_count=len(b_genes),
               holin_count=sum("holin" in c for c in selected.values()),
               holin_like_count=sum("holin_like" in c for c in selected.values()),
               lysis_protein_count=sum(bool(c) for c in selected.values()),
               partial_lysis_count=len(partial), conflicting_protein_ids=";".join(sorted(conflicts)) or NA)
    if conflicts:
        row.update(architecture="unresolved", classification_reason="conflicting_a_b_functional_assignments")
    elif len(a_genes) > 1 or len(b_genes) > 1:
        row.update(architecture="complex", classification_reason="multiple_candidate_lysis_copies")
    elif len(pairs) != 1:
        row.update(architecture="unresolved", classification_reason="incomplete_identification_not_evidence_of_absence")
    else:
        row.update(pairs[0])
        if row["intervening_cds"] == NA:
            row.update(architecture="unresolved", classification_reason="tied_circular_paths_have_different_cds_counts", confidence="unresolved")
        else:
            row.update(architecture="closely_linked" if row["intervening_cds"] <= int(threshold) else "separated",
                       classification_reason="unique_nonpartial_a_b_pair")
    return row, pairs


def _validate(genomes, genes, evidence):
    by_accession = {}
    for genome in genomes:
        accession = genome["accession"]
        if accession in by_accession or _integer(genome["length"], "genome length") <= 0:
            raise ValueError(f"Duplicate accession or invalid genome length: {accession}")
        by_accession[accession] = genome
    by_genome = defaultdict(list)
    ids = set()
    for gene in genes:
        pid, accession = gene["protein_id"], gene["accession"]
        if pid in ids or accession not in by_accession:
            raise ValueError(f"Duplicate protein or unknown genome accession: {pid}/{accession}")
        ids.add(pid)
        _parts(gene, _integer(by_accession[accession]["length"], "genome length"))
        if _integer(gene["strand"], "strand") not in {-1, 0, 1}:
            raise ValueError(f"Invalid strand: {pid}")
        by_genome[accession].append(gene)
    for accession, genome in by_accession.items():
        ordered = sorted(by_genome[accession], key=lambda g: _integer(g["cds_index"], "CDS index"))
        if [_integer(g["cds_index"], "CDS index") for g in ordered] != list(range(len(ordered))):
            raise ValueError(f"CDS indices must be unique and contiguous: {accession}")
        if [int(g["start"]) for g in ordered] != sorted(int(g["start"]) for g in ordered):
            raise ValueError(f"CDS indices must follow genomic starts: {accession}")
        if _integer(genome["cds_count"], "CDS count") != len(ordered):
            raise ValueError(f"CDS count disagrees with genome metadata: {accession}")
        by_genome[accession] = ordered
    for row in evidence:
        if row["protein_id"] not in ids:
            raise ValueError(f"Evidence references unknown protein: {row['protein_id']}")
        if row["confidence"] not in CONFIDENCE:
            raise ValueError(f"Unknown evidence confidence: {row['confidence']}")
    return by_genome


def run(cfg):
    """Stage 5: primary architectures, all pairs, and prespecified sensitivities."""
    genomes = read_tsv(ROOT / "data/processed/genomes.tsv")
    genes = read_tsv(ROOT / "data/processed/cds.tsv")
    evidence = read_tsv(ROOT / "data/processed/lysis_evidence.tsv")
    by_genome = _validate(genomes, genes, evidence)
    by_protein = {g["protein_id"]: g["accession"] for g in genes}
    ev_by_genome = defaultdict(list)
    for row in evidence:
        ev_by_genome[by_protein[row["protein_id"]]].append(row)
    threshold = _integer(cfg["architecture"]["separation_threshold"], "separation threshold")
    thresholds = sorted(set(_integer(v, "sensitivity threshold") for v in cfg["architecture"]["thresholds"]))
    if threshold < 0 or any(value < 0 for value in thresholds):
        raise ValueError("Architecture thresholds must be nonnegative")
    primary, pairs, sensitivity = [], [], []
    for genome in sorted(genomes, key=lambda g: g["accession"]):
        accession = genome["accession"]
        row, all_pairs = classify(genome, by_genome[accession], ev_by_genome[accession], threshold)
        primary.append(row)
        pairs.extend(all_pairs)
        for scope in SCOPES:
            if scope != "high_moderate":
                _, scoped_pairs = classify(genome, by_genome[accession], ev_by_genome[accession], threshold, scope)
                pairs.extend(scoped_pairs)
            for value in thresholds:
                row, _ = classify(genome, by_genome[accession], ev_by_genome[accession], value, scope)
                sensitivity.append(row)
    write_tsv(ROOT / "data/processed/architectures.tsv", primary, ARCH_FIELDS)
    write_tsv(ROOT / "data/processed/architecture_pairs.tsv", pairs, PAIR_FIELDS)
    write_tsv(ROOT / "results/tables/architectures.tsv", primary, ARCH_FIELDS)
    write_tsv(ROOT / "results/tables/architecture_pairs.tsv", pairs, PAIR_FIELDS)
    write_tsv(ROOT / "results/tables/architecture_sensitivity.tsv", sensitivity, ARCH_FIELDS)
    return {"genomes": len(primary), "eligible_pairs": len(pairs), "sensitivity_rows": len(sensitivity)}
