"""Publication figures derived only from validated study tables.

Each visual has an accompanying TSV containing its observations and aggregates.
When matched comparisons are unavailable, Figure 5 reports measured eligibility
constraints rather than drawing unsupported conservation estimates.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import MaxNLocator
import numpy as np

from .common import NA, ROOT, read_tsv, sha256, write_json, write_tsv

ARCHITECTURES = ("closely_linked", "separated", "complex", "unresolved")
COLORS = {"closely_linked": "#0072B2", "separated": "#D55E00", "complex": "#CC79A7", "unresolved": "#999999"}
SCOPE_LABELS = {"high": "High confidence", "high_moderate": "High + moderate", "annotation": "Annotation only"}
SCOPE_COLORS = {"high": "#009E73", "high_moderate": "#0072B2", "annotation": "#D55E00"}
PALETTE = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#8C6D31", "#7B6FD0")
MISSING_CLUSTERS = {"", NA, "unknown", "unresolved", "none", "na", "singleton", "unclustered", "unassigned"}
FUNCTIONS = ("lysin_a", "lysin_b", "holin", "holin_like")
RANK = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0}


def _number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def _jitter(key):
    return (int(hashlib.sha256(key.encode()).hexdigest()[:8], 16) % 1001 / 1000 - .5) * .28


def _nice(value):
    return str(value).replace("_", " ")


def _style(ax, integer=False):
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#dddddd", linewidth=.55, alpha=.75)
    if integer:
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))


def _status_legend():
    return [Patch(facecolor=COLORS[category], label=_nice(category)) for category in ARCHITECTURES]


def primary_annotations(sensitivities, cfg):
    return [row for row in sensitivities if row["scope"] == "annotation"
            and int(row["threshold"]) == int(cfg["architecture"]["separation_threshold"])]


def dataset_figure(genomes):
    hosts = sorted({row["host_genus"] for row in genomes})
    cluster_counts = Counter(row["cluster"] for row in genomes)
    clusters = sorted(cluster_counts, key=lambda c: (-cluster_counts[c], c))
    fig, axes = plt.subplots(1, 3, figsize=(12, max(4.3, len(clusters) * .25)), layout="constrained", gridspec_kw={"width_ratios": [1, 1.35, 1.25]})
    data = []
    host_colors = {host: PALETTE[index % len(PALETTE)] for index, host in enumerate(hosts)}
    counts = [sum(row["host_genus"] == host for row in genomes) for host in hosts]
    axes[0].bar(range(len(hosts)), counts, color=[host_colors[host] for host in hosts], width=.6)
    axes[0].set_xticks(range(len(hosts)), hosts, rotation=25, ha="right")
    axes[0].set(ylabel="Included genomes", title="A  Host representation")
    _style(axes[0], True)
    for index, (host, count) in enumerate(zip(hosts, counts)):
        axes[0].text(index, count, str(count), ha="center", va="bottom")
        data.append(dict(panel="host_representation", host_genus=host, n_genomes=count))
    axes[1].barh(range(len(clusters)), [cluster_counts[c] for c in clusters], color="#0072B2")
    axes[1].set_yticks(range(len(clusters)), clusters)
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Included genomes", title="B  PhagesDB clusters")
    axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    for cluster in clusters:
        data.append(dict(panel="cluster_representation", cluster=cluster, n_genomes=cluster_counts[cluster]))
    for row in sorted(genomes, key=lambda r: r["accession"]):
        x = hosts.index(row["host_genus"]) + _jitter(row["accession"])
        y = int(row["length"]) / 1000
        axes[2].scatter(x, y, color=host_colors[row["host_genus"]], s=30, edgecolor="white", linewidth=.4)
        data.append(dict(panel="genome_length", accession=row["accession"], phage=row["phage"], host_genus=row["host_genus"], cluster=row["cluster"], length_bp=int(row["length"]), length_kb=y, plot_x=x))
    axes[2].set_xticks(range(len(hosts)), hosts, rotation=25, ha="right")
    axes[2].set(ylabel="Genome length (kb)", title="C  Genome size")
    _style(axes[2])
    fig.suptitle("Dataset representation", fontsize=13)
    return fig, data, "Host and cluster counts and genome lengths describe the included cohort; cluster membership is a coarse ancestry proxy."


def architecture_figure(architectures, annotations):
    if not annotations:
        return None, [], "Annotation-only architecture rows are unavailable; the requested evidence-scope comparison cannot be plotted."
    scopes = [("high_moderate", architectures), ("annotation", annotations)]
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.6), layout="constrained")
    data = []
    bottoms = np.zeros(2)
    for category in ARCHITECTURES:
        values = [sum(row["architecture"] == category for row in rows) for _, rows in scopes]
        axes[0].bar(range(2), values, bottom=bottoms, color=COLORS[category], width=.58)
        for (scope, _), value in zip(scopes, values):
            data.append(dict(panel="architecture_status", scope=scope, architecture=category, n_genomes=value))
        bottoms += values
    axes[0].set_xticks(range(2), ["High +\nmoderate", "Annotation\nonly"])
    axes[0].set(ylabel="Genomes", title="A  Classification and unresolved calls")
    _style(axes[0], True)
    for index, (scope, rows) in enumerate(scopes):
        for row in rows:
            if row["architecture"] not in {"closely_linked", "separated"}:
                continue
            intervening, gap = _number(row.get("intervening_cds")), _number(row.get("gap_bp"))
            if intervening is None or gap is None:
                raise ValueError("Resolved architecture lacks numeric distances")
            x = index + _jitter(scope + row["accession"])
            axes[1].scatter(x, intervening, color=SCOPE_COLORS[scope], s=33, alpha=.8, edgecolor="white", linewidth=.4)
            axes[2].scatter(x, gap / 1000, color=SCOPE_COLORS[scope], s=33, alpha=.8, edgecolor="white", linewidth=.4)
            data.append(dict(panel="unique_pair_distances", scope=scope, accession=row["accession"], phage=row["phage"], architecture=row["architecture"], confidence=row["confidence"], intervening_cds=intervening, gap_bp=gap, gap_kb=gap / 1000, plot_x=x))
    for ax, title, ylabel in [(axes[1], "B  Unique A/B pairs", "Intervening CDS (gene order)"), (axes[2], "C  Unique A/B pairs", "Nearest edge gap (kb)")]:
        ax.set_xticks(range(2), ["High +\nmoderate", "Annotation\nonly"])
        ax.set(title=title, ylabel=ylabel, xlim=(-.5, 1.5), ylim=(0, None))
        _style(ax, ax == axes[1])
    axes[0].legend(handles=_status_legend(), fontsize=7, loc="upper left", frameon=False)
    fig.suptitle("Lysis architecture and functional annotation sensitivity", fontsize=13)
    return fig, data, "Distances include only unique nonpartial A/B pairs. Unresolved and complex records remain explicit in panel A; deposited annotation alone does not independently validate protein function."


def _best_labels(evidence):
    labels = defaultdict(dict)
    for row in evidence:
        category = row["category"]
        if category not in {"lysin_a", "lysin_b"}:
            continue
        pid = row["protein_id"]
        old = labels[pid].get(category)
        if old is None or RANK[row["confidence"]] > RANK[old["confidence"]]:
            labels[pid][category] = row
    result = {}
    for pid, categories in labels.items():
        best = sorted(categories.values(), key=lambda row: (-RANK[row["confidence"]], row["category"]))
        if len(best) > 1 and best[0]["confidence"] == best[1]["confidence"]:
            result[pid] = dict(category="ambiguous_a_b", confidence=best[0]["confidence"])
        else:
            result[pid] = best[0]
    return result


def genome_map_figure(genomes, cds, evidence, architectures, annotations, families):
    by_name = {row["phage"]: row for row in genomes}
    selected = [by_name[name] for name in ("D29", "Yvonnetastic") if name in by_name]
    separated = sorted([row for row in architectures if row["architecture"] == "separated" and row["confidence"] in {"high", "moderate"}], key=lambda row: (-RANK[row["confidence"]], -int(row["intervening_cds"]), row["accession"]))
    metadata = {row["accession"]: row for row in genomes}
    for row in separated:
        if row["accession"] not in {g["accession"] for g in selected}:
            selected.append(metadata[row["accession"]])
            break
    if not selected:
        # A map still requires an observed A/B candidate rather than an arbitrary
        # phage chosen simply to fill a required figure number.
        for row in annotations:
            if row["architecture"] in {"closely_linked", "separated"}:
                selected.append(metadata[row["accession"]])
                break
    if not selected:
        return None, [], "No control genome or resolved annotation-supported A/B pair is available for a meaningful comparative map."
    selected_ids = {row["accession"] for row in selected}
    gene_families = {row["protein_id"]: row["family_id"] for row in families}
    family_genomes = defaultdict(set)
    for row in cds:
        if row["accession"] in selected_ids and row["protein_id"] in gene_families:
            family_genomes[gene_families[row["protein_id"]]].add(row["accession"])
    shared = sorted([family for family, ids in family_genomes.items() if len(ids) > 1], key=lambda family: (-len(family_genomes[family]), family))[:8]
    family_colors = {family: PALETTE[index] for index, family in enumerate(shared)}
    labels = _best_labels(evidence)
    arch = {row["accession"]: row for row in architectures}
    fig, axes = plt.subplots(len(selected), 1, figsize=(12.8, 2.05 * len(selected) + 1.2), layout="constrained", squeeze=False)
    data = []
    for index, genome in enumerate(selected):
        ax = axes[index, 0]
        accession = genome["accession"]
        label_intervals = {"top": [], "bottom": []}
        data.append(dict(panel="genome_summary", accession=accession, phage=genome["phage"], length_bp=int(genome["length"]), architecture=arch[accession]["architecture"]))
        for gene in sorted([row for row in cds if row["accession"] == accession], key=lambda row: (int(row["start"]), row["protein_id"])):
            label = labels.get(gene["protein_id"], {})
            category, confidence = label.get("category", "other_cds"), label.get("confidence", NA)
            family = gene_families.get(gene["protein_id"], NA)
            fill = family_colors.get(family, "#d9d9d9")
            edge = "#D55E00" if category == "lysin_a" else "#6A3D9A" if category == "lysin_b" else "#555555"
            strand = int(gene["strand"])
            raw = gene.get("location_parts", "")
            parts = json.loads(raw) if raw and raw != NA else [[int(gene["start"]), int(gene["end"])]]
            for part_index, (start, end) in enumerate(parts):
                start, end = int(start), int(end)
                y = .14 if strand > 0 else -.14 if strand < 0 else 0
                height = .20
                head = min((end - start) * .30, int(genome["length"]) * .003)
                if strand > 0:
                    points = [(start, y - height / 2), (end - head, y - height / 2), (end, y), (end - head, y + height / 2), (start, y + height / 2)]
                elif strand < 0:
                    points = [(end, y - height / 2), (start + head, y - height / 2), (start, y), (start + head, y + height / 2), (end, y + height / 2)]
                else:
                    points = [(start, y - height / 2), (end, y - height / 2), (end, y + height / 2), (start, y + height / 2)]
                ax.add_patch(Polygon(points, closed=True, facecolor=fill, edgecolor=edge,
                                     linewidth=1.5 if category in {"lysin_a", "lysin_b"} else .35,
                                     hatch="///" if confidence == "moderate" else "..." if confidence == "low" else None))
                data.append(dict(panel="cds_parts", accession=accession, phage=genome["phage"], protein_id=gene["protein_id"], part_index=part_index, start=start, end=end, strand=strand, arrow_type="forward" if strand > 0 else "reverse" if strand < 0 else "unoriented", plot_y=y, family_id=family, family_fill=fill, outline=edge, category=category, confidence=confidence))
            if category in {"lysin_a", "lysin_b", "ambiguous_a_b"}:
                # Point to a covered base interval, including an origin-spanning
                # join whose envelope midpoint may be far from the actual CDS.
                label_start, label_end = min(parts, key=lambda interval: (-(int(interval[1]) - int(interval[0])), int(interval[0])))
                x = (int(label_start) + int(label_end)) / 2
                text = ("A" if category == "lysin_a" else "B" if category == "lysin_b" else "A/B?") + " (" + confidence + ")"
                side = "top" if category == "lysin_a" else "bottom"
                half_width = max(.03, len(text) * .0026) * int(genome["length"])
                text_interval = (x - half_width, x + half_width)
                level = next((level for level, intervals in enumerate(label_intervals[side])
                              if not any(text_interval[0] < right and text_interval[1] > left for left, right in intervals)),
                             len(label_intervals[side]))
                if level == len(label_intervals[side]):
                    label_intervals[side].append([])
                label_intervals[side][level].append(text_interval)
                label_y = (.55 + .33 * level) * (1 if side == "top" else -1)
                ax.annotate(text, (x, .25 if strand > 0 else -.25 if strand < 0 else 0), (x, label_y), ha="center", fontsize=7, color=edge,
                            arrowprops={"arrowstyle": "-", "color": edge, "linewidth": .5})
                data.append(dict(panel="lysis_labels", accession=accession, protein_id=gene["protein_id"], label=text, plot_x=x, plot_y=label_y, category=category, confidence=confidence))
        ax.axhline(0, color="#cccccc", linewidth=.5, zorder=-1)
        top_limit = max(.82, .77 + .33 * (len(label_intervals["top"]) - 1))
        bottom_limit = max(.82, .77 + .33 * (len(label_intervals["bottom"]) - 1))
        ax.set(xlim=(0, int(genome["length"])), ylim=(-bottom_limit, top_limit), yticks=[], xlabel="Original GenBank coordinate (bp)",
               title=f"{genome['phage']} | {accession} | primary: {_nice(arch[accession]['architecture'])}")
        ax.spines[["top", "left", "right"]].set_visible(False)
    legend = [Patch(facecolor=family_colors[family], label=family.removeprefix("FAM_")[:8]) for family in shared]
    legend += [Patch(facecolor="#d9d9d9", label="Other or unassigned family"), Patch(facecolor="white", edgecolor="#D55E00", linewidth=1.5, label="Lysin A outline"), Patch(facecolor="white", edgecolor="#6A3D9A", linewidth=1.5, label="Lysin B outline")]
    axes[0, 0].legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, 1.65), ncol=6, frameon=False, fontsize=6.5, title="Selected shared families (sequence similarity) and functional labels", title_fontsize=7)
    fig.suptitle("Comparative genome maps in original coordinates", fontsize=13)
    return fig, data, "Arrows retain original CDS parts and strand. Fills mark up to eight families shared among displayed genomes; other families are gray. A/B outlines and confidence labels are functional evidence, including explicitly low-confidence annotations. Family sharing does not establish identical function or conserved regulation."


def lineage_figure(genomes, architectures):
    arch = {row["accession"]: row["architecture"] for row in architectures}
    grouped = defaultdict(list)
    for genome in genomes:
        if genome["cluster"].lower() not in MISSING_CLUSTERS:
            grouped[genome["cluster"]].append(arch[genome["accession"]])
    data, weighted = [], []
    raw = Counter(arch.values())
    for cluster, values in sorted(grouped.items()):
        counts = Counter(values)
        for category in ARCHITECTURES:
            data.append(dict(panel="cluster_components", cluster=cluster, architecture=category, n_genomes=len(values), category_count=counts[category], fraction=counts[category] / len(values)))
    for category in ARCHITECTURES:
        value = sum(Counter(values)[category] / len(values) for values in grouped.values()) / len(grouped) if grouped else None
        weighted.append(value)
        data.append(dict(panel="raw_counts", architecture=category, n_genomes=raw[category], total_genomes=len(genomes)))
        if value is not None:
            data.append(dict(panel="equal_weight_clusters", architecture=category, mean_fraction=value, n_clusters=len(grouped)))
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.7), layout="constrained")
    x = np.arange(len(ARCHITECTURES))
    axes[0].bar(x, [raw[c] for c in ARCHITECTURES], color=[COLORS[c] for c in ARCHITECTURES], width=.65)
    axes[0].set(ylabel="Genomes", title=f"A  Raw cohort counts (n = {len(genomes)})")
    if grouped:
        axes[1].bar(x, weighted, color=[COLORS[c] for c in ARCHITECTURES], width=.65)
        axes[1].set(ylabel="Mean within-cluster fraction", ylim=(0, 1), title=f"B  Equal cluster weight (k = {len(grouped)})")
    else:
        axes[1].bar([0, 1], [0, len(genomes)], color=["#0072B2", "#999999"])
        axes[1].set_xticks([0, 1], ["Usable cluster", "Missing cluster"])
        axes[1].set(ylabel="Genomes", title="B  Cluster-weighting eligibility")
        data.append(dict(panel="cluster_weighting_eligibility", usable_cluster_genomes=0, missing_cluster_genomes=len(genomes)))
    for index, ax in enumerate(axes):
        if index == 0 or grouped:
            ax.set_xticks(x, [_nice(c).replace(" ", "\n") for c in ARCHITECTURES])
        _style(ax, index == 0 or not grouped)
    missing = sum(genome["cluster"].lower() in MISSING_CLUSTERS for genome in genomes)
    fig.suptitle(f"Architecture by evolutionary grouping | genomes without usable cluster: {missing}", fontsize=12)
    return fig, data, "Each usable PhagesDB cluster contributes equal weight in panel B, including its unresolved genomes. Raw counts include every genome; missing cluster labels are excluded from the weighted estimate. This comparison is descriptive and does not establish phylogenetic independence."


def _primary_pairs(pairs, cfg):
    return [row for row in pairs if float(row["family_identity"]) == float(cfg["neighborhood"]["primary_identity"])
            and int(row["window"]) == int(cfg["neighborhood"]["primary_window"])
            and int(row["architecture_threshold"]) == int(cfg["architecture"]["separation_threshold"])
            and row["confidence_set"] == "high_moderate" and row["genome_subset"] == "all" and row["metric"] == "jaccard"]


def neighborhood_figure(genomes, cds, evidence, architectures, families, pairs, cfg):
    primary = _primary_pairs(pairs, cfg)
    eligible = [row for row in primary if row["status"] == "eligible"]
    data = []
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), layout="constrained")
    if eligible:
        grouped = defaultdict(list)
        for row in sorted(eligible, key=lambda row: row["pair_id"]):
            control, target, delta = float(row["control_mean"]), float(row["target_score"]), float(row["delta"])
            if not math.isclose(target - control, delta, abs_tol=1e-9):
                raise ValueError("Figure 5 comparison delta disagrees with source scores")
            category = row["pair_architecture"]
            axes[0].scatter(control, target, color=COLORS.get(category, "#009E73"), s=35, alpha=.75)
            grouped[(row["cluster"], category)].append(delta)
            data.append(dict(panel="observed_vs_background", pair_id=row["pair_id"], cluster=row["cluster"], pair_architecture=category, target_score=target, control_mean=control, delta=delta, n_controls=int(row["n_controls"])))
        axes[0].plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=.8)
        axes[0].set(xlim=(-.03, 1.03), ylim=(-.03, 1.03), xlabel="Mean matched-background Jaccard", ylabel="Lysin B neighborhood Jaccard", title="A  Observed genome-pair comparisons")
        values = sorted(grouped.items())
        for index, ((cluster, category), deltas) in enumerate(values):
            mean = float(np.mean(deltas))
            axes[1].scatter(mean, index, color=COLORS.get(category, "#009E73"), s=45)
            data.append(dict(panel="within_cluster_contrast", cluster=cluster, pair_architecture=category, n_pairs=len(deltas), mean_delta=mean))
        axes[1].set_yticks(range(len(values)), [f"{cluster} | {_nice(category)}" for (cluster, category), _ in values])
        axes[1].axvline(0, color="#777777", linewidth=.8)
        axes[1].set(xlabel="Within-cluster mean (target - background)", title="B  Descriptive cluster contrasts", xlim=(-1.03, 1.03))
        separated_clusters = len({row["cluster"] for row in eligible if row["pair_architecture"] == "separated"})
        minimum = max(5, int(cfg["statistics"]["min_lineages"]))
        availability = f"Separated comparisons: {separated_clusters} clusters; inference requires at least {minimum}."
        caption = "Primary family identity, window, and architecture threshold are used. Controls are matched homologous nonlysis anchors in the same genome pair. Cluster contrasts are descriptive; pairs and controls are not additional independent lineages. " + availability
        fig.suptitle("Lysin B neighborhood conservation and matched background", fontsize=13)
    else:
        genes = {row["protein_id"]: row for row in cds}
        b_calls = defaultdict(set)
        for row in evidence:
            if row["category"] == "lysin_b" and row["confidence"] in {"high", "moderate"} and str(genes[row["protein_id"]].get("partial", "false")).lower() != "true":
                b_calls[genes[row["protein_id"]]["accession"]].add(row["protein_id"])
        unique = {accession: next(iter(proteins)) for accession, proteins in b_calls.items() if len(proteins) == 1}
        assignments = {row["protein_id"]: row["family_id"] for row in families}
        grouped = defaultdict(list)
        for row in genomes:
            if row["cluster"].lower() not in MISSING_CLUSTERS:
                grouped[row["cluster"]].append(row["accession"])
        possible, both, shared = 0, 0, 0
        for accessions in grouped.values():
            for index, left in enumerate(accessions):
                for right in accessions[index + 1:]:
                    possible += 1
                    if left in unique and right in unique:
                        both += 1
                        if unique[left] in assignments and unique[right] in assignments and assignments[unique[left]] == assignments[unique[right]]:
                            shared += 1
        counts = [("Within-cluster\ngenome pairs", possible), ("Unique supported B\nin both genomes", both), ("Shared B\nhomologous family", shared), ("Eligible matched\ncomparisons", 0)]
        axes[0].bar(range(len(counts)), [value for _, value in counts], color=["#999999", "#56B4E9", "#0072B2", "#D55E00"])
        axes[0].set_xticks(range(len(counts)), [label for label, _ in counts], fontsize=8)
        axes[0].set(ylabel="Genome pairs", title="A  Observed comparison eligibility")
        for index, (label, value) in enumerate(counts):
            axes[0].text(index, value, str(value), ha="center", va="bottom")
            data.append(dict(panel="pair_eligibility", criterion=label.replace("\n", " "), n_pairs=value))
        reasons = Counter(row["reason"] or "unspecified_ineligibility" for row in primary)
        if reasons:
            labels = sorted(reasons)
            axes[1].barh(range(len(labels)), [reasons[label] for label in labels], color="#999999")
            axes[1].set_yticks(range(len(labels)), [textwrap.fill(_nice(label), 35) for label in labels], fontsize=8)
            axes[1].set(xlabel="Primary genome-pair rows", title="B  Recorded exclusion reasons")
            axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
            for reason in labels:
                data.append(dict(panel="recorded_ineligibility", reason=reason, n_pairs=reasons[reason]))
        else:
            resolved = sum(row["architecture"] in {"closely_linked", "separated"} for row in architectures)
            cohort = [("Included genomes", len(genomes)), ("Unique supported Lysin B", len(unique)), ("Resolved primary A/B pair", resolved)]
            axes[1].barh(range(len(cohort)), [value for _, value in cohort], color=["#999999", "#56B4E9", "#0072B2"])
            axes[1].set_yticks(range(len(cohort)), [label for label, _ in cohort], fontsize=8)
            axes[1].set(xlabel="Genomes", title="B  Observed functional-coverage constraints")
            axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
            for label, value in cohort:
                data.append(dict(panel="functional_coverage", criterion=label, n_genomes=value))
        min_controls, min_lineages = max(5, int(cfg["statistics"]["min_controls"])), max(5, int(cfg["statistics"]["min_lineages"]))
        data.append(dict(panel="prespecified_constraints", minimum_matched_controls=min_controls, minimum_lineages=min_lineages, primary_pair_rows=len(primary), eligible_pairs=0))
        fig.suptitle("Neighborhood comparison eligibility: H2 unavailable", fontsize=13)
        caption = f"No primary genome-pair comparison meets the prespecified eligibility rules. Counts show actual supported anchors, shared families and recorded exclusions. H2 requires at least {min_controls} matched homologous control anchors per pair and {min_lineages} replicated clusters; no conservation effect or significance is estimated here."
    for ax in axes:
        _style(ax, not eligible)
    return fig, data, caption


def robustness_figure(sensitivities, evidence):
    if not sensitivities:
        return None, [], "Architecture sensitivity rows are unavailable."
    blocks = defaultdict(list)
    for row in sensitivities:
        blocks[(row["scope"], int(row["threshold"]))].append(row)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), layout="constrained")
    data = []
    for scope in ("high", "high_moderate", "annotation"):
        thresholds = sorted(threshold for sc, threshold in blocks if sc == scope)
        for category, axis in [("separated", axes[0]), ("unresolved", axes[1])]:
            values = [sum(row["architecture"] == category for row in blocks[(scope, threshold)]) for threshold in thresholds]
            axis.plot(thresholds, values, marker="o", markersize=4, linewidth=1.4, color=SCOPE_COLORS[scope], label=SCOPE_LABELS[scope])
            for threshold, value in zip(thresholds, values):
                data.append(dict(panel="threshold_sensitivity", scope=scope, threshold=threshold, architecture=category, n_genomes=value, total_genomes=len(blocks[(scope, threshold)])))
    for axis, title in [(axes[0], "A  Separated architectures"), (axes[1], "B  Unresolved architectures")]:
        axis.set(xlabel="Operational intervening-CDS threshold", ylabel="Genomes", title=title, ylim=(0, None))
        axis.xaxis.set_major_locator(MaxNLocator(integer=True))
        _style(axis, True)
    axes[0].legend(frameon=False, fontsize=7)
    best = {}
    for row in evidence:
        if row["category"] not in FUNCTIONS:
            continue
        key = row["protein_id"], row["category"]
        if key not in best or RANK[row["confidence"]] > RANK[best[key]]:
            best[key] = row["confidence"]
    bottom = np.zeros(len(FUNCTIONS))
    for confidence, color in [("high", "#009E73"), ("moderate", "#0072B2"), ("low", "#E69F00"), ("unresolved", "#999999")]:
        counts = [sum(category == function and value == confidence for (_, category), value in best.items()) for function in FUNCTIONS]
        axes[2].bar(range(len(FUNCTIONS)), counts, bottom=bottom, label=confidence, color=color, width=.65)
        bottom += counts
        for category, count in zip(FUNCTIONS, counts):
            data.append(dict(panel="functional_confidence", category=category, confidence=confidence, n_proteins=count))
    axes[2].set_xticks(range(len(FUNCTIONS)), ["Lysin A", "Lysin B", "Holin", "Holin-like"], rotation=25, ha="right")
    axes[2].set(ylabel="Unique protein-category assignments", title="C  Functional evidence")
    axes[2].legend(frameon=False, fontsize=7)
    _style(axes[2], True)
    fig.suptitle("Threshold and functional-evidence robustness", fontsize=13)
    return fig, data, "Counts use all genomes in each prespecified scope and threshold. Evidence rows are deduplicated per protein and functional category using the strongest confidence. Low confidence includes annotation-only calls; unsupported assignments do not establish gene absence."


def _save(root, identifier, title, result):
    figure, rows, caption = result
    if figure is None or not rows:
        if figure is not None:
            plt.close(figure)
        return {"figure_id": identifier, "title": title, "status": "skipped", "reason": caption, "files": [], "data_rows": 0}
    folder = root / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    table = root / f"results/tables/{identifier}_data.tsv"
    fields = list(dict.fromkeys(key for row in rows for key in row))
    write_tsv(table, [{field: row.get(field, NA) for field in fields} for row in rows], fields)
    png, pdf = folder / f"{identifier}.png", folder / f"{identifier}.pdf"
    fixed_time = datetime(2000, 1, 1, tzinfo=timezone.utc)
    try:
        figure.savefig(png, dpi=300, bbox_inches="tight", facecolor="white", metadata={"Software": "actinolysis"})
        figure.savefig(pdf, bbox_inches="tight", facecolor="white", metadata={"Title": title, "Creator": "actinolysis", "Producer": "matplotlib", "CreationDate": fixed_time, "ModDate": fixed_time})
    finally:
        plt.close(figure)
    return {"figure_id": identifier, "title": title, "status": "created", "caption": caption,
            "data_rows": len(rows), "files": [{"path": str(path.relative_to(root)).replace("\\", "/"), "sha256": sha256(path)} for path in (png, pdf, table)]}


def run(cfg):
    root = Path(ROOT)
    paths = {
        "genomes": "data/processed/genomes.tsv", "cds": "data/processed/cds.tsv",
        "architectures": "data/processed/architectures.tsv", "evidence": "data/processed/lysis_evidence.tsv",
        "sensitivities": "results/tables/architecture_sensitivity.tsv", "families": "data/processed/families.tsv",
        "pairs": "results/tables/neighborhood_pairs.tsv",
    }
    inputs = {key: read_tsv(root / path) for key, path in paths.items()}
    if not inputs["genomes"]:
        raise ValueError("Scientific figures require a nonempty observed cohort")
    ids = [row["accession"] for row in inputs["genomes"]]
    arch_ids = [row["accession"] for row in inputs["architectures"]]
    if len(set(ids)) != len(ids) or len(set(arch_ids)) != len(arch_ids) or set(ids) != set(arch_ids):
        raise ValueError("Figures require exactly one primary architecture per genome")
    annotations = primary_annotations(inputs["sensitivities"], cfg)
    annotation_ids = [row["accession"] for row in annotations]
    if annotations and (len(set(annotation_ids)) != len(annotation_ids) or set(annotation_ids) != set(ids)):
        raise ValueError("Annotation sensitivity requires exactly one row per genome at the primary threshold")
    figures = []
    settings = {"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.unicode_minus": False,
                "pdf.fonttype": 42, "ps.fonttype": 42, "svg.hashsalt": "actinolysis", "savefig.dpi": 300}
    with plt.rc_context(settings):
        jobs = [
            ("figure_01_dataset", "Dataset and evolutionary representation", lambda: dataset_figure(inputs["genomes"])),
            ("figure_02_architecture", "Genomic organization of lysis genes", lambda: architecture_figure(inputs["architectures"], annotations)),
            ("figure_03_genome_maps", "Representative comparative genome maps", lambda: genome_map_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], annotations, inputs["families"])),
            ("figure_04_lineages", "Evolutionary distribution of architecture types", lambda: lineage_figure(inputs["genomes"], inputs["architectures"])),
            ("figure_05_neighborhoods", "Neighborhood conservation or measured comparison eligibility", lambda: neighborhood_figure(inputs["genomes"], inputs["cds"], inputs["evidence"], inputs["architectures"], inputs["families"], inputs["pairs"], cfg)),
            ("figure_06_robustness", "Statistical robustness and annotation uncertainty", lambda: robustness_figure(inputs["sensitivities"], inputs["evidence"])),
        ]
        for identifier, title, job in jobs:
            figures.append(_save(root, identifier, title, job()))
    report = {"created_count": sum(row["status"] == "created" for row in figures), "figures": figures,
              "raster_dpi": 300, "pdf_metadata_date": "2000-01-01T00:00:00+00:00",
              "source_inputs": [{"path": path, "rows": len(inputs[key]), "sha256": sha256(root / path)} for key, path in paths.items()],
              "interpretation": "Figures describe observed records and declared evidence scopes. Eligibility diagnostics do not substitute for an unavailable H2 effect estimate."}
    write_json(root / "results/reports/figure_manifest.json", report)
    return report
