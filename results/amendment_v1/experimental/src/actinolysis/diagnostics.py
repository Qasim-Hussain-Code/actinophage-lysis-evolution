"""Exploratory alignments for a literature-selected unresolved control.

These local alignments never enter functional classification, architecture,
family clustering, or the neighborhood hypothesis tests. Raw alignment scores
are uncalibrated; no E-values or new functional assignments are produced.
"""
from __future__ import annotations

from collections import Counter
import math

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

from .common import NA, ROOT, read_tsv, write_json, write_tsv

CONTROL_QUERIES = (
    {"accession": "KU963248.1", "protein_accession": "AMS02591.1",
     "literature_gene": "47", "comparison_category": "lysin_a"},
    {"accession": "KU963248.1", "protein_accession": "AMS02592.1",
     "literature_gene": "48", "comparison_category": "lysin_a"},
    {"accession": "KU963248.1", "protein_accession": "AMS02705.1",
     "literature_gene": "161", "comparison_category": "lysin_b"},
)
SOURCE_URL = "https://www.ncbi.nlm.nih.gov/nuccore/KU963248.1"
LITERATURE_URL = "https://doi.org/10.1128/genomeA.00594-16"
INTERPRETATION = (
    "Exploratory local alignment of a literature-selected known control; "
    "scores are uncalibrated and have no E-values; no functional assignment "
    "or confidence is changed; absent broad-search hits do not establish gene absence."
)
FIELDS = [
    "control_phage", "accession", "protein_id", "protein_accession",
    "gene", "literature_gene", "start", "end", "strand", "product", "note",
    "comparison_category", "source_url", "literature_url", "query_selection",
    "reference_id", "reference_protein_accession", "reference_evidence_type",
    "reference_url", "query_length", "target_length", "alignment_score",
    "alignment_columns", "identical_residues", "identity",
    "query_aligned_residues", "target_aligned_residues",
    "query_coverage", "target_coverage", "status", "interpretation",
]


def local_alignment_metrics(query, target):
    """Measure one best local BLOSUM62 alignment, with gaps in identity's denominator.

    Query and target coverage count residues consumed by the local alignment,
    including residues aligned opposite gaps. Unaligned flanks are excluded.
    Identity is exact residue matches divided by all local alignment columns.
    """
    query, target = str(query).upper(), str(target).upper()
    if not query or not target:
        raise ValueError("Diagnostic alignment requires nonempty sequences")
    matrix = substitution_matrices.load("BLOSUM62")
    if (set(query) | set(target)) - set(matrix.alphabet):
        raise ValueError("Diagnostic sequence contains residues unsupported by BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(query, target)
    alignment = next(iter(alignments), None)
    if alignment is None:
        return {"query_length": len(query), "target_length": len(target),
                "alignment_score": 0.0, "alignment_columns": 0, "identical_residues": 0,
                "identity": "", "query_aligned_residues": 0, "target_aligned_residues": 0,
                "query_coverage": 0.0, "target_coverage": 0.0, "status": "no_positive_local_alignment"}
    coordinates = alignment.coordinates
    columns = identical = query_residues = target_residues = 0
    for index in range(coordinates.shape[1] - 1):
        qs, qe = (int(x) for x in coordinates[0, index:index + 2])
        ts, te = (int(x) for x in coordinates[1, index:index + 2])
        dq, dt = qe - qs, te - ts
        if dq < 0 or dt < 0 or (dq and dt and dq != dt):
            raise ValueError("Unexpected diagnostic alignment coordinate segment")
        columns += max(dq, dt)
        query_residues += dq
        target_residues += dt
        if dq and dt:
            identical += sum(left == right for left, right in zip(query[qs:qe], target[ts:te]))
    if not columns or not math.isfinite(alignment.score):
        raise ValueError("Invalid diagnostic local alignment")
    return {"query_length": len(query), "target_length": len(target),
            "alignment_score": float(alignment.score), "alignment_columns": columns,
            "identical_residues": identical, "identity": identical / columns,
            "query_aligned_residues": query_residues, "target_aligned_residues": target_residues,
            "query_coverage": query_residues / len(query), "target_coverage": target_residues / len(target),
            "status": "exploratory_alignment_only"}


def run(cfg):
    """Write targeted diagnostics without modifying any scientific analysis inputs."""
    cds = read_tsv(ROOT / "data/processed/cds.tsv")
    references = read_tsv(ROOT / "data/manifests/reference_manifest.tsv")
    sequences = {}
    for record in SeqIO.parse(ROOT / "data/processed/references.faa", "fasta"):
        if record.id in sequences:
            raise ValueError(f"Duplicated diagnostic reference FASTA identifier: {record.id}")
        sequences[record.id] = str(record.seq)
    rows = []
    for selected in CONTROL_QUERIES:
        queries = [row for row in cds if row["accession"] == selected["accession"]
                   and row["protein_accession"] == selected["protein_accession"]]
        if len(queries) > 1:
            raise ValueError(f"Known-control diagnostic query is ambiguous: {selected['protein_accession']}")
        query = queries[0] if queries else None
        for reference in references:
            if (reference["category"] != selected["comparison_category"]
                    or reference["evidence_type"] != "experimental"):
                continue
            rid = reference["reference_id"]
            if rid not in sequences:
                raise ValueError(f"Diagnostic reference sequence is missing: {rid}")
            row = {field: "" for field in FIELDS}
            row.update(selected)
            row.update(control_phage="Yvonnetastic", source_url=SOURCE_URL,
                       literature_url=LITERATURE_URL, query_selection="literature_selected_known_control",
                       reference_id=rid, reference_protein_accession=reference["resolved_protein_accession"],
                       reference_evidence_type=reference["evidence_type"],
                       reference_url=reference["reference_url"], interpretation=INTERPRETATION)
            if query is None:
                row.update(protein_id=NA, gene=NA, status="known_control_query_not_in_frozen_cohort")
            else:
                row.update({key: query.get(key, NA) for key in
                            ["protein_id", "gene", "start", "end", "strand", "product", "note"]})
                row.update(local_alignment_metrics(query["translation"], sequences[rid]))
            rows.append(row)
    write_tsv(ROOT / "results/tables/control_alignment_diagnostics.tsv", rows, FIELDS)
    report = {"method": "Biopython PairwiseAligner local BLOSUM62; gap open -10; gap extension -0.5",
              "identity_definition": "Exact residue matches divided by aligned columns including gaps",
              "coverage_definition": "Residues consumed in the local alignment divided by full sequence length",
              "query_selection": "Literature-selected known control; three exact accession-versioned proteins",
              "n_query_specifications": len(CONTROL_QUERIES), "n_reference_comparisons": len(rows),
              "statuses": dict(Counter(row["status"] for row in rows)), "interpretation": INTERPRETATION,
              "assignment_changes": False, "evalues_computed": False}
    write_json(ROOT / "results/reports/control_alignment_diagnostics.json", report)
    return report
