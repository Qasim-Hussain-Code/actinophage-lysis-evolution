"""Audit frozen files, assess independent controls, and rank observed candidates."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re

from Bio import SeqIO

from .common import NA, ROOT, read_tsv, sha256, write_json, write_tsv

EXPECTED_CONTROLS = {"Ms6": "closely_linked", "D29": "closely_linked", "Yvonnetastic": "separated"}
CANDIDATE_FIELDS = [
    "rank", "accession", "phage", "host_genus", "cluster", "subcluster", "architecture",
    "confidence", "lysin_a", "lysin_b", "intervening_cds", "gap_bp", "relative_strands",
    "partial_lysis_count", "similarity_group", "representative", "quality_basis", "candidate_status",
]


def _number(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _independent(row):
    """Require a traceable functional reference and non-contextual evidence."""
    source = row.get("source", "").lower()
    reference = row.get("reference_id", "")
    return (row.get("confidence") in {"high", "moderate"}
            and reference not in {"", NA, "none", "None"}
            and any(token in source for token in ("homology", "alignment", "hmm", "experimental"))
            and not any(token in source for token in ("neighborhood", "synteny", "context")))


def _pair_independent(architecture, by_protein):
    reasons = []
    for category, key in (("lysin_a", "lysin_a"), ("lysin_b", "lysin_b")):
        pid = architecture.get(key, NA)
        if not any(row.get("category") == category and _independent(row) for row in by_protein.get(pid, [])):
            reasons.append(f"{category}:{pid}:no_traceable_independent_high_or_moderate_evidence")
    return reasons


def audit_checksums(root):
    """Check every path/checksum row in frozen TSV and snapshot manifests."""
    root = Path(root).resolve()
    checks = []
    entries = []
    protein_entries = []
    manifests = sorted((root / "data/manifests").glob("*.tsv"))
    for manifest in manifests:
        rows = read_tsv(manifest)
        for index, row in enumerate(rows, 2):
            if "path" in row and "sha256" in row:
                entries.append((str(manifest.relative_to(root)), index, row))
            if "source_sha256" in row and "source_accession" in row:
                entries.append((str(manifest.relative_to(root)), index, {
                    "path": f"data/raw/references/{row['source_accession']}.gb", "sha256": row["source_sha256"]}))
            if "protein_sha256" in row and "reference_id" in row:
                protein_entries.append((str(manifest.relative_to(root)), index, row))
    for manifest in sorted((root / "data/manifests").glob("*.json")):
        obj = json.loads(manifest.read_text(encoding="utf-8"))
        pages = obj.get("pages", []) if isinstance(obj, dict) else []
        if isinstance(obj,dict) and obj.get('selection_snapshot'):
            entries.append((str(manifest.relative_to(root)),0,obj['selection_snapshot']))
        for index, row in enumerate(pages, 1):
            if "path" in row and "sha256" in row:
                entries.append((str(manifest.relative_to(root)), index, row))
    for manifest, index, row in entries:
        path = (root / row["path"]).resolve()
        result = {"manifest": manifest, "row": index, "path": row["path"],
                  "expected_sha256": row["sha256"], "actual_sha256": NA}
        if not path.is_relative_to(root):
            result.update(status="fail", reason="manifest_path_outside_project")
        elif not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
            result.update(status="fail", reason="invalid_sha256_format")
        elif not path.is_file():
            optional=str(row.get('optional_archive','false')).lower()=='true'
            result.update(status="unavailable_optional" if optional else "fail",
                          reason="original_archive_not_redistributed_frozen_selection_metadata_checked_separately" if optional else "frozen_file_missing")
        else:
            result["actual_sha256"] = sha256(path)
            result.update(status="pass" if result["actual_sha256"] == row["sha256"] else "fail",
                          reason="checksum_matches" if result["actual_sha256"] == row["sha256"] else "checksum_mismatch")
        checks.append(result)
    if protein_entries:
        fasta = root / "data/processed/references.faa"
        proteins = defaultdict(list)
        if fasta.is_file():
            for record in SeqIO.parse(fasta, "fasta"):
                proteins[record.id].append(str(record.seq))
        for manifest, index, row in protein_entries:
            sequences = proteins[row["reference_id"]]
            actual = hashlib.sha256(sequences[0].encode()).hexdigest() if len(sequences) == 1 else NA
            valid_expected = bool(re.fullmatch(r"[0-9a-f]{64}", row["protein_sha256"]))
            matches = valid_expected and actual == row["protein_sha256"]
            checks.append({"manifest": manifest, "row": index,
                           "path": f"data/processed/references.faa#{row['reference_id']}",
                           "expected_sha256": row["protein_sha256"], "actual_sha256": actual,
                           "status": "pass" if matches else "fail",
                           "reason": "invalid_sha256_format" if not valid_expected else "checksum_matches" if matches
                           else "reference_protein_missing_or_duplicated" if len(sequences) != 1 else "checksum_mismatch"})
    return checks


def assess_controls(controls, architectures, by_protein, annotation_architectures=None,
                    annotation_pairs=None, separation_threshold=4):
    by_name = defaultdict(list)
    for row in architectures:
        by_name[row["phage"]].append(row)
    annotation_by_name = defaultdict(list)
    for row in annotation_architectures or []:
        annotation_by_name[row["phage"]].append(row)
    pairs_by_name = defaultdict(list)
    for row in annotation_pairs or []:
        if row.get("scope") == "annotation":
            pairs_by_name[row["phage"]].append(row)
    results = []
    for control in controls:
        name, expected = control["phage"], control["expected_architecture"]
        result = {"phage": name, "expected_architecture": expected,
                  "reference_url": control.get("reference_url", NA),
                  "evidence_basis": control.get("evidence_basis", NA),
                  "availability_note": control.get("availability_note", NA)}
        annotation_rows = annotation_by_name[name]
        annotation = annotation_rows[0]["architecture"] if len(annotation_rows) == 1 else NA
        result["annotation_architecture"] = annotation
        result["annotation_recovery_status"] = (
            "expected_arrangement_recovered_from_annotation_only" if annotation == expected
            else "annotation_arrangement_contradicts_expectation" if annotation in {"closely_linked", "separated"}
            else "annotation_arrangement_not_resolved")
        component_pairs = pairs_by_name[name]
        result["annotation_component_pairs"] = [{key: pair.get(key, NA) for key in
                                                   ("lysin_a", "lysin_b", "intervening_cds", "gap_bp", "confidence")}
                                                  for pair in component_pairs]
        if annotation == "complex" and expected == "separated" and len(annotation_rows) == 1:
            expected_pairs = _number(annotation_rows[0].get("lysin_a_count")) * _number(annotation_rows[0].get("lysin_b_count"))
            distinct_pairs = {(pair.get("lysin_a"), pair.get("lysin_b")) for pair in component_pairs}
            if (expected_pairs > 0 and len(component_pairs) == len(distinct_pairs) == expected_pairs
                    and all(_number(pair.get("intervening_cds"), -1) > int(separation_threshold) for pair in component_pairs)):
                result["annotation_recovery_status"] = "annotation_only_component_separation_recovered_in_complex_architecture"
        rows = by_name[name]
        if not rows:
            result.update(status="unresolved", reason="control_not_in_included_dataset", observed_architecture=NA)
        elif len(rows) != 1:
            result.update(status="unresolved", reason="multiple_genomes_share_control_name", observed_architecture=NA)
        else:
            row = rows[0]
            result.update(accession=row["accession"], observed_architecture=row["architecture"],
                          lysin_a=row.get("lysin_a", NA), lysin_b=row.get("lysin_b", NA))
            reasons = _pair_independent(row, by_protein)
            if row["architecture"] not in {"closely_linked", "separated"}:
                reasons.insert(0, f"architecture_not_resolvable:{row.get('classification_reason', NA)}")
            if control.get("reference_url", NA) in {"", NA} or control.get("evidence_basis", NA) in {"", NA}:
                reasons.append("published_expected_architecture_reference_not_curated")
            if reasons:
                result.update(status="unresolved", reason=";".join(reasons))
            else:
                result.update(status="pass" if row["architecture"] == expected else "fail",
                              reason="independently_resolved_architecture_matches" if row["architecture"] == expected
                              else "independently_resolved_architecture_contradicts_expected_control")
        results.append(result)
    return results


def rank_candidates(genomes, architectures, by_protein, controls):
    metadata = {g["accession"]: g for g in genomes}
    known = {row["phage"] for row in controls}
    candidates = []
    for row in architectures:
        if (row["architecture"] != "separated" or row.get("confidence") not in {"high", "moderate"}
                or _number(row.get("lysin_a_count")) != 1 or _number(row.get("lysin_b_count")) != 1
                or _pair_independent(row, by_protein)):
            continue
        genome = metadata[row["accession"]]
        candidate = {key: row.get(key, NA) for key in CANDIDATE_FIELDS}
        candidate.update(similarity_group=genome.get("similarity_group", row["accession"]),
                         representative=genome.get("representative", NA),
                         quality_basis="included_complete_record_unique_nonpartial_independently_supported_a_b_pair",
                         candidate_status="previously_characterized_control_architecture" if row["phage"] in known
                         else "computationally_supported_separated_architecture_literature_status_not_assessed")
        candidates.append(candidate)
    candidates.sort(key=lambda row: (0 if row["confidence"] == "high" else 1,
                                    _number(row.get("partial_lysis_count"), 999999),
                                    0 if row["representative"] == "true" else 1,
                                    -_number(row["intervening_cds"]), -_number(row["gap_bp"]), row["accession"]))
    for index, row in enumerate(candidates, 1):
        row["rank"] = index
    return candidates


def run(cfg):
    required_manifest = ROOT / f"data/manifests/{cfg['mode']}_inputs.tsv"
    if not required_manifest.is_file():
        raise FileNotFoundError(f"Frozen manifest required for validation: {required_manifest}")
    checks = audit_checksums(ROOT)
    genomes = read_tsv(ROOT / "data/processed/genomes.tsv")
    architectures = read_tsv(ROOT / "data/processed/architectures.tsv")
    genome_ids = [row["accession"] for row in genomes]
    architecture_ids = [row["accession"] for row in architectures]
    if (len(set(genome_ids)) != len(genome_ids) or len(set(architecture_ids)) != len(architecture_ids)
            or set(genome_ids) != set(architecture_ids)):
        raise ValueError("Validation requires exactly one architecture per included genome accession")
    evidence = read_tsv(ROOT / "data/processed/lysis_evidence.tsv")
    by_protein = defaultdict(list)
    for row in evidence:
        by_protein[row["protein_id"]].append(row)
    control_path = ROOT / "references/architecture_controls.tsv"
    controls = read_tsv(control_path) if control_path.exists() else [
        {"phage": name, "expected_architecture": expected, "reference_url": NA, "evidence_basis": NA}
        for name, expected in EXPECTED_CONTROLS.items()]
    sensitivity_path = ROOT / "results/tables/architecture_sensitivity.tsv"
    annotation_architectures = [row for row in read_tsv(sensitivity_path)
                               if row.get("scope") == "annotation"
                               and _number(row.get("threshold"), -1) == int(cfg["architecture"]["separation_threshold"])] if sensitivity_path.exists() else []
    pairs_path = ROOT / "data/processed/architecture_pairs.tsv"
    annotation_pairs = [row for row in read_tsv(pairs_path) if row.get("scope") == "annotation"] if pairs_path.exists() else []
    control_results = assess_controls(controls, architectures, by_protein, annotation_architectures,
                                      annotation_pairs, cfg["architecture"]["separation_threshold"])
    candidates = rank_candidates(genomes, architectures, by_protein, controls)
    checksum_failed = not checks or any(row["status"] == "fail" for row in checks)
    control_failed = any(row["status"] == "fail" for row in control_results)
    controls_unresolved = not control_results or any(row["status"] == "unresolved" for row in control_results)
    report = {
        "status": "fail" if checksum_failed or control_failed else "incomplete" if controls_unresolved else "pass",
        "checksum_status": "fail" if checksum_failed else "pass",
        "checksum_checks": checks, "controls": control_results,
        "architecture_counts": dict(sorted(Counter(row["architecture"] for row in architectures).items())),
        "candidate_count": len(candidates),
        "candidate_ranking": "Lexicographic: high before moderate confidence; fewer additional partial lysis calls; similarity-group representative; decreasing intervening CDS; decreasing bp gap; accession tie-break.",
        "interpretation": "Candidates describe observed separation of independently supported lysis proteins. Separation, proximity, and ranking do not establish regulation, rearrangement, independent origin, or a new biological mechanism.",
        "controls_policy": "Only a uniquely resolved pair with traceable independent high/moderate functional evidence and a curated published expectation may pass. Annotation-only or missing controls are unresolved.",
        "expansion_ready": (bool(candidates) and not checksum_failed and not control_failed
            and any(r['status']=='pass' and r['expected_architecture']=='closely_linked' for r in control_results)
            and any(r['expected_architecture']=='separated' and r.get('annotation_recovery_status') in {
                'expected_arrangement_recovered_from_annotation_only',
                'annotation_only_component_separation_recovered_in_complex_architecture'} for r in control_results)),
        "expansion_reason": "At least one independently supported separated architecture is required before larger acquisition; annotation-only separation cannot satisfy this gate.",
    }
    write_tsv(ROOT / "results/tables/candidates.tsv", candidates, CANDIDATE_FIELDS)
    write_json(ROOT / "results/reports/validation.json", report)
    if checksum_failed:
        raise ValueError("Frozen-input checksum validation failed; inspect results/reports/validation.json")
    if control_failed:
        raise ValueError("Independently identified architecture contradicts a curated positive control; inspect results/reports/validation.json")
    if (ROOT/'data/processed/references.faa').exists() and (ROOT/'data/manifests/reference_manifest.tsv').exists():
        from .diagnostics import run as diagnostic_run
        diagnostic_run(cfg)
    return report
