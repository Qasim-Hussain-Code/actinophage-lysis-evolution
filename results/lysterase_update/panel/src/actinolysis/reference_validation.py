"""Held-out query recoverability and cellular-esterase specificity challenges.

Reference-source exclusion precedes interpretation of alignment hits. Search
statistics use the unchanged full target panel, so this is a cross-source
homology-recoverability assessment rather than a measured functional sensitivity.
All output and small downloaded query assets belong to a separate output folder.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import time
from urllib.parse import urlparse

from Bio import SeqIO

from .acquisition import download, validate_accession
from .common import ROOT, NA, budget_check, load_config, read_tsv, sha256, utcnow, write_json, write_tsv
from .identification import CAPPED_REFERENCE_TYPES, evaluate_hit

CONFIDENCE = {"high": 3, "moderate": 2, "low": 1, "unresolved": 0}
RAW_FIELDS = ["protein_id", "reference_id", "identity", "alignment_length", "query_coverage", "target_coverage", "evalue", "bits"]
SEARCH_SETTINGS = {"threads": 4, "sensitivity": "7.5", "kmer": "5", "alignment_mode": "3",
                   "search_evalue": "0.001", "max_seqs": 100, "search_min_identity": "0.25",
                   "search_coverage": "0.5", "coverage_mode": "0", "split_memory_limit": "1G"}
QUERY_FIELDS = ["query_id", "query_kind", "expected_category", "origin_evidence_type", "source_reference_id",
                "resolved_source_accession", "source_accession_base", "challenge_id", "source_id", "length_aa", "protein_sha256"]
THRESHOLD_FIELDS = ["minimum_identity", "minimum_query_coverage", "minimum_target_coverage", "maximum_evalue"]
HIT_FIELDS = QUERY_FIELDS + ["reference_id", "reference_category", "reference_evidence_type", "reference_source_accession",
    "reference_source_base", "identity", "alignment_length", "query_coverage", "target_coverage", "evalue", "bits",
    "score_confidence", "fold_eligible", "exclusion_reason", "accepted_call", "accepted_expected_category",
    "potential_cross_family_lysis_call"] + THRESHOLD_FIELDS
SUMMARY_FIELDS = QUERY_FIELDS + ["n_reported_hits", "n_eligible_hits", "n_excluded_hits", "n_excluded_panel_references",
    "expected_category_recovered", "accepted_categories", "competing_accepted_categories", "potential_cross_family_lysis_call",
    "status", "best_reference_id", "best_reference_category", "best_reference_evidence_type", "best_confidence",
    "best_identity", "best_query_coverage", "best_target_coverage", "best_evalue", "best_bits", "interpretation"] + THRESHOLD_FIELDS


def source_base(accession):
    """Use only the resolved nucleotide accession base; do not infer aliases."""
    validate_accession(accession, version_required=True)
    return accession.split(".")[0]


def _sequence(sequence, label):
    sequence = str(sequence).upper()
    if not sequence or not set(sequence) <= set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
        raise ValueError(f"Invalid protein sequence: {label}")
    return sequence


def read_reference_panel(root):
    root = Path(root)
    rows = read_tsv(root / "data/manifests/reference_manifest.tsv")
    refs = {row["reference_id"]: row for row in rows}
    if not rows or len(refs) != len(rows):
        raise ValueError("Validation requires a nonempty reference manifest with unique identifiers")
    sequences = {}
    for record in SeqIO.parse(root / "data/processed/references.faa", "fasta"):
        if record.id in sequences:
            raise ValueError(f"Duplicate reference FASTA identifier: {record.id}")
        sequences[record.id] = _sequence(record.seq, record.id)
    if set(sequences) != set(refs):
        raise ValueError("Reference FASTA identifiers disagree with the frozen manifest")
    for identifier, row in refs.items():
        source_base(row["resolved_source_accession"])
        if hashlib.sha256(sequences[identifier].encode()).hexdigest() != row["protein_sha256"]:
            raise ValueError(f"Reference protein checksum mismatch: {identifier}")
        if len(sequences[identifier]) != int(row["length_aa"]):
            raise ValueError(f"Reference protein length mismatch: {identifier}")
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+", identifier):
            raise ValueError(f"Unsafe reference FASTA identifier: {identifier}")
    if len(refs) > SEARCH_SETTINGS["max_seqs"]:
        raise ValueError("Panel exceeds the unchanged 100-target search candidate cap; use separately defined source-fold searches before interpreting held-out failures")
    return refs, sequences


def parse_challenge_entity(path, spec):
    matches = [record for record in SeqIO.parse(path, "fasta") if record.id.split("|")[0] == spec["source_id"]]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one PDB polymer entity {spec['source_id']}; found {len(matches)}")
    record = matches[0]
    header_parts = record.description.split("|")
    if len(header_parts) < 2 or not re.match(r"^Chains?\s+", header_parts[1]):
        raise ValueError(f"PDB entity FASTA lacks an explicit chain label: {spec['source_id']}")
    chain_text = re.sub(r"^Chains?\s+", "", header_parts[1])
    chains = set(filter(None, re.split(r"[,\s]+", chain_text)))
    expected_chains = set(spec["equivalent_author_chains"].split(";"))
    if chains != expected_chains or spec["author_chain_id"] not in chains:
        raise ValueError(f"PDB entity chain labels disagree with curation: {spec['source_id']}")
    sequence = _sequence(record.seq, spec["source_id"])
    if len(sequence) != int(spec["expected_length_aa"]):
        raise ValueError(f"PDB entity length mismatch: {spec['source_id']}")
    return sequence, record.description


def prepare_challenges(specifications, out, cfg, fetch=None):
    """Fetch at most four exact, small RCSB entities; return frozen asset data."""
    fetch = download if fetch is None else fetch
    out = Path(out)
    if not 1 <= len(specifications) <= 4:
        raise ValueError("The specificity challenge must contain between one and four proteins")
    if len({row["challenge_id"] for row in specifications}) != len(specifications):
        raise ValueError("Duplicate challenge identifiers")
    assets, sequences = [], {}
    for spec in sorted(specifications, key=lambda row: row["challenge_id"]):
        identifier = spec["challenge_id"]
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+", identifier):
            raise ValueError("Unsafe challenge identifier")
        url = urlparse(spec["sequence_url"])
        expected_path = f"/fasta/entry/{spec['pdb_id']}/display"
        if url.scheme != "https" or url.netloc != "www.rcsb.org" or url.path != expected_path or url.query:
            raise ValueError("Challenge retrieval must use its curated official RCSB entry FASTA URL")
        if spec["source_id"] != f"{spec['pdb_id']}_{spec['polymer_entity_id']}":
            raise ValueError("Challenge PDB entry and polymer-entity identifier disagree")
        if spec["analysis_use"] != "query_only_exclude_from_functional_reference_panel":
            raise ValueError("Challenge sequence is not explicitly marked query-only")
        raw = out / "raw" / f"{spec['source_id']}.faa"
        fetch(spec["sequence_url"], raw, cfg, offline=cfg.get("offline", False), maximum=64_000)
        sequence, header = parse_challenge_entity(raw, spec)
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest in {row["protein_sha256"] for row in assets}:
            raise ValueError("Duplicate challenge sequences cannot count as independent query observations")
        sequences[identifier] = sequence
        assets.append({"challenge_id": identifier, "source_id": spec["source_id"], "sequence_url": spec["sequence_url"],
                       "author_chain_id": spec["author_chain_id"], "equivalent_author_chains": spec["equivalent_author_chains"],
                       "raw_path": raw.relative_to(out).as_posix(), "raw_sha256": sha256(raw),
                       "returned_fasta_header": header, "length_aa": len(sequence), "protein_sha256": digest})
    return assets, sequences


def write_queries(out, references, reference_sequences, challenges, assets):
    out = Path(out)
    sequences, queries = {}, []
    for asset in assets:
        query_id = "CHAL__" + asset["challenge_id"]
        sequences[query_id] = challenges[asset["challenge_id"]]
        queries.append(dict(query_id=query_id, query_kind="specificity_challenge", expected_category="unassigned_cellular_esterase",
            origin_evidence_type="characterized_cellular_enzyme", source_reference_id=NA, resolved_source_accession=NA,
            source_accession_base=NA, challenge_id=asset["challenge_id"], source_id=asset["source_id"],
            length_aa=asset["length_aa"], protein_sha256=asset["protein_sha256"]))
    for identifier, ref in sorted(references.items()):
        if ref["category"] == "negative_control":
            continue
        query_id = "POS__" + identifier
        sequences[query_id] = reference_sequences[identifier]
        queries.append(dict(query_id=query_id, query_kind="positive_reference", expected_category=ref["category"],
            origin_evidence_type=ref["evidence_type"], source_reference_id=identifier,
            resolved_source_accession=ref["resolved_source_accession"], source_accession_base=source_base(ref["resolved_source_accession"]),
            challenge_id=NA, source_id=ref["resolved_protein_accession"], length_aa=len(sequences[query_id]),
            protein_sha256=hashlib.sha256(sequences[query_id].encode()).hexdigest()))
    queries.sort(key=lambda row: row["query_id"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "queries.faa").write_text("".join(f">{identifier}\n{sequences[identifier]}\n" for identifier in sorted(sequences)), encoding="ascii")
    write_tsv(out / "queries.tsv", queries, QUERY_FIELDS)
    return queries


def freeze_inputs(path, value):
    path = Path(path)
    if path.exists():
        if json.loads(path.read_text(encoding="utf-8")) != value:
            raise ValueError("Reference-validation inputs or provenance changed; use a new output directory to preserve the existing freeze")
    else:
        write_json(path, value)


def parse_hits(path, queries, references):
    query_ids = {row["query_id"] for row in queries}
    hits, seen = [], set()
    with Path(path).open(encoding="utf-8", newline="") as stream:
        for line, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if not fields:
                continue
            if len(fields) != len(RAW_FIELDS):
                raise ValueError(f"Invalid search-output field count at line {line}")
            row = dict(zip(RAW_FIELDS, fields))
            pair = row["protein_id"], row["reference_id"]
            if pair[0] not in query_ids or pair[1] not in references or pair in seen:
                raise ValueError(f"Unknown or duplicate query/reference pair at line {line}")
            seen.add(pair)
            for key in ("identity", "query_coverage", "target_coverage", "evalue", "bits"):
                value = float(row[key])
                if not math.isfinite(value) or value < 0 or (key in {"identity", "query_coverage", "target_coverage"} and value > 1):
                    raise ValueError(f"Invalid alignment {key} at line {line}")
            if int(row["alignment_length"]) <= 0:
                raise ValueError("Alignment lengths must be positive")
            hits.append(row)
    return sorted(hits, key=lambda row: (row["protein_id"], row["reference_id"]))


def _rank(hit, reference, cfg):
    confidence = evaluate_hit(hit, reference, cfg)
    capped_penalty = int(confidence == "moderate" and reference["evidence_type"] in CAPPED_REFERENCE_TYPES)
    return (-CONFIDENCE[confidence], capped_penalty, float(hit["evalue"]), -float(hit["bits"]), hit["reference_id"])


def score_queries(queries, references, hits, cfg):
    """Retain every match and every tested query, including unsuccessful folds."""
    limits = cfg["identification"]
    thresholds = dict(zip(THRESHOLD_FIELDS, [limits["min_identity"], limits["min_query_coverage"], limits["min_target_coverage"], limits["max_evalue"]]))
    by_id = {row["query_id"]: row for row in queries}
    if len(by_id) != len(queries):
        raise ValueError("Duplicate validation query identifier")
    by_query = defaultdict(list)
    alignments = []
    for hit in hits:
        query = by_id[hit["protein_id"]]
        ref = references[hit["reference_id"]]
        confidence = evaluate_hit(hit, ref, cfg)
        ref_base = source_base(ref["resolved_source_accession"])
        exclusion = NA
        if query["query_kind"] == "positive_reference":
            if hit["reference_id"] == query["source_reference_id"]:
                exclusion = "self_reference"
            elif ref_base == query["source_accession_base"]:
                exclusion = "same_resolved_source_accession_base"
            elif ref["protein_sha256"] == query["protein_sha256"]:
                exclusion = "exact_sequence_alias"
        eligible = exclusion == NA
        accepted = eligible and confidence in {"high", "moderate"}
        row = {**query, "reference_id": hit["reference_id"], "reference_category": ref["category"],
               "reference_evidence_type": ref["evidence_type"], "reference_source_accession": ref["resolved_source_accession"],
               "reference_source_base": ref_base, **{key: hit[key] for key in RAW_FIELDS[2:]},
               "score_confidence": confidence, "fold_eligible": eligible, "exclusion_reason": exclusion,
               "accepted_call": accepted, "accepted_expected_category": accepted and ref["category"] == query["expected_category"],
               "potential_cross_family_lysis_call": accepted and query["query_kind"] == "specificity_challenge" and ref["category"] in {"lysin_a", "lysin_b"}, **thresholds}
        alignments.append(row)
        by_query[query["query_id"]].append((hit, row))
    summaries = []
    for query in queries:
        observations = by_query[query["query_id"]]
        eligible = [(hit, row) for hit, row in observations if row["fold_eligible"]]
        accepted_categories = sorted({row["reference_category"] for _, row in eligible if row["accepted_call"]})
        recovered = query["expected_category"] in accepted_categories
        competing = [category for category in accepted_categories if category != query["expected_category"]]
        potential = any(row["potential_cross_family_lysis_call"] for _, row in eligible)
        expected_hits = [(hit, row) for hit, row in eligible if row["reference_category"] == query["expected_category"]]
        best = min(expected_hits or eligible, key=lambda item: _rank(item[0], references[item[0]["reference_id"]], cfg)) if eligible else None
        excluded_count = sum(ref["reference_id"] == query["source_reference_id"] or source_base(ref["resolved_source_accession"]) == query["source_accession_base"]
                             or ref["protein_sha256"] == query["protein_sha256"]
                             for ref in references.values()) if query["query_kind"] == "positive_reference" else 0
        if query["query_kind"] == "specificity_challenge":
            status = "potential_cross_family_lysis_label" if potential else "no_accepted_lysis_label"
            interpretation = "Specificity challenge; mAGP inactivity is not established by the cited characterization."
        else:
            status = "recovered_with_competing_category" if recovered and competing else "expected_category_recovered" if recovered else "expected_category_not_recovered"
            interpretation = "Leave-resolved-source-base-out homology recoverability; query origin evidence is reported separately and functional sensitivity is not estimated."
        row = {**query, "n_reported_hits": len(observations), "n_eligible_hits": len(eligible),
               "n_excluded_hits": len(observations) - len(eligible), "n_excluded_panel_references": excluded_count,
               "expected_category_recovered": recovered, "accepted_categories": ";".join(accepted_categories) or NA,
               "competing_accepted_categories": ";".join(competing) or NA, "potential_cross_family_lysis_call": potential,
               "status": status, "interpretation": interpretation, **thresholds}
        best_values = {"best_reference_id": NA, "best_reference_category": NA, "best_reference_evidence_type": NA,
                       "best_confidence": NA, "best_identity": NA, "best_query_coverage": NA,
                       "best_target_coverage": NA, "best_evalue": NA, "best_bits": NA}
        if best:
            hit, alignment = best
            best_values.update(best_reference_id=hit["reference_id"], best_reference_category=alignment["reference_category"],
                               best_reference_evidence_type=alignment["reference_evidence_type"], best_confidence=alignment["score_confidence"],
                               **{"best_" + key: hit[key] for key in ("identity", "query_coverage", "target_coverage", "evalue", "bits")})
        summaries.append({**row, **best_values})
    return alignments, summaries


def search_command(executable, query_fasta, target_fasta, output, temporary):
    return [str(executable), "easy-search", str(query_fasta), str(target_fasta), str(output), str(temporary),
            "--threads", "4", "-s", "7.5", "-k", "5", "--alignment-mode", "3", "-e", "0.001",
            "--max-seqs", "100", "--min-seq-id", "0.25", "-c", "0.5", "--cov-mode", "0",
            "--split-memory-limit", "1G", "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue,bits", "--remove-tmp-files", "1"]


def _memory():
    try:
        import resource
        scale = 1 if platform.system() == "Darwin" else 1024
        return {"peak_parent_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
                "peak_child_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * scale,
                "method": "process high-water marks; not simultaneous process-tree memory"}
    except ImportError:
        return {"peak_parent_rss_bytes": NA, "peak_child_rss_bytes": NA, "method": "resource module unavailable on this platform"}


def run(root, out, project_root=ROOT, offline=False):
    root, out, project_root = Path(root).resolve(), Path(out).resolve(), Path(project_root).resolve()
    if not out.is_relative_to(project_root):
        raise ValueError("Validation output must stay inside the budgeted project directory")
    marker = out / ".reference_validation"
    if out.exists() and any(out.iterdir()) and not marker.is_file():
        raise ValueError("Choose a new empty output directory or an existing reference-validation output directory")
    out.mkdir(parents=True, exist_ok=True)
    if marker.exists() and marker.read_text() != "actinolysis-reference-validation-v1\n":
        raise ValueError("Unexpected reference-validation directory marker")
    marker.write_text("actinolysis-reference-validation-v1\n")
    started, began = utcnow(), time.monotonic()
    cfg = load_config(root / "config/analysis.yaml")
    cfg["offline"] = offline
    budget_cfg = load_config(project_root / "config/analysis.yaml")
    budget_cfg["offline"] = offline
    refs, sequences = read_reference_panel(root)
    catalogue = project_root / "references/specificity_challenges.tsv"
    tool = cfg["tools"]["mmseqs"]
    executable = (root / tool).resolve() if "/" in tool or "\\" in tool else tool
    version = subprocess.run([str(executable), "version"], check=True, text=True, capture_output=True).stdout.strip()
    frozen_prefix = {"schema": 1, "reference_manifest_sha256": sha256(root / "data/manifests/reference_manifest.tsv"),
                     "reference_fasta_sha256": sha256(root / "data/processed/references.faa"), "challenge_catalogue_sha256": sha256(catalogue),
                     "identification_thresholds": cfg["identification"], "search_settings": SEARCH_SETTINGS,
                     "mmseqs_version": version, "source_sha256": {"reference_validation": sha256(Path(__file__)),
                     "identification": sha256(Path(evaluate_hit.__code__.co_filename))}}
    freeze_path = out / "frozen_inputs.json"
    prior = json.loads(freeze_path.read_text(encoding="utf-8")) if freeze_path.exists() else None
    if prior is not None:
        if any(prior.get(key) != value for key, value in frozen_prefix.items()):
            raise ValueError("Reference-validation inputs or provenance changed; use a new output directory before modifying the existing freeze")
        for filename, key in [("queries.faa", "query_fasta_sha256"), ("queries.tsv", "query_manifest_sha256")]:
            path = out / filename
            if path.exists() and sha256(path) != prior.get(key):
                raise ValueError(f"Frozen validation query input changed: {filename}")
    assets, challenges = prepare_challenges(read_tsv(catalogue), out, budget_cfg)
    if prior is not None and prior.get("challenge_assets") != assets:
        raise ValueError("Frozen PDB challenge asset or sequence changed; preserve the existing freeze and investigate")
    write_tsv(out / "challenge_inputs.tsv", assets)
    queries = write_queries(out, refs, sequences, challenges, assets)
    frozen = {**frozen_prefix,
              "query_fasta_sha256": sha256(out / "queries.faa"), "query_manifest_sha256": sha256(out / "queries.tsv"),
              "challenge_assets": assets}
    freeze_inputs(freeze_path, frozen)
    output = out / "raw_hits.tsv"
    temporary = out / ".tmp/mmseqs"
    temporary.mkdir(parents=True, exist_ok=True)
    command = search_command(executable, out / "queries.faa", root / "data/processed/references.faa", output, temporary)
    receipt_path = out / "search_execution.json"
    old = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    reusable = (old.get("status") == "success" and old.get("frozen_inputs_sha256") == sha256(out / "frozen_inputs.json")
                and output.exists() and old.get("raw_hits_sha256") == sha256(output))
    budget_before = budget_check(budget_cfg, 600_000_000, root=project_root)
    if not reusable:
        search_start = time.monotonic()
        receipt = {"status": "running", "started_utc": utcnow(), "command": command, "mmseqs_version": version,
                   "frozen_inputs_sha256": sha256(out / "frozen_inputs.json")}
        write_json(receipt_path, receipt)
        environment = os.environ.copy()
        environment.update(OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1")
        try:
            with (out / "search.log").open("a", encoding="utf-8") as log:
                log.write("Command arguments: " + json.dumps(command) + "\n")
                log.flush()
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT, cwd=out, env=environment)
            output.write_text("".join(line + "\n" for line in sorted(output.read_text().splitlines()) if line.strip()), encoding="utf-8")
            receipt.update(status="success", raw_hits_sha256=sha256(output))
        except Exception as error:
            receipt.update(status="failed", error=f"{type(error).__name__}: {error}")
            raise
        finally:
            receipt.update(elapsed_seconds=round(time.monotonic() - search_start, 3), memory=_memory())
            write_json(receipt_path, receipt)
    hits = parse_hits(output, queries, refs)
    alignments, summary = score_queries(queries, refs, hits, cfg)
    write_tsv(out / "alignment_evidence.tsv", alignments, HIT_FIELDS)
    write_tsv(out / "query_summary.tsv", summary, SUMMARY_FIELDS)
    groups = defaultdict(list)
    for row in summary:
        if row["query_kind"] == "positive_reference":
            groups[(row["origin_evidence_type"], row["expected_category"])].append(row)
    challenges_summary = [row for row in summary if row["query_kind"] == "specificity_challenge"]
    report = {"started_utc": started, "finished_utc": utcnow(), "elapsed_seconds": round(time.monotonic() - began, 3),
              "status": "success", "isolated_root": str(root), "output_root": str(out), "search_reused": reusable,
              "mmseqs_version": version, "search_command": command, "search_settings": SEARCH_SETTINGS,
              "identification_thresholds": cfg["identification"], "reference_count": len(refs), "query_count": len(queries),
              "query_origins": dict(sorted(Counter(row["origin_evidence_type"] for row in queries).items())),
              "positive_recoverability": [{"origin_evidence_type": evidence_type, "expected_category": category,
                  "n_tested": len(rows), "n_expected_category_recovered": sum(row["expected_category_recovered"] for row in rows),
                  "n_with_competing_accepted_categories": sum(row["competing_accepted_categories"] != NA for row in rows)}
                  for (evidence_type, category), rows in sorted(groups.items())],
              "specificity_challenges": {"n_tested": len(challenges_summary), "n_potential_cross_family_lysis_calls": sum(row["potential_cross_family_lysis_call"] for row in challenges_summary),
                  "biochemical_mAGP_inactivity_established": False},
              "holdout_policy": "Exclude self and every target with the same resolved_source_accession base before interpretation. Accession aliases are not inferred. E-values and prefiltering retain the unchanged full target panel; target count does not exceed the fixed 100-candidate cap.",
              "interpretation": "Positive results measure cross-source homology recoverability, with experimental and computational query origins separated. Challenge labels flag potential cross-family specificity problems, not proven mAGP inactivity. This small selected set does not estimate general functional sensitivity or false-discovery rate.",
              "generic_lipolytic_origin_caveat": "Generic-substrate activity supports the recorded lipolytic assay; cross-source sequence recovery does not establish native mAGP substrate specificity. These query origins remain a separate reported evidence tier.",
              "resources": {"threads": 4, "split_memory_limit": "1G", "memory": _memory(), "budget_before": budget_before,
                            "budget_after": budget_check(budget_cfg, root=project_root),
                            "output_file_bytes": sum(path.stat().st_size for path in out.rglob("*") if path.is_file())},
              "artifacts": [{"path": path.relative_to(out).as_posix(), "sha256": sha256(path)} for path in
                            [out / "frozen_inputs.json", out / "challenge_inputs.tsv", out / "queries.faa", out / "queries.tsv", output,
                             out / "alignment_evidence.tsv", out / "query_summary.tsv", receipt_path, out / "search.log"]]}
    write_json(out / "validation.json", report)
    return report
