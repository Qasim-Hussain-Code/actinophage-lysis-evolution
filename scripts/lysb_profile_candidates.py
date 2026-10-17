#!/usr/bin/env python3
"""Frozen assay-only LysB HMM diagnostics; candidates never alter primary calls."""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from actinolysis.common import budget_check, load_config, project_bytes, read_tsv, sha256, utcnow, write_json, write_tsv

SEED, FIXED_Z = 20261016, 2460
BROAD = {"D29_LysB", "Ms6_LysB", "WXIN_LysB", "Omega_LysB", "Saal_LysB", "Obama12_LysB", "SWU1_LysB", "BabyRay_LysB", "Palestino_LysB"}
NARROW = {"D29_LysB", "Ms6_LysB"}
PREDICTED = {"Suerte_LysB23", "PokyPuppy_LysB49", "EpicDab_LysB7"}
BUILDER = dict(architecture="fast", weighting="pb", effective_number="entropy", prior_scheme="alphabet",
    symfrac=.5, fragthresh=.5, wid=.62, esigma=45., eid=.62, EmL=200, EmN=200, EvL=200, EvN=200,
    EfL=100, EfN=200, Eft=.04, seed=SEED, ere=None, popen=None, pextend=None, score_matrix=None,
    window_length=None, window_beta=None)
FAMSA = dict(threads=4, guide_tree="sl", tree_heuristic=None, medoid_threshold=2000, medoid_seeds=100,
             medoid_sample=2000, medoid_evaluations=1, cluster_fraction=.1, cluster_iters=2, n_refinements=100,
             keep_duplicates=False, refine=None, scoring_matrix=None)
SEARCH = dict(Z=FIXED_Z, domZ=FIXED_Z, bias_filter=True, null2=True, seed=SEED,
              F1=.02, F2=.001, F3=1e-5, E=10., domE=10., incE=.01, incdomE=.01)
SUMMARY_FIELDS = ["model", "scan", "query_id", "query_kind", "reference_id", "expected_category", "sequence_sha256",
    "length_aa", "exact_training_alias", "reported_hits", "reported_domains", "passes_candidate_cuts",
    "potential_cross_family_flag", "computational_transfer_only", "status"]
HIT_FIELDS = ["model", "scan", "query_id", "full_bits", "full_evalue", "full_pvalue", "full_bias", "reported", "included", "Z", "domZ"]
DOMAIN_FIELDS = HIT_FIELDS + ["domain", "domain_bits", "independent_domain_evalue", "conditional_domain_evalue", "domain_bias",
    "domain_reported", "domain_included", "hmm_from", "hmm_to", "hmm_length", "target_from", "target_to", "target_length",
    "env_from", "env_to", "residue_pairs", "model_coverage", "target_coverage", "model_span_coverage", "target_span_coverage",
    "passes_candidate_cuts", "hmm_alignment", "target_alignment"]


def digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def input_hashes(paths):
    return {k: hashlib.sha256(p.read_text(encoding="utf-8").replace("\r\n", "\n").encode()).hexdigest()
            if k in {"script", "protocol"} else sha256(p) for k, p in paths.items()}


def authenticate_receipts(paths):
    receipts = json.loads(paths["stage_receipts"].read_text())
    for stage, key, relative in [("2", "cds", "data/processed/cds.tsv"), ("3", "references", "data/processed/references.faa"),
                                 ("3", "reference_manifest", "data/manifests/reference_manifest.tsv")]:
        if receipts[stage]["outputs"].get(relative) != sha256(paths[key]):
            raise ValueError(f"Frozen pipeline stage {stage} no longer authenticates {relative}")


def fasta(path):
    result = {}
    for record in SeqIO.parse(path, "fasta"):
        seq = str(record.seq).upper()
        if record.id in result or not seq or set(seq) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError(f"Duplicate or invalid protein: {record.id}")
        result[record.id] = seq
    return result


def global_metrics(a, b):
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                   open_gap_score=-10., extend_gap_score=-.5)
    alignment = next(iter(aligner.align(a, b)))
    columns = pairs = identical = 0
    for i in range(alignment.coordinates.shape[1] - 1):
        x, xx = map(int, alignment.coordinates[0, i:i+2]); y, yy = map(int, alignment.coordinates[1, i:i+2])
        columns += max(xx-x, yy-y)
        if xx > x and yy > y:
            pairs += xx-x
            identical += sum(c == d for c, d in zip(a[x:xx], b[y:yy]))
    return {"identity": identical / columns, "coverage_a": pairs / len(a), "coverage_b": pairs / len(b),
            "paired_residues": pairs, "alignment_columns": columns, "score": float(alignment.score)}


def reference_groups(ids, sequences):
    """Reference-only connected components, explicitly not independent clades."""
    parent, rows = {i: i for i in ids}, []
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for a, b in combinations(sorted(ids), 2):
        metrics = global_metrics(sequences[a], sequences[b])
        edge = metrics["identity"] >= .4 and min(metrics["coverage_a"], metrics["coverage_b"]) >= .7
        if edge:
            parent[find(b)] = find(a)
        rows.append({"reference_a": a, "reference_b": b, **metrics, "group_edge": str(edge).lower()})
    components = defaultdict(list)
    for i in sorted(ids):
        components[find(i)].append(i)
    groups = {i: "RG_" + hashlib.sha256(";".join(members).encode()).hexdigest()[:12]
              for members in components.values() for i in members}
    return groups, rows


def holdout_training(ids, held_ids, references, sequences, source_out):
    hashes = {digest(sequences[i]) for i in held_ids}
    sources = {references[i]["resolved_source_accession"].split(".")[0] for i in held_ids}
    excluded = {i for i in ids if i in held_ids or digest(sequences[i]) in hashes or
                (source_out and references[i]["resolved_source_accession"].split(".")[0] in sources)}
    training = sorted(set(ids) - excluded)
    if {digest(sequences[i]) for i in training} & hashes:
        raise ValueError("Exact held-out sequence leaked into HMM training")
    return training, sorted(excluded)


def unique_training(ids, sequences):
    aliases = defaultdict(list)
    for i in sorted(ids):
        aliases[digest(sequences[i])].append(i)
    return {members[0]: members for _, members in sorted(aliases.items())}


def coverage(hmm_alignment, target_alignment, model_length, target_length, hmm_from, hmm_to, target_from, target_to):
    if len(hmm_alignment) != len(target_alignment) or min(model_length, target_length) < 1:
        raise ValueError("Invalid profile domain alignment")
    pairs = sum(a.isalpha() and b.isalpha() for a, b in zip(hmm_alignment, target_alignment))
    if pairs > min(model_length, target_length):
        raise ValueError("Domain paired-residue count exceeds full sequence/model length")
    return dict(residue_pairs=pairs, model_coverage=pairs / model_length, target_coverage=pairs / target_length,
                model_span_coverage=(hmm_to-hmm_from+1) / model_length,
                target_span_coverage=(target_to-target_from+1) / target_length)


def candidate(full_e, domain_e, model_cov, target_cov):
    return full_e <= 1e-6 and domain_e <= 1e-6 and model_cov >= .7 and target_cov >= .7


def build_model(ids, sequences, out, label, pyfamsa, pyhmmer):
    out.mkdir(parents=True, exist_ok=True)
    aliases = unique_training(ids, sequences)
    aligner = pyfamsa.Aligner(**FAMSA)
    matrix = aligner.scoring_matrix.name
    if matrix != "PFASUM43":
        raise ValueError(f"Unexpected pinned FAMSA scoring matrix: {matrix}")
    alignment = aligner.align([pyfamsa.Sequence(id=i.encode(), sequence=sequences[i].encode()) for i in sorted(aliases)])
    text_msa = pyhmmer.easel.TextMSA(name=label, sequences=[pyhmmer.easel.TextSequence(name=r.id.decode(), sequence=r.sequence.decode()) for r in alignment])
    with (out / "alignment.sto").open("wb") as stream:
        text_msa.write(stream, "stockholm")
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet, **BUILDER)
    hmm, _, _ = builder.build_msa(text_msa.digitize(alphabet), pyhmmer.plan7.Background(alphabet))
    hmm.creation_time = None  # Normalize metadata only; scores and model parameters are unchanged.
    with (out / "model.hmm").open("wb") as stream:
        hmm.write(stream)
    write_json(out / "training.json", {"train_ids": sorted(ids), "unique_training_aliases": aliases,
        "sequence_sha256": {i: digest(sequences[i]) for i in ids}, "builder_constructor": BUILDER,
        "builder_runtime": {k: getattr(builder, k) for k in BUILDER if hasattr(builder, k) and not callable(getattr(builder, k))},
        "actual_builder_seed": builder.seed, "hmm_date_metadata_normalized": True, "famsa_constructor": FAMSA,
        "famsa_runtime": {k: getattr(aligner, k) for k in FAMSA if k != "scoring_matrix"}, "actual_scoring_matrix": matrix,
        "model_length": hmm.M, "model_sequence_count": hmm.nseq, "model_effective_sequences": hmm.nseq_effective})
    return hmm


def search_model(hmm, sequences, model, scan, pyhmmer, settings=None):
    """One worker, fixed statistical search space; preserve every returned hit/domain."""
    options = dict(SEARCH if settings is None else settings)
    alphabet = pyhmmer.easel.Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, [pyhmmer.easel.TextSequence(name=i, sequence=s).digitize(alphabet) for i, s in sorted(sequences.items())])
    hits = pyhmmer.plan7.Pipeline(alphabet, **options).search_hmm(hmm, block)
    if hits.Z != options["Z"] or hits.domZ != options["domZ"]:
        raise ValueError("HMM statistical search space was not held fixed")
    hit_rows, domain_rows = [], []
    for hit in hits:
        name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        row = dict(model=model, scan=scan, query_id=name, full_bits=hit.score, full_evalue=hit.evalue,
                   full_pvalue=hit.pvalue, full_bias=hit.bias, reported=bool(hit.reported), included=bool(hit.included), Z=hits.Z, domZ=hits.domZ)
        hit_rows.append(row)
        for number, domain in enumerate(hit.domains, 1):
            ali = domain.alignment
            ha, ta = ali.hmm_sequence, ali.target_sequence
            if isinstance(ha, bytes):
                ha, ta = ha.decode(), ta.decode()
            cov = coverage(ha, ta, hmm.M, len(sequences[name]), ali.hmm_from, ali.hmm_to, ali.target_from, ali.target_to)
            passed = bool(hit.reported and domain.reported and candidate(hit.evalue, domain.i_evalue, cov["model_coverage"], cov["target_coverage"]))
            domain_rows.append({**row, "domain": number, "domain_bits": domain.score, "independent_domain_evalue": domain.i_evalue,
                "conditional_domain_evalue": domain.c_evalue, "domain_bias": domain.bias, "domain_reported": bool(domain.reported),
                "domain_included": bool(domain.included), "hmm_from": ali.hmm_from, "hmm_to": ali.hmm_to, "hmm_length": hmm.M,
                "target_from": ali.target_from, "target_to": ali.target_to, "target_length": len(sequences[name]),
                "env_from": domain.env_from, "env_to": domain.env_to, **cov, "passes_candidate_cuts": passed,
                "hmm_alignment": ha, "target_alignment": ta})
    return hit_rows, domain_rows


def summarize_queries(model, scan, sequences, metadata, train_ids, reference_sequences, hits, domains, untestable=False):
    train_hashes = {digest(reference_sequences[i]) for i in train_ids}
    rows = []
    for i, seq in sorted(sequences.items()):
        ds, hs = [d for d in domains if d["query_id"] == i], [h for h in hits if h["query_id"] == i and h["reported"]]
        passed = any(d["passes_candidate_cuts"] for d in ds)
        kind = metadata[i]["query_kind"]
        alias = digest(seq) in train_hashes
        status = ("untestable_no_training_sequences" if untestable else "exact_training_sequence_not_independent" if alias else
                  "candidate_only_no_primary_promotion" if passed else "reported_but_candidate_cuts_failed" if hs else "no_reported_profile_hit_not_evidence_of_absence")
        rows.append(dict(model=model, scan=scan, query_id=i, query_kind=kind, reference_id=metadata[i].get("reference_id", ""),
            expected_category=metadata[i].get("expected_category", ""), sequence_sha256=digest(seq), length_aa=len(seq),
            exact_training_alias=alias, reported_hits=len(hs), reported_domains=sum(d["domain_reported"] for d in ds),
            passes_candidate_cuts=passed, potential_cross_family_flag=passed and kind in {"specificity_challenge", "non_B_curated_reference"},
            computational_transfer_only=kind == "computational_B_prediction", status=status))
    return rows


def run(root, validation_dir, out):
    started = time.perf_counter()
    import pyfamsa
    import pyhmmer
    import resource
    from importlib.metadata import version
    versions = {p: version(p) for p in ("pyfamsa", "pyhmmer", "scoring-matrices", "psutil")}
    if versions != {"pyfamsa": "0.7.0", "pyhmmer": "0.12.3", "scoring-matrices": "0.3.5", "psutil": "7.2.2"}:
        raise ValueError("Use the exact optional versions pinned in requirements-profile.txt")
    root, validation_dir, out = map(lambda p: Path(p).resolve(), (root, validation_dir, out))
    protected = [root / p for p in ("data", "references", "src", "scripts", "figures", "results/tables", "results/reports")]
    if out == root or root.is_relative_to(out) or any(out.is_relative_to(p) for p in protected) or out.is_relative_to(validation_dir) or validation_dir.is_relative_to(out):
        raise ValueError("The diagnostic needs a separate output directory")
    query_fasta = validation_dir / "queries.faa"
    if not query_fasta.exists():
        query_fasta = validation_dir / "frozen_queries.faa"
    paths = {"references": root / "data/processed/references.faa", "reference_manifest": root / "data/manifests/reference_manifest.tsv",
        "cds": root / "data/processed/cds.tsv", "queries": query_fasta, "query_manifest": validation_dir / "queries.tsv",
        "stage_receipts": root / "results/reports/pilot_stage_receipts.json", "validation_freeze": validation_dir / "frozen_inputs.json",
        "protocol": PROJECT / "docs/lysb_profile_protocol.md", "script": Path(__file__)}
    authenticate_receipts(paths)
    budget_cfg = load_config(PROJECT / "config/analysis.yaml")
    budget_check(budget_cfg, additional_bytes=50_000_000)
    freeze = {"input_sha256": input_hashes(paths), "builder": BUILDER, "famsa": FAMSA,
              "search": SEARCH, "hmm_workers": 1, "versions": versions}
    prior = out / "frozen_inputs.json"
    if prior.exists() and json.loads(prior.read_text()) != freeze:
        raise ValueError("Changed inputs/code/protocol/versions require a new diagnostic freeze")
    previous = out / "profile_report.json"
    if previous.exists():
        report = json.loads(previous.read_text())
        if not prior.exists() or any(not (out / p).is_file() or sha256(out / p) != h for p, h in report["artifacts"].items()):
            raise ValueError("Previously frozen diagnostic artifacts changed")
        return report
    reference_rows = read_tsv(paths["reference_manifest"])
    references = {r["reference_id"]: r for r in reference_rows}
    seqs, queries = fasta(paths["references"]), fasta(query_fasta)
    if len(references) != len(reference_rows) or set(references) != set(seqs):
        raise ValueError("Reference FASTA/manifest identifiers differ")
    for i, r in references.items():
        if digest(seqs[i]) != r["protein_sha256"] or len(seqs[i]) != int(r["length_aa"]):
            raise ValueError(f"Frozen reference sequence mismatch: {i}")
    if {i for i, r in references.items() if r["category"] == "lysin_b" and r["evidence_type"] in {"experimental", "experimental_generic_lipolytic"}} != BROAD:
        raise ValueError("Expected exactly the prespecified nine assay-supported B references")
    if {i for i, r in references.items() if r["category"] == "lysin_b" and r["evidence_type"] == "published_profile_prediction"} != PREDICTED:
        raise ValueError("Expected exactly three computational Gordonia B evaluation references")
    frozen = json.loads(paths["validation_freeze"].read_text())
    for name, key in [("queries", "query_fasta_sha256"), ("query_manifest", "query_manifest_sha256"), ("references", "reference_fasta_sha256"), ("reference_manifest", "reference_manifest_sha256")]:
        if freeze["input_sha256"][name] != frozen[key]:
            raise ValueError(f"Reference-validation frozen input mismatch: {name}")
    query_rows = read_tsv(paths["query_manifest"])
    if {r["query_id"] for r in query_rows} != set(queries) or len(query_rows) != len(queries):
        raise ValueError("Frozen query identifiers disagree")
    for r in query_rows:
        if digest(queries[r["query_id"]]) != r["protein_sha256"] or len(queries[r["query_id"]]) != int(r["length_aa"]):
            raise ValueError("Frozen query sequence/hash mismatch")
    challenges = [r for r in query_rows if r["query_kind"] == "specificity_challenge"]
    if len(challenges) != 3:
        raise ValueError("Exactly three frozen cellular-enzyme challenges are required")
    cds = read_tsv(paths["cds"])
    targets = {r["protein_id"]: r["translation"].upper().removesuffix("*") for r in cds}
    if len(targets) != len(cds) or len(cds) != FIXED_Z or any(not s or set(s) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO") for s in targets.values()):
        raise ValueError("The frozen 2460-CDS cohort requires unique identifiers and valid translations")
    groups, pairwise = reference_groups(BROAD, seqs)
    write_json(prior, freeze)
    write_tsv(out / "training_references.tsv", [references[i] for i in sorted(BROAD)])
    write_tsv(out / "reference_group_edges.tsv", pairwise)
    write_tsv(out / "reference_groups.tsv", [{"reference_id": i, "group": groups[i], "source": references[i]["resolved_source_accession"], "sequence_sha256": digest(seqs[i])} for i in sorted(BROAD)])
    metadata = {"REF__" + i: {"reference_id": i, "expected_category": r["category"], "query_kind":
        "assay_B_reference" if i in BROAD else "computational_B_prediction" if i in PREDICTED else "non_B_curated_reference"} for i, r in references.items()}
    validation_sequences = {"REF__" + i: s for i, s in seqs.items()}
    for r in challenges:
        validation_sequences[r["query_id"]] = queries[r["query_id"]]
        metadata[r["query_id"]] = {"query_kind": "specificity_challenge", "expected_category": r["expected_category"]}
    target_metadata = {i: {"query_kind": "cohort_CDS", "expected_category": "unknown_not_assigned"} for i in targets}
    hit_rows, domain_rows, summaries, folds = [], [], [], []
    for panel, ids in (("broad_assay9", BROAD), ("narrow_substrate2", NARROW)):
        definitions = [("full", sorted(ids), [])]
        for kind, mapping in (("source_out", {i: references[i]["resolved_source_accession"].split(".")[0] for i in ids}), ("family_out", {i: groups[i] for i in ids})):
            for group in sorted(set(mapping.values())):
                held = sorted(i for i in ids if mapping[i] == group)
                training, excluded = holdout_training(ids, held, references, seqs, kind == "source_out")
                definitions.append((kind + "_" + group, training, held))
        for fold, training, held in definitions:
            label = panel + "__" + fold
            folds.append(dict(model=label, panel=panel, fold=fold, train_ids=";".join(training), held_ids=";".join(held),
                excluded_ids=";".join(sorted(set(ids) - set(training))), n_unique_training=len(unique_training(training, seqs)),
                status="testable" if training else "untestable_no_training_sequences"))
            hmm = build_model(training, seqs, out / "models" / label, label, pyfamsa, pyhmmer) if training else None
            scans = [("reference_validation", validation_sequences, metadata), ("cohort", targets, target_metadata)] if not held else [
                ("held_out_reference", {"REF__" + i: seqs[i] for i in held}, metadata)]
            for scan, sequences, meta in scans:
                hh, dd = search_model(hmm, sequences, label, scan, pyhmmer) if hmm is not None else ([], [])
                hit_rows.extend(hh); domain_rows.extend(dd)
                summaries.extend(summarize_queries(label, scan, sequences, meta, training, seqs, hh, dd, hmm is None))
    write_tsv(out / "folds.tsv", folds)
    write_tsv(out / "full_sequence_hits.tsv", hit_rows, HIT_FIELDS)
    write_tsv(out / "domain_hits.tsv", domain_rows, DOMAIN_FIELDS)
    write_tsv(out / "query_summary.tsv", summaries, SUMMARY_FIELDS)
    write_tsv(out / "cohort_candidates.tsv", [r for r in summaries if r["scan"] == "cohort" and r["passes_candidate_cuts"]], SUMMARY_FIELDS)
    recovery = [{"model": f["model"], "n_queries": len(rows := [r for r in summaries if r["model"] == f["model"] and r["scan"] == "held_out_reference"]),
        "n_unique_query_sequences": len({r["sequence_sha256"] for r in rows}), "n_recovered": sum(r["passes_candidate_cuts"] for r in rows),
        "n_unique_sequences_recovered": len({r["sequence_sha256"] for r in rows if r["passes_candidate_cuts"]}), "status": f["status"]} for f in folds if f["held_ids"]]
    write_tsv(out / "held_out_recovery.tsv", recovery)
    if input_hashes(paths) != freeze["input_sha256"]:
        raise ValueError("A frozen input changed during HMM execution")
    budget_check(budget_cfg)
    storage_bytes = project_bytes(PROJECT, exact=True)
    report = dict(created_utc=utcnow(), diagnostic_only=True, primary_calls_modified=False, n_cohort_cds=len(cds), fixed_Z=FIXED_Z,
        fixed_domZ=FIXED_Z, hmm_workers=1, n_training_references=len(BROAD), n_reference_similarity_groups=len(set(groups.values())),
        cohort_candidate_counts={p + "__full": sum(r["model"] == p + "__full" and r["scan"] == "cohort" and r["passes_candidate_cuts"] for r in summaries) for p in ("broad_assay9", "narrow_substrate2")},
        cohort_non_training_candidate_counts={p + "__full": sum(r["model"] == p + "__full" and r["scan"] == "cohort" and r["passes_candidate_cuts"] and not r["exact_training_alias"] for r in summaries) for p in ("broad_assay9", "narrow_substrate2")},
        cross_family_flags=[r for r in summaries if r["potential_cross_family_flag"]],
        computational_transfer_matches=[r for r in summaries if r["computational_transfer_only"] and r["passes_candidate_cuts"]],
        interpretation="Candidate homology only; no biochemical activity, native substrate specificity, primary confidence promotion, gene absence or independent origins are established.",
        resources={"elapsed_seconds": time.perf_counter() - started, "peak_self_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                   "project_bytes": storage_bytes, "project_storage_measurement": "exact_file_bytes_at_completion", "hmm_workers": 1, "famsa_threads": 4},
        artifacts={p.relative_to(out).as_posix(): sha256(p) for p in sorted(out.rglob("*")) if p.is_file() and p.name != "profile_report.json"})
    write_json(previous, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "validation-dir", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.validation_dir, args.out), indent=2))


if __name__ == "__main__":
    main()
