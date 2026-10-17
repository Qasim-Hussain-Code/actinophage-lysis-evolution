"""Synthetic geometry without promoting the candidate's functional identity."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("candidate_review_tests", Path(__file__).resolve().parents[1] / "scripts/review_lysb_candidates.py")
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def test_every_supported_a_is_retained_and_training_alias_is_excluded():
    genes = [{"protein_id": f"TEST_{i}", "protein_accession": f"SYNTHETIC{i}", "accession": "TEST", "cds_index": str(i),
              "start": str(i * 100), "end": str(i * 100 + 90), "strand": "1", "partial": "false", "product": "synthetic test protein"}
             for i in range(10)]
    genomes = [{"accession": "TEST", "phage": "synthetic test genome", "host_genus": "TEST", "cluster": "TEST", "length": "1000", "topology": "linear"}]
    evidence = [{"protein_id": f"TEST_{i}", "category": "lysin_a", "confidence": "moderate", "reference_id": "SYNTHETIC_REF",
                 "source": "curated_reference_homology", "evidence_url": "https://example.invalid/test", "context_used": "false"} for i in (1, 2)]
    candidates = [{"query_id": "TEST_9", "model": m, "passes_candidate_cuts": "True", "exact_training_alias": "False"} for m in ("broad", "narrow")]
    candidates.append({"query_id": "TEST_8", "model": "broad", "passes_candidate_cuts": "True", "exact_training_alias": "True"})
    rows = review.annotate(candidates, genes, genomes, evidence)
    assert len(rows) == 2 and {r["intervening_cds"] for r in rows} == {6, 7}
    assert {r["passing_models"] for r in rows} == {"broad;narrow"}
    assert all(r["candidate_protein_id"] == "TEST_9" and r["existing_b_confidence"] == "not_detected" for r in rows)
    assert all("does not establish" in r["interpretation"] for r in rows)
    missing = review.annotate(candidates, genes, genomes, [])
    assert len(missing) == 1 and missing[0]["distance_status"] == "no_complete_independently_supported_a_detected"
