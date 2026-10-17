"""Independent analytical checks of exploratory alignment metrics."""
import pytest

from actinolysis.diagnostics import local_alignment_metrics


def test_exact_blosum62_alignment_has_known_score_and_coverage():
    # BLOSUM62 A-A is 4, so two identical alanines score exactly 8.
    result = local_alignment_metrics("AA", "AA")
    assert result["alignment_score"] == 8
    assert result["identical_residues"] == result["alignment_columns"] == 2
    assert result["identity"] == result["query_coverage"] == result["target_coverage"] == 1


def test_internal_gap_contributes_to_identity_denominator():
    # Eight W-W matches score 8*11. A four-residue gap costs 10+3*0.5.
    result = local_alignment_metrics("WWWWAAAAWWWW", "WWWWWWWW")
    assert result["alignment_score"] == pytest.approx(76.5)
    assert result["identical_residues"] == 8
    assert result["alignment_columns"] == 12
    assert result["identity"] == pytest.approx(8 / 12)
    assert result["query_coverage"] == result["target_coverage"] == 1
