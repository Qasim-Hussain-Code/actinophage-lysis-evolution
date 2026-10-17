"""Training-set boundary checks; identifiers are fixtures, not new observations."""
import importlib.util
from pathlib import Path
import pytest

spec = importlib.util.spec_from_file_location('lysb_update', Path(__file__).resolve().parents[1] / 'scripts/lysb_profile_update.py')
update = importlib.util.module_from_spec(spec)
spec.loader.exec_module(update)

def references():
    return {**{identifier: {'category': 'lysin_b', 'evidence_type': 'experimental_antimycobacterial' if identifier in update.ADDED else 'experimental'}
               for identifier in update.BROAD},
            **{identifier: {'category': 'lysin_b', 'evidence_type': 'published_profile_prediction'} for identifier in update.engine.PREDICTED}}

def test_training_whitelist_rejects_target_recruitment_and_missing_assayed_reference():
    panel = references()
    update.validate_training(panel)
    panel['SYNTHETIC_TEST_ONLY_COHORT_RECRUIT'] = {'category': 'lysin_b', 'evidence_type': 'experimental'}
    with pytest.raises(ValueError, match='prespecified'):
        update.validate_training(panel)
    del panel['SYNTHETIC_TEST_ONLY_COHORT_RECRUIT']
    del panel['Bongo_LysB_2026']
    with pytest.raises(ValueError, match='prespecified'):
        update.validate_training(panel)

def test_gordonia_predictions_cannot_enter_assay_training():
    panel = references()
    panel[next(iter(update.engine.PREDICTED))]['evidence_type'] = 'experimental'
    with pytest.raises(ValueError, match='prespecified'):
        update.validate_training(panel)

def test_metadata_normalization_preserves_numerical_model_and_ascii_parameters(tmp_path):
    import io
    pyhmmer = pytest.importorskip('pyhmmer')
    source = Path(__file__).resolve().parents[1] / 'results/amendment_v1/lysb_profile_candidates/models/broad_assay9__full/model.hmm'
    with pyhmmer.plan7.HMMFile(source) as stream: hmm = stream.read()
    before = io.BytesIO()
    hmm.write(before)
    target = tmp_path / 'normalized.hmm'
    update.normalize_metadata(hmm, target)
    omit = lambda data: [line for line in data.splitlines() if not line.startswith((b'COM ', b'DATE '))]
    assert omit(before.getvalue()) == omit(target.read_bytes())
    assert hmm.command_line is None and hmm.creation_time is None
