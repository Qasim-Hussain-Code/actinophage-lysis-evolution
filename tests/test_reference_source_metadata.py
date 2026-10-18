"""Synthetic source-host metadata; no biological observations are fixtures."""
import importlib.util
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import pytest

SPEC = importlib.util.spec_from_file_location('reference_survey', Path(__file__).resolve().parents[1] / 'scripts/survey_reference_architectures.py')
survey = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(survey)


@pytest.mark.parametrize('host,expected', [
    ('Mycobacterium synthetic_test_only', 'Mycobacterium'),
    ('Mycolicibacterium synthetic_test_only', 'Mycolicibacterium'),
    ('Gordonia synthetic_test_only', 'Gordonia'),
    ('STRAIN_TEST_ONLY', 'not_available'),
    ('', 'not_available'),
])
def test_strain_tokens_are_not_assigned_as_host_genera(host, expected):
    record = SeqRecord(Seq('ATG'))
    record.features = [SeqFeature(SimpleLocation(0, 3), type='source', qualifiers={'host':[host]})]
    preserved, genus = survey.source_host_fields(record)
    assert preserved == (host or 'not_available')
    assert genus == expected


def test_missing_source_host_remains_missing():
    record = SeqRecord(Seq('ATG'), description='Gordonia phage synthetic_test_only')
    assert survey.source_host_fields(record) == ('not_available', 'not_available')
