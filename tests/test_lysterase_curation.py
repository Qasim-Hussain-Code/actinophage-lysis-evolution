"""Synthetic source-mapping and integrity checks; no network or biological results."""
import importlib.util
from pathlib import Path
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

spec = importlib.util.spec_from_file_location('lysterase_curation', Path(__file__).resolve().parents[1] / 'scripts/curate_lysterase_update.py')
curation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(curation)

def protein(coded_by='SYN000001.1:1..12'):
    record = SeqRecord(Seq('MAA'), id='SYNPROT.1')
    record.features = [SeqFeature(FeatureLocation(0, 3), type='CDS', qualifiers={'coded_by': [coded_by]})]
    return record

def genome(sequence='ATGGCTGCTTAA'):
    record = SeqRecord(Seq(sequence), id='SYN000001.1')
    record.features = [SeqFeature(FeatureLocation(0, len(sequence), strand=1), type='CDS',
                                 qualifiers={'protein_id': ['SYNPROT.1'], 'translation': ['MAA']})]
    return record

def test_coding_source_requires_exact_protein_version_and_unambiguous_source():
    assert curation.resolve_protein(protein(), 'SYNPROT.1') == ('SYN000001.1', 'MAA')
    with pytest.raises(ValueError, match='accession-version'):
        curation.resolve_protein(protein(), 'SYNPROT.2')
    with pytest.raises(ValueError, match='ambiguous'):
        curation.resolve_protein(protein('join(SYN000001.1:1..6,SYN000002.1:7..12)'), 'SYNPROT.1')

def test_reference_requires_independent_genomic_translation_not_only_deposited_protein():
    assert curation.verify_cds(genome(), 'SYN000001.1', 'SYNPROT.1', 'MAA').type == 'CDS'
    with pytest.raises(ValueError, match='disagree'):
        curation.verify_cds(genome('ATGTGTTGTTAA'), 'SYN000001.1', 'SYNPROT.1', 'MAA')
    with pytest.raises(ValueError, match='accession-version'):
        curation.verify_cds(genome(), 'SYN000001.2', 'SYNPROT.1', 'MAA')

def test_curation_freezes_cannot_silently_replace_changed_sequence_evidence(tmp_path):
    path = tmp_path / 'synthetic_test_only.tsv'
    curation.freeze_tsv(path, [{'protein': 'SYNPROT.1', 'length': 3}])
    before = path.read_bytes()
    curation.freeze_tsv(path, [{'protein': 'SYNPROT.1', 'length': 3}])
    with pytest.raises(ValueError, match='changed'):
        curation.freeze_tsv(path, [{'protein': 'SYNPROT.1', 'length': 4}])
    assert path.read_bytes() == before
