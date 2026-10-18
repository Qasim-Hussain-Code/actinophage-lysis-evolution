from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from scripts.audit_gma6_arabinanase import analyze_record, passes_lysb_thresholds


def cds(start, end, protein_id, locus, product, aa='M'):
    return SeqFeature(
        SimpleLocation(start, end, strand=1), type='CDS',
        qualifiers={
            'protein_id': [protein_id], 'locus_tag': [locus], 'product': [product],
            'translation': [aa],
        },
    )


def test_gma6_enzyme_mapping_preserves_overlapping_neighbors():
    record = SeqRecord(Seq('A' * 100), id='KR063280.1')
    record.annotations['molecule_type'] = 'DNA'
    record.features = [
        cds(0, 20, 'UP.1', 'GMA6_44', 'hypothetical protein'),
        cds(15, 60, 'AKL88326.1', 'GMA6_45', 'putative lysin', 'M' * 15),
        cds(55, 75, 'DOWN.1', 'GMA6_46', 'putative holin'),
    ]

    row, translation = analyze_record(record, 'abc123')

    assert translation == 'M' * 15
    assert (row['start_1based'], row['end_1based'], row['length_aa']) == (16, 60, 15)
    assert row['upstream_overlap_bp'] == 5
    assert row['downstream_overlap_bp'] == 5
    assert row['downstream_product'] == 'putative holin'
    assert row['native_phage_lysis_tested'] is False
    assert 'candidate LysB replacement' in row['published_cassette_context']


def test_gma6_mapping_rejects_wrong_protein_identity():
    record = SeqRecord(Seq('A' * 100), id='KR063280.1')
    record.annotations['molecule_type'] = 'DNA'
    record.features = [
        cds(0, 20, 'UP.1', 'GMA6_44', 'hypothetical protein'),
        cds(15, 60, 'OTHER.1', 'GMA6_45', 'putative lysin'),
        cds(55, 75, 'DOWN.1', 'GMA6_46', 'putative holin'),
    ]

    try:
        analyze_record(record, 'abc123')
    except ValueError as error:
        assert 'Expected one CDS' in str(error)
    else:
        raise AssertionError('Wrong protein identifier was accepted')


def test_lysa_family_hit_cannot_be_mislabeled_as_lysb_assignment():
    assert passes_lysb_thresholds('lysin_a', 0.90, 1.0, 1.0, 0.0) is False
    assert passes_lysb_thresholds('lysin_b', 0.40, 0.70, 0.70, 1e-6) is True
