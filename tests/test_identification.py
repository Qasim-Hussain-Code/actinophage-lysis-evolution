from actinolysis.identification import annotation_category,evaluate_hit

def config():
    return {'identification':{'min_identity':.4,'min_query_coverage':.7,'min_target_coverage':.7,'max_evalue':1e-6}}

def hit(**changes):
    x={'identity':.8,'query_coverage':.9,'target_coverage':.9,'evalue':1e-20}
    x.update(changes); return x

def test_specific_labels_only():
    assert annotation_category('Lysin B')=='lysin_b'
    assert annotation_category('putative endolysin')=='lysin_a'
    assert annotation_category('holin')=='holin_like'
    for p in ['esterase','lipase','transmembrane protein','tail-associated lysozyme','hypothetical protein']:
        assert annotation_category(p)=='unresolved'

def test_independent_reference_evidence_hierarchy():
    r={'category':'lysin_b','evidence_type':'experimental'}
    assert evaluate_hit(hit(),r,config())=='moderate'
    assert evaluate_hit(hit(identity=1,query_coverage=1,target_coverage=1),r,config())=='high'
    assert evaluate_hit(hit(identity=.39),r,config())=='low'
    assert evaluate_hit(hit(target_coverage=.69),r,config())=='low'
    assert evaluate_hit(hit(evalue=.01),r,config())=='low'
    assert evaluate_hit(hit(),{**r,'evidence_type':'published_prediction'},config())=='low'
    assert evaluate_hit(hit(),{**r,'category':'negative_control'},config())=='unresolved'

def test_genbank_note_is_retained_without_independent_promotion(tmp_path,monkeypatch):
    from actinolysis import identification
    from actinolysis.common import read_tsv,write_tsv
    monkeypatch.setattr(identification,'ROOT',tmp_path)
    write_tsv(tmp_path/'data/processed/cds.tsv',[{'protein_id':'SYN000001.1__cds0000','protein_accession':'SYNPROT.1',
        'accession':'SYN000001.1','start':100,'end':112,'strand':1,'partial':'false','product':'serine hydrolase',
        'note':'Lysin B','translation':'MAA'}])
    write_tsv(tmp_path/'data/processed/genomes.tsv',[{'accession':'SYN000001.1','phage':'SYNTHETIC_TEST_ONLY',
        'host_genus':'test_only','cluster':'test_only'}])
    write_tsv(tmp_path/'data/manifests/reference_manifest.tsv',[{'reference_id':'R','category':'lysin_b','evidence_type':'experimental'}])
    (tmp_path/'data/processed/reference_hits.tsv').write_text('')
    identification.evidence(config())
    rows=read_tsv(tmp_path/'data/processed/lysis_evidence.tsv')
    assert len(rows)==1
    assert rows[0]['category']=='lysin_b' and rows[0]['confidence']=='low'
    assert rows[0]['annotation_note']=='Lysin B' and rows[0]['source']=='genbank_annotation'

def test_qualifying_full_coverage_reference_not_suppressed_by_partial_lower_evalue(tmp_path, monkeypatch):
    """A valid full-family match must survive a competing partial match."""
    from actinolysis import identification
    from actinolysis.common import read_tsv, write_tsv
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    references = [{'reference_id': name, 'category': 'lysin_b', 'evidence_type': 'experimental',
                   'reference_url': 'https://invalid.example/synthetic-test-only'} for name in ('partial', 'full')]
    write_tsv(tmp_path/'data/manifests/reference_manifest.tsv', references)
    write_tsv(tmp_path/'data/processed/cds.tsv', [{
        'protein_id': 'SYN000001.1__cds0000', 'protein_accession': 'SYNTHETIC_TEST_ONLY',
        'accession': 'SYN000001.1', 'start': 0, 'end': 900, 'strand': 1,
        'product': 'hypothetical protein', 'partial': 'false'}])
    write_tsv(tmp_path/'data/processed/genomes.tsv', [{
        'accession': 'SYN000001.1', 'phage': 'Synthetic test only',
        'host_genus': 'test_only', 'cluster': 'test_only'}])
    (tmp_path/'data/processed/reference_hits.tsv').write_text(
        'SYN000001.1__cds0000\tpartial\t0.95\t150\t0.50\t1.0\t1e-40\t400\n'
        'SYN000001.1__cds0000\tfull\t0.80\t270\t0.90\t0.90\t1e-20\t300\n')
    identification.evidence(config())
    evidence = read_tsv(tmp_path/'data/processed/lysis_evidence.tsv')
    assert len(evidence) == 1
    assert evidence[0]['reference_id'] == 'full'
    assert evidence[0]['confidence'] == 'moderate'


def _synthetic_reference(tmp_path, protein_accession, target):
    """Set up source data locally so reference tests make no network request."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from actinolysis.common import write_tsv
    source = 'AB123456'
    record = SeqRecord(Seq('ATGGCTGCTTAA'), id=source+'.1',
                       description='Synthetic test-only reference locus')
    record.annotations['molecule_type'] = 'DNA'
    record.features = [SeqFeature(FeatureLocation(0, 12, strand=1), type='CDS', qualifiers={
        'protein_id': [protein_accession], 'gene': ['2'], 'translation': ['MAA']})]
    path = tmp_path/'data/raw/references'/f'{source}.gb'
    path.parent.mkdir(parents=True)
    SeqIO.write(record, path, 'genbank')
    write_tsv(tmp_path/'references/protein_controls.tsv', [{
        'reference_id': 'synthetic_test_only', 'category': 'lysin_a', 'protein_accession': target,
        'source_accession': source, 'gene': '2', 'reference_url': 'https://invalid.example/test-only',
        'evidence_type': 'experimental', 'expected_category': 'lysin_a', 'description': 'synthetic test only'}])


def test_explicit_reference_protein_version_must_match(tmp_path, monkeypatch):
    import pytest
    from actinolysis import identification
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    _synthetic_reference(tmp_path, 'AAG100001.2', 'AAG100001.1')
    with pytest.raises(ValueError):
        identification.build_references({'offline': True})


def test_unversioned_reference_can_resolve_and_freeze_returned_version(tmp_path, monkeypatch):
    from actinolysis import identification
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    _synthetic_reference(tmp_path, 'AAG100001.2', 'AAG100001')
    rows = identification.build_references({'offline': True})
    assert rows[0]['resolved_protein_accession'] == 'AAG100001.2'
    assert rows[0]['resolved_source_accession'] == 'AB123456.1'


def test_changed_reference_evidence_cannot_reuse_old_experimental_claim(tmp_path, monkeypatch):
    """An identical sequence does not make revised scientific provenance equal."""
    import pytest
    from actinolysis import identification
    from actinolysis.common import read_tsv, write_tsv
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    _synthetic_reference(tmp_path, 'AAG100001.1', 'AAG100001')
    identification.build_references({'offline': True})
    specification = tmp_path/'references/protein_controls.tsv'
    rows = read_tsv(specification)
    rows[0]['evidence_type'] = 'published_prediction'
    write_tsv(specification, rows)
    with pytest.raises(ValueError):
        identification.build_references({'offline': True})
