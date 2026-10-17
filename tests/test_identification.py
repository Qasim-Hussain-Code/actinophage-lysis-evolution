from actinolysis.identification import annotation_category,evaluate_hit
import pytest

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


@pytest.mark.parametrize('requested_source, should_pass', [('AB123456.2', False), ('AB123456.1', True)])
def test_first_freeze_honors_explicit_source_version_even_with_exact_protein_accession(
        tmp_path, monkeypatch, requested_source, should_pass):
    from actinolysis import identification
    from actinolysis.common import read_tsv, write_tsv
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    _synthetic_reference(tmp_path, 'AAG100001.1', 'AAG100001.1')
    catalogue = tmp_path / 'references/protein_controls.tsv'
    specifications = read_tsv(catalogue)
    specifications[0]['source_accession'] = requested_source
    write_tsv(catalogue, specifications)
    raw = tmp_path / 'data/raw/references'
    (raw / f'{requested_source}.gb').write_bytes((raw / 'AB123456.gb').read_bytes())
    if should_pass:
        result = identification.build_references({'offline': True})
        assert result[0]['resolved_source_accession'] == requested_source
    else:
        with pytest.raises(ValueError, match='source accession-version mismatch'):
            identification.build_references({'offline': True})
        assert not (tmp_path / 'data/manifests/reference_manifest.tsv').exists()
        assert not (tmp_path / 'data/processed/references.faa').exists()


@pytest.mark.parametrize('replacement_type', ['published_prediction', 'published_profile_prediction', 'experimental_generic_lipolytic'])
def test_changed_reference_evidence_cannot_reuse_old_experimental_claim(tmp_path, monkeypatch, replacement_type):
    """An identical sequence does not make revised scientific provenance equal."""
    import pytest
    from actinolysis import identification
    from actinolysis.common import read_tsv, write_tsv
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    _synthetic_reference(tmp_path, 'AAG100001.1', 'AAG100001')
    identification.build_references({'offline': True})
    specification = tmp_path/'references/protein_controls.tsv'
    rows = read_tsv(specification)
    rows[0]['evidence_type'] = replacement_type
    write_tsv(specification, rows)
    with pytest.raises(ValueError):
        identification.build_references({'offline': True})


@pytest.mark.parametrize('category', ['lysin_a', 'lysin_b'])
def test_profile_supported_family_prediction_is_capped_at_moderate_even_when_exact(category):
    reference = {'category': category, 'evidence_type': 'published_profile_prediction'}
    assert evaluate_hit(hit(), reference, config()) == 'moderate'
    exact = hit(identity=1, query_coverage=1, target_coverage=1)
    assert evaluate_hit(exact, reference, config()) == 'moderate'
    assert evaluate_hit(exact, {**reference, 'evidence_type': 'experimental'}, config()) == 'high'
    assert evaluate_hit(exact, {**reference, 'evidence_type': 'published_prediction'}, config()) == 'low'


@pytest.mark.parametrize('changes', [dict(identity=.39999), dict(query_coverage=.69999),
                                    dict(target_coverage=.69999), dict(evalue=1.00001e-6)])
def test_profile_prediction_requires_every_unchanged_alignment_threshold(changes):
    reference = {'category': 'lysin_b', 'evidence_type': 'published_profile_prediction'}
    boundary = hit(identity=.4, query_coverage=.7, target_coverage=.7, evalue=1e-6)
    assert evaluate_hit(boundary, reference, config()) == 'moderate'
    assert evaluate_hit({**boundary, **changes}, reference, config()) == 'low'


@pytest.mark.parametrize('category', ['holin', 'holin_like', 'esterase', 'accessory_lysis_chaperone'])
def test_profile_prediction_does_not_promote_holins_or_generic_function_labels(category):
    reference = {'category': category, 'evidence_type': 'published_profile_prediction'}
    assert evaluate_hit(hit(identity=1, query_coverage=1, target_coverage=1), reference, config()) == 'low'
    assert evaluate_hit(hit(), {**reference, 'category': 'negative_control'}, config()) == 'unresolved'


def _synthetic_evidence_case(tmp_path, monkeypatch, references, alignments, partial='false', cfg=None,
                             guard_context_reads=False):
    """Evaluate synthetic scores without acquisition, search or biological data."""
    from actinolysis import identification
    from actinolysis.common import read_tsv, write_tsv
    monkeypatch.setattr(identification, 'ROOT', tmp_path)
    write_tsv(tmp_path/'data/manifests/reference_manifest.tsv', [
        {'reference_id': name, 'category': 'lysin_b', 'evidence_type': tier,
         'reference_url': 'https://invalid.example/synthetic-test-only'} for name, tier in references])
    write_tsv(tmp_path/'data/processed/cds.tsv', [{
        'protein_id': 'SYN000001.1__cds0000', 'protein_accession': 'SYNTHETIC_TEST_ONLY',
        'accession': 'SYN000001.1', 'start': 100, 'end': 1000, 'strand': 1,
        'product': 'hypothetical protein', 'note': '', 'partial': partial}])
    write_tsv(tmp_path/'data/processed/genomes.tsv', [{
        'accession': 'SYN000001.1', 'phage': 'Synthetic test only',
        'host_genus': 'test_only', 'cluster': 'test_only'}])
    fields = ['protein_id', 'reference_id', 'identity', 'alignment_length', 'query_coverage', 'target_coverage', 'evalue', 'bits']
    lines = []
    for reference_id, values in alignments:
        row = {'protein_id': 'SYN000001.1__cds0000', 'reference_id': reference_id,
               'alignment_length': 270, 'bits': 300, **values}
        lines.append('\t'.join(str(row[key]) for key in fields))
    (tmp_path/'data/processed/reference_hits.tsv').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    if guard_context_reads:
        allowed = {'data/manifests/reference_manifest.tsv', 'data/processed/cds.tsv', 'data/processed/genomes.tsv'}
        def independent_reads(path):
            relative = path.relative_to(tmp_path).as_posix()
            assert relative in allowed, f'Functional evidence read forbidden contextual input: {relative}'
            return read_tsv(path)
        monkeypatch.setattr(identification, 'read_tsv', independent_reads)
    identification.evidence(config() if cfg is None else cfg)
    return read_tsv(tmp_path/'data/processed/lysis_evidence.tsv')


@pytest.mark.parametrize('experimental_alignment, expected_reference, expected_confidence', [
    (hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-20), 'experimental', 'high'),
    (hit(evalue=1e-20), 'experimental', 'moderate'),
    (hit(query_coverage=.5, evalue=1e-150), 'profile', 'moderate'),
])
def test_reference_tier_ranking_keeps_valid_experimental_evidence_without_suppressing_qualifying_profile(
        tmp_path, monkeypatch, experimental_alignment, expected_reference, expected_confidence):
    references = [('prediction', 'published_prediction'), ('profile', 'published_profile_prediction'), ('experimental', 'experimental')]
    rows = _synthetic_evidence_case(tmp_path, monkeypatch, references, [
        ('prediction', hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-200, bits=2000)),
        ('profile', hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-100, bits=1000)),
        ('experimental', experimental_alignment)])
    assert len(rows) == 1
    assert rows[0]['reference_id'] == expected_reference
    assert rows[0]['confidence'] == expected_confidence


def test_profile_evidence_does_not_read_architecture_or_neighborhood_inputs(tmp_path, monkeypatch):
    import json
    class AlignmentOnlyConfig(dict):
        def __getitem__(self, key):
            assert key == 'identification', f'Functional assignment used contextual config: {key}'
            return super().__getitem__(key)
    rows = _synthetic_evidence_case(tmp_path, monkeypatch, [('profile', 'published_profile_prediction')],
        [('profile', hit(identity=1, query_coverage=1, target_coverage=1))],
        cfg=AlignmentOnlyConfig(config()), guard_context_reads=True)
    assert rows[0]['confidence'] == 'moderate'
    assert rows[0]['reference_evidence_type'] == 'published_profile_prediction'
    assert rows[0]['context_used'] == 'false' and rows[0]['source'] == 'curated_reference_homology'
    report = json.loads((tmp_path/'results/reports/functional_identification.json').read_text())
    assert report['profile_prediction_policy']['maximum_confidence'] == 'moderate'
    assert report['profile_prediction_policy']['architecture_or_neighborhood_used'] is False
    assert report['reference_evidence_types'] == {'published_profile_prediction': 1}
    assert 'expanded panel' in report['limitations'][0]


def test_partial_target_is_not_promoted_by_profile_evidence(tmp_path, monkeypatch):
    rows = _synthetic_evidence_case(tmp_path, monkeypatch, [('profile', 'published_profile_prediction')],
        [('profile', hit(identity=1, query_coverage=1, target_coverage=1))], partial='true')
    assert rows[0]['confidence'] == 'low'


def test_profile_absent_evidence_retains_original_ranking_and_report_schema(tmp_path, monkeypatch):
    import json
    references = [('weak_evalue', 'experimental'), ('strong_evalue', 'experimental'), ('old_prediction', 'published_prediction')]
    rows = _synthetic_evidence_case(tmp_path, monkeypatch, references, [
        ('weak_evalue', hit(evalue=1e-15, bits=400)),
        ('strong_evalue', hit(evalue=1e-25, bits=300)),
        ('old_prediction', hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-150, bits=1000))])
    assert rows[0]['reference_id'] == 'strong_evalue' and rows[0]['confidence'] == 'moderate'
    report = json.loads((tmp_path/'results/reports/functional_identification.json').read_text())
    assert set(report) == {'assignments', 'by_confidence', 'limitations'}
    assert report['limitations'][0] == 'Only a small experimentally characterized reference panel; divergent functions can be undetected.'


def test_generic_lipolytic_experimental_seed_is_capped_at_moderate_for_lysin_b():
    reference = {'category': 'lysin_b', 'evidence_type': 'experimental_generic_lipolytic'}
    assert evaluate_hit(hit(), reference, config()) == 'moderate'
    assert evaluate_hit(hit(identity=1, query_coverage=1, target_coverage=1), reference, config()) == 'moderate'
    assert evaluate_hit(hit(), {**reference, 'category': 'negative_control'}, config()) == 'unresolved'


@pytest.mark.parametrize('category', ['lysin_a', 'holin', 'holin_like', 'esterase'])
def test_generic_lipolytic_evidence_cannot_promote_other_categories(category):
    assert evaluate_hit(hit(identity=1, query_coverage=1, target_coverage=1),
                        {'category': category, 'evidence_type': 'experimental_generic_lipolytic'}, config()) == 'low'


@pytest.mark.parametrize('changes', [dict(identity=.39999), dict(query_coverage=.69999),
                                    dict(target_coverage=.69999), dict(evalue=1.00001e-6)])
def test_generic_lipolytic_evidence_requires_all_unchanged_alignment_cutoffs(changes):
    reference = {'category': 'lysin_b', 'evidence_type': 'experimental_generic_lipolytic'}
    boundary = hit(identity=.4, query_coverage=.7, target_coverage=.7, evalue=1e-6)
    assert evaluate_hit(boundary, reference, config()) == 'moderate'
    assert evaluate_hit({**boundary, **changes}, reference, config()) == 'low'


@pytest.mark.parametrize('experimental_alignment, expected_reference, expected_confidence', [
    (hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-20), 'experimental', 'high'),
    (hit(evalue=1e-20), 'experimental', 'moderate'),
    (hit(query_coverage=.5, evalue=1e-150), 'generic', 'moderate'),
])
def test_generic_lipolytic_seed_respects_reference_tier_ranking(
        tmp_path, monkeypatch, experimental_alignment, expected_reference, expected_confidence):
    rows = _synthetic_evidence_case(tmp_path, monkeypatch,
        [('generic', 'experimental_generic_lipolytic'), ('experimental', 'experimental')],
        [('generic', hit(identity=1, query_coverage=1, target_coverage=1, evalue=1e-100, bits=1000)),
         ('experimental', experimental_alignment)])
    assert rows[0]['reference_id'] == expected_reference
    assert rows[0]['confidence'] == expected_confidence


def test_generic_lipolytic_report_keeps_substrate_specificity_unresolved(tmp_path, monkeypatch):
    import json
    _synthetic_evidence_case(tmp_path, monkeypatch, [('generic', 'experimental_generic_lipolytic')],
        [('generic', hit(identity=1, query_coverage=1, target_coverage=1))], guard_context_reads=True)
    report = json.loads((tmp_path/'results/reports/functional_identification.json').read_text())
    assert report['reference_evidence_types'] == {'experimental_generic_lipolytic': 1}
    assert report['generic_lipolytic_policy']['maximum_confidence'] == 'moderate'
    assert report['generic_lipolytic_policy']['mycolate_substrate_specificity_established'] is False
    assert report['generic_lipolytic_policy']['architecture_or_neighborhood_used'] is False
    assert any('do not individually establish native' in limitation for limitation in report['limitations'])
