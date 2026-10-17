"""Synthetic fixtures only: evidence separation, minimum controls and weighting."""
from copy import deepcopy
import pytest
from actinolysis import annotation_sensitivity as sensitivity

def fixture():
    cfg = {'architecture': {'thresholds': [4]}, 'neighborhood': {'windows': [3]},
           'statistics': {'min_controls': 1, 'matched_controls': 20}}
    genomes, genes, evidence, memberships, architectures = [], [], [], [], []
    for acc in ('SYNTHETIC_TEST_ONLY_A', 'SYNTHETIC_TEST_ONLY_B'):
        genomes.append({'accession': acc, 'cluster': 'SYNTHETIC_TEST_ONLY', 'length': '5000',
                        'topology': 'linear', 'similarity_group': acc, 'representative': 'true'})
        for index in range(40):
            identifier = f'{acc}__cds{index:04d}'
            category = 'lysin_a' if index == 3 else 'lysin_b' if index == 20 else 'unresolved'
            genes.append({'protein_id': identifier, 'accession': acc, 'cds_index': index, 'start': index * 120,
                          'end': index * 120 + 90, 'strand': 1, 'translation': 'M' + 'A'*29, 'partial': 'false'})
            memberships.append({'family_identity': '0.4', 'protein_id': identifier, 'family_id': f'SYNFAM_{index}', 'prevalence': '1.0'})
            if category != 'unresolved': evidence.append({'protein_id': identifier, 'category': category,
                'annotation_category': category, 'confidence': 'low', 'source': 'genbank_annotation'})
        architectures.append({'accession': acc, 'scope': 'annotation', 'threshold': '4',
                              'architecture': 'separated', 'lysin_b': acc+'__cds0020'})
    return cfg, genomes, genes, evidence, memberships, architectures

@pytest.mark.parametrize('count, available', [(4, False), (5, True)])
def test_annotations_keep_low_confidence_and_cannot_lower_minimum_controls(monkeypatch, count, available):
    data = fixture()
    before = deepcopy(data)
    monkeypatch.setattr(sensitivity.nb, 'matched_controls', lambda *args: ([{'jaccard': .4, 'signed_order': .4} for _ in range(count)], count))
    pairs, _, _ = sensitivity.analyze(*data)
    assert data == before
    assert all(row['confidence'] == 'low' for row in data[3])
    assert all(row['evidence_scope'] == 'deposited_annotation_only' for row in pairs)
    assert all((row['status'] == 'descriptive_eligible') == available for row in pairs)
    assert all(row['delta'] == pytest.approx(.6) if available else row['delta'] == '' for row in pairs)

@pytest.mark.parametrize('failure', ['unresolved_architecture', 'near_identical_genomes'])
def test_annotation_assumptions_cannot_bypass_architecture_or_replication_gates(monkeypatch, failure):
    data = fixture()
    if failure == 'unresolved_architecture': data[-1][0]['architecture'] = 'unresolved'
    else: data[1][1]['similarity_group'] = data[1][0]['similarity_group']
    def forbidden(*args): raise AssertionError('Ineligible target must not reach matching')
    monkeypatch.setattr(sensitivity.nb, 'matched_controls', forbidden)
    pairs, _, _ = sensitivity.analyze(*data)
    assert all(row['status'] == 'unavailable' and row['delta'] == '' for row in pairs)

def test_descriptive_summary_weights_cluster_labels_instead_of_repeated_pairs():
    common = {'evidence_scope': 'deposited_annotation_only', 'family_identity': .4, 'window': 5,
              'architecture_threshold': 4, 'genome_subset': 'all', 'pair_architecture': 'separated',
              'metric': 'jaccard', 'status': 'descriptive_eligible'}
    rows = [{**common, 'cluster': 'SYN_A', 'delta': 1} for _ in range(3)] + [{**common, 'cluster': 'SYN_B', 'delta': 0}]
    result = sensitivity.summarize(rows)[0]
    assert result['n_available_pairs'] == 4 and result['n_cluster_labels'] == 2
    assert result['equal_cluster_mean_delta'] == .5
    assert result['formal_inference'] == 'not_performed'
