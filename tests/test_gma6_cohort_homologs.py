import pytest

from scripts.audit_gma6_cohort_homologs import independent_cluster_count, validate_fixed_cohort


def cohort_row(index, cds_count=123, inclusion='included'):
    return {'phage': f'phage-{index}', 'cds_count': str(cds_count), 'inclusion_status': inclusion}


def test_gma6_homology_audit_keeps_the_original_fixed_cohort():
    rows = [cohort_row(index) for index in range(20)]
    rows.append(cohort_row(21, inclusion='excluded'))

    assert len(validate_fixed_cohort(rows)) == 20


def test_gma6_homology_audit_rejects_cohort_expansion():
    rows = [cohort_row(index) for index in range(21)]

    with pytest.raises(ValueError, match='fixed 20-genome cohort'):
        validate_fixed_cohort(rows)


def test_gma6_homology_audit_checks_frozen_cds_count():
    rows = [cohort_row(index) for index in range(19)] + [cohort_row(20, cds_count=122)]

    with pytest.raises(ValueError, match='2,460 fixed-cohort CDS'):
        validate_fixed_cohort(rows)


def test_two_homologs_in_one_cluster_count_as_one_replication():
    hits = [{'phage': 'Jollymon', 'cluster': 'CP'}, {'phage': 'Stillion', 'cluster': 'CP'}]

    assert independent_cluster_count(hits) == 1
