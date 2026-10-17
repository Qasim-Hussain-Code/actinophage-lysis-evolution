"""Standalone supplement tests with synthetic scores, sequences and receipts only."""
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

from actinolysis.common import NA, read_tsv, sha256, write_tsv

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/archive_cohort_hits.py'
SPEC = importlib.util.spec_from_file_location('archive_cohort_hits', SCRIPT)
archive_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(archive_module)


def prepared(tmp_path):
    outer = tmp_path / 'project'; root = outer / 'isolated'; out = outer / 'supplement'
    cfg = dict(mode='pilot', seed=1, sample_sizes={'pilot': 1}, sampling={},
               identification={'min_identity': .4, 'min_query_coverage': .7, 'min_target_coverage': .7, 'max_evalue': 1e-6},
               architecture={}, neighborhood={}, statistics={},
               resources={'max_project_bytes': 10_000_000_000, 'reserve_bytes': 0})
    for base in [root, outer]:
        (base / 'config').mkdir(parents=True, exist_ok=True)
        (base / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg))
    ids = [f'SYNTHETIC_TEST_ONLY:cds{i:04}' for i in range(6)]
    cds = [dict(protein_id=pid, protein_accession=f'SYNTHETIC_TEST_ONLY_{index}.1', accession='SYNTHETIC_TEST_ONLY.1',
                cds_index=index, start=index * 30, end=index * 30 + 9, strand=1,
                location_parts=json.dumps([[index * 30, index * 30 + 9]]),
                partial='true' if index in {1, 4} else 'false',
                product='Lysin B' if index == 3 else 'synthetic test protein', note='')
           for index, pid in enumerate(ids)]
    write_tsv(root / 'data/processed/cds.tsv', cds)
    write_tsv(root / 'data/processed/genomes.tsv', [{'accession': 'SYNTHETIC_TEST_ONLY.1', 'cds_count': 6}])
    (root / 'data/processed/proteins.faa').write_text(''.join(f'>{pid}\nMAA\n' for index, pid in enumerate(ids) if index != 4))
    refs = [dict(reference_id=rid, category=category, evidence_type=evidence,
                 resolved_source_accession=f'SYNTHETIC_{rid}.1', resolved_protein_accession=f'SYNTHETIC_{rid}_PROTEIN.1',
                 reference_url='https://example.invalid/synthetic_test_only')
            for rid, category, evidence in [('exp', 'lysin_b', 'experimental'),
                                            ('profile', 'lysin_b', 'published_profile_prediction'),
                                            ('generic', 'lysin_b', 'experimental_generic_lipolytic'),
                                            ('prediction', 'lysin_a', 'published_prediction'),
                                            ('neg', 'negative_control', 'experimental')]]
    write_tsv(root / 'data/manifests/reference_manifest.tsv', refs)
    (root / 'data/processed/references.faa').write_text(''.join(f'>{row["reference_id"]}\nMAA\n' for row in refs))
    hits = [(ids[0], 'exp', '.8'), (ids[0], 'profile', '.8'), (ids[0], 'neg', '1'),
            (ids[1], 'generic', '1'), (ids[2], 'neg', '1'), (ids[5], 'prediction', '1')]
    (root / 'data/processed/reference_hits.tsv').write_text(''.join(f'{pid}\t{rid}\t{identity}\t3\t1\t1\t1e-20\t100\n' for pid, rid, identity in hits))
    evidence = [dict(protein_id=pid, reference_id=rid, category=category, source=source, confidence=confidence)
                for pid, rid, category, source, confidence in [
                    (ids[0], 'exp', 'lysin_b', 'curated_reference_homology', 'moderate'),
                    (ids[1], 'generic', 'lysin_b', 'curated_reference_homology', 'low'),
                    (ids[3], NA, 'lysin_b', 'genbank_annotation', 'low'),
                    (ids[5], 'prediction', 'lysin_a', 'curated_reference_homology', 'low')]]
    write_tsv(root / 'data/processed/lysis_evidence.tsv', evidence)
    source_paths = ['src/actinolysis/common.py', 'scripts/run_all.py', 'src/actinolysis/quality.py',
                    'src/actinolysis/identification.py', 'references/protein_controls.tsv']
    for relative in source_paths:
        path = root / relative; path.parent.mkdir(parents=True, exist_ok=True)
        if relative.endswith('identification.py'):
            path.write_bytes(Path(archive_module.identification.__file__).read_bytes())
        else:
            path.write_text('SYNTHETIC TEST ONLY SOURCE RECEIPT\n')
    receipts = {}
    for number, outputs in archive_module.REQUIRED_OUTPUTS.items():
        sources = [name for name in source_paths if name != ('src/actinolysis/identification.py' if number == '2' else 'src/actinolysis/quality.py')]
        if number == '2':
            sources.remove('references/protein_controls.tsv')
        receipts[number] = dict(outputs={name: sha256(root / name) for name in outputs},
                                source_sha256={name: sha256(root / name) for name in sources},
                                scientific_config_sha256=archive_module.scientific_fingerprint(cfg), seconds=1)
    receipt_path = root / 'results/reports/pilot_stage_receipts.json'
    receipt_path.parent.mkdir(parents=True)
    receipt_path.write_text(json.dumps(receipts))
    (root / 'results/reports/reference_search.json').write_text(json.dumps(dict(
        command=['SYNTHETIC_TEST_ONLY_MMSEQS_NOT_EXECUTED'],
        hit_file_sha256=sha256(root / 'data/processed/reference_hits.tsv'), reference_count=len(refs))))
    return outer, root, out, ids


def test_archive_preserves_all_hits_and_every_cds_without_promoting_calls(tmp_path):
    outer, root, out, ids = prepared(tmp_path)
    original_ledger = (root / 'data/processed/lysis_evidence.tsv').read_bytes()
    report = archive_module.archive(root, out, project_root=outer)
    hits = read_tsv(out / 'all_cohort_reference_hits.tsv')
    statuses = {row['protein_id']: row for row in read_tsv(out / 'protein_detection_status.tsv')}
    assert len(hits) == 6 and len(statuses) == 6
    assert (root / 'data/processed/lysis_evidence.tsv').read_bytes() == original_ledger
    p1 = {row['reference_id']: row for row in hits if row['protein_id'] == ids[0]}
    assert p1['exp']['selected_for_category'] == 'true'
    assert p1['profile']['accepted_under_frozen_rules'] == 'true' and p1['profile']['selected_for_category'] == 'false'
    assert p1['neg']['negative_control'] == 'true' and p1['neg']['accepted_under_frozen_rules'] == 'false'
    partial = next(row for row in hits if row['protein_id'] == ids[1])
    assert partial['confidence_before_partial_cap'] == 'moderate' and partial['score_confidence'] == 'low'
    assert partial['accepted_under_frozen_rules'] == 'false'
    assert statuses[ids[0]]['accepted_hit_count'] == '2' and statuses[ids[0]]['accepted_categories'] == 'lysin_b'
    assert statuses[ids[0]]['accepted_category_count'] == '1'
    assert statuses[ids[2]]['detection_status'] == 'unresolved_only_negative_control_hits'
    assert statuses[ids[3]]['detection_status'] == 'unresolved_no_reported_reference_hit'
    assert statuses[ids[3]]['selected_categories'] == 'lysin_b' and statuses[ids[3]]['supported_selected_categories'] == NA
    assert statuses[ids[4]]['detection_status'] == 'not_searched_no_eligible_translation'
    assert statuses[ids[5]]['detection_status'] == 'unresolved_reported_hits_below_criteria_or_evidence_tier'
    assert report['counts']['accepted_reference_hits'] == 2 and report['counts']['searched_proteins'] == 5
    assert report['archive_code_sha256'] == sha256(SCRIPT)
    assert all(sha256(out / name) == digest for name, digest in report['outputs_sha256'].items())
    before = {name: (out / name).read_bytes() for name in report['outputs_sha256']}
    archive_module.archive(root, out, project_root=outer)
    assert before == {name: (out / name).read_bytes() for name in before}


@pytest.mark.parametrize('relative', ['data/processed/reference_hits.tsv', 'data/processed/lysis_evidence.tsv',
                                      'data/processed/cds.tsv', 'src/actinolysis/identification.py'])
def test_changed_frozen_input_rejected_before_creating_output(tmp_path, relative):
    outer, root, out, _ = prepared(tmp_path)
    with (root / relative).open('a') as stream:
        stream.write('TAMPERED SYNTHETIC TEST ONLY\n')
    with pytest.raises(ValueError, match='frozen input changed'):
        archive_module.archive(root, out, project_root=outer)
    assert not out.exists()


def test_changed_input_preserves_already_published_archive(tmp_path):
    outer, root, out, _ = prepared(tmp_path)
    archive_module.archive(root, out, project_root=outer)
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    with (root / 'data/processed/reference_hits.tsv').open('a') as stream:
        stream.write('TAMPERED SYNTHETIC TEST ONLY\n')
    with pytest.raises(ValueError, match='frozen input changed'):
        archive_module.archive(root, out, project_root=outer)
    assert before == {path.name: path.read_bytes() for path in out.iterdir()}


def test_configuration_change_is_not_silently_rescored(tmp_path):
    outer, root, out, _ = prepared(tmp_path)
    path = root / 'config/analysis.yaml'; cfg = yaml.safe_load(path.read_text())
    cfg['identification']['min_identity'] = .2; path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match='scientific configuration changed'):
        archive_module.archive(root, out, project_root=outer)
    assert not out.exists()


def test_selected_assignment_must_match_reported_hit_and_recomputed_confidence(tmp_path):
    outer, root, out, _ = prepared(tmp_path)
    path = root / 'data/processed/lysis_evidence.tsv'
    rows = read_tsv(path); rows[0]['confidence'] = 'high'; write_tsv(path, rows)
    receipt_path = root / 'results/reports/pilot_stage_receipts.json'
    receipts = json.loads(receipt_path.read_text()); receipts['4']['outputs']['data/processed/lysis_evidence.tsv'] = sha256(path)
    receipt_path.write_text(json.dumps(receipts))
    with pytest.raises(ValueError, match='disagrees with frozen selected evidence'):
        archive_module.archive(root, out, project_root=outer)
    assert not out.exists()


def test_changed_published_table_is_not_silently_replaced(tmp_path):
    outer, root, out, _ = prepared(tmp_path)
    archive_module.archive(root, out, project_root=outer)
    path = out / 'protein_detection_status.tsv'
    path.write_text('TAMPERED SYNTHETIC TEST ONLY\n')
    before = {path.name: path.read_bytes() for path in out.iterdir()}
    with pytest.raises(ValueError, match='Existing supplement output changed'):
        archive_module.archive(root, out, project_root=outer)
    assert before == {path.name: path.read_bytes() for path in out.iterdir()}


def test_search_execution_metadata_must_match_the_receipted_hit_file(tmp_path):
    outer, root, out, _ = prepared(tmp_path)
    path = root / 'results/reports/reference_search.json'; report = json.loads(path.read_text())
    report['hit_file_sha256'] = '0' * 64; path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match='Search execution report disagrees'):
        archive_module.archive(root, out, project_root=outer)
    assert not out.exists()
