"""Synthetic-only replay-integrity, alias and scientific-export tests."""
import hashlib
import json

import pytest

from actinolysis.common import NA, sha256, write_json, write_tsv
from scripts import replay_exploratory_screen as replay
from scripts import run_exploratory_screen as runner


def test_replay_compares_scientific_bytes_and_rejects_any_difference(tmp_path):
    original, repeated = tmp_path / 'original', tmp_path / 'repeated'
    for root in (original, repeated):
        for name in replay.SCIENCE:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'SYNTHETIC TEST ONLY scientific artifact\n')
    assert len(replay.compare_science(original, repeated)) == len(replay.SCIENCE)
    (repeated / 'data/processed/cds.tsv').write_bytes(b'SYNTHETIC CHANGED CDS')
    with pytest.raises(ValueError, match='replay differs'):
        replay.compare_science(original, repeated)


@pytest.fixture
def receipted(tmp_path):
    write_json(tmp_path / 'execution_inputs.json', {'synthetic_test_only': True})
    for paths in runner.OUTPUTS.values():
        for name in paths:
            path = tmp_path / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'SYNTHETIC TEST ONLY\n')
    for name in ('data/manifests/eligibility_snapshot.tsv', 'data/manifests/phagesdb_snapshot.json',
                 'data/manifests/acquisition_ledger.tsv', 'data/raw/genbank/SYN000000.1.gb',
                 'data/raw/references/SYN000001.1.gb'):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'SYNTHETIC TEST ONLY\n')
    write_tsv(tmp_path / 'data/manifests/intermediate_inputs.tsv', [
        {'path': 'data/raw/genbank/SYN000000.1.gb'}])
    write_json(tmp_path / 'results/reports/exploratory_execution.json', {
        'status': 'exploratory_complete_awaiting_replay', 'completed_stages': [1, 2, 3, 4, 5],
        'h2_analysis_performed': False})
    receipts = {str(number): {'status': 'exploratory_stage_complete',
        'execution_inputs_sha256': sha256(tmp_path / 'execution_inputs.json'),
        'outputs': {name: sha256(tmp_path / name) for name in runner.stage_outputs(tmp_path, number)}}
        for number in range(1, 6)}
    write_json(tmp_path / 'results/reports/exploratory_stage_receipts.json', receipts)
    return tmp_path, receipts


def test_replay_checks_every_receipted_output(receipted):
    root, _ = receipted
    assert 'data/processed/reference_hits.tsv' in replay.verify_stage_receipts(root)
    (root / 'data/processed/reference_hits.tsv').write_bytes(b'SYNTHETIC ALTERED HITS')
    with pytest.raises(ValueError, match='output checksum mismatch'):
        replay.verify_stage_receipts(root)


def test_replay_requires_all_five_receipts(receipted):
    root, receipts = receipted
    del receipts['3']
    write_json(root / 'results/reports/exploratory_stage_receipts.json', receipts)
    with pytest.raises(ValueError, match='Exactly five'):
        replay.verify_stage_receipts(root)


def test_replay_rejects_receipt_path_escape_before_reading_external_file(receipted):
    root, receipts = receipted
    receipts['3']['outputs']['../outside'] = 'a' * 64
    write_json(root / 'results/reports/exploratory_stage_receipts.json', receipts)
    with pytest.raises(ValueError, match='Invalid exploratory stage 3'):
        replay.verify_stage_receipts(root)


def test_exact_alias_flags_use_full_sequence_hashes_not_similarity_scores():
    cds = {'exact': {'accession': 'SYN000001.2', 'translation': 'MAA'},
           'different': {'accession': 'SYN000002.1', 'translation': 'MAB'},
           'missing': {'accession': 'SYN000003.1', 'translation': ''}}
    refs = {'SYNTHETIC_REF': {'resolved_source_accession': 'SYN000001.1',
                             'protein_sha256': hashlib.sha256(b'MAA').hexdigest()}}
    flags = replay.alias_flags(cds, refs)
    assert flags['exact']['reference_source_accession_match'] is True
    assert flags['exact']['exact_reference_sequence_alias'] is True
    assert flags['different']['exact_reference_sequence_alias'] is False
    assert flags['missing']['protein_sha256'] == NA


def test_public_science_keeps_no_hit_denominators_and_discards_private_fields(tmp_path):
    accession = 'SYN000001.1'
    private = 'PRIVATE_TEST_ONLY_DO_NOT_EXPORT'
    manifest = {'accession': accession, 'phage': 'SYNTHETIC_TEST_ONLY', 'host_genus': 'Gordonia',
        'cluster': 'SYNTHETIC_ONLY', 'subcluster': NA, 'source_url': 'https://example.invalid/source',
        'sha256': 'a' * 64, 'host': private, 'path': private}
    write_tsv(tmp_path / 'data/manifests/intermediate_inputs.tsv', [manifest])
    write_tsv(tmp_path / 'data/processed/genomes.tsv', [{**manifest, 'inclusion_status': 'included',
        'exclusion_reason': 'not_applicable', 'source_host': private, 'cds_count': 4}])
    sequences = ['MAA', 'MCC', 'MNN', '']
    cds = [{'protein_id': f'{accession}__cds{i}', 'protein_accession': f'SYNTHETIC_ONLY_{i}.1',
        'accession': accession, 'cds_index': i, 'start': i * 9, 'end': i * 9 + 9,
        'location_parts': json.dumps([[i * 9, i * 9 + 9]]), 'strand': 1,
        'partial': 'false' if seq else 'true', 'product': private, 'note': private,
        'translation': seq, 'translation_status': 'verified' if seq else 'missing'}
        for i, seq in enumerate(sequences)]
    write_tsv(tmp_path / 'data/processed/cds.tsv', cds)
    (tmp_path / 'data/processed/proteins.faa').write_text(''.join(
        f'>{row["protein_id"]}\n{row["translation"]}\n' for row in cds if row['translation']))
    refs = [{'reference_id': rid, 'category': category, 'evidence_type': 'experimental',
        'resolved_source_accession': source, 'resolved_protein_accession': f'SYNTHETIC_{rid}.1',
        'source_url': 'https://example.invalid/source', 'reference_url': 'https://example.invalid/assay',
        'source_sha256': 'b' * 64, 'protein_sha256': hashlib.sha256(seq.encode()).hexdigest(),
        'length_aa': len(seq), 'description': private}
        for rid, category, source, seq in [('B', 'lysin_b', accession, 'MAA'),
                                         ('negative', 'negative_control', 'SYN000099.1', 'MNN')]]
    write_tsv(tmp_path / 'data/manifests/reference_manifest.tsv', refs)
    (tmp_path / 'data/processed/reference_hits.tsv').write_text(
        f'{cds[0]["protein_id"]}\tB\t1\t3\t1\t1\t0\t100\n'
        f'{cds[2]["protein_id"]}\tnegative\t1\t3\t1\t1\t0\t100\n')
    evidence = [{'protein_id': cds[0]['protein_id'], 'accession': accession, 'category': 'lysin_b',
        'reference_id': 'B', 'source': 'curated_reference_homology', 'confidence': 'high',
        'reference_evidence_type': 'experimental', 'annotation_note': private, 'product': private}]
    write_tsv(tmp_path / 'data/processed/lysis_evidence.tsv', evidence)
    architecture = {**manifest, 'architecture': 'unresolved', 'scope': 'high_moderate', 'threshold': 4}
    for name in ('architectures.tsv', 'architecture_sensitivity.tsv'):
        write_tsv(tmp_path / 'results/tables' / name, [architecture])
    for name in ('exclusions.tsv', 'near_identical_pairs.tsv', 'architecture_pairs.tsv', 'annotation_coverage.tsv'):
        write_tsv(tmp_path / 'results/tables' / name, [], ['synthetic_test_only'])
    cfg = {'identification': {'min_identity': .4, 'min_query_coverage': .7,
                             'min_target_coverage': .7, 'max_evalue': 1e-6}}
    tables, summary = replay.public_tables(tmp_path, [{'accession': accession, 'pilot': 'true'}], cfg)
    statuses = {row['protein_id']: row for row in tables['protein_detection_status.tsv'][0]}
    assert len(statuses) == 4
    assert statuses[cds[1]['protein_id']]['detection_status'] == 'unresolved_no_reported_reference_hit'
    assert statuses[cds[3]['protein_id']]['detection_status'] == 'not_searched_no_eligible_translation'
    assert statuses[cds[0]['protein_id']]['exact_reference_sequence_alias'] == 'true'
    assert summary['counts']['cds'] == 4 and summary['counts']['searched_proteins'] == 3
    assert summary['counts']['supported_B_genomes'] == 1
    assert summary['h2_analysis_performed'] is False
    serialized = json.dumps({'tables': tables, 'summary': summary})
    assert private not in serialized
    for _, fields in tables.values():
        assert not {'translation', 'note', 'annotation_note', 'product', 'host', 'source_host', 'path'} & set(fields)
