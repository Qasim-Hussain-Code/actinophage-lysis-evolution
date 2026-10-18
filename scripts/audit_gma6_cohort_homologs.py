#!/usr/bin/env python3
"""Search the fixed pilot CDS set for GMA6 PhageGH4185 homologs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from actinolysis.common import write_json, write_tsv

MANIFEST = ROOT / 'results/tables/dataset_manifest.tsv'
LYSIS_EVIDENCE = ROOT / 'results/lysterase_update/panel/results/tables/lysis_evidence.tsv'
GMA6_SOURCE = ROOT / '.experiments/gma6/raw/KR063280.1.gb'
OUT = ROOT / 'results/gma6_arabinanase'
WORK = ROOT / '.experiments/gma6/cohort_search'
MMSEQS = ROOT / '.tools/mmseqs/bin/mmseqs'
QUERY_PROTEIN = 'AKL88326.1'
QUERY_LOCUS = 'GMA6_45'
QUERY_HEADER = f'{QUERY_LOCUS}|{QUERY_PROTEIN}|KR063280.1'
FORMAT = 'query,target,fident,alnlen,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,evalue,bits'
SEARCH_PRESET = {
    'minimum_sequence_identity': 0.25, 'minimum_query_and_target_coverage': 0.50,
    'maximum_evalue': 0.001, 'sensitivity': 7.5, 'alignment_mode': 3,
}
FIELDS = [
    'query_protein', 'phage', 'accession', 'host_genus', 'cluster', 'locus_tag', 'protein_id',
    'deposited_product', 'cds_index', 'start_1based', 'end_1based', 'downstream_locus_tag',
    'downstream_product', 'identity_to_GMA6', 'alignment_length', 'query_start_1based',
    'query_end_1based', 'target_start_1based', 'target_end_1based', 'query_length', 'target_length',
    'query_coverage', 'target_coverage', 'evalue', 'bit_score', 'current_pipeline_category',
    'current_pipeline_confidence', 'current_reference_id', 'current_reference_evidence_type',
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_fixed_cohort(manifest_rows):
    cohort = [row for row in manifest_rows if row['inclusion_status'] == 'included']
    if len(cohort) != 20:
        raise ValueError(f'Expected the fixed 20-genome cohort, found {len(cohort)}')
    total_cds = sum(int(row['cds_count']) for row in cohort)
    if total_cds != 2460:
        raise ValueError(f'Expected 2,460 fixed-cohort CDS, found {total_cds}')
    return cohort


def independent_cluster_count(hits) -> int:
    return len({row['cluster'] for row in hits})


def load_inputs(manifest_path: Path, gma6_path: Path):
    manifest_rows = list(csv.DictReader(manifest_path.open(encoding='utf-8'), delimiter='\t'))
    cohort = validate_fixed_cohort(manifest_rows)

    query_record = SeqIO.read(gma6_path, 'genbank')
    if query_record.id != 'KR063280.1':
        raise ValueError(f'Unexpected GMA6 record: {query_record.id}')
    query_features = [feature for feature in query_record.features if feature.type == 'CDS']
    query_matches = [feature for feature in query_features
                     if feature.qualifiers.get('protein_id', [''])[0] == QUERY_PROTEIN]
    if len(query_matches) != 1 or query_matches[0].qualifiers.get('locus_tag', [''])[0] != QUERY_LOCUS:
        raise ValueError('The pinned PhageGH4185 CDS could not be resolved uniquely')
    query_translation = query_matches[0].qualifiers.get('translation', [''])[0]
    if len(query_translation) != 1087 or '*' in query_translation:
        raise ValueError('Unexpected PhageGH4185 translation length or internal stop')

    sequences = []
    protein_features = {}
    source_hashes = {}
    for row in cohort:
        source = ROOT / row['path']
        source_hash = sha256(source)
        if source_hash != row['sha256']:
            raise ValueError(f'Frozen source checksum mismatch: {row["accession"]}')
        record = SeqIO.read(source, 'genbank')
        if record.id != row['accession']:
            raise ValueError(f'Frozen source version mismatch: {row["accession"]}')
        cds = [feature for feature in record.features if feature.type == 'CDS']
        if len(cds) != int(row['cds_count']):
            raise ValueError(f'Frozen CDS count mismatch: {row["phage"]}')
        source_hashes[row['accession']] = source_hash
        for index, feature in enumerate(cds, 1):
            translation = feature.qualifiers.get('translation', [''])[0]
            if not translation or '*' in translation:
                raise ValueError(f'Invalid deposited CDS translation: {row["phage"]}/{index}')
            locus = feature.qualifiers.get('locus_tag', [f'cds{index}'])[0]
            protein_id = feature.qualifiers.get('protein_id', ['unresolved'])[0]
            header = f'{row["phage"]}|{row["accession"]}|{locus}|{protein_id}'
            sequences.append(f'>{header}\n{translation}\n')
            protein_features[header] = (row, index, feature, cds)
    if len(protein_features) != 2460:
        raise ValueError(f'Expected 2,460 cohort CDS, found {len(protein_features)}')
    return cohort, source_hashes, query_translation, sequences, protein_features


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--lysis-evidence', type=Path, default=LYSIS_EVIDENCE)
    parser.add_argument('--gma6-source', type=Path, default=GMA6_SOURCE)
    parser.add_argument('--mmseqs', type=Path, default=MMSEQS)
    parser.add_argument('--out-dir', type=Path, default=OUT)
    parser.add_argument('--work-dir', type=Path, default=WORK)
    parser.add_argument('--precomputed-hits', type=Path,
                        help='raw MMseqs TSV for split Windows/WSL execution')
    parser.add_argument('--replay-hits', type=Path,
                        help='optional independent raw search result to verify byte-identical replay')
    parser.add_argument('--mmseqs-version', default='not-recorded')
    parser.add_argument('--extract-only', action='store_true')
    args = parser.parse_args()
    for path in (args.manifest, args.lysis_evidence, args.gma6_source):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.precomputed_hits is None and not args.mmseqs.is_file():
        raise FileNotFoundError(args.mmseqs)

    cohort, source_hashes, query_sequence, sequences, feature_map = load_inputs(args.manifest, args.gma6_source)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    query_fasta = args.work_dir / 'PhageGH4185.faa'
    cohort_fasta = args.work_dir / 'cohort_CDS.faa'
    query_fasta.write_text(f'>{QUERY_HEADER}\n{query_sequence}\n', encoding='ascii')
    cohort_fasta.write_text(''.join(sequences), encoding='ascii')
    if args.extract_only:
        return

    raw_hits = args.work_dir / 'cohort_hits.tsv'
    tmp = args.work_dir / 'mmseqs_tmp'
    command = [str(args.mmseqs), 'easy-search', str(query_fasta), str(cohort_fasta), str(raw_hits), str(tmp),
               '--threads', '4', '-s', '7.5', '-k', '5', '--alignment-mode', '3', '--max-seqs', '100',
               '--min-seq-id', '0.25', '-c', '0.5', '--cov-mode', '0', '--split-memory-limit', '1G',
               '-e', '0.001', '--format-output', FORMAT]
    if args.precomputed_hits is not None:
        raw_hits = args.precomputed_hits
        if not raw_hits.is_file():
            raise FileNotFoundError(raw_hits)
        version = args.mmseqs_version
    else:
        search = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        if search.returncode:
            raise RuntimeError(f'MMseqs failed ({search.returncode}): {search.stderr[-2000:]}')
        version = subprocess.run([str(args.mmseqs), 'version'], cwd=ROOT, capture_output=True,
                                 text=True, check=True).stdout.strip()

    evidence_by_accession = {}
    with args.lysis_evidence.open(encoding='utf-8', newline='') as handle:
        for evidence in csv.DictReader(handle, delimiter='\t'):
            evidence_by_accession[evidence['protein_accession']] = evidence

    output_rows = []
    with raw_hits.open(encoding='utf-8') as handle:
        for hit in csv.DictReader(handle, fieldnames=FORMAT.split(','), delimiter='\t'):
            target_header = hit['target']
            if target_header not in feature_map:
                raise ValueError(f'Hit target is absent from frozen cohort FASTA: {target_header}')
            row, cds_index, feature, cds = feature_map[target_header]
            protein_id = feature.qualifiers.get('protein_id', [''])[0]
            evidence = evidence_by_accession.get(protein_id, {})
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            downstream = cds[cds_index] if cds_index < len(cds) else None
            values = [hit['fident'], hit['qcov'], hit['tcov'], hit['evalue'], hit['bits']]
            identity, qcov, tcov, evalue, bits = map(float, values)
            output_rows.append({
                'query_protein': QUERY_PROTEIN, 'phage': row['phage'], 'accession': row['accession'],
                'host_genus': row['host_genus'], 'cluster': row['cluster'],
                'locus_tag': feature.qualifiers.get('locus_tag', [f'cds{cds_index}'])[0],
                'protein_id': protein_id, 'deposited_product': feature.qualifiers.get('product', [''])[0],
                'cds_index': cds_index, 'start_1based': start, 'end_1based': end,
                'downstream_locus_tag': (downstream.qualifiers.get('locus_tag', [''])[0] if downstream else ''),
                'downstream_product': (downstream.qualifiers.get('product', [''])[0] if downstream else ''),
                'identity_to_GMA6': identity, 'alignment_length': int(hit['alnlen']),
                'query_start_1based': int(hit['qstart']), 'query_end_1based': int(hit['qend']),
                'target_start_1based': int(hit['tstart']), 'target_end_1based': int(hit['tend']),
                'query_length': int(hit['qlen']), 'target_length': int(hit['tlen']),
                'query_coverage': qcov, 'target_coverage': tcov, 'evalue': hit['evalue'], 'bit_score': bits,
                'current_pipeline_category': evidence.get('category', ''),
                'current_pipeline_confidence': evidence.get('confidence', ''),
                'current_reference_id': evidence.get('reference_id', ''),
                'current_reference_evidence_type': evidence.get('reference_evidence_type', ''),
            })
    output_rows.sort(key=lambda item: (float(item['evalue']), -item['bit_score'], item['accession']))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.out_dir / 'cohort_homologs.tsv'
    write_tsv(table_path, output_rows, FIELDS)

    clusters = sorted({row['cluster'] for row in output_rows})
    report = {
        'status': 'pass', 'fixed_cohort_genomes': len(cohort),
        'fixed_cohort_CDS_translations': len(feature_map), 'source_record_sha256': source_hashes,
        'cohort_manifest_sha256': sha256(args.manifest), 'gma6_source_sha256': sha256(args.gma6_source),
        'lysis_evidence_sha256': sha256(args.lysis_evidence),
        'query_fasta_sha256': sha256(query_fasta), 'cohort_CDS_fasta_sha256': sha256(cohort_fasta),
        'query_protein': QUERY_PROTEIN, 'query_locus_tag': QUERY_LOCUS,
        'MMseqs_search_hits': len(output_rows), 'hit_phages': sorted({row['phage'] for row in output_rows}),
        'hit_clusters': clusters, 'independent_cluster_replication': independent_cluster_count(output_rows),
        'search_preset': SEARCH_PRESET,
        'mmseqs_version': version,
        'mmseqs_command': ['mmseqs', 'easy-search', '<query_fasta>', '<cohort_CDS_fasta>',
                           '<output_tsv>', '<temporary_directory>'] + command[6:],
        'raw_mmseqs_hits_sha256': sha256(raw_hits), 'cohort_homolog_table_sha256': sha256(table_path),
        'script_sha256': sha256(Path(__file__).resolve()),
        'interpretation': (
            'The fixed pilot contains near-full-length homologs in the same Gordonia cluster. Their existing LysA '
            'assignments derive from profile-prediction homology; direct arabinanase activity is established '
            'only for GMA6 PhageGH4185. This prioritizes a same-cluster substrate-specific assay but supplies '
            'no independent lineage replication, validated LysA/B pair, or native lysis result.'),
    }
    write_json(args.out_dir / 'cohort_homolog_report.json', report)
    if args.replay_hits is not None:
        if not args.replay_hits.is_file():
            raise FileNotFoundError(args.replay_hits)
        replay = {
            'status': 'pass' if sha256(raw_hits) == sha256(args.replay_hits) else 'fail',
            'raw_result_sha256': sha256(raw_hits), 'independent_replay_sha256': sha256(args.replay_hits),
            'byte_identical': raw_hits.read_bytes() == args.replay_hits.read_bytes(),
            'scope': 'Independent MMseqs2 run against the same 20-genome cohort FASTA using a fresh output and temporary directory.',
        }
        write_json(args.out_dir / 'cohort_replay_validation.json', replay)
        if replay['status'] != 'pass':
            raise ValueError('Independent cohort MMseqs2 replay differs from the recorded raw result')


if __name__ == '__main__':
    main()
