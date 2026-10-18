#!/usr/bin/env python3
"""Compare displaced LysB candidates with the frozen assay-reference panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
ROOT = Path(__file__).resolve().parents[1]

OUT = ROOT / 'results/displaced_lysB_candidates'
REFERENCE_FASTA = ROOT / 'results/lysterase_update/reference_validation/queries.faa'
MMSEQS = ROOT / '.tools/mmseqs/bin/mmseqs'
CANDIDATE_FASTA = ROOT / '.experiments/rhodococcus_followup/candidates.faa'
SOURCE_REPORT = OUT / 'audit_report.json'
SEARCH_WORK = ROOT / '.experiments/rhodococcus_followup/reference_search_work'
FORMAT = 'query,target,fident,alnlen,qcov,tcov,evalue,bits'
HIT_FIELDS = ['phage', 'source_accession', 'candidate_protein', 'reference_id', 'identity', 'alignment_length',
              'query_coverage', 'reference_coverage', 'evalue', 'bit_score',
              'passes_existing_identity_coverage_evalue_thresholds']
MIN_IDENTITY = 0.40
MIN_QUERY_COVERAGE = 0.70
MIN_TARGET_COVERAGE = 0.70
MAX_EVALUE = 1e-6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate-fasta', type=Path, default=CANDIDATE_FASTA)
    parser.add_argument('--source-report', type=Path, default=SOURCE_REPORT)
    parser.add_argument('--reference-fasta', type=Path, default=REFERENCE_FASTA)
    parser.add_argument('--out-dir', type=Path, default=OUT)
    parser.add_argument('--work-dir', type=Path, default=SEARCH_WORK)
    parser.add_argument('--mmseqs', type=Path, default=MMSEQS)
    args = parser.parse_args()
    if not args.candidate_fasta.is_file() or not args.source_report.is_file() or not args.reference_fasta.is_file() or not args.mmseqs.is_file():
        raise FileNotFoundError('Candidate FASTA, source audit report, frozen reference FASTA and MMseqs executable are required')

    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)
    config_path = ROOT / 'config/analysis.yaml'
    config_text = config_path.read_text(encoding='utf-8')
    configured = {}
    for key in ('min_identity', 'min_query_coverage', 'min_target_coverage', 'max_evalue'):
        match = re.search(rf'(?m)^  {key}:\s*([0-9.eE+-]+)\s*$', config_text)
        if not match:
            raise ValueError(f'Could not read {key} from the identification configuration')
        configured[key] = float(match.group(1))
    expected = {
        'min_identity': MIN_IDENTITY, 'min_query_coverage': MIN_QUERY_COVERAGE,
        'min_target_coverage': MIN_TARGET_COVERAGE, 'max_evalue': MAX_EVALUE,
    }
    if configured != expected:
        raise ValueError(f'Configured identification thresholds changed: {configured}')
    raw_hits = work / 'all_reference_hits.tsv'
    temporary = work / 'mmseqs_tmp'
    command = [str(args.mmseqs), 'easy-search', str(args.candidate_fasta), str(args.reference_fasta),
               str(raw_hits), str(temporary), '--threads', '4',
               '-s', '7.5', '-k', '5', '--alignment-mode', '3', '--max-seqs', '100',
               '--min-seq-id', '0.25', '-c', '0.5', '--cov-mode', '0',
               '--split-memory-limit', '1G', '-e', '0.001',
               '--format-output', FORMAT]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f'MMseqs failed ({result.returncode}): {result.stderr[-2000:]}')

    hits = []
    with raw_hits.open(encoding='utf-8') as handle:
        for record in csv.DictReader(handle, fieldnames=FORMAT.split(','), delimiter='\t'):
            if 'lysb' not in record['target'].lower():
                continue
            identity = float(record['fident'])
            qcov = float(record['qcov'])
            tcov = float(record['tcov'])
            evalue = float(record['evalue'])
            query_parts = record['query'].split('|')
            hits.append({
                'phage': query_parts[0], 'source_accession': query_parts[1], 'candidate_protein': query_parts[2],
                'reference_id': record['target'], 'identity': identity, 'alignment_length': int(record['alnlen']),
                'query_coverage': qcov, 'reference_coverage': tcov, 'evalue': record['evalue'],
                'bit_score': float(record['bits']),
                'passes_existing_identity_coverage_evalue_thresholds': (
                    identity >= MIN_IDENTITY and qcov >= MIN_QUERY_COVERAGE and
                    tcov >= MIN_TARGET_COVERAGE and evalue <= MAX_EVALUE),
            })
    hits.sort(key=lambda row: (row['phage'], float(row['evalue']), row['reference_id']))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.out_dir / 'reference_panel_hits.tsv'
    with table_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(hits)
    by_phage = {}
    for row in hits:
        by_phage.setdefault(row['phage'], []).append(row)
    best = {}
    for phage, rows in by_phage.items():
        row = min(rows, key=lambda hit: (float(hit['evalue']), -hit['bit_score']))
        best[phage] = {key: row[key] for key in
                       ('reference_id', 'identity', 'query_coverage', 'reference_coverage', 'evalue',
                        'passes_existing_identity_coverage_evalue_thresholds')}

    version = subprocess.run([str(args.mmseqs), 'version'], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    report = {
        'status': 'pass', 'candidate_sequences_searched': 3,
        'candidate_sequences_with_LysB_reference_hits': len(by_phage),
        'candidate_sequences_meeting_existing_LysB_assignment_thresholds': sum(
            any(row['passes_existing_identity_coverage_evalue_thresholds'] for row in rows)
            for rows in by_phage.values()),
        'existing_thresholds': {
            'minimum_identity': MIN_IDENTITY,
            'minimum_query_coverage': MIN_QUERY_COVERAGE,
            'minimum_target_coverage': MIN_TARGET_COVERAGE,
            'maximum_evalue': MAX_EVALUE,
        },
        'analysis_config_sha256': sha256(config_path),
        'best_LysB_reference_hit_per_candidate': best,
        'mmseqs_version': version,
        'mmseqs_command': ['mmseqs', 'easy-search', '<candidate_fasta>', '<reference_fasta>',
                           '<output_tsv>', '<temporary_directory>'] + command[6:],
        'source_records_sha256': json.loads(args.source_report.read_text(encoding='utf-8'))['source_records_sha256'],
        'reference_fasta_sha256': sha256(args.reference_fasta),
        'candidate_fasta_sha256': sha256(args.candidate_fasta),
        'script_sha256': sha256(Path(__file__).resolve()), 'table_sha256': sha256(table_path),
        'interpretation': 'Sub-threshold similarity does not upgrade candidate function; these matches are exploratory comparisons against the frozen assayed LysB panel.',
    }
    (args.out_dir / 'reference_hit_report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
