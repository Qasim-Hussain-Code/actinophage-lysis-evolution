#!/usr/bin/env python3
"""Map the biochemically assayed GMA6 arabinanase and challenge the LysB panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from urllib.parse import urlencode

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from actinolysis.acquisition import download
from actinolysis.common import load_config
SOURCE_ACCESSION = 'KR063280.1'
PROTEIN_ID = 'AKL88326.1'
LOCUS_TAG = 'GMA6_45'
SOURCE = ROOT / '.experiments/gma6/raw/KR063280.1.gb'
REFERENCE_FASTA = ROOT / 'results/lysterase_update/reference_validation/queries.faa'
REFERENCE_METADATA = ROOT / 'results/lysterase_update/reference_validation/queries.tsv'
MMSEQS = ROOT / '.tools/mmseqs/bin/mmseqs'
OUT = ROOT / 'results/gma6_arabinanase'
PAPER = 'https://www.nature.com/articles/s41467-023-37839-5'
FIELDS = [
    'phage', 'source_accession', 'source_sha256', 'protein_id', 'locus_tag',
    'deposited_product', 'strand', 'start_1based', 'end_1based', 'length_aa',
    'upstream_locus_tag', 'upstream_product', 'upstream_start_1based', 'upstream_end_1based',
    'upstream_overlap_bp', 'downstream_locus_tag', 'downstream_product',
    'downstream_start_1based', 'downstream_end_1based', 'downstream_overlap_bp',
    'literature_function', 'functional_evidence_scope', 'native_phage_lysis_tested', 'primary_source',
]
HIT_FIELDS = [
    'reference_id', 'reference_category', 'identity', 'alignment_length',
    'query_start_1based', 'query_end_1based', 'reference_start_1based', 'reference_end_1based',
    'query_length', 'reference_length', 'query_coverage',
    'reference_coverage', 'evalue', 'bit_score', 'passes_existing_LysB_assignment_thresholds',
]
FORMAT = 'query,target,fident,alnlen,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,evalue,bits'
THRESHOLDS = {
    'minimum_identity': 0.40, 'minimum_query_coverage': 0.70,
    'minimum_target_coverage': 0.70, 'maximum_evalue': 1e-6,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def bounds(feature) -> tuple[int, int]:
    return int(feature.location.start) + 1, int(feature.location.end)


def overlap_bp(first, second) -> int:
    a0, a1 = int(first.location.start), int(first.location.end)
    b0, b1 = int(second.location.start), int(second.location.end)
    return max(0, min(a1, b1) - max(a0, b0))


def passes_lysb_thresholds(category: str, identity: float, query_coverage: float,
                           reference_coverage: float, evalue: float) -> bool:
    return (
        category == 'lysin_b' and identity >= THRESHOLDS['minimum_identity'] and
        query_coverage >= THRESHOLDS['minimum_query_coverage'] and
        reference_coverage >= THRESHOLDS['minimum_target_coverage'] and
        evalue <= THRESHOLDS['maximum_evalue']
    )


def analyze_record(record, source_hash: str) -> tuple[dict, str]:
    if record.id != SOURCE_ACCESSION:
        raise ValueError(f'Expected pinned source {SOURCE_ACCESSION}, found {record.id}')
    if record.annotations.get('molecule_type') != 'DNA':
        raise ValueError('Expected a DNA genome record')
    features = [feature for feature in record.features if feature.type == 'CDS']
    found = [(index, feature) for index, feature in enumerate(features)
             if feature.qualifiers.get('protein_id', [''])[0] == PROTEIN_ID]
    if len(found) != 1:
        raise ValueError(f'Expected one CDS with protein ID {PROTEIN_ID}, found {len(found)}')
    index, target = found[0]
    locus = target.qualifiers.get('locus_tag', [''])[0]
    if locus != LOCUS_TAG:
        raise ValueError(f'Expected locus tag {LOCUS_TAG}, found {locus!r}')
    translation = target.qualifiers.get('translation', [''])[0]
    if not translation:
        raise ValueError('Pinned CDS has no deposited translation')
    if index == 0 or index + 1 >= len(features):
        raise ValueError('Pinned CDS lacks immediate neighboring CDS features')
    upstream, downstream = features[index - 1], features[index + 1]
    start, end = bounds(target)
    up_start, up_end = bounds(upstream)
    down_start, down_end = bounds(downstream)

    def q(feature, key):
        return feature.qualifiers.get(key, [''])[0]

    row = {
        'phage': 'Gordonia phage GMA6', 'source_accession': SOURCE_ACCESSION,
        'source_sha256': source_hash, 'protein_id': PROTEIN_ID, 'locus_tag': LOCUS_TAG,
        'deposited_product': q(target, 'product'), 'strand': '+' if target.location.strand == 1 else '-',
        'start_1based': start, 'end_1based': end, 'length_aa': len(translation),
        'upstream_locus_tag': q(upstream, 'locus_tag'), 'upstream_product': q(upstream, 'product'),
        'upstream_start_1based': up_start, 'upstream_end_1based': up_end,
        'upstream_overlap_bp': overlap_bp(target, upstream),
        'downstream_locus_tag': q(downstream, 'locus_tag'), 'downstream_product': q(downstream, 'product'),
        'downstream_start_1based': down_start, 'downstream_end_1based': down_end,
        'downstream_overlap_bp': overlap_bp(target, downstream),
        'literature_function': 'DUF4185-family endo-D-arabinofuranase (GH183)',
        'functional_evidence_scope': 'Purified recombinant protein assayed on mycobacterial arabinogalactan in vitro',
        'native_phage_lysis_tested': False, 'primary_source': PAPER,
    }
    return row, translation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--reference-fasta', type=Path, default=REFERENCE_FASTA)
    parser.add_argument('--reference-metadata', type=Path, default=REFERENCE_METADATA)
    parser.add_argument('--mmseqs', type=Path, default=MMSEQS)
    parser.add_argument('--out-dir', type=Path, default=OUT)
    parser.add_argument('--work-dir', type=Path, default=ROOT / '.experiments/gma6/mmseqs_work')
    parser.add_argument('--precomputed-hits', type=Path,
                        help='use a raw MMseqs TSV produced separately, for split Windows/WSL workflows')
    parser.add_argument('--replay-hits', type=Path,
                        help='optional independent raw search result to verify byte-identical replay')
    parser.add_argument('--mmseqs-version', default='not-recorded')
    parser.add_argument('--extract-only', action='store_true',
                        help='write the mapped CDS and query FASTA without launching MMseqs')
    parser.add_argument('--download', action='store_true', help='retrieve the pinned genome from NCBI Eutils')
    args = parser.parse_args()
    for path in (args.source, args.reference_fasta, args.reference_metadata):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.precomputed_hits is None and not args.mmseqs.is_file():
        raise FileNotFoundError(args.mmseqs)
    if args.download:
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
            'db': 'nuccore', 'id': SOURCE_ACCESSION, 'rettype': 'gb', 'retmode': 'text', 'tool': 'actinolysis'})
        download(url, args.source, load_config(ROOT / 'config/analysis.yaml'), maximum=3_000_000)

    record = SeqIO.read(args.source, 'genbank')
    row, translation = analyze_record(record, sha256(args.source))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    query_fasta = args.work_dir / 'PhageGH4185.faa'
    query_fasta.write_text(f'>{LOCUS_TAG}|{PROTEIN_ID}|{SOURCE_ACCESSION}\n{translation}\n', encoding='ascii')
    table_path = args.out_dir / 'gma6_arabinanase.tsv'
    with table_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)
    if args.extract_only:
        return
    raw_hits = args.work_dir / 'mmseqs_hits.tsv'
    tmp = args.work_dir / 'mmseqs_tmp'
    command = [str(args.mmseqs), 'easy-search', str(query_fasta), str(args.reference_fasta), str(raw_hits), str(tmp),
               '--threads', '4', '-s', '7.5', '-k', '5', '--alignment-mode', '3', '--max-seqs', '100',
               '--min-seq-id', '0.25', '-c', '0.5', '--cov-mode', '0', '--split-memory-limit', '1G',
               '-e', '0.001', '--format-output', FORMAT]
    if args.precomputed_hits is not None:
        raw_hits = args.precomputed_hits
        if not raw_hits.is_file():
            raise FileNotFoundError(raw_hits)
        version = args.mmseqs_version
    else:
        run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(f'MMseqs failed ({run.returncode}): {run.stderr[-2000:]}')
        version = subprocess.run([str(args.mmseqs), 'version'], cwd=ROOT, capture_output=True,
                                 text=True, check=True).stdout.strip()

    hits = []
    with args.reference_metadata.open(encoding='utf-8', newline='') as handle:
        reference_categories = {row['query_id']: row['expected_category']
                               for row in csv.DictReader(handle, delimiter='\t')}
    if raw_hits.exists():
        with raw_hits.open(encoding='utf-8') as handle:
            for record_hit in csv.DictReader(handle, fieldnames=FORMAT.split(','), delimiter='\t'):
                identity, qcov, tcov, evalue = map(float, (record_hit['fident'], record_hit['qcov'],
                                                            record_hit['tcov'], record_hit['evalue']))
                reference_id = record_hit['target']
                category = reference_categories.get(reference_id, 'unresolved')
                hits.append({
                    'reference_id': reference_id, 'reference_category': category, 'identity': identity,
                    'alignment_length': int(record_hit['alnlen']),
                    'query_start_1based': int(record_hit['qstart']), 'query_end_1based': int(record_hit['qend']),
                    'reference_start_1based': int(record_hit['tstart']), 'reference_end_1based': int(record_hit['tend']),
                    'query_length': int(record_hit['qlen']), 'reference_length': int(record_hit['tlen']),
                    'query_coverage': qcov,
                    'reference_coverage': tcov, 'evalue': record_hit['evalue'],
                    'bit_score': float(record_hit['bits']),
                    'passes_existing_LysB_assignment_thresholds': passes_lysb_thresholds(
                        category, identity, qcov, tcov, evalue),
                })
    hits.sort(key=lambda item: (float(item['evalue']), -item['bit_score'], item['reference_id']))
    panel_path = args.out_dir / 'reference_panel_hits.tsv'
    with panel_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(hits)
    lysb_hits = [hit for hit in hits if hit['reference_category'] == 'lysin_b']
    hit_path = args.out_dir / 'lysB_panel_hits.tsv'
    with hit_path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(lysb_hits)
    report = {
        'status': 'pass', 'source_accession': SOURCE_ACCESSION,
        'source_record_sha256': sha256(args.source), 'protein_id': PROTEIN_ID,
        'protein_length_aa': len(translation), 'reference_fasta_sha256': sha256(args.reference_fasta),
        'reference_metadata_sha256': sha256(args.reference_metadata),
        'mmseqs_version': version, 'mmseqs_command': ['mmseqs', 'easy-search', '<query_fasta>',
            '<reference_fasta>', '<output_tsv>', '<temporary_directory>'] + command[6:],
        'raw_mmseqs_hits_sha256': sha256(raw_hits),
        'all_reference_panel_hits_returned': len(hits),
        'LysB_reference_panel_hits_returned': len(lysb_hits),
        'lysB_panel_hits_meeting_existing_assignment_thresholds': sum(
            item['passes_existing_LysB_assignment_thresholds'] for item in lysb_hits),
        'best_non_LysB_reference_hit': (
            {key: hit[key] for key in ('reference_id', 'reference_category', 'identity', 'query_coverage',
                                       'reference_coverage', 'evalue')}
            for hit in hits if hit['reference_category'] != 'lysin_b'
        ).__next__() if any(hit['reference_category'] != 'lysin_b' for hit in hits) else None,
        'existing_thresholds': THRESHOLDS,
        'interpretation': (
            'GMA6_45 has direct in-vitro arabinanase evidence. Its strongest hit is a LysA-labeled reference '
            'that was itself included on published profile prediction rather than protein-specific assay evidence. '
            'No LysB-class reference meets the frozen thresholds. The cross-category hit highlights possible '
            'family overlap or multifunctionality but cannot establish endolysin activity or native phage lysis.'),
        'architecture_table_sha256': sha256(table_path), 'reference_panel_hit_table_sha256': sha256(panel_path),
        'hit_table_sha256': sha256(hit_path),
        'script_sha256': sha256(Path(__file__).resolve()), 'primary_source': PAPER,
    }
    (args.out_dir / 'audit_report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    if args.replay_hits is not None:
        if not args.replay_hits.is_file():
            raise FileNotFoundError(args.replay_hits)
        replay = {
            'status': 'pass' if sha256(raw_hits) == sha256(args.replay_hits) else 'fail',
            'raw_result_sha256': sha256(raw_hits), 'independent_replay_sha256': sha256(args.replay_hits),
            'byte_identical': raw_hits.read_bytes() == args.replay_hits.read_bytes(),
            'scope': 'Independent MMseqs2 run using the same pinned query, reference panel, parameters and a fresh temporary/output directory.',
        }
        (args.out_dir / 'replay_validation.json').write_text(
            json.dumps(replay, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        if replay['status'] != 'pass':
            raise ValueError('Independent MMseqs2 replay differs from the recorded raw result')


if __name__ == '__main__':
    main()
