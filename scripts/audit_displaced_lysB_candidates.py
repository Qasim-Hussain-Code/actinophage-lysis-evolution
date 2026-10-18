#!/usr/bin/env python3
"""Map literature-proposed displaced LysB candidates to pinned genome CDSs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from Bio import SeqIO
from actinolysis.acquisition import download
from actinolysis.common import load_config, write_json, write_tsv

RAW = ROOT / '.experiments/rhodococcus_followup/raw'
OUT = ROOT / 'results/displaced_lysB_candidates'
CITATION = 'https://journals.asm.org/doi/10.1128/aem.01952-10'

SOURCES = {
    'GU580941.1': {
        'phage': 'ReqiPepy6', 'A': ('ADD80901.1', 'LysA'), 'B': ('ADD80924.1', 'LysB'),
        'B_evidence': 'putative LysB; sequence similarity reported; no protein-specific assay',
        'assay_scope': 'Whole-phage infection and host-range experiments; no LysB-specific functional attribution',
    },
    'GU580942.1': {
        'phage': 'ReqiPoco6', 'A': ('ADD81010.1', 'LysA'), 'B': ('ADD81032.1', 'LysB'),
        'B_evidence': 'putative LysB; sequence similarity reported; no protein-specific assay',
        'assay_scope': 'Whole-phage infection and host-range experiments; no LysB-specific functional attribution',
    },
    'GU580940.1': {
        'phage': 'ReqiDocB7', 'A': ('ADD80813.1', 'LysA'), 'B': ('ADD80824.1', 'LysB'),
        'B_evidence': 'possible LysB; weak cutinase and phage LysB sequence similarity reported',
        'assay_scope': 'Whole-phage host killing in soil; no LysB-specific functional attribution',
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def feature_bounds(feature):
    return int(feature.location.start) + 1, int(feature.location.end)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=Path, default=RAW)
    parser.add_argument('--out-dir', type=Path, default=OUT)
    parser.add_argument('--download', action='store_true', help='retrieve pinned GenBank versions from NCBI Eutils')
    args = parser.parse_args()
    if args.download:
        cfg = load_config(ROOT / 'config/analysis.yaml')
        for accession in SOURCES:
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
                'db': 'nuccore', 'id': accession, 'rettype': 'gb', 'retmode': 'text', 'tool': 'actinolysis'})
            download(url, args.input_dir / f'{accession}.gb', cfg, maximum=3_000_000)

    rows = []
    source_hashes = {}
    for accession, spec in SOURCES.items():
        path = args.input_dir / f'{accession}.gb'
        if not path.is_file():
            raise FileNotFoundError(path)
        record = SeqIO.read(path, 'genbank')
        if record.id != accession:
            raise ValueError(f'{accession}: version mismatch, found {record.id}')
        if record.annotations.get('molecule_type') != 'DNA':
            raise ValueError(f'{accession}: expected DNA genome')
        cds = [feature for feature in record.features if feature.type == 'CDS']
        by_protein = {}
        for index, feature in enumerate(cds):
            proteins = feature.qualifiers.get('protein_id', [])
            if proteins:
                if proteins[0] in by_protein:
                    raise ValueError(f'{accession}: duplicate CDS protein {proteins[0]}')
                by_protein[proteins[0]] = (index, feature)
        if record.annotations.get('topology') == 'circular':
            raise ValueError(f'{accession}: circular genome requires circular distance handling')
        (a_id, a_label), (b_id, b_label) = spec['A'], spec['B']
        if a_id not in by_protein or b_id not in by_protein:
            raise ValueError(f'{accession}: pinned LysA/B CDS was not found')
        ai, a = by_protein[a_id]
        bi, b = by_protein[b_id]
        a_product = a.qualifiers.get('product', [''])[0]
        b_product = b.qualifiers.get('product', [''])[0]
        if a_label.lower() not in a_product.lower() or b_label.lower() not in b_product.lower():
            raise ValueError(f'{accession}: expected current feature labels, found {a_product!r}, {b_product!r}')
        if bi <= ai:
            raise ValueError(f'{accession}: expected LysB after LysA in linear CDS order')
        a_start, a_end = feature_bounds(a)
        b_start, b_end = feature_bounds(b)
        intervening = bi - ai - 1
        gap_bp = b_start - a_end - 1
        if gap_bp < 0:
            raise ValueError(f'{accession}: pinned candidate genes overlap; review interpretation')
        source_hashes[accession] = sha256(path)
        rows.append({
            'phage': spec['phage'], 'host': 'Rhodococcus equi', 'source_accession': accession,
            'source_sha256': source_hashes[accession], 'lysA_protein': a_id,
            'lysA_gene_label': a.qualifiers.get('locus_tag', [''])[0],
            'lysA_start_1based': a_start, 'lysA_end_1based': a_end,
            'lysA_candidate_evidence': 'sequence/domain-based annotation in the genome study; no protein-specific assay reported',
            'lysB_candidate_protein': b_id, 'lysB_candidate_gene_label': b.qualifiers.get('locus_tag', [''])[0],
            'lysB_start_1based': b_start, 'lysB_end_1based': b_end,
            'lysB_candidate_evidence': spec['B_evidence'], 'intervening_CDS': intervening,
            'intervening_gap_bp': gap_bp,
            'genome_article_reports_protein_specific_LysB_assay': False,
            'whole_phage_assay_scope': spec['assay_scope'], 'primary_source': CITATION,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.out_dir / 'displaced_candidate_architectures.tsv'
    write_tsv(table_path, rows)
    report = {
        'status': 'pass', 'sources_checked': len(rows), 'candidate_pairs_mapped': len(rows),
        'protein_specific_LysA_assays_reported': 0, 'protein_specific_LysB_assays_reported': 0,
        'interpretation': 'Literature-proposed candidates are coordinate-mapped but remain unvalidated for LysB activity; whole-phage phenotypes do not assign function to these proteins.',
        'source_records_sha256': source_hashes,
        'script_sha256': sha256(Path(__file__).resolve()),
        'table_sha256': sha256(table_path),
        'source_article': CITATION,
    }
    write_json(args.out_dir / 'audit_report.json', report)


if __name__ == '__main__':
    main()
