#!/usr/bin/env python3
"""Verify added experimental LysA/B allele mappings against pinned GenBank records."""
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
from Bio.Seq import Seq
from actinolysis.acquisition import download
from actinolysis.common import load_config, sha256, write_json, write_tsv

RAW = ROOT / '.experiments/assayed_reference_followup'
OUT = ROOT / 'results/assayed_reference_followup'

SOURCES = {
    'EU203571.2': {
        'phage': 'Giles', 'A': ('ABW88426.1', '31', 402), 'B': ('ABW88427.1', '32', 411),
        'A_evidence': 'gene-deletion/complementation lysis phenotype',
        'B_evidence': 'gene-deletion delays host lysis; D29 LysB complementation',
        'citation_A': 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003957',
        'citation_B': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC2774421/',
        'mapping': 'gene and protein accession in deposited genome; native lysB deletion study',
    },
    'KR029087.1': {
        'phage': 'PDRPxv', 'A': ('AKF12448.1', '51', 429), 'B': ('AKF12449.1', '52', 451),
        'A_evidence': 'purified protein; M. smegmatis activity assay',
        'B_evidence': 'purified protein; esterase and M. smegmatis activity assays',
        'citation_A': 'https://link.springer.com/article/10.1007/s11274-020-02858-7',
        'citation_B': 'https://link.springer.com/article/10.1007/s11274-020-02858-7',
        'mapping': 'paper names source genome and gp51/gp52; B source CDS is 1 nt longer than the reported 1355-nt gene; construct sequence not reported',
    },
    'KC172839.1': {
        'phage': 'BTCU-1', 'A': ('AGI61688.1', '7', 521), 'B': ('AGI61689.1', '8', 246),
        'A_evidence': 'purified protein; M. smegmatis activity assay',
        'B_evidence': 'purified protein; M. smegmatis activity assay',
        'citation_A': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6332426/',
        'citation_B': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6332426/',
        'mapping': 'reported ORF8 amplicon length and reverse primer agree; forward primer has one synonymous mismatch',
    },
    'MZ343158.1': {
        'phage': 'RitSun', 'A': ('QXG07407.1', '37', 454), 'B': ('QXG07408.1', '38', 356),
        'A_evidence': 'conference abstract reports recombinant LysA lysozyme and antimycobacterial assays; exact protein accession not stated',
        'B_evidence': 'purified protein; esterase, lysis-from-without, growth and biofilm assays',
        'citation_A': 'https://acbrdu.edu/pdf/Abstract_Book_FBR-2022_25102022.pdf',
        'citation_B': 'https://www.nature.com/articles/s41598-024-68636-9',
        'mapping': 'reported B primers map to the deposited ORF38 CDS after removal of restriction tails',
        'mapping_A': 'conference abstract names RitSun LysA activity but does not specify the tested protein accession or gene; QXG07407.1 is the adjacent source annotation only',
    },
}

def position(feature):
    return int(feature.location.start) + 1, int(feature.location.end)

def primer_matches(feature, record, forward, reverse, minimum_forward, minimum_reverse):
    cds = str(feature.extract(record.seq)).upper()
    forward_bases = max((i for i in range(1, len(forward) + 1) if forward[-i:] == cds[:i]), default=0)
    reverse_bases = max((i for i in range(1, len(reverse) + 1)
                         if str(Seq(reverse[-i:]).reverse_complement()) == cds[-i:]), default=0)
    if forward_bases < minimum_forward or reverse_bases < minimum_reverse:
        raise ValueError('Published primer binding sequence fails source-CDS boundary check')
    return forward_bases, reverse_bases

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=Path, default=RAW)
    parser.add_argument('--out-dir', type=Path, default=OUT)
    parser.add_argument('--download', action='store_true', help='retrieve the four version-pinned GenBank records from NCBI Eutils')
    args = parser.parse_args()
    if args.download:
        cfg = load_config(ROOT / 'config/analysis.yaml')
        for accession in SOURCES:
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
                'db': 'nuccore', 'id': accession, 'rettype': 'gb', 'retmode': 'text', 'tool': 'actinolysis'})
            download(url, args.input_dir / f'{accession}.gb', cfg, maximum=3_000_000)
    rows, summary = [], []
    for accession, spec in SOURCES.items():
        source = args.input_dir / f'{accession}.gb'
        if not source.is_file(): raise FileNotFoundError(source)
        record = SeqIO.read(source, 'genbank')
        if record.id != accession: raise ValueError(f'{accession}: version mismatch: {record.id}')
        if record.annotations.get('molecule_type') != 'DNA': raise ValueError(f'{accession}: expected DNA')
        cds = {f.qualifiers.get('protein_id', [''])[0]: f for f in record.features if f.type == 'CDS'}
        selected = {}
        for category in ('A', 'B'):
            protein, gene, length = spec[category]
            feature = cds.get(protein)
            if feature is None: raise ValueError(f'{accession}: missing protein {protein}')
            translation = feature.qualifiers.get('translation', [''])[0]
            if len(translation) != length: raise ValueError(f'{accession}: unexpected {category} translation length')
            selected[category] = feature
            start, end = position(feature)
            rows.append({
                'phage': spec['phage'], 'source_accession': accession, 'role': category,
                'gene_number': gene, 'protein_accession': protein, 'length_aa': length,
                'source_cds_length_nt': len(feature.extract(record.seq)),
                'start_1based': start, 'end_1based': end, 'assay_evidence': spec[f'{category}_evidence'],
                'mapping_basis': spec.get(f'mapping_{category}', spec['mapping']), 'citation': spec[f'citation_{category}'],
                'source_record_url': f'https://www.ncbi.nlm.nih.gov/nuccore/{accession}',
                'source_sha256': sha256(source),
            })
        a, b = selected['A'], selected['B']
        a_start, a_end = position(a); b_start, b_end = position(b)
        intervening = sum(1 for f in record.features if f.type == 'CDS' and
                          min(a_end, b_end) < int(f.location.start) + 1 < max(a_start, b_start))
        if intervening != 0: raise ValueError(f'{accession}: expected adjacent A/B source loci')
        summary.append({
            'phage': spec['phage'], 'source_accession': accession,
            'A_protein_accession': spec['A'][0], 'B_protein_accession': spec['B'][0],
            'A_length_aa': spec['A'][2], 'B_length_aa': spec['B'][2],
            'A_source_cds_length_nt': len(a.extract(record.seq)),
            'B_source_cds_length_nt': len(b.extract(record.seq)),
            'A_B_intervening_CDS': intervening, 'A_B_overlap_bp': max(0, min(a_end,b_end)-max(a_start,b_start)+1),
            'A_B_gap_bp': max(0, max(a_start,b_start)-min(a_end,b_end)-1),
            'both_enzymes_functionally_tested': spec['A_evidence'] != 'LysA functional status not assessed in cited study',
            'A_evidence_level': 'conference_abstract_sequence_unresolved' if spec['phage']=='RitSun' else 'full_primary_study',
            'evidence_scope': f"A: {spec['A_evidence']}; B: {spec['B_evidence']}",
            'source_sha256': sha256(source), 'citation_A': spec['citation_A'], 'citation_B': spec['citation_B'],
        })

    # BTCU-1 primer-derived insert: identical length and translated protein despite
    # the single base change introduced by the published forward primer.
    btc = SeqIO.read(args.input_dir / 'KC172839.1.gb', 'genbank')
    btc_b = next(f for f in btc.features if f.type == 'CDS' and f.qualifiers.get('protein_id') == ['AGI61689.1'])
    btc_cds = str(btc_b.extract(btc.seq)).upper()
    forward = 'GTGAGCGACCGCTGGCT'
    reverse = 'TCATGTCAAGTCGCGTAGAAACT'
    reverse_site = str(Seq(reverse).reverse_complement())
    mismatch = [i for i, (x, y) in enumerate(zip(forward, btc_cds), 1) if x != y]
    if len(btc_cds) != 741 or len(forward) != 17 or mismatch != [16] or not btc_cds.endswith(reverse_site):
        raise ValueError('BTCU-1 reported primers do not map as expected')
    primer_insert = forward + btc_cds[len(forward):]
    if str(Seq(primer_insert).translate(table=11, cds=True)) != btc_b.qualifiers['translation'][0]:
        raise ValueError('BTCU-1 primer-derived insert changes the translated protein')

    # PDRPxv and RitSun PCR primers include cloning tails. Check their 3' binding
    # segments against the relevant deposited CDS boundaries; do not claim that
    # the unsequenced recombinant constructs are byte-for-byte source alleles.
    primer_specs = {
        'KR029087.1': {
            'AKF12448.1': ('ATGCCGGGATCCATGCCCATGCGC', 'GCGAAGCTTCTACGGTTTGAGGGCG', 6, 16),
            'AKF12449.1': ('AGGGGATCCAGCATGACACAG', 'CACACTCGAGCGCTACGCGGGCAG', 9, 12),
        },
        'MZ343158.1': {
            'QXG07408.1': ('CCGGAATTCATGAAACTCAACGGGGTATATG', 'CCCAAGCTTTCATGCGACACCCCTCATC', 22, 19),
        },
    }
    primer_audit = []
    for accession, targets in primer_specs.items():
        record = SeqIO.read(args.input_dir / f'{accession}.gb', 'genbank')
        cds = {f.qualifiers.get('protein_id', [''])[0]: f for f in record.features if f.type == 'CDS'}
        for protein, (forward, reverse, expected_f, expected_r) in targets.items():
            fwd, rev = primer_matches(cds[protein], record, forward, reverse, expected_f, expected_r)
            primer_audit.append({'source_accession': accession, 'protein_accession': protein,
                                 'forward_primer_3prime_match_nt': fwd, 'reverse_primer_3prime_match_nt': rev,
                                 'construct_sequence_reported': False})
    primer_audit.append({'source_accession': 'KC172839.1', 'protein_accession': 'AGI61689.1',
                         'forward_primer_3prime_match_nt': 1, 'reverse_primer_3prime_match_nt': 23,
                         'construct_sequence_reported': False})
    pdr = SeqIO.read(args.input_dir / 'KR029087.1.gb', 'genbank')
    pdr_b = next(f for f in pdr.features if f.type == 'CDS' and f.qualifiers.get('protein_id') == ['AKF12449.1'])
    if len(pdr_b.extract(pdr.seq)) != 1356:
        raise ValueError('PDRPxv LysB source CDS length differs from the audited accession')

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir / 'source_allele_evidence.tsv', rows)
    write_tsv(args.out_dir / 'source_pair_architectures.tsv', summary)
    write_tsv(args.out_dir / 'primer_mapping.tsv', primer_audit)
    artifacts = {}
    for path in sorted(args.out_dir.glob('*.tsv')):
        artifacts[path.name] = sha256(path)
    report = {
        'status': 'pass', 'sources_checked': len(summary), 'mapped_proteins': len(rows),
        'sources_with_any_primary_report_for_both_A_and_B': sum(r['both_enzymes_functionally_tested'] for r in summary),
        'sources_with_full_article_evidence_for_both_A_and_B': sum(r['both_enzymes_functionally_tested'] and r['A_evidence_level']=='full_primary_study' for r in summary),
        'sources_with_conference_abstract_A_evidence': sum(r['A_evidence_level']=='conference_abstract_sequence_unresolved' for r in summary),
        'sources_with_zero_intervening_CDS': sum(r['A_B_intervening_CDS'] == 0 for r in summary),
        'BTCU1_primer_mismatch_position_1based': mismatch[0],
        'BTCU1_primer_derived_B_protein_matches_deposited_translation': True,
        'PDRPxv_B_source_CDS_length_nt': len(pdr_b.extract(pdr.seq)),
        'PDRPxv_B_paper_reported_gene_length_nt': 1355,
        'PDRPxv_B_exact_construct_sequence_resolved': False,
        'source_records_sha256': {accession: sha256(args.input_dir / f'{accession}.gb') for accession in SOURCES},
        'script_sha256': sha256(Path(__file__)),
        'protocol_sha256': sha256(ROOT / 'docs/assayed_reference_allele_protocol.md'),
        'interpretation': 'Added assay sources are linked controls; this audit does not change the locked cohort or infer population frequencies.',
        'artifacts': artifacts,
    }
    write_json(args.out_dir / 'audit_report.json', report)
    print(json.dumps(report, indent=2))

if __name__ == '__main__': main()
