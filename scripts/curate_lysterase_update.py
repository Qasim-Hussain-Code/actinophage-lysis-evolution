#!/usr/bin/env python3
"""Resolve the six Table 1 proteins and verify their deposited CDS translations."""
from __future__ import annotations
import argparse
import hashlib
import io
from pathlib import Path
import re
import sys
import time
from urllib.parse import urlencode

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
from actinolysis.acquisition import download
from actinolysis.common import budget_check, load_config, read_tsv, sha256, write_json, write_tsv
from actinolysis.quality import translation_check

DOI = 'https://doi.org/10.1128/jb.00410-26'
TARGETS = [('D29', 'AAC18452.1', '12'), ('MissWhite', 'ASZ75177.1', '11'),
           ('Bxb1', 'AAG59714.1', '9'), ('Constella', 'AXF51531.1', '41'),
           ('Bongo', 'YP_009604896.1', '38'), ('phiGD57-1', 'QPO17264.1', '25')]

def url(database, accession, kind):
    return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode({
        'db': database, 'id': accession, 'rettype': kind, 'retmode': 'text', 'tool': 'actinolysis'})

def resolve_protein(record, expected):
    if record.id != expected:
        raise ValueError('Protein accession-version differs from the curated Table 1 link')
    sequences = {str(record.seq)}
    coding = [f.qualifiers['coded_by'][0] for f in record.features
              if f.type == 'CDS' and len(f.qualifiers.get('coded_by', [])) == 1]
    if len(coding) != 1:
        raise ValueError('Protein does not have exactly one explicit coding source')
    sources = set(re.findall(r'([A-Z][A-Z0-9_]*\.\d+):', coding[0]))
    if len(sources) != 1 or len(sequences) != 1:
        raise ValueError('Protein coding source is ambiguous')
    return sources.pop(), str(record.seq)

def verify_cds(record, source, protein, sequence):
    if record.id != source:
        raise ValueError('Nucleotide accession-version differs from the protein coding source')
    features = [f for f in record.features if f.type == 'CDS' and f.qualifiers.get('protein_id') == [protein]]
    if len(features) != 1:
        raise ValueError('Expected one exactly accessioned CDS in the nucleotide source')
    translated, status = translation_check(features[0], record)
    if translated != sequence or status != 'verified':
        raise ValueError('Protein sequence and independently translated genomic CDS disagree')
    return features[0]

def freeze_tsv(path, rows):
    normalized = [{key: str(value) for key, value in row.items()} for row in rows]
    if path.exists() and read_tsv(path) != normalized:
        raise ValueError('Curated sequence evidence changed; preserve the existing freeze and use a new analysis version')
    if not path.exists(): write_tsv(path, rows)

def freeze_json(path, value):
    import json
    if path.exists() and json.loads(path.read_text()) != value:
        raise ValueError('Curation provenance changed; preserve the existing freeze and use a new analysis version')
    if not path.exists(): write_json(path, value)

def run(offline=False):
    cfg = load_config(ROOT / 'config/analysis.yaml')
    budget_check(cfg, 25_000_000)
    raw = ROOT / 'data/raw/lysterase_update'
    protein_path = raw / 'table1_proteins.gb'
    protein_url = url('protein', ','.join(p for _,p,_ in TARGETS), 'gp')
    download(protein_url, protein_path, cfg, offline, maximum=2_000_000)
    proteins = {r.id: r for r in SeqIO.parse(protein_path, 'genbank')}
    if set(proteins) != {p for _,p,_ in TARGETS}:
        raise ValueError('Retrieved protein set differs from Table 1')
    existing = {r['reference_id']: r for r in read_tsv(ROOT / 'data/manifests/reference_manifest.tsv')}
    ledger, aliases, sources = [], [], []
    for phage, protein, gene in TARGETS:
        source, sequence = resolve_protein(proteins[protein], protein)
        source_url = url('nuccore', source, 'gbwithparts')
        source_path = raw / f'{source}.gb'
        if not offline and not source_path.exists(): time.sleep(.4)
        download(source_url, source_path, cfg, offline, maximum=3_000_000)
        record = SeqIO.read(source_path, 'genbank')
        feature = verify_cds(record, source, protein, sequence)
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        sources.append({'phage': phage, 'protein_accession': protein, 'source_accession': source,
            'source_url': source_url, 'source_sha256': sha256(source_path), 'protein_sha256': digest,
            'length_aa': len(sequence), 'translation_status': 'verified',
            'start': int(feature.location.start), 'end': int(feature.location.end),
            'strand': feature.location.strand})
        if phage == 'D29':
            if digest != existing['D29_LysB']['protein_sha256']:
                raise ValueError('Published D29 protein is not the existing reference sequence')
            aliases.append({'reference_id': 'D29_LysB', 'protein_accession': protein,
                            'existing_protein_accession': existing['D29_LysB']['resolved_protein_accession'],
                            'protein_sha256': digest, 'action': 'existing_exact_sequence_alias'})
            continue
        ledger.append({'reference_id': phage.replace('-', '_') + '_LysB_2026', 'category': 'lysin_b',
            'phage': phage, 'protein_accession': protein, 'source_accession': source, 'gene': gene,
            'length_aa': len(sequence), 'reference_url': DOI, 'evidence_type': 'experimental_antimycobacterial',
            'expected_category': 'lysin_b', 'assay': 'Purified-protein antimycobacterial assay',
            'assay_outcome': 'positive', 'panel_action': 'add', 'evidence_locator': 'Table 1; Figure 2',
            'source_record_url': f'https://www.ncbi.nlm.nih.gov/nuccore/{source}',
            'caveat': 'Moderate maximum; native mycolate-substrate cleavage and phage lysis requirement unresolved'})
    out = ROOT / 'results/lysterase_update/curation'
    freeze_tsv(ROOT / 'references/lysterase_2026_expansion.tsv', ledger)
    freeze_tsv(out / 'source_sequences.tsv', sources)
    freeze_tsv(out / 'reference_aliases.tsv', aliases)
    freeze_json(out / 'source_manifest.json', {'paper': DOI, 'protein_record_url': protein_url,
        'protein_records_sha256': sha256(protein_path), 'protein_records': len(proteins),
        'added_references': len(ledger), 'aliases': len(aliases),
        'protocol_sha256': sha256(ROOT / 'docs/lysterase_update_protocol.md'),
        'curation_script_sha256': sha256(__file__)})
    print('Verified six source translations; five new references and one existing sequence alias.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true')
    run(parser.parse_args().offline)
