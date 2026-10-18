#!/usr/bin/env python3
"""Survey complete assay-reference source genomes with the unchanged classifier."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from Bio import SeqIO
import yaml
from actinolysis.common import NA, budget_check, load_config, read_tsv, sha256, write_json, write_tsv
from actinolysis.quality import extract_record

BASE = ROOT / 'results/lysterase_update/panel'
PUBLIC = ROOT / 'results/reference_architecture_controls'
PROTOCOL = ROOT / 'docs/reference_architecture_protocol.md'

def source_host_fields(record):
    feature = next((f for f in record.features if f.type == 'source'), None)
    host = feature.qualifiers.get('host', feature.qualifiers.get('lab_host',[NA]))[0] if feature else NA
    first = host.split()[0] if host else NA
    genus = first if first in {'Mycobacterium', 'Mycolicibacterium', 'Gordonia'} else NA
    return host or NA, genus

def cached_source(row):
    relative = Path('data/raw/references') / (row['source_accession'] + '.gb')
    prefixes = ['', '.experiments/lysterase_2026', '.experiments/lysterase_2026_replay',
                '.experiments/amendment_profile_v1', '.experiments/amendment_experimental_v1']
    candidates = [ROOT / prefix / relative for prefix in prefixes]
    candidates.append(ROOT / 'data/raw/lysterase_update' / (row['resolved_source_accession'] + '.gb'))
    for candidate in candidates:
        if not candidate.is_file(): continue
        if sha256(candidate) != row['source_sha256']: raise ValueError('Cached source differs from the frozen reference manifest')
        return candidate
    raise FileNotFoundError('Acquire the version-pinned reference records with the 39-reference cohort replay first')

def prepare(runtime):
    if runtime.exists() and any(runtime.iterdir()): raise ValueError('Preserve the existing workspace; resume it or select a new name')
    budget_check(load_config(ROOT / 'config/analysis.yaml'), additional_bytes=100_000_000)
    training_path = ROOT / 'results/lysb_profile_update/training_references.tsv'
    profile_report = json.loads((ROOT / 'results/lysb_profile_update/profile_report.json').read_text())
    if sha256(training_path) != profile_report['artifacts']['training_references.tsv']: raise ValueError('Assay-reference set changed')
    training = read_tsv(training_path)
    if len(training) != 14 or any(row['category'] != 'lysin_b' for row in training): raise ValueError('Expected fourteen assay-supported B references')
    references = read_tsv(BASE / 'data/manifests/reference_manifest.tsv')
    if len(references) != 39: raise ValueError('Expected the unchanged 39-reference catalogue')
    for directory in ('src','references'):
        shutil.copytree(ROOT / directory, runtime / directory, ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    (runtime / 'scripts').mkdir()
    shutil.copyfile(ROOT / 'scripts/run_all.py', runtime / 'scripts/run_all.py')
    shutil.copyfile(BASE / 'references/protein_controls.tsv', runtime / 'references/protein_controls.tsv')
    (runtime / 'data/manifests').mkdir(parents=True)
    shutil.copyfile(BASE / 'data/manifests/reference_manifest.tsv', runtime / 'data/manifests/reference_manifest.tsv')
    for row in references:
        target = runtime / 'data/raw/references' / (row['source_accession'] + '.gb')
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(cached_source(row), target)
        elif sha256(target) != row['source_sha256']: raise ValueError('Conflicting reference source checksums')
    by_source = defaultdict(list)
    for row in training: by_source[row['resolved_source_accession']].append(row)
    cohort, selection = [], []
    for accession, aliases in sorted(by_source.items()):
        reference = aliases[0]
        source = cached_source(reference)
        record = SeqIO.read(source, 'genbank')
        if record.id != accession: raise ValueError('Source accession-version changed')
        if len({row['source_sha256'] for row in aliases}) != 1: raise ValueError('Source aliases have conflicting bytes')
        phage = reference['reference_id'].split('_LysB')[0].replace('phiGD57_1','phiGD57-1')
        host, genus = source_host_fields(record)
        relative = f'data/raw/genbank/{accession}.gb'
        metadata = {'phage':phage, 'host':host, 'host_genus':genus, 'cluster':NA, 'subcluster':NA,
            'metadata_source':'GenBank_assay_reference_source', 'topology_evidence':record.annotations.get('topology',NA),
            'accession':accession, 'retrieved_utc':NA, 'sha256':sha256(source), 'path':relative,
            'source_url':reference['source_url'], 'inclusion_status':'pending_qc', 'exclusion_reason':'not_applicable'}
        try:
            extract_record(record, metadata)
            status, reason = 'eligible_complete_source', ''
            cohort.append(metadata)
            target = runtime / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        except ValueError as error:
            status, reason = 'excluded_source', str(error)
        selection.append({'source_accession':accession, 'phage':phage, 'assay_reference_ids':';'.join(sorted(row['reference_id'] for row in aliases)),
            'source_sha256':sha256(source), 'source_url':reference['source_url'], 'source_host':host, 'host_genus':genus,
            'length_bp':len(record.seq), 'source_description':record.description, 'input_status':status, 'exclusion_reason':reason})
    if not cohort or len({row['phage'] for row in cohort}) != len(cohort): raise ValueError('Empty or ambiguous source-control cohort')
    write_tsv(runtime / 'data/manifests/pilot_inputs.tsv', cohort)
    cfg = load_config(ROOT / 'config/analysis.yaml')
    cfg['sample_sizes'] = {'pilot':len(cohort)}
    cfg['sampling'] = {'host_genera':sorted({row['host_genus'] for row in cohort}), 'controls':[row['phage'] for row in cohort]}
    cfg['tools']['mmseqs'] = str((ROOT / cfg['tools']['mmseqs']).resolve())
    (runtime / 'config').mkdir()
    (runtime / 'config/analysis.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False), encoding='utf-8')
    frame = runtime / 'data/manifests/eligibility_snapshot.tsv'
    write_tsv(frame, [{**row,'seq_finished':True,'length':'not_available'} for row in cohort])
    write_json(runtime / 'data/manifests/phagesdb_snapshot.json', {'pages':[], 'rows':len(cohort),
        'method':'complete_assay_reference_source_survey', 'source_database':'GenBank',
        'selection_snapshot':{'path':'data/manifests/eligibility_snapshot.tsv','sha256':sha256(frame),'rows':len(cohort),
            'basis':'Every eligible complete source supplying a frozen assayed B reference; no population sampling'}})
    write_tsv(runtime / 'source_selection.tsv', selection)
    snapshot = {p.relative_to(runtime).as_posix():sha256(p) for directory in ('src','scripts','references')
                for p in sorted((runtime / directory).rglob('*')) if p.is_file()}
    write_json(runtime / 'survey_inputs.json', {'protocol_sha256':sha256(PROTOCOL), 'driver_sha256':sha256(__file__),
        'reference_training_set_sha256':sha256(training_path), 'reference_manifest_sha256':sha256(BASE / 'data/manifests/reference_manifest.tsv'),
        'source_selection_sha256':sha256(runtime / 'source_selection.tsv'), 'snapshot_sha256':snapshot,
        'selection':'All unique resolved assay-B source records passing unchanged whole-genome QC'})

def summarize(runtime, target):
    if target.exists(): raise ValueError('Preserve the existing survey output; use --no-export for a replay')
    inputs = json.loads((runtime / 'survey_inputs.json').read_text())
    if sha256(runtime / 'data/manifests/reference_manifest.tsv') != inputs['reference_manifest_sha256']:
        raise ValueError('The functional reference catalogue changed during the survey')
    execution = json.loads((runtime / 'results/reports/pilot_execution.json').read_text())
    receipts = json.loads((runtime / 'results/reports/pilot_stage_receipts.json').read_text())
    if execution['status'] != 'success' or set(receipts) != {str(i) for i in range(1,6)}: raise ValueError('Stages 1 through 5 must complete')
    for receipt in receipts.values():
        for name,digest in receipt['outputs'].items():
            if sha256(runtime / name) != digest: raise ValueError('Survey stage output changed')
    target.mkdir(parents=True)
    shutil.copytree(runtime / 'results/tables', target / 'tables')
    (target / 'manifests').mkdir()
    for name in ('pilot_inputs.tsv','reference_manifest.tsv'):
        shutil.copyfile(runtime / 'data/manifests' / name, target / 'manifests' / name)
    shutil.copyfile(runtime / 'source_selection.tsv', target / 'source_selection.tsv')
    genomes = read_tsv(runtime / 'data/processed/genomes.tsv')
    cds = read_tsv(runtime / 'data/processed/cds.tsv')
    evidence = read_tsv(runtime / 'data/processed/lysis_evidence.tsv')
    training = read_tsv(ROOT / 'results/lysb_profile_update/training_references.tsv')
    architectures = read_tsv(runtime / 'results/tables/architectures.tsv')
    summary, controls = [], []
    for genome in genomes:
        acc = genome['accession']
        arc = next(row for row in architectures if row['accession'] == acc)
        local_refs = [row for row in training if row['resolved_source_accession'] == acc]
        assay_cds = []
        for ref in local_refs:
            genes = [row for row in cds if row['accession'] == acc and row['protein_accession'] == ref['resolved_protein_accession']]
            if len(genes) != 1: raise ValueError('Assayed source protein does not resolve uniquely after source QC')
            gene = genes[0]
            if hashlib.sha256(gene['translation'].encode()).hexdigest() != ref['protein_sha256']: raise ValueError('Assayed source sequence changed')
            calls = [row for row in evidence if row['protein_id'] == gene['protein_id'] and row['category'] == 'lysin_b']
            if len(calls) != 1: raise ValueError('Expected one selected B-family evidence row for the assayed source protein')
            call = calls[0]
            controls.append({'phage':genome['phage'],'source_accession':acc,'assayed_reference_id':ref['reference_id'],
                'protein_accession':gene['protein_accession'],'protein_id':gene['protein_id'],'confidence':call['confidence'],
                'selected_reference_id':call['reference_id'],'reference_evidence_type':call['reference_evidence_type'],
                'protein_sha256':ref['protein_sha256'],'source_self_match_expected':'true'})
            assay_cds.append(gene['protein_id'])
        summary.append({'phage':genome['phage'],'source_accession':acc,'host_genus':genome['host_genus'],
            'length_bp':genome['length'],'cds':genome['cds_count'],'assayed_b_protein_ids':';'.join(assay_cds),
            **{key:arc[key] for key in ('lysin_a_count','lysin_b_count','architecture','classification_reason','confidence',
                'lysin_a','lysin_b','intervening_cds','gap_bp','relative_strands')}})
    write_tsv(target / 'architecture_summary.tsv', summary)
    write_tsv(target / 'assayed_source_recovery.tsv', controls)
    raw_fields = ['protein_id','reference_id','identity','alignment_length','query_coverage','target_coverage','evalue','bits']
    import csv
    with (runtime / 'data/processed/reference_hits.tsv').open() as stream:
        hits = list(csv.DictReader(stream, fieldnames=raw_fields, delimiter='\t'))
    write_tsv(target / 'reported_reference_alignments.tsv', hits, raw_fields)
    write_json(target / 'survey_report.json', {'analysis':'assay_reference_source_architecture_survey','completed_stages':[1,2,3,4,5],
        'n_assayed_b_references':14,'n_source_records':len(read_tsv(runtime / 'source_selection.tsv')),
        'n_complete_sources_selected':len(read_tsv(runtime / 'data/manifests/pilot_inputs.tsv')),
        'n_genomes_passing_pipeline_qc':len(genomes),'n_cds':len(cds),'architecture_counts':dict(Counter(row['architecture'] for row in summary)),
        'reference_count':39,'input_provenance':inputs,
        'interpretation':'Control-selected source genomes; self-matches are expected. No population frequency, independent functional sensitivity or neighborhood inference is estimated.',
        'primary_cohort_modified':False,'expansion_gate_modified':False})
    artifact_hashes = {p.relative_to(target).as_posix():sha256(p) for p in sorted(target.rglob('*')) if p.is_file()}
    write_json(target / 'artifact_manifest.json', {'artifacts':artifact_hashes})
    print(json.dumps({'included_sources':len(genomes),'cds':len(cds),'architectures':dict(Counter(row['architecture'] for row in summary))}))

def run(args):
    runtime = ROOT / '.experiments' / args.workspace
    if not runtime.resolve().is_relative_to((ROOT / '.experiments').resolve()): raise ValueError('Workspace is outside .experiments')
    if not args.resume: prepare(runtime)
    frozen = json.loads((runtime / 'survey_inputs.json').read_text())
    if frozen['driver_sha256'] != sha256(__file__) or frozen['protocol_sha256'] != sha256(PROTOCOL): raise ValueError('Survey code or protocol changed')
    for name,digest in frozen['snapshot_sha256'].items():
        if sha256(runtime / name) != digest: raise ValueError('Isolated survey source snapshot changed')
    command = [sys.executable,'-B',str(runtime / 'scripts/run_all.py'),'--threads','4','--offline','--stop-after','5']
    if args.resume: command += ['--resume-from',str(args.resume)]
    subprocess.run(command,cwd=runtime,check=True)
    if not args.no_export: summarize(runtime,PUBLIC)
    elif PUBLIC.exists():
        expected = json.loads((PUBLIC / 'artifact_manifest.json').read_text())['artifacts']
        names = [name for name in expected if name.startswith('tables/')]
        for name in names:
            if sha256(runtime / 'results' / name) != expected[name]: raise ValueError('Source-survey scientific table differs on replay')
        print(json.dumps({'reproduced_scientific_tables':len(names),'status':'pass'}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',default='reference_source_architectures_v2')
    parser.add_argument('--resume',type=int,choices=range(1,6))
    parser.add_argument('--no-export',action='store_true')
    run(parser.parse_args())
