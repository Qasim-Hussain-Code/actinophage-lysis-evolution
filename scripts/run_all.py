#!/usr/bin/env python3
"""Deterministic, fail-fast pipeline with hashed stage receipts."""
from __future__ import annotations
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
sys.dont_write_bytecode=True

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.cache/matplotlib'))
os.environ.setdefault('OMP_NUM_THREADS','4')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from actinolysis.common import budget_check, load_config, project_bytes, read_tsv, record_environment, sha256, utcnow, write_json

STAGES={1:('acquisition','run'),2:('quality','run'),3:('identification','run'),4:('identification','evidence'),
        5:('architecture','run'),6:('neighborhoods','run'),7:('statistics','run'),8:('validation','run'),9:('figures','run')}
OUTPUTS={1:[],2:['data/processed/genomes.tsv','data/processed/cds.tsv','data/processed/proteins.faa'],
  3:['data/processed/reference_hits.tsv','data/manifests/reference_manifest.tsv','data/processed/references.faa'],4:['data/processed/lysis_evidence.tsv'],
  5:['data/processed/architectures.tsv','data/processed/architecture_pairs.tsv','results/tables/architecture_sensitivity.tsv'],
  6:['data/processed/families.tsv','results/tables/neighborhood_pairs.tsv','results/tables/neighborhood_controls.tsv','results/tables/neighborhoods.tsv'],
  7:['results/tables/evolutionary_tests.tsv'],8:['results/reports/validation.json','results/tables/candidates.tsv','results/tables/control_alignment_diagnostics.tsv'],9:['results/reports/figure_manifest.json']}

def scientific_fingerprint(cfg):
    fields=['mode','seed','sample_sizes','sampling','identification','architecture','neighborhood','statistics']
    return hashlib.sha256(json.dumps({k:cfg[k] for k in fields},sort_keys=True).encode()).hexdigest()

def source_fingerprint(number):
    names=[ROOT/'src/actinolysis'/f'{STAGES[number][0]}.py',ROOT/'src/actinolysis/common.py',ROOT/'scripts/run_all.py']
    if number in {3,4}: names.append(ROOT/'references/protein_controls.tsv')
    if number==8: names.append(ROOT/'references/architecture_controls.tsv')
    if number==8: names.append(ROOT/'src/actinolysis/diagnostics.py')
    return {str(p.relative_to(ROOT)):sha256(p) for p in names}

def main(argv=None,stage=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default='config/analysis.yaml')
    parser.add_argument('--mode',choices=['pilot','intermediate','full'])
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--resume-from',type=int,choices=range(1,10),default=1)
    parser.add_argument('--stop-after',type=int,choices=range(1,10),default=9)
    parser.add_argument('--validation-only',action='store_true')
    args=parser.parse_args(argv)
    cfg=load_config(ROOT/args.config)
    if args.mode: cfg['mode']=args.mode
    if not 1<=args.threads<=4: parser.error('--threads must be between1 and4')
    cfg.update(threads=args.threads,offline=args.offline)
    os.chdir(ROOT)
    if cfg['mode']!='pilot':
        gate=ROOT/'results/reports/validation.json'
        pilot=ROOT/'results/reports/pilot_stage_receipts.json'
        if (not gate.exists() or not pilot.exists() or not all(str(s) in json.loads(pilot.read_text()) for s in range(1,10))
                or not json.loads(gate.read_text()).get('expansion_ready',False)):
            raise RuntimeError('Expansion blocked: complete the pilot and validate independently supported separated architectures before larger acquisition')
    baseline=ROOT/'results/reports/resource_baseline.json'
    if baseline.exists(): cfg['resources']['effective_budget_bytes']=json.loads(baseline.read_text())['effective_budget_bytes']
    else:
        free=shutil.disk_usage(ROOT).free
        cfg['resources']['effective_budget_bytes']=min(cfg['resources']['max_project_bytes'],free-cfg['resources']['reserve_bytes'])
        write_json(baseline,{'initial_free_bytes':free,'effective_budget_bytes':cfg['resources']['effective_budget_bytes'],'reserve_bytes':cfg['resources']['reserve_bytes']})
    start=time.monotonic()
    report={'started_utc':utcnow(),'mode':cfg['mode'],'config':cfg,'environment':record_environment(),'stages':[]}
    receipts=ROOT/f"results/reports/{cfg['mode']}_stage_receipts.json"
    previous=json.loads(receipts.read_text()) if receipts.exists() else {}
    requested=[8] if args.validation_only else ([stage] if stage else list(range(args.resume_from,args.stop_after+1)))
    # A resumed dependency must still match its recorded checksum.
    if min(requested)>1:
        for number in range(1,min(requested)):
            receipt=previous.get(str(number))
            if receipt is None: raise RuntimeError(f'Stage{number} has no receipt; run it first')
            if receipt.get('scientific_config_sha256')!=scientific_fingerprint(cfg) or receipt.get('source_sha256')!=source_fingerprint(number):
                raise RuntimeError(f'Stage{number} config or source changed; rerun from stage1 to regenerate dependent outputs')
            for path,digest in receipt['outputs'].items():
                if not (ROOT/path).exists() or sha256(ROOT/path)!=digest: raise ValueError(f'Resume input changed: {path}')
    try:
        for number in requested:
            previous={k:v for k,v in previous.items() if int(k)<number}
            write_json(receipts,previous)
            budget_check(cfg,0); began=time.monotonic()
            print(f'Stage{number}: {STAGES[number][0]}',flush=True)
            module=importlib.import_module('actinolysis.'+STAGES[number][0])
            getattr(module,STAGES[number][1])(cfg)
            paths=OUTPUTS[number] if number!=1 else [f"data/manifests/{cfg['mode']}_inputs.tsv"]
            if number==1:
                paths=paths+[r['path'] for r in read_tsv(ROOT/paths[0])]
                paths+=['data/manifests/phagesdb_snapshot.json','data/manifests/eligibility_snapshot.tsv']
                paths += [r['path'] for r in json.loads((ROOT/'data/manifests/phagesdb_snapshot.json').read_text())['pages'] if (ROOT/r['path']).exists()]
            if number==3:
                paths=paths+[str(p.relative_to(ROOT)).replace('\\','/') for p in (ROOT/'data/raw/references').glob('*.gb')]
            if any(not (ROOT/path).exists() for path in paths): raise ValueError(f'Stage{number} missing required outputs')
            receipt={'outputs':{path:sha256(ROOT/path) for path in paths},'seconds':round(time.monotonic()-began,3),
                'scientific_config_sha256':scientific_fingerprint(cfg),'source_sha256':source_fingerprint(number)}
            previous[str(number)]=receipt
            report['stages'].append({'stage':number,**receipt})
            write_json(receipts,previous)
            budget_check(cfg,0)
        report['status']='success'
    except Exception as error:
        report.update(status='failed',error=f'{type(error).__name__}: {error}')
        raise
    finally:
        report.update(finished_utc=utcnow(),elapsed_seconds=round(time.monotonic()-start,3),project_bytes=project_bytes(),
            free_bytes=shutil.disk_usage(ROOT).free,
            peak_self_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            peak_child_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024)
        write_json(ROOT/f"results/reports/{cfg['mode']}_execution.json",report)

if __name__=='__main__': main()
