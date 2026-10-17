"""Snapshot or compare frozen scientific tables and figures by SHA-256."""
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from actinolysis.common import sha256,write_json
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--snapshot',action='store_true')
args=parser.parse_args()
paths=sorted([*ROOT.glob('results/tables/*.tsv'),*ROOT.glob('figures/*.png'),*ROOT.glob('figures/*.pdf')])
current={str(p.relative_to(ROOT)).replace('\\','/'):sha256(p) for p in paths}
baseline=ROOT/'.tmp/reproduction_baseline.json'
if args.snapshot:
    write_json(baseline,current)
    print(f'Snapshotted {len(current)} scientific artifacts')
else:
    if not baseline.exists(): raise SystemExit('Run --snapshot before the frozen offline rerun')
    old=json.loads(baseline.read_text())
    differences=[p for p in sorted(set(old)|set(current)) if old.get(p)!=current.get(p)]
    report={'status':'pass' if not differences else 'fail','compared_artifacts':len(current),'differences':differences,
        'scope':'Every results/tables TSV and figures PNG/PDF on this machine and pinned environment; execution timing, logs and resource/cache reports are excluded.'}
    write_json(ROOT/'results/reports/reproduction.json',report)
    print(json.dumps(report))
    if differences: raise SystemExit(1)
