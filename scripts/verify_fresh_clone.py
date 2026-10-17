#!/usr/bin/env python3
"""Compare committed pilot scientific outputs against a clean local clone replay."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from actinolysis.common import read_tsv, sha256, utcnow, write_json
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--clone', default='.tmp/fresh_clone')
args = parser.parse_args()
clone = (ROOT / args.clone).resolve()
paths = sorted([*ROOT.glob('results/tables/*.tsv'), *ROOT.glob('figures/*.png'), *ROOT.glob('figures/*.pdf')])
differences, artifacts = [], []
for path in paths:
    relative = path.relative_to(ROOT)
    actual = sha256(clone / relative) if (clone / relative).is_file() else None
    expected = sha256(path)
    artifacts.append({'path': str(relative).replace('\\', '/'), 'expected_sha256': expected, 'actual_sha256': actual})
    if actual != expected:
        differences.append(str(relative))
source_checks = []
for row in read_tsv(clone / 'data/manifests/pilot_inputs.tsv'):
    source_checks.append({'path': row['path'], 'pass': sha256(clone / row['path']) == row['sha256']})
for row in read_tsv(clone / 'data/manifests/reference_manifest.tsv'):
    source = clone / f"data/raw/references/{row['source_accession']}.gb"
    source_checks.append({'path': str(source.relative_to(clone)).replace('\\', '/'), 'pass': sha256(source) == row['source_sha256']})
execution = json.loads((clone / 'results/reports/pilot_execution.json').read_text())
commit = subprocess.run(['git', '-C', str(clone), 'rev-parse', 'HEAD'], check=True, capture_output=True, text=True).stdout.strip()
passed = not differences and all(r['pass'] for r in source_checks) and execution['status'] == 'success' and len(execution['stages']) == 9
report = {'status': 'pass' if passed else 'fail', 'verified_utc': utcnow(), 'clone_commit': commit,
    'method': 'Local clean Git clone; no copied raw or processed inputs; downloaded version-pinned source records; all nine stages rerun with uncached family clustering.',
    'software': 'Existing pinned WSL Python environment and MMseqs2 executable reused by absolute path; software was not independently reinstalled.',
    'compared_artifacts': len(artifacts), 'differences': differences, 'artifacts': artifacts, 'source_checks': source_checks,
    'execution_seconds': execution['elapsed_seconds'], 'peak_child_rss_bytes': execution['peak_child_rss_bytes']}
write_json(ROOT / 'results/reports/fresh_clone_reproduction.json', report)
print(json.dumps({k: v for k, v in report.items() if k not in {'artifacts', 'source_checks'}}))
if not passed:
    raise SystemExit(1)
