#!/usr/bin/env python3
"""Authenticate two profile runs and compare their complete frozen artifact sets."""
import argparse
import hashlib
import json
from pathlib import Path

def authenticate(root):
    report = json.loads((root / 'profile_report.json').read_text())
    for name, digest in report['artifacts'].items():
        target = (root / name).resolve()
        if not target.is_relative_to(root.resolve()): raise ValueError('Artifact path escapes the diagnostic directory')
        if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Profile artifact integrity failure: {name}')
    return report['artifacts']

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected', required=True, type=Path)
    parser.add_argument('--observed', required=True, type=Path)
    args = parser.parse_args()
    expected, observed = authenticate(args.expected), authenticate(args.observed)
    differences = sorted(name for name in set(expected) | set(observed) if expected.get(name) != observed.get(name))
    print(json.dumps({'status': 'pass' if not differences else 'fail', 'compared_artifacts': len(expected), 'differences': differences}))
    if differences: raise SystemExit(1)
