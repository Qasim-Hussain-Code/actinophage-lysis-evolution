#!/usr/bin/env python3
"""Run query-only source-held-out and cellular-esterase reference validation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from actinolysis.reference_validation import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="Read-only isolated analysis root containing its config and frozen reference panel")
    parser.add_argument("--out", required=True, type=Path, help="New validation output folder inside the project")
    parser.add_argument("--offline", action="store_true", help="Require previously downloaded challenge assets")
    args = parser.parse_args()
    report = run(args.root, args.out, project_root=ROOT, offline=args.offline)
    print(json.dumps({"status": report["status"], "queries": report["query_count"], "reference_count": report["reference_count"], "output": str(args.out)}))


if __name__ == "__main__":
    main()
