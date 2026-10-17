#!/bin/sh
set -eu
# Run in project root; no system packages or global Python changes.
mkdir -p results/reports
python3 - <<'PY'
from pathlib import Path
import json,shutil
free=shutil.disk_usage('.').free
cap=min(10_000_000_000,free-4_000_000_000)
if cap<1_000_000_000: raise SystemExit('Refusing install: insufficient bounded storage')
Path('.tools').mkdir(exist_ok=True)
Path('results/reports').mkdir(parents=True,exist_ok=True)
Path('results/reports/resource_baseline.json').write_text(json.dumps({'initial_free_bytes':free,'effective_budget_bytes':cap,'install_estimate_bytes':700_000_000,'genome_and_analysis_estimate_bytes':600_000_000,'reserve_bytes':4_000_000_000},indent=2)+'\n')
PY
python3 -m venv --without-pip .venv
curl -fL --retry 2 --max-time 120 https://bootstrap.pypa.io/get-pip.py -o .tools/get-pip.py
.venv/bin/python .tools/get-pip.py --no-cache-dir
requirements_file=requirements.txt
if [ -f requirements-lock.txt ]; then requirements_file=requirements-lock.txt; fi
.venv/bin/python -m pip install --no-cache-dir --no-compile -r "$requirements_file"
curl -fL --retry 2 --max-time 120 https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-avx2.tar.gz -o .tools/mmseqs-linux-avx2.tar.gz
tar -xzf .tools/mmseqs-linux-avx2.tar.gz -C .tools
.tools/mmseqs/bin/mmseqs version
.venv/bin/python scripts/software_inventory.py
