"""Run and retain the actual test result in the active Python environment."""
import json
from pathlib import Path
import re
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
start=time.monotonic()
result=subprocess.run([sys.executable,'-B','-m','pytest','-q'],cwd=ROOT,capture_output=True,text=True)
log=ROOT/'results/logs/tests.txt';log.parent.mkdir(parents=True,exist_ok=True)
log.write_text(result.stdout+result.stderr,encoding='utf-8')
match=re.search(r'(\d+) passed',result.stdout)
report={'status':'pass' if result.returncode==0 else 'fail','passed':int(match.group(1)) if match else None,
        'python':sys.version.split()[0],'elapsed_seconds':round(time.monotonic()-start,3),'log':'results/logs/tests.txt'}
(ROOT/'results/reports/tests.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(report))
raise SystemExit(result.returncode)
