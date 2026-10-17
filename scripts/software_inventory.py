"""Measure installed immutable software once; refresh after changing software."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
result={'immutable_directories':{name:sum(p.stat().st_size for p in (ROOT/name).rglob('*') if p.is_file()) for name in ['.venv','.tools']},
        'policy':'Refresh this inventory after software changes. Pipeline disables bytecode writes; downloads and outputs are measured separately at every guard.'}
path=ROOT/'results/reports/software_storage.json'; path.parent.mkdir(parents=True,exist_ok=True)
path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(result),flush=True)
