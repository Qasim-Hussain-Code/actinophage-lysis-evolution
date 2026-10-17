"""Record a separate filesystem measurement and the latest measured run."""
import json
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from actinolysis.common import project_bytes,utcnow,write_json
execution=json.loads((ROOT/'results/reports/pilot_execution.json').read_text())
report={'measured_utc':utcnow(),'logical_project_file_bytes':project_bytes(exact=True),'free_bytes':shutil.disk_usage(ROOT).free,
        'effective_budget_bytes':execution['config']['resources']['effective_budget_bytes'],
        'reserve_bytes':execution['config']['resources']['reserve_bytes'],
        'peak_parent_rss_bytes':execution['peak_self_rss_bytes'],'peak_child_rss_bytes':execution['peak_child_rss_bytes'],
        'execution_seconds':execution['elapsed_seconds'],
        'measurement_method':'Native filesystem sum of regular file sizes, including environment, tools, caches, raw inputs, outputs and current Git files. Existing external system Python referents are reused. Filesystem allocation overhead is not included.',
        'memory_limitations':'Linux process high-water marks, not a sampled concurrent tree total; this report describes the linked full pilot replay.'}
write_json(ROOT/'results/reports/resource_usage.json',report)
print(json.dumps(report))
