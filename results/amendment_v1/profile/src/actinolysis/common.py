"""Small, explicit I/O and resource guards shared by pipeline stages."""
from __future__ import annotations
import csv
import hashlib
import json
import platform
import shutil
from datetime import datetime, timezone
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]
NA = "not_available"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))

def write_tsv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError(f"Cannot infer schema for empty table: {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)

def write_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")

def load_config(path):
    cfg = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if cfg["mode"] not in {"pilot", "intermediate", "full"}:
        raise ValueError("mode must be pilot, intermediate or full")
    return cfg

def project_bytes(root=ROOT, exact=False):
    root=Path(root)
    inventory=root/'results/reports/software_storage.json'
    frozen={} if exact or not inventory.exists() else json.loads(inventory.read_text())['immutable_directories']
    size=sum(frozen.values())
    # Installed software is immutable during analysis. Reuse its measured size
    # to avoid thousands of remote-filesystem stat RPCs at every stage guard.
    for child in root.iterdir():
        if child.name in frozen: continue
        if child.is_file(): size+=child.stat().st_size
        elif child.is_dir(): size+=sum(p.stat().st_size for p in child.rglob('*') if p.is_file())
    return size

def budget_check(cfg, additional_bytes=0, root=ROOT):
    used = project_bytes(root)
    free = shutil.disk_usage(root).free
    cap = min(int(cfg["resources"]["max_project_bytes"]), int(cfg["resources"].get("effective_budget_bytes", 10_000_000_000)))
    reserve = int(cfg["resources"]["reserve_bytes"])
    if additional_bytes < 0 or used + additional_bytes > cap or free - additional_bytes < reserve:
        raise RuntimeError(f"Storage guard: used={used}, proposed={additional_bytes}, free={free}, cap={cap}, reserve={reserve}")
    return {"project_bytes": used, "free_bytes": free, "cap_bytes": cap}

def utcnow():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")

def record_environment():
    from importlib.metadata import version
    return {"python": platform.python_version(), "platform": platform.platform(),
            "packages": {p:version(p) for p in ["biopython", "numpy", "pandas", "scipy", "statsmodels", "matplotlib", "pyyaml", "pytest"]}}
