"""Bounded public downloads, frozen sampling and accession-specific records."""
from __future__ import annotations
import hashlib
import csv
import io
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from Bio import SeqIO
from .common import ROOT, NA, budget_check, read_tsv, sha256, utcnow, write_json, write_tsv

ACCESSION = re.compile(r"^[A-Z]{1,6}_?\d{5,10}(?:\.\d+)?$")
UA = "actinolysis/0.1 (public comparative genomics; bounded academic retrieval)"

def validate_accession(accession, version_required=False):
    if not ACCESSION.fullmatch(accession) or (version_required and "." not in accession):
        raise ValueError(f"Invalid accession{('-version' if version_required else '')}: {accession}")
    return accession

def download(url, path, cfg, offline=False, maximum=None):
    path = Path(path)
    if path.exists():
        return path
    if offline:
        raise FileNotFoundError(f"Offline input missing: {path}")
    limit = int(maximum or cfg["resources"]["max_download_bytes"])
    budget_check(cfg, limit)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":UA})
            with urllib.request.urlopen(req, timeout=60) as res, tmp.open("wb") as stream:
                size = 0
                while block := res.read(65536):
                    size += len(block)
                    if size > limit: raise RuntimeError(f"Download exceeds {limit} bytes: {url}")
                    stream.write(block)
            tmp.replace(path)
            return path
        except Exception:
            tmp.unlink(missing_ok=True)
            if attempt == 2: raise
            time.sleep(2 ** attempt)

def value_name(value, key):
    return str(value.get(key) or NA) if isinstance(value,dict) else str(value or NA)

def normalize_metadata(row):
    host = row.get("isolation_host") or row.get("host_strain") or row.get("host") or NA
    if isinstance(host,dict):
        host = (' '.join(str(host[k]) for k in ['genus','species','strain_name'] if host.get(k))
                or host.get('name') or host.get('host_name') or NA)
    genus = row.get("host_genus") or str(host).split()[0]
    if isinstance(genus,dict): genus = genus.get("genus_name") or genus.get("genus") or genus.get("name") or NA
    accession = str(row.get("genbank_accession") or row.get("accession") or "").strip()
    return {"phage":row["phage_name"],"accession":accession,"host":str(host),"host_genus":str(genus),
        "cluster":value_name(row.get("pcluster"),"cluster"),"subcluster":value_name(row.get("psubcluster"),"subcluster"),
        "length":row.get("genome_length") or NA,"seq_finished":row.get("seq_finished") is True,
        "metadata_source":"PhagesDB_API","topology_evidence":str(row.get("end_type") or NA)}

def select_sample(rows, cfg):
    """Deterministic round-robin clusters, two per lineage before extra genomes.

    Fixed controls are selected first. Rank uses a seed/name hash, never outcomes.
    """
    eligible = [r for r in rows if r["host_genus"] in cfg["sampling"]["host_genera"] and r["seq_finished"] and ACCESSION.fullmatch(r["accession"])]
    by_name = {r["phage"]:r for r in eligible}
    chosen = [by_name[n] for n in cfg["sampling"]["controls"] if n in by_name]
    missing = [n for n in cfg["sampling"]["controls"] if n not in by_name]
    if missing: raise ValueError(f"Required published controls missing from eligible metadata: {missing}")
    count = int(cfg["sample_sizes"][cfg["mode"]])
    seed = cfg["seed"]
    rank = lambda s: hashlib.sha256(f"{seed}:{s}".encode()).hexdigest()
    genera = cfg["sampling"]["host_genera"]
    groups = {}
    for r in eligible:
        # Missing cluster is never one common lineage.
        lineage = r["cluster"] if r["cluster"] not in {NA,"Singleton","Unassigned","None"} else f"singleton:{r['phage']}"
        groups.setdefault((r["host_genus"],lineage),[]).append(r)
    for group in groups.values(): group.sort(key=lambda r:rank(r["phage"]))
    # Mode expansion nests the same deterministic lineage order. Pilot prioritizes
    # paired lineages, since one genome per cluster cannot address conservation.
    group_order = {g:sorted([k for k in groups if k[0]==g],key=lambda k:(-min(len(groups[k]),2),rank(k[1]))) for g in genera}
    selected_names = {r["phage"] for r in chosen}
    seen_acc = {r["accession"].split('.')[0] for r in chosen}
    for index in range(max(len(v) for v in group_order.values())):
        for genus in genera:
            if index >= len(group_order[genus]): continue
            key = group_order[genus][index]
            for r in groups[key][:2]:
                acc = r["accession"].split('.')[0]
                if r["phage"] not in selected_names and acc not in seen_acc and len(chosen)<count:
                    chosen.append(r); selected_names.add(r["phage"]); seen_acc.add(acc)
    return chosen[:count]

def run(cfg):
    offline = cfg.get("offline",False)
    raw = ROOT / "data/raw"
    metadata_file = raw / "phagesdb/eligible_metadata.json"
    frozen = ROOT / f"data/manifests/{cfg['mode']}_inputs.tsv"
    selection_snapshot=ROOT/'data/manifests/eligibility_snapshot.tsv'
    provenance_path=ROOT/'data/manifests/phagesdb_snapshot.json'
    if provenance_path.exists() and json.loads(provenance_path.read_text()).get('selection_snapshot') and not selection_snapshot.exists():
        raise FileNotFoundError('Frozen eligibility snapshot missing; restore it from version control')
    if not metadata_file.exists() and not selection_snapshot.exists():
        if offline: raise FileNotFoundError("Frozen PhagesDB metadata is missing")
        # The API silently ignores seq_finished=True. The official sequenced-only
        # TSV avoids fetching detailed metadata for all unsequenced isolates.
        url = 'https://phagesdb.org/data/?set=seq&type=full'
        page = download(url,raw/'phagesdb/sequenced_metadata.tsv',cfg,maximum=12_000_000)
        with page.open(encoding='utf-8-sig',newline='') as stream:
            all_rows=list(csv.DictReader(stream,delimiter='\t'))
        required={'Phage Name','Host','Cluster','Subcluster','Accession #','Finished Sequence?'}
        if not all_rows or not required <= all_rows[0].keys(): raise ValueError('Unexpected PhagesDB sequenced TSV schema')
        normalized=[{'phage':r['Phage Name'],'accession':r['Accession #'].strip(),'host':r['Host'],
            'host_genus':r['Host'].split()[0] if r['Host'].strip() else NA,
            'cluster':r['Cluster'] if r['Cluster'] not in {'','None'} else NA,
            'subcluster':r['Subcluster'] if r['Subcluster'] not in {'','None'} else NA,
            'length':r.get('Genome Length(bp)',NA),'seq_finished':r['Finished Sequence?']=='True',
            'metadata_source':'PhagesDB_sequenced_TSV','topology_evidence':r.get('End Type',NA)} for r in all_rows]
        pages=[{'url':url,'path':str(page.relative_to(ROOT)).replace('\\','/'),'sha256':sha256(page)}]
        write_json(metadata_file,normalized)
        write_json(ROOT/"data/manifests/phagesdb_snapshot.json",{"retrieved_utc":utcnow(),"pages":pages,"rows":len(normalized),"url":url})
    rows = read_tsv(selection_snapshot) if selection_snapshot.exists() else json.loads(metadata_file.read_text(encoding="utf-8"))
    for row in rows: row['seq_finished']=str(row['seq_finished']).lower()=='true'
    snapshot=ROOT/'data/manifests/phagesdb_snapshot.json'
    provenance=json.loads(snapshot.read_text(encoding='utf-8'))
    if any('\\' in p['path'] for p in provenance['pages']):
        for p in provenance['pages']: p['path']=p['path'].replace('\\','/')
        write_json(snapshot,provenance)
    # The current public export is not a historical archive. Freeze only the
    # analysis-relevant metadata, so a future clone never substitutes live
    # classifications for the cohort's original selection frame.
    if selection_snapshot.exists():
        expected_hash=provenance.get('selection_snapshot',{}).get('sha256')
        if not expected_hash or sha256(selection_snapshot)!=expected_hash:
            raise ValueError('Frozen eligibility metadata checksum mismatch')
    else:
        rows=[r for r in rows if r['host_genus'] in cfg['sampling']['host_genera']]
        write_tsv(selection_snapshot,rows)
        provenance['selection_snapshot']={'path':'data/manifests/eligibility_snapshot.tsv','sha256':sha256(selection_snapshot),
            'rows':len(rows),'basis':'Normalized analysis fields for target host genera from the original export; no finder or contact fields'}
        for page in provenance['pages']: page['optional_archive']=True
        write_json(snapshot,provenance)
    for page in provenance['pages']:
        archive=ROOT/page['path']
        if archive.exists() and sha256(archive)!=page['sha256']:
            raise ValueError('Original PhagesDB archive checksum mismatch')
    if frozen.exists():
        manifest = read_tsv(frozen)
        expected=select_sample(rows,cfg)
        if [(r['phage'],r['accession'].split('.')[0]) for r in manifest] != [(r['phage'],r['accession'].split('.')[0]) for r in expected]:
            raise ValueError('Frozen cohort conflicts with requested sampling design; preserve the old freeze and use a new dataset workspace for a different cohort')
        for r in manifest:
            validate_accession(r["accession"],True)
            path = ROOT/r["path"]
            download(r["source_url"],path,cfg,offline,maximum=3_000_000)
            if sha256(path)!=r["sha256"]: raise ValueError(f"Checksum mismatch {path}")
        return
    if offline: raise FileNotFoundError(f"Frozen manifest missing: {frozen}")
    selected = select_sample(rows,cfg)
    manifest = []
    for row in selected:
        acc = validate_accession(row["accession"])
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db":"nuccore","id":acc,"rettype":"gbwithparts","retmode":"text","tool":"actinolysis"})
        path = download(url,raw/f"genbank/{acc}.gb",cfg,maximum=3_000_000)
        record = SeqIO.read(path,"genbank")
        validate_accession(record.id,True)
        if record.id.split('.')[0] != acc.split('.')[0]: raise ValueError(f"Accession mismatch: {acc} vs {record.id}")
        version_path = raw/f"genbank/{record.id}.gb"
        if path != version_path:
            if version_path.exists():
                if sha256(path)!=sha256(version_path):
                    raise ValueError(f'New retrieval conflicts with frozen version-specific annotation: {record.id}')
                path.unlink()
            else: path.replace(version_path)
        version_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({"db":"nuccore","id":record.id,"rettype":"gbwithparts","retmode":"text","tool":"actinolysis"})
        manifest.append({**{k:row[k] for k in ["phage","host","host_genus","cluster","subcluster","metadata_source","topology_evidence"]},
            "accession":record.id,"retrieved_utc":utcnow(),"sha256":sha256(version_path),"path":str(version_path.relative_to(ROOT)).replace('\\','/'),
            "source_url":version_url,"inclusion_status":"pending_qc","exclusion_reason":"not_applicable"})
        print(f"Acquired {row['phage']} {record.id}",flush=True)
        time.sleep(0.40)
    write_tsv(frozen,manifest)
