"""Frozen metadata replay must never replace history with a live database."""
import hashlib
import json
import pytest
from actinolysis import acquisition,validation
from actinolysis.common import sha256,write_json,write_tsv

def prepared(tmp_path,monkeypatch):
    monkeypatch.setattr(acquisition,'ROOT',tmp_path)
    cfg={'mode':'pilot','seed':1,'offline':False,'sample_sizes':{'pilot':2},
        'sampling':{'host_genera':['Gordonia'],'controls':['SyntheticA']}}
    rows=[{'phage':name,'accession':acc,'host':'Gordonia synthetic_test_only','host_genus':'Gordonia',
        'cluster':'TEST_ONLY','subcluster':'not_available','length':'12000','seq_finished':True,
        'metadata_source':'synthetic_test_only','topology_evidence':'linear'}
        for name,acc in [('SyntheticA','SYN000001'),('SyntheticB','SYN000002')]]
    selected=acquisition.select_sample(rows,cfg)
    snapshot=tmp_path/'data/manifests/eligibility_snapshot.tsv'
    write_tsv(snapshot,rows)
    provenance=tmp_path/'data/manifests/phagesdb_snapshot.json'
    write_json(provenance,{'selection_snapshot':{'path':str(snapshot.relative_to(tmp_path)),'sha256':sha256(snapshot)},
        'pages':[{'path':'data/raw/metadata.tsv','sha256':'a'*64,'optional_archive':True}]})
    payload=b'SYNTHETIC TEST ONLY record payload for acquisition checksum tests'
    manifest=tmp_path/'data/manifests/pilot_inputs.tsv'
    write_tsv(manifest,[{'phage':r['phage'],'accession':r['accession']+'.1','path':f"data/raw/{r['accession']}.gb",
        'source_url':f"https://example.invalid/{r['accession']}.1",'sha256':hashlib.sha256(payload).hexdigest()} for r in selected])
    requests=[]
    def fake_download(url,path,cfg,offline=False,maximum=None):
        assert 'phagesdb' not in url.lower()
        requests.append(url);path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(payload);return path
    monkeypatch.setattr(acquisition,'download',fake_download)
    return cfg,snapshot,provenance,manifest,requests

def test_fresh_replay_uses_frozen_frame_without_live_metadata(tmp_path,monkeypatch):
    cfg,snapshot,provenance,manifest,requests=prepared(tmp_path,monkeypatch)
    before=(provenance.read_bytes(),manifest.read_bytes())
    acquisition.run(cfg)
    assert len(requests)==2
    assert before==(provenance.read_bytes(),manifest.read_bytes())
    assert not (tmp_path/'data/raw/phagesdb/eligible_metadata.json').exists()
    checks=validation.audit_checksums(tmp_path)
    assert sum(c['status']=='pass' for c in checks)==3
    assert sum(c['status']=='unavailable_optional' for c in checks)==1

def test_frozen_metadata_change_is_rejected_before_download(tmp_path,monkeypatch):
    cfg,snapshot,provenance,manifest,requests=prepared(tmp_path,monkeypatch)
    snapshot.write_text(snapshot.read_text().replace('TEST_ONLY','MODIFIED'))
    with pytest.raises(ValueError,match='eligibility metadata checksum'):
        acquisition.run(cfg)
    assert not requests

def test_missing_frozen_frame_cannot_be_refreshed_silently(tmp_path,monkeypatch):
    cfg,snapshot,provenance,manifest,requests=prepared(tmp_path,monkeypatch)
    snapshot.unlink()
    with pytest.raises(FileNotFoundError,match='Frozen eligibility'):
        acquisition.run(cfg)
    assert not requests
