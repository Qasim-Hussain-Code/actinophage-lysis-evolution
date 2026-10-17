from types import SimpleNamespace
import pytest
from actinolysis import common

def test_budget_reserve_and_cap(monkeypatch,tmp_path):
    (tmp_path/'input').write_bytes(b'x'*40)
    cfg={'resources':{'max_project_bytes':100,'effective_budget_bytes':90,'reserve_bytes':20}}
    monkeypatch.setattr(common.shutil,'disk_usage',lambda p:SimpleNamespace(free=100))
    assert common.budget_check(cfg,50,tmp_path)['project_bytes']==40
    with pytest.raises(RuntimeError): common.budget_check(cfg,51,tmp_path)
    monkeypatch.setattr(common.shutil,'disk_usage',lambda p:SimpleNamespace(free=30))
    with pytest.raises(RuntimeError): common.budget_check(cfg,11,tmp_path)
