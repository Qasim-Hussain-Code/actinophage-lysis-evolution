"""Selection regression tests; synthetic records are explicitly test-only."""
import copy
from pathlib import Path

import pytest

from actinolysis import acquisition
from actinolysis.common import NA, load_config, read_tsv


def synthetic_row(index, cluster="TEST_ONLY"):
    return {"phage": f"SyntheticTestOnly{index}", "accession": f"SYN{index:06d}.1",
            "host_genus": "Gordonia", "cluster": cluster, "seq_finished": True}


def config(count, controls=()):
    return {"mode": "pilot", "seed": 1701, "sample_sizes": {"pilot": count},
            "sampling": {"host_genera": ["Gordonia"], "controls": list(controls)}}


def identities(rows):
    return [(r["phage"], r["accession"].split(".")[0]) for r in rows]


def test_expansion_fills_beyond_two_per_lineage_and_preserves_prefix():
    rows = [synthetic_row(i, f"TEST_ONLY_{i % 2}") for i in range(1, 13)]
    small = acquisition.select_sample(rows, config(4))
    large = acquisition.select_sample(rows, config(10))
    assert len(large) == 10
    assert large[:4] == small
    assert large == acquisition.select_sample(list(reversed(rows)), config(10))
    assert {r["cluster"] for r in large[:4]} == {"TEST_ONLY_0", "TEST_ONLY_1"}


def test_accession_aliases_do_not_exhaust_first_two_lineage_slots():
    rows = [synthetic_row(i) for i in range(1, 7)]
    rows[1]["accession"] = rows[0]["accession"].replace(".1", ".2")
    result = acquisition.select_sample(rows, config(5))
    assert len(result) == len({r["accession"].split(".")[0] for r in result}) == 5
    assert result == acquisition.select_sample(list(reversed(rows)), config(5))


def test_insufficient_unique_pool_fails_instead_of_silently_underfilling():
    rows = [synthetic_row(1), synthetic_row(2)]
    rows[1]["accession"] = rows[0]["accession"].replace(".1", ".2")
    with pytest.raises(ValueError, match="Requested 2 genomes but only 1"):
        acquisition.select_sample(rows, config(2))


def test_controls_require_distinct_accessions():
    rows = [synthetic_row(1), synthetic_row(2)]
    rows[1]["accession"] = rows[0]["accession"].replace(".1", ".2")
    with pytest.raises(ValueError, match="controls must have distinct"):
        acquisition.select_sample(rows, config(2, [r["phage"] for r in rows]))


def test_sample_size_must_accommodate_every_required_control():
    rows = [synthetic_row(1), synthetic_row(2)]
    with pytest.raises(ValueError, match="accommodate all required controls"):
        acquisition.select_sample(rows, config(1, [r["phage"] for r in rows]))


@pytest.mark.parametrize("missing", [None, "", " ", "None", "Unassigned", "Singleton"])
def test_missing_cluster_labels_are_individual_lineages(missing):
    rows = [synthetic_row(i, NA if i <= 4 else "TEST_ONLY_PAIRED") for i in range(1, 7)]
    baseline = acquisition.select_sample(rows, config(3))
    modified = copy.deepcopy(rows)
    for row in modified[:4]:
        row["cluster"] = missing
    assert acquisition.select_sample(modified, config(3)) == [
        next(r for r in modified if r["phage"] == selected["phage"]) for selected in baseline]


def test_frozen_pilot_identity_order_and_nested_requested_sizes():
    root = Path(__file__).resolve().parents[1]
    snapshot = root / "data/manifests/eligibility_snapshot.tsv"
    manifest = root / "data/manifests/pilot_inputs.tsv"
    if not snapshot.exists() or not manifest.exists():
        pytest.skip("Frozen source selection frame and pilot manifest are unavailable in this checkout")
    cfg = load_config(root / "config/analysis.yaml")
    rows = read_tsv(snapshot)
    for row in rows:
        row["seq_finished"] = row["seq_finished"].lower() == "true"
    frozen = read_tsv(manifest)
    previous = identities(frozen)
    for mode in ("pilot", "intermediate", "full"):
        cfg["mode"] = mode
        selected = acquisition.select_sample(rows, cfg)
        assert len(selected) == cfg["sample_sizes"][mode]
        current = identities(selected)
        assert current[:len(previous)] == previous
        assert len({acc for _, acc in current}) == len(current)
        assert current == identities(acquisition.select_sample(list(reversed(rows)), cfg))
        previous = current
