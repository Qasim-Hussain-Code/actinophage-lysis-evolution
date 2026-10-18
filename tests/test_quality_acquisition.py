"""Analytical acquisition/QC tests; every sequence here is synthetic test-only."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from actinolysis import acquisition, quality
from actinolysis.common import NA, read_tsv, sha256, write_tsv


FIXTURES = Path(__file__).parent / "fixtures"
spec = importlib.util.spec_from_file_location("build_test_fixture", FIXTURES / "build_fixture.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def test_fixture_rebuild_is_deterministic_and_explicitly_synthetic(tmp_path):
    path = tmp_path / "rebuilt.gb"
    with path.open('w',encoding='ascii',newline='\n') as stream:
        SeqIO.write(builder.make_record(), stream, "genbank")
    assert path.read_bytes() == (FIXTURES / "synthetic_test_only.gb").read_bytes()
    record = SeqIO.read(path, "genbank")
    assert len(record) == 12000
    assert "Synthetic test-only" in record.description
    assert record.annotations["organism"] == "synthetic construct"


def test_nested_metadata_and_explicit_missing_values():
    normalized = acquisition.normalize_metadata({
        "phage_name": "Synthetic", "host_strain": {"name": "Gordonia terrae 3612"},
        "host_genus": {"genus": "Gordonia"}, "genbank_accession": "AB123456.1",
        "pcluster": {"cluster": "CZ"}, "psubcluster": {"subcluster": "CZ1"},
        "seq_finished": True, "genome_length": 12000})
    assert normalized["host"] == "Gordonia terrae 3612"
    assert normalized["host_genus"] == "Gordonia"
    assert normalized["cluster"] == "CZ" and normalized["subcluster"] == "CZ1"
    assert normalized["seq_finished"] is True
    missing = acquisition.normalize_metadata({"phage_name": "Synthetic", "host": "Mycobacterium smegmatis"})
    assert missing["host_genus"] == "Mycobacterium"
    assert missing["cluster"] == NA and missing["subcluster"] == NA
    assert missing["seq_finished"] is False


@pytest.mark.parametrize("accession", ["NC_001416.1", "AF454123.2", "KJ000001", "OP123456.1"])
def test_accession_validation_accepts_nucleotide_forms(accession):
    assert acquisition.validate_accession(accession) == accession


@pytest.mark.parametrize("accession", ["", "not an accession", "../AF123456.1", "AF1234.1", "af123456.1", "AF123456.x", "AF123456.1;ls"])
def test_accession_validation_rejects_malformed_values(accession):
    with pytest.raises(ValueError, match="Invalid accession"):
        acquisition.validate_accession(accession)


def test_frozen_accession_requires_version():
    with pytest.raises(ValueError, match="accession-version"):
        acquisition.validate_accession("AF123456", version_required=True)


def test_extract_coordinate_parts_strands_offset_and_partial():
    record = SeqIO.read(FIXTURES / "synthetic_test_only.gb", "genbank")
    metadata, rows = quality.extract_record(record, {"phage": "Synthetic", "accession": record.id})
    by_tag = {r["locus_tag"].removeprefix("SYNTHETIC_TEST_ONLY_"): r for r in rows}
    assert metadata["length"] == 12000 and metadata["cds_count"] == 5
    assert metadata["topology"] == "circular"
    assert by_tag["alt_start"]["translation"] == "MAA"
    assert by_tag["alt_start"]["translation_status"] == "verified"
    assert by_tag["reverse"]["strand"] == -1 and by_tag["reverse"]["translation"] == "MPG"
    assert (by_tag["reverse"]["start"], by_tag["reverse"]["end"]) == (400, 412)
    assert by_tag["codon_offset"]["translation_status"] == "verified"
    assert by_tag["partial"]["partial"] == "true"
    assert by_tag["origin_join"]["cds_index"] == 0
    assert json.loads(by_tag["origin_join"]["location_parts"]) == [[11994, 12000], [0, 6]]
    assert by_tag["origin_join"]["translation_status"] == "verified"
    assert all(r["partial"] == "false" for name, r in by_tag.items() if name != "partial")


def _translation(sequence, protein, **qualifiers):
    record = SeqRecord(Seq(sequence))
    feature = SeqFeature(FeatureLocation(0, len(sequence), strand=1), type="CDS",
                         qualifiers={"translation": [protein], **qualifiers})
    return quality.translation_check(feature, record)


def test_translation_missing_mismatch_pseudogene_and_exception():
    assert _translation("ATGGCTTAA", "")[1] == "missing_translation"
    assert _translation("ATGGCTTAA", "MP")[1] == "translation_mismatch"
    assert _translation("ATGGCTTAA", "MA", pseudo=[""])[1] == "pseudogene"
    assert _translation("ATGGCTTAA", "MP", transl_except=["(pos:4..6,aa:OTHER)"])[1] == "translation_exception_unresolved"
    assert _translation("ATGGCTTA", "MA")[1] == "incomplete_codon"


def test_arbitrary_first_codon_is_not_valid_alternative_initiation():
    assert _translation("GTGGCTGCTTAA", "MAA")[1] == "verified"
    assert _translation("TTTGCTGCTTAA", "MAA")[1] == "translation_mismatch"


def test_multiple_terminal_stops_are_not_silently_removed():
    assert _translation("ATGTAATAA", "M")[1] == "translation_mismatch"

def test_biological_translation_exception_is_unresolved_not_genome_mismatch():
    """A recoded stop can be legitimate and needs separate handling."""
    assert _translation("ATGTGATAA", "MU",
                        transl_except=["(pos:4..6,aa:Sec)"])[1] == "translation_exception_unresolved"


def test_canonical_duplicate_hash_and_sketch_reverse_complement():
    sequence = "AACCGTTAGGCTACGTGAACTTGCAAAGC"
    reverse = str(Seq(sequence).reverse_complement())
    assert quality.canonical_hash(sequence) == quality.canonical_hash(reverse)
    assert quality.canonical_hash(sequence) != quality.canonical_hash(sequence[:-1] + "T")
    assert quality.sketch(sequence, k=7, size=100) == quality.sketch(reverse, k=7, size=100)
    assert len(quality.sketch(sequence, k=7, size=4)) <= 4


def test_quality_rejects_noncomplete_record_and_excess_ambiguity():
    record = builder.make_record()
    record.description = "Synthetic test-only phage partial genome"
    with pytest.raises(ValueError, match="not_explicitly_complete"):
        quality.extract_record(record, {})
    record = builder.make_record()
    record.seq = Seq("N" * 121 + str(record.seq)[121:])
    with pytest.raises(ValueError, match="more_than_one_percent"):
        quality.extract_record(record, {})


@pytest.mark.parametrize("organism,division", [
    ("Synthetic test-only bacteriophage", "VRL"),
    ("Synthetic test-only virus", "PHG"),
])
def test_virus_title_can_use_explicit_phage_record_metadata(organism, division):
    record = builder.make_record()
    record.description = "Synthetic test-only virus, complete genome"
    record.annotations.update(organism=organism, data_file_division=division)
    metadata, _ = quality.extract_record(record, {})
    assert metadata["inclusion_status"] == "included"
    record.description = "Synthetic test-only virus, partial genome"
    with pytest.raises(ValueError, match="not_explicitly_complete"):
        quality.extract_record(record, {})


def test_virus_title_without_phage_classification_is_rejected():
    record = builder.make_record()
    record.description = "Synthetic test-only virus, complete genome"
    record.annotations.update(organism="Synthetic test-only virus", data_file_division="VRL")
    with pytest.raises(ValueError, match="not_phage_record"):
        quality.extract_record(record, {})


def _sample_row(name, accession, genus, cluster):
    return dict(phage=name, accession=accession, host_genus=genus, cluster=cluster, seq_finished=True)


def _sample_cfg():
    return {"mode": "pilot", "seed": 1701, "sample_sizes": {"pilot": 8},
            "sampling": {"host_genera": ["Mycobacterium", "Gordonia"],
                         "controls": ["ControlM", "ControlG"], "max_per_cluster": 4}}


def test_sampling_deterministic_controls_genera_and_accession_deduplication():
    rows = [_sample_row("ControlM", "AF100001.1", "Mycobacterium", "A"),
            _sample_row("ControlG", "AF100002.1", "Gordonia", "C")]
    for index in range(3, 18):
        genus = "Mycobacterium" if index % 2 else "Gordonia"
        rows.append(_sample_row(f"Synthetic{index}", f"AF{100000 + index}.1", genus, f"C{index // 2}"))
    rows += [_sample_row("DuplicateVersion", "AF100001.2", "Mycobacterium", "A"),
             _sample_row("WrongHost", "AF199999.1", "Escherichia", "X")]
    cfg = _sample_cfg()
    result = acquisition.select_sample(rows, cfg)
    assert result == acquisition.select_sample(list(reversed(rows)), cfg)
    assert len(result) == 8
    assert [r["phage"] for r in result[:2]] == ["ControlM", "ControlG"]
    assert {r["host_genus"] for r in result} == {"Mycobacterium", "Gordonia"}
    assert len({r["accession"].split(".")[0] for r in result}) == len(result)
    assert "WrongHost" not in {r["phage"] for r in result}


def test_sampling_missing_control_is_explicit_failure():
    with pytest.raises(ValueError, match="controls missing"):
        acquisition.select_sample([], _sample_cfg())


def test_download_offline_uses_existing_input_and_requires_missing_input(tmp_path):
    existing = tmp_path / "frozen.gb"
    existing.write_text("synthetic test only")
    assert acquisition.download("https://invalid.example/unused", existing, {}, offline=True) == existing
    with pytest.raises(FileNotFoundError, match="Offline input missing"):
        acquisition.download("https://invalid.example/unused", tmp_path / "missing.gb", {}, offline=True)


def test_quality_stage_excludes_duplicate_and_fasta_excludes_partial(tmp_path, monkeypatch):
    original = builder.make_record()
    duplicate = copy.deepcopy(original)
    duplicate.id = "SYN000002.1"
    rows = []
    for record in (original, duplicate):
        path = tmp_path / "data/raw" / f"{record.id}.gb"
        path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(record, path, "genbank")
        rows.append(dict(accession=record.id, phage="Synthetic test only", path=str(path.relative_to(tmp_path)),
                         sha256=sha256(path), inclusion_status="pending_qc", exclusion_reason="not_applicable"))
    write_tsv(tmp_path / "data/manifests/pilot_inputs.tsv", rows)
    monkeypatch.setattr(quality, "ROOT", tmp_path)
    quality.run({"mode": "pilot"})
    included = read_tsv(tmp_path / "data/processed/genomes.tsv")
    excluded = read_tsv(tmp_path / "results/tables/exclusions.tsv")
    assert len(included) == 1 and len(excluded) == 1
    assert excluded[0]["exclusion_reason"] == "exact_sequence_duplicate_of:SYN000001.1"
    fasta = list(SeqIO.parse(tmp_path / "data/processed/proteins.faa", "fasta"))
    assert len(fasta) == 4
    assert not any("SYN000002.1" in record.id for record in fasta)


def test_validation_never_counts_annotation_recovery_as_independent_pass():
    from actinolysis.validation import assess_controls
    control = {"phage": "Synthetic", "expected_architecture": "separated",
               "reference_url": "https://invalid.example/test-only", "evidence_basis": "test-only expectation"}
    architecture = {"phage": "Synthetic", "accession": "SYN000001.1", "architecture": "unresolved",
                    "classification_reason": "incomplete_identification_not_evidence_of_absence"}
    annotation = {**architecture, "architecture": "separated"}
    result = assess_controls([control], [architecture], {}, [annotation])[0]
    assert result["status"] == "unresolved"
    assert result["annotation_recovery_status"] == "expected_arrangement_recovered_from_annotation_only"
    assert "no_traceable_independent" in result["reason"]


def test_complex_annotation_component_recovery_requires_every_pair_and_stays_unresolved():
    from actinolysis.validation import assess_controls
    control = {"phage": "Synthetic", "expected_architecture": "separated",
               "reference_url": "https://invalid.example/test-only", "evidence_basis": "test-only expectation"}
    architecture = {"phage": "Synthetic", "accession": "SYN000001.1", "architecture": "unresolved",
                    "classification_reason": "incomplete_identification_not_evidence_of_absence"}
    annotation = {**architecture, "architecture": "complex", "lysin_a_count": 2, "lysin_b_count": 1}
    pairs = [dict(phage="Synthetic", scope="annotation", lysin_a=a, lysin_b="b", intervening_cds=distance,
                  gap_bp=50000 + distance, confidence="annotation_only") for a, distance in [("a1", 10), ("a2", 12)]]
    recovered = assess_controls([control], [architecture], {}, [annotation], pairs)[0]
    assert recovered["status"] == "unresolved"
    assert recovered["annotation_architecture"] == "complex"
    assert recovered["annotation_recovery_status"] == "annotation_only_component_separation_recovered_in_complex_architecture"
    assert len(recovered["annotation_component_pairs"]) == 2
    incomplete = assess_controls([control], [architecture], {}, [annotation], pairs[:1])[0]
    assert incomplete["annotation_recovery_status"] == "annotation_arrangement_not_resolved"
    contradicted = assess_controls([control], [architecture], {}, [annotation], [pairs[0], {**pairs[1], "intervening_cds": 2}])[0]
    assert contradicted["annotation_recovery_status"] == "annotation_arrangement_not_resolved"


def test_validation_checksum_audit_detects_mutation_and_path_escape(tmp_path):
    from actinolysis.validation import audit_checksums
    frozen = tmp_path / "data/raw/test_only.txt"
    frozen.parent.mkdir(parents=True)
    frozen.write_text("synthetic test only")
    write_tsv(tmp_path / "data/manifests/pilot_inputs.tsv", [{"path": "data/raw/test_only.txt", "sha256": sha256(frozen)}])
    assert audit_checksums(tmp_path)[0]["status"] == "pass"
    frozen.write_text("synthetic test only altered")
    assert audit_checksums(tmp_path)[0]["reason"] == "checksum_mismatch"
    write_tsv(tmp_path / "data/manifests/pilot_inputs.tsv", [{"path": "../outside.txt", "sha256": "0" * 64}])
    assert audit_checksums(tmp_path)[0]["reason"] == "manifest_path_outside_project"


def test_validation_reference_source_and_protein_checksums(tmp_path):
    import hashlib
    from actinolysis.validation import audit_checksums
    source = tmp_path / "data/raw/references/SYN000001.1.gb"
    source.parent.mkdir(parents=True)
    source.write_text("synthetic test only source")
    fasta = tmp_path / "data/processed/references.faa"
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">synthetic_test_only\nMAA\n")
    write_tsv(tmp_path / "data/manifests/reference_manifest.tsv", [{
        "reference_id": "synthetic_test_only", "source_accession": "SYN000001.1", "source_sha256": sha256(source),
        "protein_sha256": hashlib.sha256(b"MAA").hexdigest()}])
    checks = audit_checksums(tmp_path)
    assert len(checks) == 2 and all(row["status"] == "pass" for row in checks)
    fasta.write_text(">synthetic_test_only\nMAL\n")
    assert audit_checksums(tmp_path)[1]["reason"] == "checksum_mismatch"


def test_validation_independent_control_can_pass_or_contradict_expectation():
    from actinolysis.validation import assess_controls
    control = {"phage": "Synthetic", "expected_architecture": "separated",
               "reference_url": "https://invalid.example/test-only", "evidence_basis": "test-only expectation"}
    architecture = {"phage": "Synthetic", "accession": "SYN000001.1", "architecture": "separated", "lysin_a": "a", "lysin_b": "b"}
    evidence = {pid: [{"category": category, "confidence": "moderate", "source": "curated_reference_homology", "reference_id": "synthetic_test_reference"}]
                for pid, category in [("a", "lysin_a"), ("b", "lysin_b")]}
    assert assess_controls([control], [architecture], evidence)[0]["status"] == "pass"
    architecture["architecture"] = "closely_linked"
    assert assess_controls([control], [architecture], evidence)[0]["status"] == "fail"


def test_validation_stage_preserves_unresolved_controls_and_audits_frozen_inputs(tmp_path, monkeypatch):
    from actinolysis import validation
    monkeypatch.setattr(validation, "ROOT", tmp_path)
    frozen = tmp_path / "data/raw/synthetic_test_only.txt"
    frozen.parent.mkdir(parents=True)
    frozen.write_text("synthetic test only")
    write_tsv(tmp_path / "data/manifests/pilot_inputs.tsv", [{"path": "data/raw/synthetic_test_only.txt", "sha256": sha256(frozen)}])
    metadata = {"accession": "SYN000001.1", "phage": "Synthetic test only", "host_genus": "test_only",
                "cluster": "test_only", "similarity_group": "SYN000001.1", "representative": "true"}
    architecture = {**metadata, "architecture": "separated", "confidence": "high", "lysin_a_count": 1,
                    "lysin_b_count": 1, "lysin_a": "a", "lysin_b": "b", "intervening_cds": 10,
                    "gap_bp": 1000, "partial_lysis_count": 0}
    evidence = [{"protein_id": pid, "category": category, "confidence": "high", "source": "curated_reference_homology",
                 "reference_id": "synthetic_test_reference"} for pid, category in [("a", "lysin_a"), ("b", "lysin_b")]]
    write_tsv(tmp_path / "data/processed/genomes.tsv", [metadata])
    write_tsv(tmp_path / "data/processed/architectures.tsv", [architecture])
    write_tsv(tmp_path / "data/processed/lysis_evidence.tsv", evidence)
    cfg = {"mode": "pilot", "architecture": {"separation_threshold": 4}}
    report = validation.run(cfg)
    assert report["status"] == "incomplete" and report["checksum_status"] == "pass"
    assert all(row["status"] == "unresolved" for row in report["controls"])
    assert report["candidate_count"] == 1
    assert read_tsv(tmp_path / "results/tables/candidates.tsv")[0]["rank"] == "1"
    frozen.write_text("synthetic test only altered")
    with pytest.raises(ValueError, match="checksum validation failed"):
        validation.run(cfg)
    failed = json.loads((tmp_path / "results/reports/validation.json").read_text())
    assert failed["status"] == "fail"
