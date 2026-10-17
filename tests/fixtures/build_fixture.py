"""Build a deterministic SYNTHETIC TEST-ONLY GenBank coordinate fixture.

This artificial sequence and all functional labels are software test inputs.
They are not downloaded phage data and must never enter scientific results.
Run ``python tests/fixtures/build_fixture.py`` to regenerate the committed .gb.
"""
from pathlib import Path
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def make_record():
    rng = random.Random(1701)
    bases = list("".join(rng.choice("ACGT") for _ in range(12000)))

    def insert(start, sequence):
        bases[start:start + len(sequence)] = sequence

    insert(100, "GTGGCTGCTTAA")  # Alternative bacterial initiation, MAA.
    insert(400, str(Seq("ATGCCCGGGTAA").reverse_complement()))  # Reverse MPG.
    insert(700, "AATGAAACCCTAA")  # codon_start 2, MKP.
    insert(1000, "ATGCCGCTGTAA")  # Fuzzy genomic start, MPL.
    insert(11994, "ATGAAA")
    insert(0, "CCCTAA")  # Joined across the circular origin, MKP.
    record = SeqRecord(Seq("".join(bases)), id="SYN000001.1", name="SYNTESTONLY",
                       description="Synthetic test-only phage complete genome; not a biological observation")
    record.annotations = {"molecule_type": "DNA", "topology": "circular", "date": "01-JAN-2000",
                          "organism": "synthetic construct", "source": "synthetic test-only data"}
    record.features = [SeqFeature(FeatureLocation(0, 12000), type="source",
                                  qualifiers={"organism": ["synthetic construct"], "host": ["test_only"]})]
    cases = [
        (FeatureLocation(100, 112, strand=1), "alt_start", "MAA", {}),
        (FeatureLocation(400, 412, strand=-1), "reverse", "MPG", {}),
        (FeatureLocation(700, 713, strand=1), "codon_offset", "MKP", {"codon_start": ["2"]}),
        (FeatureLocation(BeforePosition(1000), 1012, strand=1), "partial", "MPL", {}),
        (CompoundLocation([FeatureLocation(11994, 12000, strand=1), FeatureLocation(0, 6, strand=1)]),
         "origin_join", "MKP", {}),
    ]
    for location, tag, translation, extra in cases:
        record.features.append(SeqFeature(location, type="CDS", qualifiers={
            "locus_tag": [f"SYNTHETIC_TEST_ONLY_{tag}"], "translation": [translation],
            "product": ["synthetic test-only protein"], "transl_table": ["11"], **extra}))
    return record


def main():
    destination = Path(__file__).with_name("synthetic_test_only.gb")
    with destination.open('w',encoding='ascii',newline='\n') as stream:
        SeqIO.write(make_record(), stream, "genbank")
    return destination


if __name__ == "__main__":
    print(main())
