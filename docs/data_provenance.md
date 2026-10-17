# Data sources and processing

## Genome cohort

The analysis uses 20 complete GenBank genomes. The [cohort manifest](../data/manifests/pilot_inputs.tsv) records the accession version, retrieval date, source URL and SHA-256 checksum for each genome. Record versions and translated CDS are verified before analysis. Deposited products and feature notes are retained separately from homology-based functional assignments.

PhagesDB supplied host and cluster labels for stratification. Selection included D29 and Yvonnetastic as published architecture controls, followed by a seeded, genus-interleaved selection of paired clusters. The original eligibility frame contained 3,760 entries from the two host genera. Its retrieval record and checksum are retained in [the source manifest](../data/manifests/phagesdb_snapshot.json). Bulk metadata exports are retrieved from the source rather than included in this repository. The published 20-genome manifest fixes the cohort for downstream reproduction; replaying that cohort does not repeat selection from the complete eligibility frame.

Cluster labels describe the source classification at the acquisition date. They are ancestry proxies rather than independently reconstructed clades. PhagesDB end-type text is distinct from topology in the frozen GenBank record; circular permutation does not imply adjacency across record ends.

## Functional references

The original panel contains characterized Ms6 and D29 proteins, WXIN LysB and a D29 structural-protein specificity control. Ms6 AF319619 is a partial lysis locus and is used as a reference rather than as a complete cohort genome. Yvonnetastic provides an annotation-level architecture comparison.

The reference-expansion analyses use 21 assay-anchored references and 34 references after including published profile predictions. Generic ester-hydrolysis assays and computational predictions remain distinct evidence tiers. Accessions, protein versions, assay descriptions and primary citations are recorded in the reference catalogues. A missing annotation or failed match does not establish gene absence.

## Data availability

GenBank sequences are retrieved by accession version through NCBI EFetch. Protein controls and the three RCSB specificity challenges have source URLs and sequence checksums. Source records, local software and intermediate caches are downloaded or generated during execution. Original software is available under the [MIT license](../LICENSE); source datasets retain their original terms.

Analysis records include software versions, scientific settings and stage output hashes. Reproduction compares the scientific tables and figures by checksum. [Methods](methods.md), [validation](validation.md) and [reference-expansion results](amendment_results.md) describe the analytical scope and limitations.

PhagesDB: Russell and Hatfull (2017), [doi:10.1093/bioinformatics/btw711](https://doi.org/10.1093/bioinformatics/btw711). Genome and experimental sources are listed in the [literature ledger](../references/literature_audit.tsv).
