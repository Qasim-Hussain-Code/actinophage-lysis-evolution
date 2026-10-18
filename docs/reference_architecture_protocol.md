# Architecture survey of assay-reference source genomes

Design recorded on 18 October 2026 before extracting complete-source lysis architectures. The survey evaluates the complete genomes supplying the fourteen assay-supported B references. It uses source records already acquired for reference curation and is analysed separately from the frozen 20-genome cohort.

## Source selection

Include every unique resolved nucleotide source supplying a B protein in the fourteen-reference assay panel. Retain all source aliases. A genome must pass the unchanged sequence, coordinate and translation QC and be explicitly described as a complete phage genome or complete sequence. Partial reference loci, including the Ms6 lysis locus, remain reference proteins but are excluded as whole-genome controls. Report every exclusion without replacement by a more favorable accession or annotation.

Source versions and checksums come from the existing 39-reference manifest. Hosts are read from the source record; missing host or ancestry metadata remain explicit missing values. Selection does not use lysis-gene positions or the resulting architecture. This control-selected collection does not estimate population architecture frequencies.

## Fixed identification and architecture

Search all eligible CDS against the unchanged 39-reference catalogue using the same MMseqs2 parameters and 40% identity, 70% bidirectional coverage and E-value 1e-6 assignment criteria. Retain the existing evidence-tier caps. Run the existing architecture classifier at separation cutoffs 2, 4, 6 and 10. Preserve all component pairs and complex or unresolved classifications.

For each source, record recovery of the assayed B accession, supported A and B counts, architecture category, CDS-order separation, base-pair gap and evidence scope. Identify any separated arrangement supported under the existing high/moderate operational rules, with the assay-substrate limitations retained. Source self-matches are expected positive controls, not independent sensitivity estimates. The survey cannot validate native mycolate cleavage or infection requirements beyond the underlying assays.

## Interpretation and integrity

The survey tests whether the reference collection supplies additional inspectable architecture controls and whether incomplete A identification persists even among assay-source genomes. It performs no neighborhood inference, family clustering, phylogenetic reconstruction or population hypothesis test. It does not change the original cohort, functional calls, primary H2 result or expansion gate.

Raw records, code and derived outputs are checksum-authenticated. Only acquisition, QC, reference identification, evidence extraction and architecture stages are required. A source-control survey is distinct from completing the nine-stage biological study. All intermediate and full population-cohort modes remain gated.

## Input-check amendment, 18 October 2026

Inspection of the first source-selection ledger exposed a title-only phage check that excluded D29 RefSeq NC_001900.2: its definition says “Mycobacterium virus D29,” whereas its organism explicitly says “Mycobacterium phage D29” and its division is PHG. The record check now accepts an explicit phage designation in either the definition or organism, or the PHG division. [NCBI defines PHG as bacteriophage sequences](https://www.ncbi.nlm.nih.gov/genbank/release/238/). Completeness, size, coordinates and translation criteria are unchanged. The original twenty-genome cohort remains frozen. This correction was specified from the input ledger before interpreting source-genome architectures.

Raw host strings are preserved. A source-host token is assigned to the host-genus field only when it explicitly names Mycobacterium, Mycolicibacterium or Gordonia. A strain identifier such as GD35 remains in the source-host field with genus unavailable; it is not treated as a genus. No genus is inferred from the viral name. The first ledger and its code snapshot are retained separately; all fourteen source records are reconsidered under the amended input check without architecture-based replacement.
