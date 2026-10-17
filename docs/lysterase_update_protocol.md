# Assay-reference update

This exploratory follow-up adds experimentally tested reference sequences to the frozen 20-genome analysis. The design was recorded on 18 October 2026 before searching the new references against the cohort. The original pilot, two expanded panels and nine-reference profile diagnostic remain separately frozen.

## Reference selection

Table 1 of [Maneekul et al. (2026)](https://doi.org/10.1128/jb.00410-26), published on 13 October 2026, identifies six purified phage proteins tested for antimycobacterial activity. The D29 entry is checked against the existing D29 LysB sequence and retained as an alias rather than another independent reference. The other five entries are MissWhite ASZ75177.1, Bxb1 AAG59714.1, Constella AXF51531.1, Bongo YP_009604896.1 and phiGD57-1 QPO17264.1.

Each protein record must resolve to an explicitly versioned nucleotide source. The corresponding deposited CDS must have the same protein accession and amino-acid sequence, and its nucleotide translation must verify. Ambiguous source mapping or conflicting bytes stop the update. Raw records are retrieved separately; source URLs and checksums accompany the derived reference ledger.

## Functional evidence

The new evidence type, `experimental_antimycobacterial`, records activity of the purified reference protein against mycobacteria. This assay supports at most moderate LysB-family homology, including an exact target match. It does not independently establish cleavage of native mycolyl-arabinogalactan linkages, endogenous expression or a requirement for phage lysis. It cannot promote LysA, holins or other categories.

The numerical filters remain 40% amino-acid identity, 70% coverage of both proteins and E-value at most 1e-6. Partial CDS remain low confidence. At equal moderate confidence, the original experimentally characterized tier retains priority over capped evidence. No target annotation, coordinates, neighborhood score or architecture result enters functional assignment.

## Analysis and validation

The combined 34-reference panel gains only the five new proteins. All 20 genomes and 2,460 CDS are retained. All nine pipeline stages are rerun in an isolated workspace, with the same family thresholds, separation cutoffs, matching rules and inference requirements. Changes in supported assignments, component pairs, architecture categories and comparison eligibility are reported against the combined-panel baseline. Existing primary outputs are preserved.

Validation includes exact reference extraction, the D29 alias check, source-held-out query searches with exact-sequence aliases excluded, and the three existing cellular-enzyme challenges. These selected controls cannot estimate population sensitivity, specificity or a false-discovery rate. Non-detection remains unresolved rather than evidence of absence. Increasing functional coverage does not by itself establish a displaced architecture or satisfy the minimum lineage requirement for H2.

## Reproduction

The update records input, code and output checksums, software versions and numerical settings. Public results contain the scientific evidence and derived outputs. Local execution inventories are excluded.
