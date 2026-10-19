# Validation

This document reports validation of the archived 20-genome pilot with its original 10-reference panel. Later reference-panel checks are reported with their respective analyses; the [50-genome exploratory results](exploratory_screen_results.md) provide current detection, architecture and reproduction outcomes.

Validation distinguishes source integrity, software behavior and biological-control recovery.

The frozen offline replay completed stages 1 through 9 and generated six figures. All 39 audited scientific artifacts (tables and figure files) reproduced with identical SHA-256 hashes. Source integrity passed while biological validation remained incomplete.

## Source and sequence integrity

The observed validation report contains 42 passing checksum rows covering the 20 cohort GenBank records, repeated reference-source checks, the metadata export and ten reference protein sequences. These are audit rows, not 42 distinct genomes or files. All 20 selected genomes passed implemented QC and all 2,460 CDS translations verified. [Validation report](../results/reports/validation.json), [QC report](../results/reports/quality_control.json).

Reference extraction resolved the following source identities without guessing version suffixes:

| Reference | Resolved source | Resolved protein |
|---|---|---|
| D29 LysA | NC_001900.2 | NP_046825.1 |
| D29 LysB | NC_001900.2 | NP_046827.1 |
| D29 holin | NC_001900.2 | NP_046826.1 |
| Ms6 Gp1 | AF319619.1 | AAG48317.1 |
| Ms6 LysA | AF319619.1 | AAG48318.1 |
| Ms6 LysB | AF319619.1 | AAG48319.1 |
| Ms6 Gp4 | AF319619.1 | AAG48320.1 |
| Ms6 Gp5 | AF319619.1 | AAG48321.1 |
| WXIN LysB | OR813930.2 | WPH57662.1 |
| D29 structural negative control | NC_001900.2 | NP_046832.1 |

The reference manifest records both source and protein hashes. Explicitly versioned targets must match exactly. A corrected scientific evidence tier or category cannot silently reuse an older manifest merely because its sequence is unchanged. [Reference manifest](../data/manifests/reference_manifest.tsv).

## Biological controls

| Control | Independent validation | Explanation |
|---|---|---|
| D29 linked lysis arrangement | Pass | Exact supported LysA, holin and LysB sequences recovered in the cohort record. |
| Ms6 linked arrangement | Unresolved in complete-genome cohort | AF319619 is a reference locus; an accessioned complete Ms6 genome was unavailable in the audited metadata. |
| Yvonnetastic displaced arrangement | Unresolved independently | The locus includes two deposited LysA components and a LysB note, but the independent reference panel does not support those functions under the primary criteria. |

Yvonnetastic's gene 161, AMS02705.1, has a serine hydrolase product and a Lysin B note. The original product-only parsing omitted that note. Its preservation is an annotation fix, not a new homology or biochemical result. Annotation-only component-pair separation is reported separately for this complex architecture. It cannot satisfy independent positive-control validation.

A targeted diagnostic for genes 47, 48 and 161 uses local BLOSUM62 alignments with gap-open -10 and extension -0.5. Two independent analytical tests check an exact AA alignment score of eight and a gap-containing example with score 76.5 and identity 8/12. The diagnostic is selected because this is a known unresolved control; raw scores are uncalibrated, no E-values are supplied and no assignment is promoted. [Diagnostic report](../results/reports/control_alignment_diagnostics.json).

The current biological validation status is incomplete, despite passing source checks. Zero candidates in the primary ranking means no qualifying supported separated architecture was available. It does not mean that no displaced lysis architecture exists.

## Software checks

Tests cover coordinates, overlaps, reverse strands, circular boundaries, translations, ambiguity, sampling, protein families, matching, exact null calibration, multiplicity correction, source integrity, resumption and figure-data consistency. Synthetic fixtures are test-only and are excluded from biological result tables. [Test report](../results/reports/tests.json).

## Conditions for stronger validation

A defensible expanded study requires independently supported recovery of divergent separated controls, enough eligible separated genome pairs across informative lineages, suitable homologous backgrounds and evidence that the ancestry and sign-symmetry assumptions are reasonable. The present pilot has completed checksum and output reconciliation, software tests and visual inspection of its six figures; these checks must be repeated after scientific changes.

