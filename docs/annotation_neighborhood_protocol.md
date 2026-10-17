# Annotation-conditioned neighborhood sensitivity

Design recorded on 18 October 2026 before computing annotation-conditioned conservation outcomes. The independent high/moderate analysis has no eligible H2 contrast. This descriptive sensitivity measures neighborhood conservation conditional on the deposited labels and quantifies how evidence scope affects comparison eligibility. It is separate from the primary functional analysis.

## Inputs and anchor definitions

Use the unchanged 20-genome cohort, verified CDS, 39-reference evidence catalogue and existing protein families. Annotation-only architectures come from the already generated architecture-sensitivity table. Deposited product names and notes select putative A/B components; they never become independent functional evidence. Unique, complete annotated B anchors must be in the same single-copy sequence family in both genomes. Unresolved or multiply annotated A/B architectures remain ineligible.

No functional confidence field or primary architecture is changed. The output evidence scope is `deposited_annotation_only`, including when a selected locus also has independent support. Side-by-side primary comparisons are copied from the existing high/moderate results.

## Fixed comparisons

Use the original family thresholds 0.30, 0.40 and 0.50, windows of 3, 5 and 10 flanking CDS, separation cutoffs 2, 4, 6 and 10, and all-genome/representative subsets. The base configuration is 0.40 identity, five flanking CDS and cutoff four. Linear records do not wrap; only explicitly circular records may wrap.

The scoring and matching functions remain unchanged: flanking-family Jaccard similarity and signed order are compared with single-copy homologous background anchors in the same genome pair. Match boundary completeness and linear-position bins exactly, both protein lengths within a factor of two, and family prevalence within 0.20. Exclude the target family and families containing any lysis-related evidence or deposited lysis label. Select up to 20 controls using covariates, never conservation scores; retain the original minimum of five. Exclude pairs in the same near-identical screening group.

Report every pair and exclusion, selected controls, target scores, background means and conditional contrasts. Contrasts require at least five controls. Summaries average pairs within cluster labels and then weight labels equally, with separated and linked architectures kept distinct. Sensitivities are displayed together rather than selecting a favorable configuration.

## Interpretation limits

These are descriptive contrasts for putatively annotated loci. Deposited labels can incorporate genomic context or shared annotation practices, so they cannot independently validate a functional conservation hypothesis. There are no formal p-values, bootstrap confidence intervals or claims of biochemical function, selection, rearrangement mechanism or independent origin. A conditional contrast cannot replace the unavailable primary H2 estimate. The comparison determines how much can be measured under annotation assumptions and which additional functional validation would change eligibility.
