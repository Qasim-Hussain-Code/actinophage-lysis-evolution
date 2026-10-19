# Pilot results

This is the archived 20-genome pilot report. The later [50-genome exploratory screen](exploratory_screen_results.md) is reported separately.

Reporting date: 17 October 2026. The pilot has exposed a functional-identification feasibility limit: only D29 has an independently supported Lysin A/Lysin B pair. There are no supported separated pairs and no contributing lineage contrasts for the primary neighborhood test. This is an untestable comparison, not evidence against neighborhood conservation.

Corrected frozen offline replays completed stages 1 through 9 and generated six figures on 17 October 2026. Counts below were reconciled against the generated tables. The latest execution report records the exact run times and resource measurements. Software execution succeeded; biological validation remains incomplete and expansion is gated. Neither intermediate nor full expansion was executed.

## Cohort and quality control

Twenty accession-versioned genomes were acquired and all 20 passed the implemented record and translation checks. They contain 1,493,450 total bases and 2,460 CDS translations, all verified. Genome lengths range from 42,392 to 145,208 bp; GC content ranges from 59.748716% to 67.831099%. No genome was excluded at this QC stage.

| Cohort property | Observed count |
|---|---:|
| Mycobacterium-host phages | 11 |
| Gordonia-host phages | 9 |
| Named clusters represented | 11 |
| Clusters with two sampled genomes | 9 |
| Single sampled control members | 2 |
| Near-duplicate screening groups | 19 |

The current metadata places D29 in cluster A and Yvonnetastic in DD. Historical publications treated Yvonnetastic as a singleton; that historical description is not silently substituted for the current frozen label. One pair, MK864265.1 and MZ005680.1, has canonical 21-mer sketch Jaccard 0.998001998 and is grouped by the 0.95 screening threshold. The screen is not ANI or an inferred clade. [Dataset manifest](../results/tables/dataset_manifest.tsv), [exclusions](../results/tables/exclusions.tsv), [screened pairs](../results/tables/near_identical_pairs.tsv).

## Functional support is the limiting measurement

The panel contains nine experimental functional references and one experimentally identified virion structural protein as a specificity control. The final evidence ledger has 54 rows: three high confidence, zero moderate and 51 low confidence. The three high-confidence assignments are all exact-reference matches in D29: LysA AAC18450.1, holin AAC18451.1 and LysB AAC18452.1. The cohort's GenBank protein accessions differ from the corresponding RefSeq identifiers in the reference panel, but the aligned protein sequences match. Evidence-row counts are not counts of independent genomes.

At the prespecified 40% identity, 70% query coverage, 70% target coverage and E-value at most 1e-6 requirements, none of the other genomes has a supported A/B pair. Deposited labels and weaker reference matches remain low confidence. These are measurements of identification coverage under a small reference set, not measurements of the prevalence of lysis genes in nature. The D29 capsid negative control matched the expected capsid sequence and did not propagate a lysis assignment; that limited control does not establish specificity against all esterases or hydrolases. [Functional evidence](../results/tables/lysis_evidence.tsv), [coverage](../results/tables/annotation_coverage.tsv), [negative-control hits](../results/tables/negative_control_hits.tsv).

## Architecture depends strongly on evidence scope

At the primary cutoff of more than four intervening CDS:

| Evidence scope | Closely linked | Separated | Complex | Unresolved |
|---|---:|---:|---:|---:|
| Independent high/moderate assignments | 1 | 0 | 0 | 19 |
| Deposited annotation only | 9 | 4 | 3 | 4 |

The four annotation-only separated calls are Phrappuccino, Settecandela, VanLee and Colossa. They are not independently supported separated candidates. Multiple putative LysA components also occur, so complex cases are not forced into a unique A/B pair. Counts at alternative thresholds, and all individual pairs, remain available for review. [Primary architectures](../results/tables/architectures.tsv), [all-scope pairs](../results/tables/architecture_pairs.tsv), [sensitivity table](../results/tables/architecture_sensitivity.tsv).

D29 contains a holin that overlaps the adjacent LysA and LysB boundaries. Its principal CDS-order count is one intervening gene, while the supplemental complete-gap containment count is zero; the physical edge gap remains 418 bp. Treating the strict gap-containment count as the principal gene-separation measure would miss that intervening holin. This correction changes the measurement definition without promoting any functional assignment.

## Why the published displaced control remains unresolved

The exact frozen Yvonnetastic record, KU963248.1, contains:

| Gene | Protein | Zero-based half-open interval | Deposited evidence |
|---|---|---|---|
| 47 | AMS02591.1 | [33145, 33718) | Product and note identify LysA |
| 48 | AMS02592.1 | [33724, 34624) | Product and note identify LysA |
| 161 | AMS02705.1 | [84593, 85232) | Product says serine hydrolase; note says Lysin B |

An earlier product-only parser missed the LysB note. Correcting that omission preserves the note as low-confidence deposited evidence. It does not establish experimentally characterized LysB activity or create an independently supported pair. The two predicted LysA components require a complex architecture with explicit component pairs. Gene 47 to gene 161 has 113 intervening CDS and a 50,875-bp physical gap; gene 48 to gene 161 has 112 intervening CDS and a 49,969-bp gap. Annotation-only component separation is recovered without forcing those components into a single gene. [Frozen primary sequence record](https://www.ncbi.nlm.nih.gov/nuccore/KU963248.1), [published prediction](https://doi.org/10.1128/genomeA.00594-16).

None of these three proteins had a row in the emitted broad reference-search hits. The lack of an emitted hit does not identify which search or prefilter condition caused non-detection and cannot prove absence of a lysis function. A separate targeted diagnostic locally aligns these literature-selected proteins against same-category experimental references using BLOSUM62, gap open -10 and extension -0.5. The seven uncalibrated local alignments have identities from 19.8% to 26.8%, with variable coverage; no E-values or functional assignments are produced. These exploratory scores document sequence divergence from this panel and do not validate or refute the published functional prediction. [Diagnostic table](../results/tables/control_alignment_diagnostics.tsv).

## Neighborhood analysis produced no inferential sample

All 2,460 proteins entered study-specific clustering. The observed family counts are 1,411 at identity 0.30, 1,460 at 0.40 and 1,496 at 0.50, under coverage 0.80 and E-value 0.001. Family construction does not depend on a lysis-neighborhood label. [Memberships](../results/tables/family_memberships.tsv).

The configured windows, thresholds, confidence sets and representative subsets produce 2,592 comparison rows, but none is eligible for H2 inference. The data contain zero supported within-cluster separated pairs, zero contributing cluster effects and zero conditional H2 tests. No effect size, uncertainty interval, p-value or q-value is estimable for the primary question. Empty inferential outputs are retained rather than populated with zeros. [Comparison rows](../results/tables/neighborhood_pairs.tsv), [test ledger](../results/tables/evolutionary_tests.tsv).

H1 remains a descriptive architecture summary conditional on detectability. H3 independent reorganization remains untestable. H4 demonstrates an identification and annotation-interpretation problem, including the Yvonnetastic note omission, but cannot separate every undetected function from true gene absence. No result establishes an evolutionary mechanism.

## Feasibility decision

A larger cohort alone would not remedy the observed failure to recover a published separated control with independent functional evidence. The implemented gate records `expansion_ready: false` and prevents intermediate and full acquisition while this identification limitation is unresolved. The current study is a documented pilot feasibility result, not a completed 300-genome test of H2.

The next scientific improvement would require independently justified, provenance-tracked functional evidence with measured performance on divergent controls, then renewed evaluation of whether enough separated lineages and suitable matched backgrounds exist. The observed thresholds have not been relaxed to obtain a desired result. Any revised evidence strategy would constitute a documented methodological amendment followed by regeneration and validation.

Six generated figures describe the dataset, architecture evidence, maps, lineages, neighborhood availability and robustness. Their successful generation is a software output, not successful validation of the biological hypotheses. [Figure manifest](../results/reports/figure_manifest.json).
