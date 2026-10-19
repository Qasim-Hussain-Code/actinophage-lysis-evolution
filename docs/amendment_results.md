# Reference-expansion results

This archived analysis uses the original 20-genome cohort and the 10-, 21- and 34-reference panels. The later [50-genome exploratory screen](exploratory_screen_results.md) uses the 39-reference panel and a separate stages 1–5 protocol.

The post-pilot amendment kept the same 20 genomes, 11 PhagesDB cluster labels, 2,460 verified translations and numerical homology thresholds. All nine WSL stages completed for both expanded panels. The original scientific outputs remain separately preserved.

| Panel | References | Genomes with supported A | Genomes with supported B | Genomes with both | Supported proteins |
|---|---:|---:|---:|---:|---:|
| Original | 10 | 1 | 1 | 1 | 3 |
| Expanded assay panel | 21 | 5 | 1 | 1 | 7 |
| Assay plus published profile predictions | 34 | 10 | 1 | 1 | 12 |

Supported means the operational high/moderate homology tiers, not a calibrated probability of activity. Generic-lipolytic and published-prediction references are capped at moderate, including exact matches. The five additional prediction-supported assignments are Lysin A: Yvonnetastic, Colossa, VanLee, Jollymon and Stillion. The four assay-supported additions are Lysin A in Balomoji, Contagion, ZenTime222 and Bernal13. No original supported call was lost.

The larger panel partly resolves the functional-coverage limitation: A coverage increases from 1/20 to 10/20 genomes. It does not resolve the B limitation. D29 remains the sole supported A/B pair. Every panel therefore retains one linked and 19 unresolved primary architectures, with zero supported separated genomes and zero eligible primary H2 lineage contrasts. Effect, confidence interval and p value remain unestimable. Expansion of the full nine-stage conservation analysis remains gated; the later 50-genome detection screen does not bypass that gate.

## Independent query checks

Source-held-out searches exclude each reference and all references from the same resolved nucleotide-accession base before interpretation. Under the expanded combined panel, 4/7 experimentally anchored A queries, 3/3 original B queries and 6/6 generic-lipolytic B queries recover their category from another source. Among published predictions, 1/10 A and 0/3 B queries do so. Holin and chaperone queries do not recover across sources. These are cross-source sequence-transfer observations in a small selected reference set, not population functional sensitivity estimates; closely related sequences are not independent validation families.

None of three related cellular esterases receives a qualifying lysis-family assignment. The challenges include two close structural relatives from the original D29 study. Three selected challenges cannot establish a false-discovery rate, and their native mycolate inactivity has not been demonstrated by the cited sources. [Evidence review](amendment_evidence_review.md), [challenge provenance](challenge_evidence.md), [complete query ledgers](../results/amendment_v1/profile/reference_validation/query_summary.tsv).

## Matching and multipart architecture

The secondary component analysis preserves all independently supported A components and requires B to be distal to every detected A to qualify as separated. It does not force complex loci into a unique pair. No within-cluster pair in any current panel has supported B anchors in both genomes, so the availability funnel stops at identification before background matching. This is a measured feasibility result rather than a statement that B genes are absent. The original annotation-only candidate arrangements remain uncertain.

The full reported cohort-reference alignments and a detection-status row for every CDS are retained in each panel's `search_evidence/` folder. Partial proteins and negative controls cannot gain acceptance through that supplement. The inherited architecture and statistical settings are unchanged. All amended analyses are exploratory; none replaces the original primary result.

[Panel comparison](../results/amendment_v1/comparison/panel_summary.tsv), [changed assignments](../results/amendment_v1/comparison/supported_protein_changes.tsv), [comparison figure](../results/amendment_v1/comparison/reference_expansion.png), [secondary-method protocol](component_amendment.md).
