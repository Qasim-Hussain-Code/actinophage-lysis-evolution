# Statistical assumptions and estimands

The unit of biological replication is a phage lineage represented operationally by a usable PhagesDB cluster label. This is a coarse approximation. Genome pairs and their matched control neighborhoods reuse genomes, and are not independent observations. The implemented summaries keep that dependence visible and do not interpret large pair counts as large effective sample sizes.

## Hypotheses supported by the implementation

| Hypothesis | Implemented result | Conditions and limits |
|---|---|---|
| H1: architecture varies among evolutionary groups | Cluster-specific frequencies, resolved denominators and threshold sensitivity | Descriptive. No validated phylogeny or exchangeable architecture null is fitted. |
| H2: separated LysB neighborhoods are unusually conserved | LysB similarity minus matched homologous background, summarized equally across clusters | Conditional inference only with enough replicated clusters and credible independence and sign symmetry. |
| H3: similar separation reflects recurrent reorganization | Explicit untestable status | Similar placement alone cannot determine independent origins, direction or mechanism. Gene and whole-genome histories are required. |
| H4: annotation limitations explain apparent incompleteness | Confidence and annotation-only summaries | Descriptive. Functional non-detection cannot establish true gene loss or quantify all missed functions. |

A failure to test a hypothesis is recorded as such. It is not assigned a nonsignificant p-value or interpreted as confirmation of a null mechanism.

## Prespecified primary comparison

The primary H2 configuration uses family identity 0.40, a window of five genes on each side, separation defined by more than four intervening CDS, high/moderate functional assignments, Jaccard similarity, and all otherwise eligible genomes. The primary comparison includes only pairs where both architectures are separated. Same-near-duplicate-group pairs are excluded. Representative-only comparisons are sensitivity analyses.

For genome pair `p`, let `T_p` be LysB-neighborhood similarity and `C_p` the mean over its retained matched non-lysis controls. The pair effect is `delta_p = T_p - C_p`, in [-1, 1]. The cluster effect is the mean of eligible pair effects in that cluster. The reported study effect is the unweighted mean of those cluster effects. Consequently the estimand concerns contributing clusters in the selected cohort, not the frequency of a phenotype across all sequenced or environmental phages.

Controls share the same pair of genomes. Their anchor families must be single-copy homologs, exclude families with lysis-associated evidence, and satisfy boundary, position, length and prevalence matching. At most 20 controls are retained, and fewer than five makes a pair untestable. Control conservation is never used to select controls. Nevertheless the target family is not randomly assigned to a genomic position, so matching does not guarantee exchangeability. Residual functional, regulatory and regional constraints can remain.

## Sign-flip null and uncertainty

The two-sided sign-flip procedure treats the cluster contrasts as independent and sign-symmetric under the null. It tests whether their mean is unusually far from zero under changing the sign of whole cluster effects. This is a conditional observational null, not a randomization supplied by the sampling design. Matching alone does not justify it. Shared genes across clusters, mosaic genomes and uncertain cluster boundaries can violate independence. Skewed contrasts can violate sign symmetry.

The implementation requires at least five contributing clusters for conditional inference. It enumerates all sign combinations for at most 20 nonzero effects. Zero effects are omitted from enumeration because they change no statistic. With more than 20 nonzero effects it uses 9,999 random sign vectors and `(exceedances + 1)/(draws + 1)`. Exact enumeration includes ties. With only five nonzero effects the smallest attainable two-sided p-value is 0.0625; the eligibility gate is not a guarantee that significance is attainable.

A 95% percentile interval comes from 5,000 bootstrap samples of whole cluster effects with replacement. Clusters receive equal weight in each sample. The bootstrap assumes those effects reasonably represent a population of independent lineages. It can be unstable with few lineages, and a narrow or degenerate interval is not proof that biological uncertainty is small. Intervals produced below the inference gate are descriptive.

The primary report includes leave-one-cluster-out effects and conditional tests when the reduced sample still passes the gate. These measure sensitivity to individual contributing lineages. They do not repair a failed independence assumption or establish causality.

## Multiple comparisons and robustness

The single prespecified primary test has its own correction family. Its BH-adjusted q-value equals its p-value when it is the only finite test in that family. All finite exploratory H2 sensitivity p-values are adjusted together using Benjamini-Hochberg. Missing or untestable p-values remain missing, rather than entering the correction as 1.

Sensitivity settings include identities 0.30/0.40/0.50, windows 3/5/10, architecture cutoffs 2/4/6/10, high/moderate versus high-only evidence, all eligible versus representative-only genomes, and Jaccard versus signed-order similarity. Closely linked and mixed-architecture comparisons are descriptive. Configurations reuse the same sequences and are correlated. BH interpretation retains its usual dependence assumptions; adjusted exploratory q-values do not independently validate the underlying null.

Every configured result is retained. A favorable threshold cannot replace the primary specification after seeing the data. A reversal across settings should be described as sensitivity, not hidden by selecting the setting that produces the desired conclusion. High-only analysis is especially restrictive because the present high-confidence rule requires an effectively exact reference match.

## Nulls, absence and generalization

A mean effect near zero means the chosen score did not distinguish target neighborhoods from the selected backgrounds in the contributing cohort. It does not establish that lysis biology, regulation or evolutionary constraint is absent. An untestable or low-information result must be distinguished from evidence of a negligible effect. No universal power claim or smallest detectable biological effect is implemented.

Missing cluster labels and singletons do not enter within-cluster H2 inference. Resolved architectures require detected supported A and B proteins, so inference is conditional on annotation detectability. The characterized reference set favors certain lineages and can induce nonrandom inclusion. Results cannot be generalized to unresolved phages merely by treating them as gene-negative.

The sketch screen is not ANI or a phylogenetic correction. Equal cluster weighting reduces domination by many sampled pairs but does not remove all ancestry confounding. Without gene-history, alignment and rearrangement evidence, the analysis cannot infer independent origins, selection, adaptive displacement or shared regulation.
