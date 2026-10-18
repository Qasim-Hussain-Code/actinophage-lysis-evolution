# Genomic architecture and evolutionary conservation of actinobacteriophage lysis machinery

This repository asks whether physically separated lysis genes retain conserved genomic neighborhoods after accounting for shared ancestry, background gene order and uncertain functional annotation. The fixed 20-genome cohort contains 2,460 verified CDS translations. Reference expansion improves supported Lysin A coverage from 1/20 to 10/20 genomes. The latest [assay-reference follow-up](docs/lysterase_update_results.md) increases supported Lysin B coverage from 1/20 to 3/20 with a 39-reference panel. D29 remains the only complete supported A/B pair, so the primary conservation test remains untestable. The original pilot and earlier reference panels are preserved separately.

The original remote-profile diagnostic identifies four candidate CDS representing three sequences. The later assay update independently supports the sequence shared by its two M1 candidates, which exactly matches assayed Bongo LysB. The two E candidates retain low primary confidence. The [original diagnostic](docs/lysb_profile_results.md) and [updated assignment ledger](results/lysterase_update/comparison/confidence_changes.tsv) distinguish profile candidates from assay-reference homology.

An [expanded fourteen-reference profile analysis](docs/lysb_profile_update_results.md) recovers 13/14 held-out assay-reference queries. Its remaining failure has strong domain similarity but insufficient full-target coverage. No additional Gordonia query meets every candidate criterion, and primary functional calls remain unchanged.

A separate [sensitivity to deposited annotations](docs/annotation_neighborhood_results.md) yields four descriptive matched comparisons. The two pairs with separated annotations have base Jaccard contrasts of 0.000 and −0.413 relative to background, with window-dependent variation. These conditional measurements retain uncertain functional labels and do not supply a primary H2 estimate.

A [survey of assay-reference source genomes](docs/reference_architecture_results.md) adds 13 complete source records with 1,649 verified CDS. Five sources have supported closely linked A/B pairs; eight remain unresolved because A identification is incomplete. This control-selected collection supplies no supported displaced arrangement and does not change the representative cohort or expansion gate.

A targeted [assayed-allele audit](docs/assayed_reference_allele_results.md) maps four additional published source genomes to their LysA/B genes and assay scope. Three full articles report assays for both enzymes; a conference abstract reports both for RitSun, but does not identify its tested LysA allele. All four source pairs are adjacent or overlapping. A one-nucleotide PDRPxv CDS/report discrepancy and a synonymous BTCU-1 primer change are retained explicitly. These sources add linked controls without resolving the separated-gene comparison; pilot classifications remain fixed.

A separate [displaced-candidate audit](docs/displaced_lysB_candidate_audit.md) maps three literature-proposed LysB-like genes in *Rhodococcus* phages. Two occupy virion structural regions and one is ten annotated CDS from LysA, but none has protein-specific functional validation in the cited study. All three have sub-threshold matches to the assay-reference LysB panel and none meets the project assignment cutoff. Whole-phage activity is kept distinct from attribution to these candidate genes.

The [GMA6 arabinanase audit](docs/gma6_arabinanase_audit.md) adds a biochemically characterized phage cell-wall hydrolase as an orthogonal specificity challenge. Its strongest panel hit is a LysA-labeled profile-prediction sequence; no LysB-class hit is returned. This cross-category match motivates a family/function caveat but supplies no native lysis result or new eligible LysA/B pair.

## Research question and prior work

Lysin A degrades peptidoglycan, whereas characterized Lysin B proteins can cleave mycolate-containing envelope linkages. Holins participate in membrane permeabilization, but a membrane protein annotation alone cannot establish holin function. These functional categories are heterogeneous and cannot be identified reliably from proximity alone.

Displaced Lysin B genes and diverse lysis-cassette arrangements have been reported in [Yvonnetastic](https://doi.org/10.1128/genomeA.00594-16), [Gordonia comparative genomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC5559632/) and the [2022 Gordonia lysis analysis](https://doi.org/10.1371/journal.pone.0276603). This analysis evaluates whether conservation around separated Lysin B exceeds conservation around matched homologous non-lysis anchors in the same genome pairs. The [literature review](docs/literature_and_novelty.md) distinguishes experimental findings from computational predictions.

| Hypothesis | Planned comparison and null | Current status |
|---|---|---|
| H1: proximity differs among evolutionary groups | Describe architecture within clusters and compare lineage-balanced frequencies; formal group effects require usable, replicated functional measurements | Descriptive only |
| H2: separated Lysin B retains excess neighborhood conservation | Equal-weight cluster mean of Lysin B minus matched homologous background conservation; null mean contrast zero under independent, sign-symmetric cluster effects | Untestable: no contributing clusters |
| H3: similar displacement reflects recurrent reorganization | Require sufficiently resolved gene and genome histories to distinguish inheritance, recombination and independent changes | Untestable from cluster names alone |
| H4: some apparent incompleteness reflects annotation limitations | Compare deposited labels, independently supported assignments and confidence sensitivities; non-detection does not establish gene loss | Descriptive annotation audit |

## Pilot data and observations

The frozen cohort contains **20 complete genomes, 11 named clusters and 2,460 verified CDS translations**. Eleven genomes infect Mycobacterium and nine infect Gordonia. All selected records passed the implemented QC checks; none was excluded. A canonical 21-mer sketch screen groups one near-identical pair, leaving 19 screening groups. These groups are similarity proxies, not a phylogeny or ANI estimate.

Selection is deterministic and independent of observed architecture. Published controls are included first, followed by two accessioned representatives per eligible cluster in a seeded, genus-interleaved order. The control allocation and cluster sampling explain the 11/9 host split; numerical balance was not imposed afterward. Current metadata, accession versions, source URLs, retrieval times and SHA-256 values are retained in the [input manifest](data/manifests/pilot_inputs.tsv). [Dataset QC](results/reports/quality_control.json) and the [result manifest](results/tables/dataset_manifest.tsv) provide the measured cohort details.

The original 10-reference panel yields three high-confidence functional assignments, all in D29: Lysin A, holin and Lysin B. It yields no moderate-confidence assignments. Deposited product labels and feature notes identify 51 candidate proteins; weaker homology and deposited labels remain low confidence. The following pilot tables describe that preserved baseline; the expanded 21- and 34-reference results are reported separately above and in the amendment results.

At the operational separation cutoff of more than four intervening CDS:

| Evidence scope | Closely linked | Separated | Complex | Unresolved |
|---|---:|---:|---:|---:|
| Independently supported high/moderate assignments | 1 | 0 | 0 | 19 |
| Deposited annotation only | 9 | 4 | 3 | 4 |

These are conditional classification counts, not estimates of biological gene absence. Annotation-only separation in Phrappuccino, Settecandela, VanLee and Colossa does not establish independently supported lysis function. See [functional evidence](results/tables/lysis_evidence.tsv), [architectures](results/tables/architectures.tsv), [all component pairs](results/tables/architecture_pairs.tsv) and [confidence/threshold sensitivities](results/tables/architecture_sensitivity.tsv).

D29 has one intervening CDS, its overlapping holin, and a 418 bp edge gap between Lysin A and Lysin B. Yvonnetastic illustrates a different uncertainty: its record contains two putative Lysin A components and a serine hydrolase whose feature note says Lysin B. Both deposited component pairs are separated, with 112 or 113 intervening CDS and gaps of 49,969 or 50,875 bp. The annotation-level displaced arrangement is recovered while its unique independently supported architecture remains unresolved. Targeted local alignments are recorded as exploratory diagnostics and never upgrade these functions.

Study-specific clustering produces 1,460 protein families at the primary 40% identity threshold, with 1,411 and 1,496 at the alternative 30% and 50% thresholds. However, none of the nine within-cluster genome pairs has the supported anchors needed for the primary matched comparison. Consequently, the primary effect size, uncertainty interval, p-value and q-value are **not estimable**. This is insufficient evidence, not a negative conservation result. No candidate recurrent rearrangement, independent origin, horizontal transfer or gene loss is established. [Detailed results](docs/results.md) retain these limitations in the interpretation.

## Reproducing the analysis

The workflow uses Python 3.14.4 and MMseqs2 release 18-8cc5c.

From the repository root inside WSL or compatible Linux, reproduce the frozen cohort:

```bash
sh scripts/bootstrap_wsl.sh
. .venv/bin/activate
python scripts/replay_cohort.py --panel original --threads 4
python -B -m pytest -q
```

The bootstrap checks available storage, creates a local environment, installs pinned dependencies without a package cache or bytecode compilation, and obtains CPU-only MMseqs2 release 18-8cc5c. It uses `requirements-lock.txt` when present. The AVX2 binary requires an AVX2-capable CPU; a compatible official CPU build can be selected by editing the explicit tool path. The downloaded archive and installed tool are local, untracked software.

The three reference panels can be reproduced separately:

```bash
python scripts/replay_cohort.py --panel original --threads 4
python scripts/replay_cohort.py --panel experimental --threads 4
python scripts/replay_cohort.py --panel profile --threads 4
```

Each command creates a fresh workspace under `.experiments/`, retrieves version-pinned source records and compares the regenerated scientific artifacts with the corresponding published panel. Existing replay workspaces are preserved; use a fresh checkout for another complete replay. Offline execution accepts `--offline --cache-root /path/to/acquired-records`. The reconstructed metadata frame contains only the frozen study cohort; it does not repeat selection from the original eligibility population. Changed source annotations fail checksum verification, even when the sequence version is unchanged. See [cohort reproduction](docs/cohort_reproduction.md).

The coordinating runner verifies upstream hashes, scientific configuration and source fingerprints before resumption. Earlier-stage reruns invalidate downstream receipts. For example, resume a generated original-panel workspace with `python .experiments/public_replay_original/scripts/run_all.py --offline --resume-from 6`. Missing, changed or invalid critical inputs stop execution.

Intermediate and full modes have configured limits of 50 and 300 genomes:

```bash
python scripts/run_all.py --config config/analysis.yaml --mode intermediate --threads 4
python scripts/run_all.py --config config/analysis.yaml --mode full --threads 4
```

These commands currently stop at the scientific expansion gate before downloading a larger cohort. They require a completed pilot, recovered controls and independently supported separated candidates. They have not been executed as completed biological analyses. More genomes alone would not remedy the present reference-coverage failure.

The [analysis versions](docs/milestones.md) summarize the pilot and exploratory follow-up analyses.

## Pipeline and statistical safeguards

| Stage | Entry point | Output |
|---|---|---|
| 1 | `scripts/01_fetch_data.py` | Bounded acquisition and frozen manifest |
| 2 | `scripts/02_quality_control.py` | Record, coordinate, translation and duplicate checks |
| 3 | `scripts/03_identify_lysis.py` | Curated-reference MMseqs2 alignments |
| 4 | `scripts/04_annotation_evidence.py` | Independent evidence and deposited-label catalogue |
| 5 | `scripts/05_architecture.py` | CDS-order separation and all eligible component pairs |
| 6 | `scripts/06_neighborhoods.py` | Study-specific families and matched comparisons |
| 7 | `scripts/07_evolutionary_tests.py` | Cluster effects, testability and sensitivity ledger |
| 8 | `scripts/08_validate_candidates.py` | Input hashes, controls, candidates and targeted diagnostics |
| 9 | `scripts/09_generate_figures.py` | Figures and supporting data/checksum manifest |

Functional support requires a characterized reference, at least 40% identity, at least 70% coverage of both proteins and E-value at most 1e-6 for moderate confidence. High confidence requires a full-length exact match to the characterized sequence. Context never upgrades an assignment. Primary assignments use pairwise homology. A separate exploratory profile diagnostic evaluates remote LysB candidates; independent gene calling was not performed.

Family clustering uses MMseqs2 single-step set cover, explicit alignment identity, 80% bidirectional coverage and E-value at most 0.001. The five-residue prefilter seed and 1 GB prefilter allowance fit this small analysis; the prefilter remains heuristic. Representative-member edges do not imply identical function or all-pairs identity within a family.

Neighborhood windows contain three, five or ten CDS on either side of Lysin B. The primary metric is flanking-family Jaccard similarity; signed order/orientation is secondary. Controls are single-copy homologous non-lysis anchors in the same genome pair, matched without using the conservation score for boundary completeness, linear position, protein length and family prevalence. At least five matched anchors are required. Genome pairs are averaged within clusters before cluster bootstrap or two-sided cluster sign-flip inference. At least five contributing clusters are required, but that threshold cannot establish their phylogenetic independence. Exact enumeration and simulation tests check the statistical implementation; the unavailable biological comparisons do not validate exchangeability of lysis and background anchors.

Only one H2 configuration is prespecified as primary. Benjamini-Hochberg correction applies separately to exploratory sensitivities. High-only assignments, representative genomes, separation cutoffs, family thresholds, neighborhood widths and leave-one-cluster-out analyses are implemented. On this pilot, none supplies an inferential sample. Full definitions and assumptions are in [methods](docs/methods.md) and [statistical assumptions](docs/statistical_assumptions.md).

## Figures

The six figures summarize cohort composition, architecture, genome maps, lineage distributions, comparison eligibility and sensitivity. Captions and source tables are listed in the [figure manifest](results/reports/figure_manifest.json).

| Figure | Content | PNG | PDF |
|---|---|---|---|
| 1 | Cohort and lineage representation | [PNG](figures/figure_01_dataset.png) | [PDF](figures/figure_01_dataset.pdf) |
| 2 | Architecture and separation by evidence scope | [PNG](figures/figure_02_architecture.png) | [PDF](figures/figure_02_architecture.pdf) |
| 3 | D29 and Yvonnetastic maps in original coordinates | [PNG](figures/figure_03_genome_maps.png) | [PDF](figures/figure_03_genome_maps.pdf) |
| 4 | Raw and equal-cluster architecture proportions | [PNG](figures/figure_04_lineages.png) | [PDF](figures/figure_04_lineages.pdf) |
| 5 | Matched-comparison eligibility | [PNG](figures/figure_05_neighborhoods.png) | [PDF](figures/figure_05_neighborhoods.pdf) |
| 6 | Confidence and cutoff sensitivity | [PNG](figures/figure_06_robustness.png) | [PDF](figures/figure_06_robustness.pdf) |

## Citation and references

Citation metadata are available in [CITATION.cff](CITATION.cff). The original software is distributed under the [MIT license](LICENSE).

Primary datasets, experimental evidence and methods are recorded in the [literature ledger](references/literature_audit.tsv), [protein catalogue](references/protein_controls.tsv) and [data sources](docs/data_provenance.md).
