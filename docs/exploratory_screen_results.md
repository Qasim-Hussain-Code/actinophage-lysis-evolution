# Expanded exploratory screen

The fixed 50-record screen contains 25 Mycobacterium-host and 25 Gordonia-host phages across 25 reported cluster labels. All records pass the implemented QC checks, and all 5,915 annotated CDS translations are verified and searched. The original 20 accession versions and source checksums are preserved. The unchanged 39-reference classifier recovers five unique linked A/B pairs and no unique separated pair.

![Category detection and architecture](../results/exploratory50_v1/figures/exploratory_detection_architecture.png)

## Detection and architecture

| Host genus | Genomes | Verified CDS | A-family | B-family | Both categories | Unique linked | Complex | Unresolved |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Mycobacterium | 25 | 2,976 | 9 | 13 | 5 | 5 | 0 | 20 |
| Gordonia | 25 | 2,939 | 12 | 2 | 1 | 0 | 2 | 23 |
| Total | 50 | 5,915 | 21 | 15 | 6 | 5 | 2 | 43 |

A/B detection columns overlap. Architecture columns partition the cohort; unique separated architectures number zero in both genera. Multiple A candidates can produce a complex classification even when no B candidate is supported. These counts describe reference-supported assignments under the stated evidence caps, rather than biochemical validation or population prevalence. [Coverage by group](../results/exploratory50_v1/coverage_by_group.tsv), [evidence tiers](../results/exploratory50_v1/coverage_by_evidence_tier.tsv).

| Phage | Accession | Reported cluster | Intervening CDS | A/B overlap | Pair confidence |
|---|---|---|---:|---:|---|
| D29 | AF022214.2 | A | 1 | 0 bp | High |
| Dewey | OP867018.1 | Q | 0 | 1 bp | Moderate |
| Hail | MZ005672.1 | Q | 0 | 1 bp | Moderate |
| Typha | MK494099.1 | Y | 0 | 1 bp | Moderate |
| AuntGwenStacy | QB092740.1 | Y | 0 | 1 bp | Moderate |

![Original gene order around resolved pairs](../results/exploratory50_v1/linked_contexts/linked_architectures.png)

D29 recovers the linked positive control. The four additional uniquely resolved pairs occur in two reported clusters and remain homology-based predictions. John316 contains two A-family candidates and one B-family candidate; Sticker17 contains two A-family candidates and no supported B-family candidate. Their multiplicity is retained explicitly. [Architectures](../results/exploratory50_v1/architectures.tsv), [component pairs](../results/exploratory50_v1/architecture_pairs.tsv).

The primary partition of five linked, two complex and 43 unresolved genomes is identical at separation cutoffs of 2, 4, 6 and 10 intervening CDS. A 21-mer sketch screen groups the 50 genomes into 47 screening groups. These groups are sequence-similarity proxies, not ANI estimates or independently reconstructed lineages. [Threshold sensitivity](../results/exploratory50_v1/architecture_sensitivity.tsv), [near-identical pairs](../results/exploratory50_v1/near_identical_pairs.tsv), [genome QC](../results/exploratory50_v1/genome_qc.tsv).

### Complex separated candidate

John316 (MN945905.1, subcluster CS2) has two A-family candidates upstream of one B-family candidate. Both component pairs meet the separation criterion:

| A-family CDS | B-family CDS | Intervening CDS | Nucleotide gap |
|---|---|---:|---:|
| MN945905.1__cds0036 | MN945905.1__cds0051 | 14 | 5,962 bp |
| MN945905.1__cds0037 | MN945905.1__cds0051 | 13 | 5,129 bp |

The best reference matches are PokyPuppy_LysA33, PokyPuppy_LysA34 and PokyPuppy_LysB49, respectively, at 88.0%, 86.3% and 71.0% identity. All three references derive from [published profile predictions](https://doi.org/10.1371/journal.pone.0276603), giving moderate category assignments rather than independent biochemical support. Multiple A candidates retain the pre-specified complex classification; these separated components are not counted as unique separated architectures. John316 is the only B-positive genome in reported cluster CS in this cohort, so even a component-based exploratory eligibility rule would supply no within-cluster replicated B anchor. These observations identify a candidate for functional study without making the primary conservation comparison estimable.

## Evidence and interpretation

The screen retains all 237 reported reference alignments and a status row for every CDS. Seven query proteins are exact sequence aliases of panel references. Exact aliases and repeated within-cluster sequences cannot supply independent functional validation. Accession-base flags can miss biological aliases represented by different GenBank and RefSeq records. [All hits](../results/exploratory50_v1/all_reference_hits.tsv), [per-CDS status and aliases](../results/exploratory50_v1/protein_detection_status.tsv).

Identity, both coverage cutoffs, E-value and confidence caps are unchanged from the [recorded protocol](exploratory_screen_protocol.md). Computationally predicted references support at most moderate family assignments. Generic lipolytic or antimycobacterial assays do not independently establish native envelope-substrate cleavage or an endogenous role in phage lysis. The [GMA6 audit](gma6_arabinanase_audit.md) illustrates why related hydrolase sequences and native substrate function must remain distinct.

The previously measured reference-panel challenge is retained and authenticated against the same reference sequences, search settings and identification source code. It is not a new challenge experiment on this expanded cohort and does not estimate population specificity. D29's supported linked arrangement and Yvonnetastic's annotation-defined complex arrangement both reproduce.

Stages 1 to 5 are complete. Neighborhood clustering and the H2 test are not performed for this exploratory cohort. No unique supported separated pair is recovered, so this screen supplies no primary conservation-effect estimate. The 300-genome analysis has not been run. Descriptive detection and architecture are the supported outcomes; gene loss, recurrent rearrangement and regulatory mechanisms remain unresolved.

## Reproduction

The source accessions, original retrieval URLs and record checksums are in [selected_genomes.tsv](../results/exploratory50_v1/selected_genomes.tsv). The 39-reference catalogue and sequence manifest remain in the [frozen panel](../results/lysterase_update/panel/). Scientific tables are listed in the [artifact manifest](../results/exploratory50_v1/artifact_manifest.json).

An isolated cache-only replay repeats all five stages, including an uncached reference search. All **22 scientific artifacts** match byte for byte. Its [comparison report](../results/exploratory50_v1/replay_validation.json) distinguishes reused source records from a fresh acquisition and retains configuration, source-code and panel checksums.

The public accession-reproduction tool uses the published cohort directly, without the original sampling frame or preparation files:

```bash
.venv/bin/python -B scripts/reproduce_exploratory_screen.py --workspace public_reproduction
```

It retrieves exact versions of 50 genome records and 27 unique reference-source records, checks every source checksum, reruns acquisition through architecture classification, and compares public scientific tables and summaries with the published hashes. The [fresh-source run](../results/exploratory50_v1/public_reproduction.json) matches all 16 public scientific artifacts with the same pinned software environment. Software is reused for this check; independent software reinstallation is not claimed.

A separate [fresh public GitHub-clone replay](../results/exploratory50_v1/github_clone_reproduction.json) matches those 16 scientific artifacts and all eight expanded-figure files: two PNGs, two vector PDFs, two plotted source tables and two caption/checksum manifests. It uses the 77 checksum-verified source records as a cache and reuses the same pinned software. The report identifies the tested public commit. The software suite passes 313 tests in WSL; software correctness and deterministic reproduction do not establish biochemical function.

For a strictly cache-only reproduction, supply a directory containing the verified `genbank/` and `references/` source files:

```bash
.venv/bin/python -B scripts/reproduce_exploratory_screen.py --workspace cache_reproduction --source-cache SOURCE_CACHE
```

Missing records, a changed accession version, altered source bytes or mismatched scientific outputs stop the run. Failed sources are retained for inspection and never replaced with a different genome.

Both expanded figures have a 300 dpi PNG, vector PDF, plotted source TSV and caption/checksum manifest: [coverage and architecture](../results/exploratory50_v1/figures/figure_manifest.json) and [local gene order](../results/exploratory50_v1/linked_contexts/figure_manifest.json). Regenerate them from the public tables in new directories:

```bash
.venv/bin/python -B scripts/plot_exploratory_screen.py --output .experiments/figure_reproduction/coverage
.venv/bin/python -B scripts/plot_linked_contexts.py --output .experiments/figure_reproduction/contexts
```
