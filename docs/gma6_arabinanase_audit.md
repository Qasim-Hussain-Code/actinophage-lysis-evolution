# GMA6 arabinanase and lysis-family audit

## Question

Can a biochemically characterized phage cell-wall hydrolase serve as a specificity challenge for the project’s LysB reference panel, and what does it reveal about sequence-based LysA assignments?

## Source mapping

The 2023 study expressed the GMA6 DUF4185 protein (PhageGH4185) and reported endo-D-arabinanase activity on mycobacterial arabinogalactan. The deposited GMA6 record maps this protein to CDS GMA6_45, protein AKL88326.1, in accession KR063280.1. The CDS is 1,087 aa at positions 38,130–41,393 on the plus strand. Its immediate downstream CDS is GMA6_46, annotated “putative holin,” at 41,390–41,899; the two deposited intervals overlap by four bases. GMA6_44 overlaps the upstream boundary by eight bases. These are annotation coordinates, not experimentally resolved translation boundaries.

The enzyme assay is direct evidence for arabinan hydrolysis by the recombinant protein. The paper does not test GMA6_45 in native phage infection, demonstrate host-cell lysis by this protein, or establish a LysB-type mycolate-esterase activity. Its placement beside a putative holin is genomic context, not proof of its role during infection. The paper's supplementary Fig. S9 proposes that DUF4185 may replace LysB in some phage lysis cassettes, either as a separate enzyme or fused to LysA. This is an architecture-level hypothesis, not an experimentally demonstrated infection mechanism; the available text does not resolve which specific arrangement applies to GMA6_45. [Primary study and supplement](https://www.nature.com/articles/s41467-023-37839-5); [deposited GMA6 genome](https://www.ncbi.nlm.nih.gov/nuccore/KR063280.1).

## Frozen-panel search

The translated GMA6_45 sequence was searched against the frozen 14-reference assay panel using the same MMseqs2 sensitivity and reporting settings used for the displaced-candidate audit. Category and assignment thresholds were read from the existing reference metadata and configuration. A LysB assignment requires at least 40% identity, 70% query and target coverage, and E-value at most 1e-6.

| Reference | Category and evidence tier | Identity | Alignment | Query coverage | Reference coverage | LysB threshold |
|---|---|---:|---:|---:|---:|---|
| Clawz LysA65 (QKY79977.1) | LysA-labeled; published profile prediction, no protein-specific assay | 53.2% | GMA6 1–1,032 aa; Clawz 1–856 aa | 94.9% | 92.7% | Not applicable: LysA category |

No LysB-class reference was returned by the search. The best match is therefore not a qualifying LysB hit. Clawz LysA65 is itself a computationally assigned multidomain protein in the cited Gordonia analysis; this alignment does not upgrade it to an experimentally established endolysin.

The close, nearly full-length match links an arabinanase-tested protein to a LysA-labeled sequence family. This is consistent with shared ancestry, a broader substrate range, or an annotation-family overlap; the present evidence cannot distinguish these possibilities. The two pilot proteins assigned to this Clawz reference are Jollymon CDS 64 and Stillion CDS 63. Their deposited annotations and sequence similarity remain homology evidence, not direct arabinanase assays. [Clawz source study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0276603).

In the fixed cohort, Jollymon XSB45410.1 and Stillion YGX17269.1 match Clawz LysA65 at 90.8% and 90.9% identity, respectively, over full query and reference coverage. Both genes are followed by deposited holin annotations. These two calls are from the same Gordonia cluster and do not constitute independent lineage replications. Their relationship to the tested GMA6 arabinanase makes substrate-specific follow-up worthwhile, while leaving their own activities unknown. The underlying rows are in the [assay-reference evidence table](../results/lysterase_update/search_evidence/all_cohort_reference_hits.tsv) and [lysis evidence table](../results/lysterase_update/panel/results/tables/lysis_evidence.tsv).

A direct search of all 2,460 translated CDS in the fixed cohort returned exactly these two near-full-length homologs to GMA6_45. Candidate reporting used MMseqs2 sensitivity 7.5, minimum sequence identity 25%, minimum query and target coverage 50%, and E-value at most 0.001. Jollymon XSB45410.1 and Stillion YGX17269.1 each show 52.8% identity over 1,032 aa, covering 94.9% of PhageGH4185 and 92.7% of each target. The existing pipeline calls both moderate LysA based on the profile-predicted Clawz reference; each deposited LysA is followed by an annotated holin. Both belong to cluster CP, so this does not add independent lineage support. Their strong similarity to the experimentally arabinanase-active GMA6 protein makes them specific candidates for biochemical substrate testing, but does not establish that they hydrolyze arabinogalactan or contribute to phage lysis.

## Interpretation

This is an orthogonal biochemical specificity challenge, not a LysB positive control and not a validated negative control for all lysis activity. The published replacement proposal makes DUF4185 a plausible third cell-wall enzyme class to consider when a LysB assignment is missing, while the cross-category match shows why substrate-specific evidence and broad “lysin A” profile labels must remain separate. A phage enzyme can hydrolyze a cell-wall arabinan substrate in vitro while its contribution to native lysis, any additional activity, and its relationship to canonical LysA function remain unresolved.

The result does not add a functionally supported LysA/B pair, recover an independently validated displaced pair, or change the eligibility of the primary neighborhood-conservation test. It does add a characterized phage hydrolase and a cross-category sequence result to the specificity audit.

## Reproduction files

- [Coordinate and evidence table](../results/gma6_arabinanase/gma6_arabinanase.tsv)
- [All returned reference-panel hits](../results/gma6_arabinanase/reference_panel_hits.tsv)
- [LysB-class hits](../results/gma6_arabinanase/lysB_panel_hits.tsv)
- [Machine-readable report and checksums](../results/gma6_arabinanase/audit_report.json)
- [Independent raw-search replay validation](../results/gma6_arabinanase/replay_validation.json)
- [Audit script](../scripts/audit_gma6_arabinanase.py)
- [Fixed-cohort homolog table](../results/gma6_arabinanase/cohort_homologs.tsv)
- [Fixed-cohort search report](../results/gma6_arabinanase/cohort_homolog_report.json)
- [Independent cohort-search replay](../results/gma6_arabinanase/cohort_replay_validation.json)
- [Cohort-search script](../scripts/audit_gma6_cohort_homologs.py)

Recreate the source from accession KR063280.1 and the frozen reference panel. In Linux or WSL, run `python scripts/audit_gma6_arabinanase.py --download --extract-only`, then the MMseqs2 command below from the repository root. Finalize with the audit script using `--precomputed-hits .experiments/gma6/mmseqs_work/mmseqs_hits.tsv --replay-hits .experiments/gma6/replay_hits.tsv --mmseqs-version 8cc5ce367b5638c4306c2d7cfc652dd099a4643f`. Produce the independent replay file by repeating the same MMseqs2 command with a new output and temporary directory.

For the pilot homolog scan, first run `python scripts/audit_gma6_cohort_homologs.py --extract-only` to verify the pinned manifest and source records and construct a FASTA containing exactly 2,460 CDS translations. Run the equivalent MMseqs2 settings with this query and cohort FASTA, then finalize with `--precomputed-hits` and `--replay-hits`. The scan does not acquire new genomes or publish the cohort sequences.

```sh
.tools/mmseqs/bin/mmseqs easy-search \
  .experiments/gma6/mmseqs_work/PhageGH4185.faa \
  results/lysterase_update/reference_validation/queries.faa \
  .experiments/gma6/mmseqs_work/mmseqs_hits.tsv \
  .experiments/gma6/mmseqs_work/mmseqs_tmp \
  --threads 4 -s 7.5 -k 5 --alignment-mode 3 --max-seqs 100 \
  --min-seq-id 0.25 -c 0.5 --cov-mode 0 --split-memory-limit 1G -e 0.001 \
  --format-output query,target,fident,alnlen,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,evalue,bits
```
