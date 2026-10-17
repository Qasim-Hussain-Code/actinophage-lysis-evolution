# Candidate-only Lysin B profile results

The frozen diagnostic completed in WSL on 2026-10-17. It searched all 2460 CDS translations from the unchanged 20-genome pilot with the prespecified nine-assay-reference model and the D29/Ms6 substrate-specific sensitivity model. Nineteen separately built models covered the full models and their source-held-out and reference-similarity-group-held-out folds. Both statistical search-space parameters remained 2460 throughout. Model construction, candidate thresholds and residue-pair coverage follow the [frozen protocol](lysb_profile_protocol.md); the [execution report](../results/amendment_v1/lysb_profile_candidates/profile_report.json) authenticates the saved artifacts and records resource usage.

## Cohort candidates

The broad model yielded five passing CDS, including the exact D29 training-sequence alias. Excluding that in-training observation leaves four candidate CDS, all from Mycobacterium phages in clusters E and M. The narrow model retained the two cluster M candidates plus the same D29 alias. The four non-training CDS represent three distinct amino-acid sequences: IPhane7 and Glaske16 encode identical candidate proteins. These counts are not independent discoveries or experimentally validated functions. All four proteins already have deposited Lysin B annotations and low-confidence B assignments in the frozen amended evidence. No assignment was upgraded. The [candidate ledger](../results/amendment_v1/lysb_profile_candidates/cohort_candidates.tsv), [domain measurements](../results/amendment_v1/lysb_profile_candidates/domain_hits.tsv) and [descriptive annotation review](../results/amendment_v1/lysb_profile_review/candidate_review.tsv) retain the identifiers and supporting measurements.

| Phage and cluster | GenBank protein; deposited gene | Passing models | Existing supported A geometry |
| --- | --- | --- | --- |
| Contagion, E | AGT13166.1; 34 | Broad | 2 intervening CDS; 714 bp edge gap |
| Balomoji, E | XEN17273.1; 34 | Broad | 2 intervening CDS; 714 bp edge gap |
| IPhane7, M | AXQ52679.1; 38 | Broad and narrow | No complete independently supported A detected |
| Glaske16, M | WNN95621.1; 39 | Broad and narrow | No complete independently supported A detected |

The Contagion and Balomoji A assignments have moderate confidence from homology to an experimental reference, rather than experimental characterization of these phage proteins. Both candidate-to-A pairs are on the same strand. The cluster M entries do not establish absence of an A gene. These are candidate-CDS coordinate comparisons, not newly classified Lysin A/B architectures; every eligible A component is retained in the review table.

For the broad model, paired model/target coverage was 88.43%/83.16% in both cluster E candidates and 91.04%/75.54% in the cluster M candidates. The E proteins also produced small E-values with the narrow model, but its paired model coverage was 66.28%, below the frozen 70% cut, while target coverage was 80%. Their broad-only status therefore depends on the model's coverage representation and must not be read as absence of narrower profile similarity. The cluster M proteins passed both models, with narrow-model paired coverage of 87.50%/93.19%. No threshold was changed after observing these results.

## Reference recovery and specificity challenges

| Model | Source-held-out recovery | Similarity-group-held-out recovery |
| --- | --- | --- |
| Broad assay panel | 9/9 distinct reference sequences | 7/9 sequences across four proxy groups |
| D29/Ms6 sensitivity | 2/2 reference sequences | 2/2 reference sequences |

All folds rebuilt their training-only MSA and HMM. The broad model failed to recover Obama12 and Palestino when their entire two-reference similarity group was excluded. Both still had small full-sequence and independent-domain E-values, but paired target coverage was 66.87%, below the frozen cut; model coverage was 81.34%. This loss is preserved in [held-out recovery](../results/amendment_v1/lysb_profile_candidates/held_out_recovery.tsv). Source exclusion alone consequently gives a more favorable impression of recoverability than removal of these related references. The four groups are reference-only similarity proxies, not independent clades. The narrow results involve only two references and do not establish broad functional sensitivity.

Neither full model produced a candidate-cut cross-family flag among the same three frozen PDB cellular-enzyme challenges or among the 22 non-B curated reference proteins. These small, selected and reused challenge sets do not estimate population specificity or a false-discovery rate. Reported/included flags remain available in the [full-sequence ledger](../results/amendment_v1/lysb_profile_candidates/full_sequence_hits.tsv) and domain ledger; no-hit query statuses remain in the query summary.

None of the three held-out Gordonia computational B predictions passed the complete candidate criteria. This is not a blanket failure to detect profile similarity: EpicDab LysB7 had a reported broad-model domain with independent E = 2.22e-9, but paired model/target coverage was only 44.03%/41.55%. PokyPuppy also had a reported broad-model hit that failed the cuts; Suerte had no reported hit. These outcomes limit transfer to divergent Gordonia candidates and do not establish gene absence or justify retrospectively reducing coverage requirements. The complete [query summary](../results/amendment_v1/lysb_profile_candidates/query_summary.tsv) preserves these negative and partial results.

## Interpretation and reproduction

This diagnostic adds inspectable candidate homology evidence in two Mycobacterium clusters. It does not resolve a displaced Gordonia lysis system, establish native mycolate-substrate specificity, demonstrate recurrent reorganization, or make the primary H2 comparison testable. The primary functional catalogue, architecture classifications, hypothesis tests and expansion gate were left unchanged.

The diagnostic uses one HMM search worker and four FAMSA alignment threads.

After preparing the unchanged isolated profile-panel inputs, install the pinned optional dependencies in the WSL environment:

```sh
.venv/bin/python -m pip install -r requirements-profile.txt
```

The original command authenticates and reuses the saved completed output without rescoring:

```sh
PYTHONPATH=src .venv/bin/python -B scripts/lysb_profile_candidates.py \
  --root .experiments/amendment_profile_v1 \
  --validation-dir results/amendment_v1/profile/reference_validation \
  --out results/amendment_v1/lysb_profile_candidates
```

The frozen input set includes execution receipts. A fresh pipeline replay creates new receipts and therefore correctly fails authentication against the original output, even when the scientific inputs match. To perform an actual rescan after replaying the isolated profile-panel inputs and reference validation, use a new ignored output directory:

```sh
PYTHONPATH=src .venv/bin/python -B scripts/lysb_profile_candidates.py \
  --root .experiments/amendment_profile_v1 \
  --validation-dir results/amendment_v1/profile/reference_validation \
  --out .tmp/lysb_profile_replay

PYTHONPATH=src .venv/bin/python -B scripts/review_lysb_candidates.py \
  --root .experiments/amendment_profile_v1 \
  --diagnostic .tmp/lysb_profile_replay \
  --out .tmp/lysb_profile_replay_review
```

Changed inputs, code, protocol or library versions also require a new output freeze. The separate review of the original saved diagnostic is reproduced with:

```sh
PYTHONPATH=src .venv/bin/python -B scripts/review_lysb_candidates.py \
  --root .experiments/amendment_profile_v1 \
  --diagnostic results/amendment_v1/lysb_profile_candidates \
  --out results/amendment_v1/lysb_profile_review
```

The review authenticates the completed HMM artifacts and the original cohort, CDS, genome-metadata and existing-evidence receipts before measuring coordinates. Its [provenance record](../results/amendment_v1/lysb_profile_review/review_provenance.json) links the results to those unchanged inputs. Eleven diagnostic tests and one descriptive-review test passed in WSL, including real-library checks of training isolation, exact-alias exclusion, single-reference folds, fixed-search-space E-value scaling, deterministic model export, residue-pair coverage and receipt authentication.
