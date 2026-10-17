# Fourteen-reference LysB profile follow-up

Design recorded on 18 October 2026 before inspecting the updated profile scores. The assay-reference follow-up supports B-family assignments in three cohort genomes but leaves the complete A/B architecture count at one. This secondary diagnostic evaluates whether the five added assay proteins change remote-homology recovery and candidate coverage.

## Fixed model design

The broad model uses the original nine assay-supported B references plus exactly MissWhite_LysB_2026, Bxb1_LysB_2026, Constella_LysB_2026, Bongo_LysB_2026 and phiGD57_1_LysB_2026. No cohort protein or Gordonia computational prediction enters training. All aliases are retained, and exact sequences are represented once. A separate two-reference D29/Ms6 sensitivity model retains the original substrate-anchored scope.

Alignment, HMM construction, search settings, candidate criteria and residue-pair coverage use the unchanged functions from `scripts/lysb_profile_candidates.py` and the [original protocol](lysb_profile_protocol.md). The full model and every fold rebuild their own training-only alignment and HMM. There is no target recruitment, alignment trimming or threshold tuning. Both Z and domZ remain 2460; full-sequence and independent-domain E-values must be at most 1e-6, and one domain must pair residues covering at least 70% of both the complete model and target.

## Validation and comparisons

Source-held-out folds exclude the resolved nucleotide-source base and exact-sequence aliases. Similarity-group-held-out folds use connected components defined only among the fourteen assay references by the original global-alignment thresholds. The same group mapping is restricted to the narrow model's reference subset; an empty training set is explicitly untestable. Similarity groups remain proxies rather than independent evolutionary clades.

Evaluation includes all 2460 cohort CDS, the three frozen cellular-enzyme challenges, non-B curated references and the three Gordonia computational B predictions. Raw sequence and domain scores, alignment strings, failures and no-hit statuses are retained. A full-model training alias is excluded from non-training candidate counts. In particular, both M1 CDS match the newly added Bongo training sequence and cannot be counted as new independent full-model candidates.

Report broad and narrow candidate counts, unique sequence counts, held-out recovery, cross-family flags and Gordonia-query outcomes alongside the original nine-reference diagnostic. No profile result changes the 39-reference primary functional calls, architecture classifications or expansion gate. The small selected validation set cannot estimate population sensitivity, specificity or a false-discovery rate.

## Integrity and reproduction

Input hashes bind the scientific CDS and reference files, query inputs, both protocols, the original scoring engine, the new driver and pinned library versions. Stage receipts authenticate those inputs but their execution-specific bytes are not model inputs. This permits an actual rescan against a fresh scientific replay when the scientific input bytes, code and software are identical. Changed scientific inputs or model artifacts require a separate output freeze. Public artifacts contain scientific measurements, model settings and provenance without local execution inventories.

HMM export clears creation time and command-line metadata without changing emission probabilities, transition probabilities, calibration parameters or scores. The numerical model settings and code hashes remain in the scientific freeze. These metadata properties follow the pinned [PyHMMER HMM interface](https://pyhmmer.readthedocs.io/en/stable/api/plan7/hmms.html#pyhmmer.plan7.HMM.command_line).
