# Reproducing the frozen genome cohort

This document covers the archived 20-genome pilot. Reproduction of the expanded 50-genome screen is documented in [expanded results](exploratory_screen_results.md#reproduction).

The pilot cohort contains 20 accession-versioned genomes. `scripts/replay_cohort.py` reconstructs a metadata frame containing those genomes, executes all nine pipeline stages in a fresh workspace and compares the resulting tables and figures with the selected reference panel.

This procedure reproduces the downstream analysis of the frozen cohort. It does not repeat sampling from the original 3,760-entry eligibility frame or update historical cluster labels from a current database export. Source versions, reference sequences and checksums remain fixed.

## Execution

From a fresh checkout in WSL or Linux:

```sh
sh scripts/bootstrap_wsl.sh
. .venv/bin/activate
python scripts/replay_cohort.py --panel original --threads 4
python scripts/replay_cohort.py --panel experimental --threads 4
python scripts/replay_cohort.py --panel profile --threads 4
```

Each panel has a separate workspace beneath `.experiments/`. An existing workspace is preserved rather than overwritten. MMseqs is read from the local installation; `--mmseqs /path/to/mmseqs` selects another executable. Offline execution uses `--offline --cache-root /path/to/acquired-records`; every reused source file must match its frozen SHA-256 checksum.

The original panel has 10 functional references, the assay-expanded panel has 21, and the combined panel has 34. All use the same 20 genomes, homology cutoffs, architecture definitions and neighborhood settings. The three specificity challenges and the remote-LysB model evaluation are documented separately in [reference results](amendment_results.md) and [profile results](lysb_profile_results.md).

## Verification

Each panel generates 27 scientific tables and six figure pairs, for 39 compared artifacts. The replay report records the complete artifact set, expected and observed checksums, and successful completion of all nine stages. New workspaces rebuild protein-family clustering without a pre-existing family cache.

The optional profile libraries are required for the full profile test suite:

```sh
python -m pip install -r requirements-profile.txt
python -B -m pytest -q
```

Reproduction checks validate the software execution and its frozen inputs. They do not establish biochemical activity, independent evolutionary origins or a testable conservation comparison.

## Measured reproduction

All three nine-stage replays completed. All 117 scientific artifacts matched their published counterparts byte-for-byte. The complete test suite passed 214 tests with the optional profile dependencies installed.

| Panel | Scientific artifacts matched |
|---|---:|
| Original | 39/39 |
| Assay-expanded | 39/39 |
| Combined | 39/39 |

The replay used checksum-verified source records and pinned Python/MMseqs software. It rebuilt family clustering in each workspace. This check did not independently reinstall the software environment.

Machine-readable reports: [original](../results/reports/public_cohort_replay_original.json), [assay-expanded](../results/reports/public_cohort_replay_experimental.json), [combined](../results/reports/public_cohort_replay_profile.json), and [tests](../results/reports/tests.json).
