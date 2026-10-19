# Exploratory 50-genome detection screen

Plan recorded on 19 October 2026 before acquisition of any additional genome for this screen. The screen asks whether the existing evidence strategy detects supported Lysin A/B architectures beyond the fixed 20-genome pilot. Fifty is the planned number of selected records, not an observed number of acquired or QC-passing genomes. No new detection counts or biological results are asserted here.

## Selection and acquisition

Use the checksum-verified frozen eligibility frame, seed `20261016`, Mycobacterium and Gordonia host genera, D29 and Yvonnetastic controls, and the existing `acquisition.select_sample` algorithm with a target of 50. This orders candidates independently of lysis outcomes, alternates genera, prioritizes groups with paired representatives and deduplicates base accessions. Missing cluster labels remain separate rather than defining one shared lineage.

Before retrieval, assert that the ordered first 20 selected base accessions reproduce the frozen pilot in order. Preserve the exact accession versions and source checksums of those 20 records; retrieve only additional selected records. A nesting mismatch stops the screen for investigation. Freeze the selection ledger, sampling settings, eligibility-frame checksum and reference-panel checksum before searching new proteins. Do not replace unavailable records, QC failures or near duplicates with more favorable genomes. Record every selected record and its acquisition or exclusion outcome.

Retrieve accessioned GenBank records with bounded streaming, at most three attempts and the existing per-record size limit. Verify accession identity, freeze each successful returned version and source checksum immediately, and reject changed frozen records. Availability failure after bounded retries leaves an incomplete screen and stops acquisition; a complete result requires all 50 selected source records to be acquired. QC failures are subsequently excluded without replacement. Keep raw sources, bulk eligibility metadata and intermediate caches in an isolated, untracked workspace. Reuse only checksum-verified cached inputs.

## Unchanged QC and evidence

Apply the current whole-genome and CDS checks: 10,000–500,000 bases, permitted IUPAC symbols, at most 1% ambiguous bases, explicit phage and whole-genome completeness evidence, CDS annotations, valid coordinates and verified translations. Retain current handling of partial features and unresolved translations. There is no additional 30-kb length filter. Retain current duplicate-screen outcomes without replacing excluded or grouped records.

Use the frozen 39-reference catalogue at `results/lysterase_update/panel/references/protein_controls.tsv` and its matching reference manifest. Preserve the existing evidence-tier caps and best-reference rules. Moderate assignments require at least 40% amino-acid identity, at least 70% coverage of both proteins and E-value at most `1e-6`; high confidence remains restricted to eligible full-length exact experimental-reference matches. Published profile predictions and the newer antimycobacterial assay tier retain their moderate-confidence ceilings. Generic lipolytic evidence retains its existing substrate limitations.

Neither deposited labels nor proximity, gene order or conservation can upgrade function. Keep annotation-only, weak, partial and independently supported evidence distinct. Preserve complex architectures and all component pairs. Use the existing separation thresholds 2, 4, 6 and 10 intervening CDS, with greater than four as the inherited main configuration. This configuration does not make the new screen confirmatory. Any GMA6-like arabinanase homology is a separate exploratory candidate measurement and cannot automatically confer a LysA or LysB assignment.

## Execution, controls and outputs

Run acquisition, QC, reference identification, evidence extraction and architecture stages only, in isolation from the locked pilot. Reuse pinned software, at most four threads, the existing MMseqs search parameters and `1G` split-memory setting. Apply the project-wide effective storage cap, the smaller of 10 GB or initial free space minus a 4 GB reserve, including the isolated workspace. A source-integrity failure, contradictory validated control or resource-guard failure stops execution; retain the failure record for diagnosis.

D29 must retain its exact supported linked arrangement. Yvonnetastic remains an annotation-level displaced comparison unless independent evidence changes in a separately documented amendment. Retain the existing structural and cellular-enzyme specificity checks under unchanged settings. Report reference-source matches and exact training-sequence aliases separately; they do not estimate independent recovery performance. The selected controls cannot establish population sensitivity, specificity or a false-discovery rate.

Minimum scientific outputs are the frozen protocol and configuration hashes; selection and nesting ledger; accession-version/source manifest; acquisition and QC exclusions; genome, CDS and duplicate-screen counts; every reported reference alignment; per-CDS detection status including no reported hit; functional-evidence tables; architecture and component-pair tables; and summaries by host genus, cluster and evidence tier. Identify pilot versus additional records and exact-sequence aliases. Report candidate counts with their denominators, supporting references and substrate uncertainties. Record code, input and output checksums and an isolated replay comparison before publishing a completed result.

## Scope and progression

Complete the prespecified selected set without stopping early for favorable findings. The screen estimates detection yield conditional on this sampling frame and reference strategy. Non-detection does not establish gene absence. Paired clusters containing supported separated architectures may be listed as candidates for later feasibility review, but they do not establish matched-background availability or H2 eligibility. No family clustering, neighborhood conservation inference, population prevalence estimate, phylogenetic reconstruction or evolutionary-mechanism claim is part of this screen.

The existing main-pipeline expansion gate remains unchanged. Its `expansion_ready` flag requires at least one independently supported separated candidate, passing source checksums, no independently contradicted architecture control, a passing linked control, and annotation-level recovery of an expected displaced control, including recovery of all separated component pairs in a complex architecture. It does not require the displaced control itself to pass independent functional validation. This exploratory screen neither sets that flag nor runs the normal intermediate or full cohort mode. Completion does not authorize automatic progression to H2 testing or 300 genomes. Any later analysis requires a separately recorded design and the relevant functional, matching, lineage and assumption checks.

## Data availability

Publish selected accession versions, phage and host-genus labels, cluster classifications, lengths, QC outcomes, source URLs, checksums, evidence tiers, measurements and methods. Source records and bulk metadata exports are obtained from their originating providers under the applicable source terms.

An accession-based replay reproduces the frozen selected cohort. It does not independently reconstruct selection from the unpublished full eligibility frame; report that distinction explicitly. The original pilot, its outputs, its H2 status and earlier analysis freezes remain unchanged.
