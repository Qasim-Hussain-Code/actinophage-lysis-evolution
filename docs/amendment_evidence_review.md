# Evidence review for the identification amendment

Reviewed 17 October 2026. This amendment addresses the limited diversity of the original experimental reference panel. The original pilot and its unresolved assignments remain a separate recorded analysis. Adding independently characterized sequences changes the reference panel, so any rerun needs its own configuration, input manifest and results. This review does not claim that identification or hypothesis testing has already improved.

The machine-readable reviews are [experimental_lysA_expansion.tsv](../references/experimental_lysA_expansion.tsv) and [experimental_lysB_expansion.tsv](../references/experimental_lysB_expansion.tsv). `panel_action=add` identifies usable additions; `existing`, `exclude` and `hold_sequence_mapping` must not be admitted automatically. Accession versions below were read from primary deposited records or explicit database cross-references, rather than inferred from adjacent identifiers. The LysA table deliberately includes tested constructs without positive results so that their exclusion is auditable.

## LysA evidence and accession mapping

The 2012 study tested twelve candidate LysA constructs. Positive mycobacterial ATP-release and culture-lysis results support D29 gp10, Kostya gp33 and L5 gp10. Positive peptidoglycan zymograms support Bxz1 gp236, Che8 gp32 and Corndog gp69 through the 2009 experiments and 2012 follow-up. These assays do not establish every predicted catalytic bond. [2012 study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034052), [2009 study, Figure 1A](https://pmc.ncbi.nlm.nih.gov/articles/PMC2774421/).

| Addition | Exact deposited mapping | Assay anchor |
|---|---|---|
| Kostya LysA | EU816591.1, gene 33, ACF34200.1, 477 aa | 2012 Figure 9 |
| L5 LysA | Z18946.1, gene 10, CAA79386.1, 292 aa | 2012 Figure 9 |
| Bxz1 LysA | AY129337.1, gene 236, AAN16868.1, 476 aa | 2009 Figure 1A |
| Che8 LysA | AY129330.2, gene 32, AAN12430.1, 424 aa | 2009 Figure 1A |
| Giles LysA | EU203571.2, gene 31, ABW88426.1, 402 aa | 2008 BRED Figure 3 |

Primary sequence records: [Kostya](https://www.ncbi.nlm.nih.gov/nuccore/EU816591.1), [L5](https://www.ncbi.nlm.nih.gov/nuccore/Z18946.1), [Bxz1](https://www.ncbi.nlm.nih.gov/nuccore/AY129337.1), [Che8](https://www.ncbi.nlm.nih.gov/nuccore/AY129330.2), [Giles](https://www.ncbi.nlm.nih.gov/nuccore/EU203571.2).

Giles provides a distinct genetic anchor. Its gene-31 deletion mutant propagated on a strain expressing Corndog lysA, required induction, and failed on the control strain. This supports a lysis role and functional substitution; it is not purified Giles enzymology. The 2012 reference to complementation as unpublished work is therefore not the evidence anchor used here. [Marinelli and colleagues, 2008, Figure 3](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003957), [Dedrick and colleagues, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3641587/).

Corndog is held out pending start-site reconciliation: the 2012 prose describes 439 residues, whereas its Table S1 and AAN02001.1 give 406. The figure-S3 label Che9 gp32 conflicts with Che8 in the main text and Table S1; it must not create another reference. Barnyard, Brujita, Bxz2, PLot and Rosebush had no positive lysis result in that screen. They are not validated non-lysins. TM4 and Fruitloop comparisons do not themselves supply new positive assays. [2012 supporting files](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034052#sec007).

Assay conditions matter. The 2023 Hammy study reproduced strong D29-mediated lysis, but Kostya caused a plate growth defect without conspicuous liquid-culture lysis in that system. Hammy growth inhibition and Waterfoul's lack of a phenotype do not establish their catalytic function. This supports retaining assay-specific descriptions rather than treating host growth inhibition as universal proof of endolysin activity. [Primary 2023 study, Figure 4](https://academic.oup.com/g3journal/article/13/12/jkad240/7338749).

## LysB evidence reviewed separately from predicted substrate specificity

The comparative structural paper belongs to the January 2020 issue and was published online on 27 December 2019. It modeled thirty proteins, but reports successful expression and enzymatic testing for seven: D29, Omega, Saal, Obama12, SWU1, BabyRay and Palestino. Its Table 3 provides positive para-nitrophenyl laurate hydrolysis for all seven. Palestino had no detected palmitate activity, and BabyRay's very small palmitate estimate has substantial uncertainty. The thirty modeled proteins are not thirty experimentally validated references. [Primary article, sections 3.3 and 3.4 and Table 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC7022511/), [publisher DOI](https://doi.org/10.3390/biom10010045).

The follow-up produced four soluble active proteins, D29, Omega, Saal and Obama12, from eight attempted constructs. It measures ester and Tween hydrolysis. Although its abstract discusses mycolylarabinogalactan-peptidoglycan activity collectively, the displayed substrate-specific LC-MS experiment and Figures 4 and 5 explicitly identify D29. This review does not extend that individual native-substrate demonstration to the other three proteins. Failed expression of Enkosi, Echild, DS6A or Pumpkin is not negative functional evidence. [Primary follow-up](https://pmc.ncbi.nlm.nih.gov/articles/PMC7246905/), [publisher DOI](https://doi.org/10.3390/ijms21093176).

Generic ester hydrolysis establishes biochemical activity under the reported assay conditions. It is less specific evidence for a native mycolylarabinogalactan esterase role than direct substrate cleavage or a native phage deletion and rescue experiment. Such additions can supply experimentally anchored homology references while retaining the assay limitation in every propagated call. The amended catalogue uses `experimental_generic_lipolytic` and caps this evidence at moderate, including an exact sequence match. This does not retroactively change the original pilot reference policy. Neither a generic serine-hydrolase match nor local gene placement independently proves LysB function.

| Addition | Paper Table S2 identifier | Deposited protein and genome | Length |
|---|---|---|---|
| Omega LysB | Q854L1 | AAN12697.1, AY129338.1, gene 53 | 290 aa |
| Saal LysB | W6AV51 | AHI61324.1, KJ025956.1, gene 31 | 333 aa |
| Obama12 LysB | W0LJG8 | AHG24142.1, KF959563.1, gene 10 | 326 aa |
| SWU1 LysB | I1V1E2 | YP_006382930.1, NC_017973.1 | 254 aa |
| BabyRay LysB | A0A1D8EW72 | AOT25464.1, KX683423.1, locus SEA_BABYRAY_10 | 244 aa |
| Palestino LysB | A0A1J0MAW0 | APD17532.1, KY083061.1, gene 10 | 326 aa |

The six identifiers and lengths come from the primary supporting Table S2. Current UniProt records explicitly connect [Omega Q854L1](https://www.uniprot.org/uniprotkb/Q854L1/entry), [Saal W6AV51](https://www.uniprot.org/uniprotkb/W6AV51/entry), [SWU1 I1V1E2](https://www.uniprot.org/uniprotkb/I1V1E2/entry) and [BabyRay A0A1D8EW72](https://www.uniprot.org/uniprotkb/A0A1D8EW72/entry) to the deposited proteins. The paper identifiers for Obama12 and Palestino are inactive in the current UniProt service, which records the deletion reason as not being part of a reference proteome. Their historical version-1 records explicitly map [Obama12 W0LJG8 to AHG24142.1](https://rest.uniprot.org/unisave/W0LJG8?format=txt&versions=1) and [Palestino A0A1J0MAW0 to APD17532.1](https://rest.uniprot.org/unisave/A0A1J0MAW0?format=txt&versions=1). Independent official sequence-archive records group the historical identifiers and corresponding EMBL protein sequence versions: [Obama12 UniParc UPI000232F0FA](https://www.uniprot.org/uniparc/UPI000232F0FA) has sequence MD5 `952FFF60CFD8E6EB8D5BBCC088F993CF`, and [Palestino UniParc UPI00078BB698](https://www.uniprot.org/uniparc/UPI00078BB698) has MD5 `4D7BA84E898BFD645DF4154F9503D14F`. Both are therefore eligible additions after archival mapping verification, rather than identification from matching names or lengths alone. The original [Obama12 GenBank record](https://www.ncbi.nlm.nih.gov/nuccore/KF959563.1) and [Palestino GenBank record](https://www.ncbi.nlm.nih.gov/nuccore/KY083061.1) provide the frozen genomic sources. Ordinary homolog calls from all six additions should retain the generic-substrate limitation.

## Conditions on the amended analysis

The numerical sequence thresholds should stay fixed during this reference expansion. Reference diversity, coverage and function-specific negative controls should be evaluated separately. Whole-protein thresholds can still miss modular LysA proteins whose component domains have different histories. Domain or profile evidence therefore needs a declared separate tier, with model provenance and independent validation. Published predictions must retain their prediction label.

Exact experimental sequences used to construct the panel are development controls. Detection of these same sequences is a pipeline check, not a held-out sensitivity estimate. Leave-one-reference-out or leave-one-family-out evaluation is preferable for assessing transfer, and distant generic hydrolases are useful specificity challenges. The absence of a broad-search hit remains unresolved detection rather than gene absence.

The amended panel must be selected by the documented assays and sequence mappings, not by whether it recovers a desired separated architecture. Yvonnetastic's deposited annotations and literature-predicted separation can check coordinate and parsing behavior; they must not train an independent functional call based on that same separation. Any change in the testable sample, architecture counts or inferential eligibility must be reported as a new result beside the original pilot.

## Access and inspection record

The PLOS article and small supporting PDFs were read from the primary publisher. GenBank records supplied exact protein versions, genes and translated lengths. The LysB full texts were inspected through the Europe PMC full-text service when live PMC or publisher pages were restricted. The Biomolecules supporting archive is about 27 MB; only its small accession-table PDF was read by HTTP byte ranges, about 540 KB per inspection, instead of retaining the entire archive. The table gives the UniProt identifiers used in the experimental study. No domain database or structure database was installed for this review.
