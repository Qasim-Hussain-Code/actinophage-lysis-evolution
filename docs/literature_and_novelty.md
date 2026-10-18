# Literature audit and scope of contribution

Audit date: 17 October 2026. This is a targeted primary-source audit, not a systematic review or proof that no other relevant publication exists. The machine-readable evidence ledger is [literature_audit.tsv](../references/literature_audit.tsv). No claim of novelty follows from this audit.

## What is already established

Variable lysis-gene placement is established prior art. The 2017 study by Pope and colleagues compared 79 Gordonia phages and described frequent separation of predicted LysB from other lysis genes, as well as lysis functions embedded among minor tail genes and split endolysin domains. These are comparative genomic observations and functional predictions. They do not establish the biochemical roles of all annotated proteins or the direction of a rearrangement. [Primary article](https://journals.asm.org/doi/10.1128/mbio.01069-17).

Yvonnetastic is an explicit published displaced-architecture control: its 2016 genome announcement reports a putative LysB over 50 kb from the other predicted lysis functions in the 98,136-bp KU963248 genome. Reproducing that separation is a validation check, not a new biological finding. [Primary announcement](https://journals.asm.org/doi/full/10.1128/genomea.00594-16).

Pollenz, Bland and Pope (2022) studied 16 representative cassette organizations spanning 77 G. rubripertincta phages. Their analysis already compares LysA flanking phams, identifies distal LysB candidates and examines membrane proteins. The proteins called holin-like were inferred using topology, homology and context; this is not experimental holin validation. Their work overlaps substantially with a simple survey of cassette diversity. [Primary article](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0276603).

Additional prior art extends beyond those starting references. DoobyDoo's 2025 announcement identifies predicted LysA gene 35 and LysB gene 86 in separate genome regions (PQ362665). NoShow's 2023 announcement describes a cluster AB genome whose lysis proteins resemble cluster H2 proteins. Such cross-cluster sharing cautions against equating whole-genome cluster membership with the history of a particular locus. [DoobyDoo](https://www.micropublication.org/journals/biology/micropub-biology-001471/), [NoShow](https://journals.asm.org/doi/10.1128/mra.00649-23).

LysA modularity also predates this project. Payne and Hatfull's 2012 study presents 25 domain organizations and selected functional assays. Some domains have different inferred histories, so a single full-length LysA tree may be misleading. Domain-specific comparisons and explicit alignment diagnostics are preferable when domain architectures differ. [Primary study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0034052).

## Experimental anchors and the limits of transfer

A small reference set should distinguish direct experimental support for the exact protein from inferred function in a homolog. A characterized reference supports a family assignment only to the extent supported by alignment coverage, specificity and conserved functional features. It does not make every matching protein experimentally characterized.

| Reference | Literature-supported identifiers and record-resolved versions | Direct evidence and appropriate use |
|---|---|---|
| D29 LysB | gp12; NP_046827.1; genomic accession AF022214; structure 3HC7 | Purified enzyme and structural characterization. Function-specific LysB reference. The Giles deletion experiment in the same paper must be attributed to Giles, not D29. |
| D29 LysA | gp10; NP_046825.1 in NC_001900.2 | Lytic expression assays in Payne and Hatfull (2012); version resolved in the frozen reference manifest. |
| D29 holin | gp11; NP_046826.1 in NC_001900.2 | Functional expression, aggregation and mutagenesis studies; version resolved in the frozen reference manifest. |
| Ms6 accessory chaperone | gp1; AAG48317; lysis locus AF319619 | Protein-interaction and endolysin-delivery evidence; do not use gene proximity alone to propagate this label. |
| Ms6 LysA | gp2; AAG48318; lysis locus AF319619 | Genetics and expression establish two translation products from the same reading frame. They are not two independent CDS. |
| Ms6 LysB | gp3; AAG48319; lysis locus AF319619 | Lipolytic assays, specific envelope-substrate assays and a native deletion phenotype. |
| Ms6 holin-like proteins | gp4: AAG48320.1; gp5: AAG48321.1; locus AF319619.1 | Native deletion/timing and interaction evidence; exact versions resolved from the frozen locus. |
| WXIN LysB | WPH57662; genomic accession OR813930 | Purified-protein assays and mutagenesis in 2026; useful divergent biochemical reference. Generic ester-substrate results should not be represented as a direct native-lifecycle knockout or a complete substrate-specificity determination. |
| Yvonnetastic / DoobyDoo | KU963248 / PQ362665 | Published predicted architecture controls; do not put their predicted LysB proteins into the experimentally characterized reference tier. |

Gordonia phage GMA6 provides a separate substrate-specific enzyme anchor: recombinant PhageGH4185 (GMA6_45; AKL88326.1) hydrolyzes mycobacterial arabinogalactan in vitro. Its sequence aligns across nearly its full length to Clawz LysA65, a profile-predicted reference in the project panel, but yields no LysB-class hit. Supplementary Fig. S9 in the enzyme study proposes DUF4185 as a LysB replacement in some phage cassettes, either as a separate protein or fused with LysA. This supports an alternative architecture-level hypothesis, while leaving native lysis and the specific GMA6 cassette mechanism untested. Clawz LysA65 also remains uncharacterized. See the [GMA6 audit](gma6_arabinanase_audit.md). [Primary enzyme study and supplement](https://www.nature.com/articles/s41467-023-37839-5); [Clawz source study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0276603).

The D29 LysB protein accession is supported by the [NCBI protein record](https://www.ncbi.nlm.nih.gov/protein/NP_046827.1) and independent reuse in the [WXIN experimental methods](https://link.springer.com/article/10.1186/s12866-025-04683-0). Sequence-version suffixes and coordinates must come from downloaded records, never from incrementing neighboring accession numbers. Literature gene numbering can differ after reannotation.

The core experimental evidence comprises [Payne et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2774421/), [Gil et al. (2008)](https://pubmed.ncbi.nlm.nih.gov/18451045/), [Gil et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/20093291/), [Catalão et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/20545844/), [Ms6 holin experiments (2011)](https://journals.asm.org/doi/full/10.1128/jb.01519-10), [Ms6 nested translation study (2011)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0020515), and [Gp1 binding-domain experiments (2011)](https://journals.asm.org/doi/10.1128/jb.00380-11).

Ms6 is especially informative because its LysB-deletion derivative can remain viable while releasing progeny inefficiently. The word “efficient” in the 2017 paper's title is essential to its meaning; do not recast this as absolute essentiality. The historical AF319619 entry is a lysis-locus reference, not evidence that a complete Ms6 genome is eligible for the genome cohort. [Deletion and cryo-electron tomography study](https://www.mdpi.com/1999-4915/9/11/343).

For classifier validation, reference sequences reused as search targets are development controls, not independent held-out sensitivity estimates. Use held-out characterized proteins or leave-one-reference-family-out evaluation where possible. D29 gp17, resolved as NP_046832.1, is experimentally identified by amino-terminal sequencing of a virion band in Ford and colleagues' genome paper. It supplies a structural specificity control, not proof of broad specificity against hydrolases or esterases or absence of every secondary function. [Primary paper](https://seaphages.org/media/forums/attachments/3b1d0020-e87c-4045-900d-00502e606468.pdf).

## Recent literature changes the interpretation

Bouras and colleagues' Current Protocols article is verified as Current Protocols 6(7):e70405, first published 7 July 2026, DOI 10.1002/cpz1.70405. It describes sequence, structure and synteny-based annotation tools. Phold and Phynteny are not used for primary functional assignments in this analysis. Its plotting discussion also warns that displaying a genome as a circle does not establish circular topology. [Publisher page](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.70405), [primary full-text archive](https://pmc.ncbi.nlm.nih.gov/articles/PMC13338908/).

The inspected 2025 Phynteny source is a preprint. Because the method uses genomic context to predict broad functions, a Phynteny “lysis” label cannot independently demonstrate conservation of a lysis neighborhood. If used at all, such predictions should be segregated as context-supported exploratory annotations. [Version 2 primary preprint](https://www.biorxiv.org/content/10.1101/2025.07.28.667340v2).

Li and colleagues' WXIN study was published 16 January 2026, despite “2025” occurring in its DOI. It provides a divergent LysB candidate with experimental esterase and bactericidal assays, including catalytic-residue mutagenesis. Bactericidal results depend on the reported surfactant conditions. [Primary article](https://link.springer.com/article/10.1186/s12866-025-04683-0).

James and colleagues' D29 study appeared online 13 July 2026 (mBio 17(8):e0073226; issue date 12 August). It reports concentration- and infection-stage-dependent LysB effects, including impaired progeny release after deletion and inhibition following exogenous treatment. The authors restrict their conclusions to D29 and leave open whether the effect is evolved regulation or an off-target consequence of enzyme accumulation. Architectural separation must not be used to infer expression timing or this phenotype in other phages. [Primary study](https://journals.asm.org/doi/10.1128/mbio.00732-26).

## The comparison this project can attempt

A defensible proposed contribution is a reproducible, uncertainty-aware comparison of LysB neighborhood conservation against matched background synteny after accounting for relatedness and uneven sampling. This is a research target, not a finding and not an assertion that nobody has done it before.

The analysis should answer whether separated LysB neighborhoods are more conserved than appropriate non-lysis neighborhoods among comparable genomes, not merely whether two related phages share adjacent genes. Controls should match gene density, local family prevalence, boundary availability and whole-genome relatedness. Cluster membership is a useful operational stratum, not a fully resolved phylogeny. Results from sparse strata must be marked untestable or descriptive.

Evidence suggesting repeated reorganization would require homologous flanking anchors, informative whole-genome relationships, and distinctions from genome rotation, reverse-complement orientation, annotation changes and broader rearrangement. Similar distances in unrelated genomes alone cannot establish recurrence, direction, selection or a mechanism.

Absence of a confident LysB call means undetected or unresolved unless sequence- and model-specific sensitivity supports a stronger conclusion. Preserve predicted functions, independent homology evidence, ambiguous alternatives and exclusions as separate fields. A negative or underpowered result is a valid outcome.

This study makes no priority claim. Functional and evolutionary conclusions are restricted to the evidence summarized in the ledger.

## Literature coverage

Searches covered primary studies of Gordonia LysB displacement and Yvonnetastic; Ms6 and D29 lysis genetics, holins and LysB; actinophage lysis synteny/mosaicism; and 2025/2026 studies. Primary publisher pages, deposited author articles and PubMed records were used. Nonprimary search hits were discovery aids only; they are not the evidence anchors in the ledger.


The ledger is deliberately explicit about “primary indexed text,” “primary full text,” and “abstract/metadata.” Rows based only on an abstract must be read at that evidence depth; detailed methods and uninspected supplements remain unresolved until obtained.
