from __future__ import annotations

import argparse
import re
from pathlib import Path
from textwrap import dedent


TITLE = "# Shared transcriptomic alterations between osteoarthritis and ovarian cancer reveal transcriptional divergence across distinct cellular contexts"


MANUSCRIPT = r"""
**Article type:** Original Research  
**Running title:** Transcriptomic divergence across OA and OC  
**Authors:** [Author names to be added]  
**Affiliations:** [Affiliations to be added]  
**Corresponding author:** [Name, postal address, and email to be added]  
**Target journal:** [To be selected before submission]

## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) differ in tissue of origin, natural history, and clinical phenotype, yet both involve extracellular-matrix remodeling, inflammatory signaling, and cellular stress. We asked whether their transcriptomes contain reproducible overlap and, critically, whether shared membership corresponds to shared direction and cellular context.

### Methods

Public human bulk-transcriptomic cohorts were analyzed on separate OA and OC tracks. Differentially expressed genes (DEGs) were compared using prespecified thresholds, followed by Gene Ontology (GO), Kyoto Encyclopedia of Genes and Genomes (KEGG), and Hallmark analyses. Fold-change and normalized-enrichment-score signs were compared explicitly. Disease-specific weighted gene co-expression network analysis (WGCNA), a transparent ten-gene evidence summary, and an auxiliary STRING association graph were used for characterization rather than mechanistic inference. Five single-cell datasets were processed with dataset-specific adapters; representative count-level OA and OC atlases were used to localize candidate expression with exact source-defined cell labels. Exploratory classification was retained as a supplementary reproducibility assessment.

### Results

At false-discovery rate (FDR) <0.05 and absolute log2 fold change >=1, 2,008 OA and 2,310 OC DEGs were identified, including 286 shared genes. Of these, 146 (51.0%) were directionally concordant and 140 (49.0%) discordant. Shared genes were enriched for extracellular-matrix, immune/cytokine, stress/hypoxia, and cell-cycle themes; cell cycle was the only KEGG term passing FDR <0.05. Ten Hallmark pathways were significant in both diseases, but six had opposite enrichment directions. The fixed candidate set (SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2) likewise showed heterogeneous disease-specific effects. Single-cell analyses localized these signals to different source-defined cellular contexts in OA and OC without treating labels as homologous states. Exploratory classification showed cohort-dependent discrimination.

### Conclusions

OA and OC exhibit partially shared transcriptomic alterations, but overlap frequently masks opposite gene and pathway states and distinct cellular localization. The results support context-dependent transcriptional divergence rather than a conserved cross-disease mechanism, diagnostic signature, or shared therapeutic target set.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; directionality; single-cell RNA sequencing; cellular context

## Introduction

Osteoarthritis is a whole-joint disease characterized by cartilage degeneration, extracellular-matrix remodeling, subchondral bone changes, and altered inflammatory and chondrocyte states [1,3,14-16]. Ovarian cancer, particularly high-grade serous disease, is a malignant process marked by genomic instability, dissemination, and a heterogeneous stromal and immune microenvironment [2,4,5,17,18]. These diseases therefore differ fundamentally in tissue architecture, comparator samples, disease kinetics, and clinical consequences.

Cross-disease transcriptomic comparisons often begin with an overlap of DEGs. An overlap can be informative, but membership alone does not establish that a gene changes in the same direction, participates in the same pathway state, or is expressed by comparable cells. This distinction is especially important when a chronic degenerative cartilage contrast is compared with a malignant tumor-reference contrast. Without direction-aware and cell-aware analysis, recurring stress or remodeling transcripts can be overinterpreted as evidence for a shared pathogenic mechanism.

We addressed three linked questions. First, which transcriptional alterations are shared between OA and OC? Second, do shared genes and pathways move in the same direction? Third, which exact source-defined cell populations carry the selected signals within representative single-cell atlases? OA and OC were analyzed separately and connected only through explicit comparisons of gene membership, direction, functional annotation, and cellular localization. WGCNA and STRING supplied secondary association context, while exploratory classification was deliberately moved to the supplement. This structure tests a focused proposition: partially shared transcriptomic alterations can coexist with disease-specific molecular states.

## Methods

### Study design and reproducibility

The study followed separate OA and OC analysis tracks (Figure 1). The primary evidence chain was discovery-cohort differential expression, shared-gene identification, functional characterization, direction analysis, transparent candidate characterization, and single-cell localization. Disease-specific WGCNA and an auxiliary STRING graph were included as association layers. Classification was treated as a supplementary reproducibility check and not as development of a clinical model. Other exploratory modules that did not directly answer the revised discovery-direction-localization question were outside the V3.2 submission scope; their historical outputs remain preserved in the earlier V3.1 project directory.

All analysis settings, data-source manifests, source data for figures, random seeds, and scope decisions are stored in the reproducible project. Randomized procedures used seed 20260726. The V3.2 presentation layer is regenerated by `run_submission_v32.ps1`; upstream data-processing and quality-control records remain part of the same project. This separation preserves provenance while preventing superseded exploratory modules from dominating the submitted narrative [22].

### Bulk cohorts and preprocessing

The OA discovery cohort was knee cartilage GSE114007 (20 OA and 18 non-OA samples) [3]. The OC discovery cohort was GSE18520 (53 advanced high-grade serous tumors and 10 normal ovarian surface-epithelium samples) [4]. External cohorts comprised GSE117999 (10 OA and 10 controls) and GSE82107 (10 OA and 7 controls) for OA, and GSE54388 (16 tumors and 6 normal samples) and GSE12470 (43 serous carcinomas and 10 normal peritoneal samples) for OC [5]. Cohort definitions, platforms, preprocessing decisions, and limitations are recorded in Table S1.

Expression matrices were checked for identifier agreement, duplicate features, missing or non-finite values, and phenotype balance. Repository metadata, rather than unsupervised clustering, defined disease and reference groups. OA and OC cohorts were normalized and modeled separately. As an unsupervised quality audit, principal-component analysis and sample correlations were calculated within each discovery cohort using the 1,000 genes with highest sample variance. No sample was removed solely because of its position in the principal-component or correlation space (Figure S6; Table S31).

### Differential expression and shared-gene direction

Differential expression was estimated separately with limma [6]. The primary threshold was Benjamini-Hochberg FDR <0.05 and absolute log2 fold change >=1 [7]. A gene was considered shared only if it met the primary rule in both diseases. Shared genes were classified into four quadrants: higher in both, lower in both, higher in OA/lower in OC, or lower in OA/higher in OC. Same-sign quadrants were designated concordant and opposite-sign quadrants discordant. Sensitivity analysis crossed FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0, and 1.5 (Figure S1; Table S3).

### Functional enrichment and pathway direction

GO and KEGG over-representation analyses used the 286 shared DEGs as the input set. Multiple testing was controlled by the Benjamini-Hochberg method. The main GO panel shows ten significant, non-redundant terms selected to represent cell-cycle, matrix, immune, and stress categories; the complete result is retained in Table S11b. All FDR-significant KEGG terms are shown.

Ranked Hallmark enrichment was estimated independently from complete OA and OC gene rankings using the same gene-set definitions and size limits [8,9]. Normalized enrichment scores (NES) were paired by Hallmark identifier. A pathway was called jointly significant when FDR <0.05 in both diseases. Jointly significant pathways were classified as concordant when NES signs matched and discordant when signs differed. This analysis compares transcriptional states; it does not establish pathway activity at the protein level or regulatory mechanism.

### Co-expression and candidate-gene characterization

Signed WGCNA networks were fitted separately in the OA and OC discovery cohorts [10]. Primary disease-associated modules were selected from module-trait correlations. Stability was assessed by soft-power perturbation, 2,000 fixed-module sample bootstraps, and leave-one-sample-out estimates (Figure S1; Table S4). Because external module preservation was unavailable, WGCNA was used as supporting co-expression context rather than proof of a conserved network.

The fixed ten-gene set comprised SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2. All ten met the shared-DEG rule. The set was retained as an interpretable evidence summary rather than an optimized predictive signature. Table S16 makes each evidence layer explicit: DEG status and direction, OA and OC WGCNA membership, original LASSO/random-forest votes, strict nested selection frequency, STRING connectivity, and single-cell detection. Deterministic top-5, top-10, and top-15 sensitivity analyses evaluated dependence on panel size without selecting a new endpoint (Figure S2; Tables S29-S30) [11-13].

### Auxiliary STRING association analysis

The 286 shared DEGs were submitted to STRING v12.0 for *Homo sapiens* identifier mapping without adding neighboring proteins [21]. The main auxiliary graph retained high-confidence physical associations (score >=0.700). Mapped isolates remained in denominator-based summaries. Degree and connected components were descriptive; graph topology was not used to select the ten candidates. STRING edges represent database associations and do not demonstrate tissue-specific binding or a shared OA-OC mechanism.

### Exploratory external-cohort classification

The fixed-direction molecular summary was evaluated in two external cohorts per disease using unsupervised structure, standardized score contrasts, label permutation, leave-one-out sensitivity, effect sizes, and receiver-operating-characteristic curves. These analyses were prespecified as retrospective cohort-separation checks, not as a locked clinical probability model. Because OA and OC contrasts differ in tissue, comparator, and signal amplitude, AUC estimates were interpreted within cohort and retained only in Supplementary Figure 3 and Tables S6 and S21-S22.

### Single-cell datasets, quality control, and localization

Five public single-cell datasets were processed with dataset-specific adapters: OA GSE104782, GSE169454, and GSE255460, and OC GSE154600 and GSE180661 [14-19]. Count-level datasets underwent cell- and gene-level quality control, mitochondrial-content assessment, doublet handling when count-compatible data and sample structure were available, normalization, dimensional reduction, and use of exact source annotations. TPM-only data were not treated as raw counts. Dataset-specific incompatibilities and permitted analysis layers are reported in Table S9.

The main single-cell figure was restricted to representative count-level atlases GSE255460 (OA) and GSE154600 (OC). OA and OC were not integrated into one latent space. UMAPs used deterministic display subsampling and exact source labels. For each candidate gene, the fraction of cells with detected expression was summarized within eligible source-defined labels. Labels with fewer than 100 cells and labels explicitly designated Unassigned, Ambiguous, or Other were excluded from ranking. The six leading labels per atlas are displayed. Detection fractions were interpreted within atlas and were not compared numerically between OA and OC. All five embeddings and additional sample-aware UCell and pseudobulk summaries are retained in supplementary materials (Figures S5 and Tables S9, S10, S19, S24, and S30) [20].

### Statistical reporting

Two-sided tests were used unless otherwise specified. Multiple-testing procedures and effect-direction definitions were fixed before interpretation. Exact sample or cell units are reported for each analysis layer. Confidence intervals accompany key resampling estimates. Null, unstable, and unfavorable results were retained. No causal, diagnostic, prognostic, or therapeutic inference was assigned to observational transcriptomic, database-association, or retrospective classification results.

## Results

### OA and OC shared 286 primary transcriptomic alterations

Separate discovery analyses identified 2,008 OA and 2,310 OC DEGs at FDR <0.05 and absolute log2 fold change >=1 (Figure 2A-C; Table S2). Their intersection contained 286 genes. Threshold sensitivity changed the absolute overlap, as expected, but retained a substantial shared set across all six prespecified combinations (Figure S1; Table S3). Unsupervised bulk audits showed phenotype-associated structure together with within-group heterogeneity; no outcome-informed sample exclusion was performed (Figure S6; Table S31).

The 286 shared genes did not constitute one uniform transcriptional state. Direction ordering separated 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC (Figure 2D). Thus, 146 genes (51.0%) were concordant and 140 (49.0%) discordant.

### Shared genes recurred in matrix, immune, stress, and cell-cycle themes

GO analysis identified significant terms spanning nuclear chromosome segregation and cell-cycle control, extracellular-matrix organization and structure, immune-effector and cytokine responses, integrated stress signaling, hypoxia, and oxygen-level responses (Figure 3A; Table S11b). These categories describe recurring biological themes rather than one conserved disease mechanism. No GO aging term passed the reported enrichment threshold; aging is therefore considered only as a possible contextual hypothesis in the Discussion.

KEGG analysis yielded one FDR-significant term, Cell cycle (Figure 3B; Table S11c). Hallmark analysis identified 22 significant pathways in OA and 16 in OC, including ten significant in both. Their paired directions were heterogeneous: four had matching NES signs and six had opposite signs (Figure 3C; Table S18).

### Direction analysis exposed gene- and pathway-level divergence

The near-balanced concordant and discordant classes were visible across the full shared set (Figure 4A). Seven of the ten fixed candidates were discordant. For example, SOX9 was lower in OA (log2 fold change -2.1) and higher in OC (+2.6); EFEMP1 was higher in OA (+3.1) and lower in OC (-2.8). JUNB, AKAP12, and DIRAS3 changed in the same direction across diseases, but concordant sign alone did not imply an identical cell source or function (Figure 4B; Table S16).

Pathway direction showed the same boundary. Epithelial-mesenchymal transition, coagulation, KRAS signaling up, complement, MTORC1 signaling, and glycolysis were significant in both diseases but had opposite NES signs, whereas E2F targets, G2M checkpoint, mitotic spindle, and apoptosis were concordant (Figure 4C; Table S18). Shared statistical significance therefore did not mean a shared pathway state.

### Candidate genes had heterogeneous co-expression and STRING support

The primary disease-associated WGCNA module-trait correlations were -0.951 in OA and -0.879 in OC. Both retained stable negative signs across 2,000 bootstrap resamples, although the smaller OA cohort and absence of independent module preservation limit network interpretation (Figure 5A; Figure S1; Table S4).

STRING mapped 275 of 286 shared DEGs. The high-confidence physical graph contained 62 edges among 46 connected products; mapped isolates were retained in the audit. Only JUNB among the ten fixed candidates was connected in this graph, while the other nine were mapped isolates (Figure 5B; Tables S25a-e). Accordingly, STRING topology neither generated nor validated the candidate set.

Candidate fold changes and the evidence matrix demonstrate why the ten genes should be read as a transparent descriptive summary (Figure 5C-D). WGCNA membership, original model votes, nested frequencies, PPI connectivity, and single-cell detection were not uniform across genes. Panel-size sensitivity preserved a predominance of discordant membership in deterministic top-5, top-10, and top-15 sets, while small-panel pathway rankings were not treated as stable mechanistic results (Figure S2; Tables S29-S30).

### Representative single-cell atlases localized candidates to distinct cellular contexts

Across five adapters, 1,187,436 cells were audited and 1,025,361 passed dataset-specific quality-control rules (Table S9). The main analysis used count-level GSE255460 and GSE154600 to retain a consistent localization question while avoiding cross-disease integration (Figure 6A-B). Exact source labels were preserved. In OA, leading labels included EC, HomC, HTC, preHTC, ProC, and RepC. In OC, leading labels included B.cell, Endothelial.cell, Fibroblast, Myeloid.cell, Ovarian.cancer.cell, and T.cell (Figure 6C).

Candidate localization differed across these atlases. SOX9 was broadly detected across several OA chondrocyte labels but was sparse across the displayed OC labels. JUNB was broadly detected in both atlases but occupied different source-defined cellular ecosystems. AKAP12 and CFI showed label-dependent detection within the OC atlas, whereas EFEMP1 detection was more prominent in selected OA labels. These are within-atlas localization patterns, not quantitative OA-versus-OC cell-type comparisons. In particular, Fibroblast was not relabeled as a cancer-associated fibroblast and Ovarian.cancer.cell was not relabeled as tumor epithelium.

All five atlas embeddings and panel-size localization sensitivity are shown in the supplement (Figures S2 and S5; Tables S24 and S30). They reinforce the central observation: the same gene name can be embedded in different cellular compositions and should not be assigned one cross-disease biological meaning.

### Exploratory classification was cohort dependent

Exploratory classification analysis showed cohort-dependent discrimination performance, consistent with biological heterogeneity among datasets (Figure S3; Tables S6 and S21-S22). Strong separation in OC tumor-reference cohorts and modest separation in OA cohorts were not interpreted as competing estimates of one diagnostic model. This module was retained only as a secondary reproducibility assessment.

## Discussion

### Shared genes are not necessarily shared programs

This study identified a reproducible overlap of 286 DEGs between OA and OC, but approximately half changed in opposite directions. The overlap was enriched for matrix remodeling, immune/cytokine signaling, stress/hypoxia responses, and cell-cycle processes. These themes are biologically plausible in both chronic tissue degeneration and malignancy, yet their recurrence does not require one common disease program. A shared gene list can instead capture generic responses to injury, altered tissue composition, proliferation, metabolic pressure, or extracellular-matrix turnover.

The Hallmark results make this distinction quantitative. Six of ten pathways significant in both diseases had opposite NES signs, and seven of ten candidates had discordant fold-change signs. Consequently, the data do not support the stronger title or conclusion that OA and OC share a conserved molecular program. The more defensible interpretation is partial overlap accompanied by transcriptional divergence.

### Distinct tissues and cell populations can give the same gene different meanings

OA cartilage and ovarian tumor-reference cohorts differ in tissue architecture, cellular composition, disease duration, and signal amplitude. Bulk overlap can therefore arise from different mixtures of responding cells. The single-cell analysis addressed this issue without forcing OA and OC into one integrated embedding or renaming source annotations. Candidate detection was distributed among OA cartilage labels such as HomC, HTC, preHTC, ProC, and RepC, and among OC labels such as Fibroblast, Ovarian.cancer.cell, endothelial, myeloid, B, and T cells.

These observations support context-dependent localization but do not prove that cell context causes the bulk direction differences. Definitive testing will require spatial transcriptomics, matched protein measurements, lineage-resolved perturbation, and models that preserve tissue architecture. Organoid or co-culture systems could test whether candidate responses change with inflammatory, matrix, or malignant microenvironments.

### Candidate, network, and classification results should remain secondary

The ten genes improve interpretability because their evidence can be inspected across DEG direction, WGCNA, original model votes, STRING, and single-cell localization. They should not be presented as a diagnostic signature. Only JUNB was connected in the high-confidence physical STRING graph, emphasizing that PPI topology did not select or validate the panel. Likewise, WGCNA describes within-cohort co-expression and does not establish preservation between diseases.

Exploratory classification was deliberately moved to the supplement. OA and OC classification tasks are biologically different: OA compares chronic degeneration within cartilage, whereas OC cohorts compare malignant tissue with ovarian or peritoneal reference samples. Cohort-dependent discrimination is therefore compatible with the main conclusion and should not be repaired with a nomogram or decision-curve analysis in a study that does not claim clinical prediction.

### Strengths and limitations

Strengths include separate disease-specific discovery, prespecified DEG thresholds, explicit gene and pathway directions, complete reporting of the shared set, WGCNA perturbation and bootstrap checks, transparent candidate accounting, auxiliary STRING mapping with isolates retained, five dataset-specific single-cell adapters, exact cell labels, figure source data, and a one-command V3.2 presentation build. Unfavorable and null observations were retained rather than hidden.

The study is observational and retrospective. Tissues, platforms, comparator definitions, disease stages, and clinical composition differed across cohorts. Complete technical and clinical covariates were unavailable for consistent adjustment. DEG membership depended on statistical thresholds; WGCNA lacked independent module-preservation testing; the candidate set retained post-discovery ranking inputs; and single-cell datasets differed in annotation granularity, controls, and sample replication. Detection fractions do not measure protein abundance or cell function. STRING is database derived. Exploratory classification was not a locked clinical model. Genetic causal inference was beyond the scope of this focused transcriptomic comparison. Spatial, protein-level, prospective, perturbational, and multi-cancer specificity validation remain necessary.

## Conclusions

OA and OC share 286 primary transcriptomic alterations, but the overlap is almost evenly divided between concordant and discordant genes, and most jointly significant Hallmark pathways change in opposite directions. Candidate genes occupy different source-defined cellular contexts, while co-expression, STRING, and classification provide secondary association evidence only. These findings define context-dependent transcriptional divergence across two distinct diseases; they do not establish a conserved mechanism, diagnostic signature, or common therapeutic target.

## Data and code availability

All GEO accession identifiers, cohort definitions, parameters, source-data tables, and audit records are provided in the supplementary materials and reproducible project. The V3.2 presentation layer can be regenerated with `run_submission_v32.ps1`. A public repository URL and archival DOI will be added before submission: [repository URL to be added]; [Zenodo DOI to be added]. Controlled or licensed source data remain subject to their original access terms.

## Ethics approval and consent to participate

This study reanalyzed publicly available, de-identified datasets. No new participant recruitment or intervention was performed. The authors should confirm the target journal's wording and the ethics statements of each source study before submission.

## Consent for publication

Not applicable to the present secondary analysis; author confirmation is required before submission.

## Competing interests

[To be completed by all authors before submission.]

## Funding

[To be completed before submission.]

## Author contributions

[To be completed using the target journal's contributor taxonomy.]

## Acknowledgements

[To be completed before submission.]

## References

1. Glyn-Jones S, Palmer AJR, Agricola R, et al. Osteoarthritis. *Lancet*. 2015;386:376-387. doi:10.1016/S0140-6736(14)60802-3.
2. Lheureux S, Gourley C, Vergote I, Oza AM. Epithelial ovarian cancer. *Lancet*. 2019;393:1240-1253. doi:10.1016/S0140-6736(18)32552-2.
3. Fisch KM, Gamini R, Alvarez-Garcia O, et al. Identification of transcription factors responsible for dysregulated networks in human osteoarthritis cartilage by global gene expression analysis. *Osteoarthritis Cartilage*. 2018;26:1531-1538. doi:10.1016/j.joca.2018.07.012.
4. Mok SC, Bonome T, Vathipadiekal V, et al. A gene signature predictive for outcome in advanced ovarian cancer identifies a novel survival factor: microfibril-associated glycoprotein 2. *Cancer Cell*. 2009;16:521-532. doi:10.1016/j.ccr.2009.10.018.
5. Yoshihara K, Tajima A, Komata D, et al. Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis. *Cancer Sci*. 2009;100:1421-1428. doi:10.1111/j.1349-7006.2009.01204.x.
6. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43:e47. doi:10.1093/nar/gkv007.
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B*. 1995;57:289-300.
8. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA*. 2005;102:15545-15550. doi:10.1073/pnas.0506580102.
9. Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database Hallmark gene set collection. *Cell Syst*. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004.
10. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559. doi:10.1186/1471-2105-9-559.
11. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33:1-22. doi:10.18637/jss.v033.i01.
12. Breiman L. Random forests. *Mach Learn*. 2001;45:5-32. doi:10.1023/A:1010933404324.
13. Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*. 2006;7:91. doi:10.1186/1471-2105-7-91.
14. Ji Q, Zheng Y, Zhang G, et al. Single-cell RNA-seq analysis reveals the progression of human osteoarthritis. *Ann Rheum Dis*. 2019;78:100-110. doi:10.1136/annrheumdis-2017-212863.
15. Fu W, Hettinghouse A, Chen Y, et al. 14-3-3 epsilon is an intracellular component of TNFR2 receptor complex and its activation protects against osteoarthritis. *Ann Rheum Dis*. 2021;80:1615-1627. doi:10.1136/annrheumdis-2021-220000.
16. Fan Y, Bian X, Sun S, et al. Unveiling inflammatory and prehypertrophic cell populations as key contributors to knee cartilage degeneration in osteoarthritis using multi-omics data integration. *Ann Rheum Dis*. 2024;83:776-790. doi:10.1136/ard-2023-224420.
17. Geistlinger L, Oh S, Ramos M, et al. Multiomic analysis of subtype evolution and heterogeneity in high-grade serous ovarian carcinoma. *Cancer Res*. 2020;80:4335-4345. doi:10.1158/0008-5472.CAN-20-0521.
18. Vazquez-Garcia I, Uhlitz F, Ceglia N, et al. Ovarian cancer mutational processes drive site-specific immune evasion. *Nature*. 2022;612:778-786. doi:10.1038/s41586-022-05496-1.
19. Germain PL, Lun A, Garcia Meixide C, Macnair W, Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Research*. 2021;10:979. doi:10.12688/f1000research.73600.2.
20. Andreatta M, Carmona SJ. UCell: robust and scalable single-cell gene signature scoring. *Comput Struct Biotechnol J*. 2021;19:3796-3798. doi:10.1016/j.csbj.2021.06.043.
21. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Res*. 2023;51:D638-D646. doi:10.1093/nar/gkac1000.
22. Wilkinson MD, Dumontier M, Aalbersberg IJJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data*. 2016;3:160018. doi:10.1038/sdata.2016.18.
"""


SUPPLEMENTARY_INDEX = r"""
## Supplementary table index

- **Table S1:** Data sources, cohorts, platforms, and analysis roles.
- **Table S2:** The 286 shared DEGs with OA and OC effect estimates and direction classes.
- **Table S3a-b:** DEG-threshold sensitivity and gene membership.
- **Table S4:** WGCNA soft-power, bootstrap, and leave-one-out stability.
- **Table S5:** Secondary nested candidate-ranking stability.
- **Table S6:** External exploratory classification metrics.
- **Table S9:** Single-cell QC, compatibility, and permitted analysis layers.
- **Table S10:** Candidate single-cell evidence and eligible pseudobulk summaries.
- **Table S11a-c:** Hallmark, GO, and KEGG enrichment results.
- **Table S16:** Transparent candidate-gene evidence matrix.
- **Table S18:** Complete Hallmark pathway-direction matrix.
- **Table S19:** Gene-cell-function context matrix using exact source labels.
- **Tables S21-S22:** External cohort separability, effect sizes, and reliability details.
- **Tables S24a-c:** Cell-context specificity and sample-aware UCell summaries.
- **Tables S25a-e:** STRING mapping, edges, topology, and label-permutation audits.
- **Tables S29-S30:** Candidate-set size and localization sensitivity.
- **Tables S31a-b:** Bulk discovery PCA and sample-correlation QC.
"""


RESPONSE = r"""
# V3.2 targeted revision response matrix

This revision implements the strategic recommendation to replace analysis accumulation with a focused discovery-direction-localization narrative. The submitted V3.2 package is distinct from the preserved V3.1 historical output.

| Recommendation | V3.2 action | Location | Scientific boundary |
|---|---|---|---|
| Reframe the paper around shared alterations rather than a shared mechanism | Rewrote the title, abstract, Introduction, Results, Discussion, and conclusion around partial overlap plus transcriptional divergence | Title; Abstract; Introduction; Figures 1 and 7 | The suggested phrase “conserved molecular programs” was not used because 140/286 genes and 6/10 jointly significant Hallmarks were directionally discordant |
| Stop adding analyses | Fixed a seven-main-figure structure and excluded MR, CellChat/NicheNet, TF-miRNA, TCGA/immune/HPA, DCA, and nomogram from V3.2 submission scope | V32 scope-decision table; Figure 1 | Historical outputs remain in V3.1 for provenance; exclusion does not imply deletion or invalidity |
| Prioritize bulk discovery | Made DEG overlap and direction-ordered shared genes the first result | Figure 2; Tables S2-S3 | Shared membership is not interpreted as shared causation |
| Strengthen functional interpretation | Consolidated representative significant GO terms, all significant KEGG results, and paired Hallmark directions | Figure 3; Tables S11 and S18 | GO terms summarize recurring themes; no aging enrichment is claimed |
| Make direction heterogeneity central | Added explicit gene quadrants, candidate fold-change heatmap, and OA-versus-OC Hallmark direction comparison | Figure 4 | Opposite NES signs are descriptive and do not prove regulatory decoupling |
| Retain WGCNA and PPI only if useful | Kept WGCNA as disease-specific co-expression support and STRING as an auxiliary association panel | Figure 5; Figure S1; Tables S4 and S25 | PPI topology did not select candidates; only JUNB was connected among the ten candidates |
| Clarify why ten genes were retained | Added a transparent evidence map and explicit Methods sentence that the set is an interpretable summary rather than an optimized predictive signature | Figure 5D; Table S16 | No diagnostic or therapeutic claim is made |
| Upgrade single-cell interpretation without forced integration | Focused the main text on representative count-level GSE255460 and GSE154600, exact source labels, and within-atlas candidate localization; retained all five atlases in the supplement | Figure 6; Figure S5; Tables S9, S10, S19, S24, and S30 | OA and OC labels are not treated as homologous; Fibroblast and Ovarian.cancer.cell are not relabeled |
| Demote AUC and machine learning | Moved exploratory classification to Supplementary Figure 3 and shortened the main Results to one bounded paragraph | Figure S3; Tables S5-S6 and S21-S22 | Classification is a retrospective reproducibility assessment, not a clinical model |
| Use conventional article figure styling | Rebuilt all main figures at 185-mm width with white backgrounds, Arial typography, consistent panel letters, restrained line weights, and blue-orange color-blind-aware encoding with shape redundancy | Figures 1-7; Figures S1-S7 | Styling follows general SCI conventions; no journal-specific compliance is claimed before a target journal is chosen |
| Simplify the Discussion | Reorganized around why overlap occurs, why the diseases remain different, and what experiments are needed | Discussion | Spatial, organoid/co-culture, protein, and perturbation work is presented as future validation, not completed evidence |
"""


def build(legends_path: Path) -> str:
    legends = legends_path.read_text(encoding="utf-8").strip()
    body = TITLE + "\n\n" + dedent(MANUSCRIPT).strip()
    body += "\n\n## Figure legends\n\n" + legends
    body += "\n\n" + dedent(SUPPLEMENTARY_INDEX).strip() + "\n"
    body = re.sub(r"\n{3,}", "\n\n", body)
    return body


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-legends", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()
    manuscript = build(args.figure_legends)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(manuscript, encoding="utf-8")
    args.response_output.write_text(dedent(RESPONSE).strip() + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.response_output}")


if __name__ == "__main__":
    main()
