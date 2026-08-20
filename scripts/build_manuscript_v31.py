from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import build_manuscript_v22 as v22


TITLE = "# Context-dependent transcriptomic convergence between osteoarthritis and ovarian cancer: multi-layer evidence from bulk transcriptomics, single-cell ecosystems, and genetic-liability analyses"


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in different tissues and have different clinical trajectories. We asked whether their transcriptomes nevertheless show partial molecular convergence, and whether any overlap is directionally concordant, cellularly comparable, reproducible across cohorts, or supported by disease-to-disease inherited liability.

### Methods

Public human bulk and single-cell transcriptomic datasets were analyzed on separate OA and OC tracks. Direction-aware differential-expression overlap and paired Hallmark enrichment were integrated with disease-specific co-expression evidence and a transparent ten-gene evidence panel. Candidate selection was distinguished from strict nested stability estimation and retrospective cross-cohort molecular separability. Panel-size sensitivity compared deterministic top-5, top-10, and top-15 extensions of the original model-vote ranking. Five single-cell atlases were evaluated using exact source labels, detection-based cell-context specificity, sample-aware UCell scoring, and eligible pseudobulk contrasts. Unsupervised discovery-cohort principal-component and sample-correlation audits, supplementary protein-association and communication analyses, contextual tissue resources, and bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes, with 286 shared: 146 (51.0%) were concordant and 140 (49.0%) discordant. Six of 10 Hallmark sets significant in both diseases had opposite normalized-enrichment-score signs. Discordant membership remained predominant in the top-5, top-10, and top-15 evidence panels (60.0%, 70.0%, and 66.7%), although small-panel Hallmark over-representation did not survive false-discovery-rate correction. Across 1,025,361 quality-control-pass single cells, candidate localization remained atlas specific; OA detection profiles consistently favored HomC in one count-level atlas, whereas the highest OC label shifted from Ovarian.cancer.cell at top 5 to Fibroblast at top 10 and top 15. Retrospective separability was modest in OA cohorts and strong in OC tumor-reference contrasts. Bidirectional MR found no evidence that genetic liability to OA affected OC risk, or that genetic liability to OC affected OA risk, under the selected datasets, instruments, and assumptions.

### Conclusions

OA and OC show partial context-dependent transcriptomic convergence, but shared membership does not imply shared direction, pathway state, cellular meaning, universal separability, or disease-to-disease inherited causality. The evidence panel is an interpretable summary rather than an optimized diagnostic signature, and the analyses define testable contexts rather than a shared disease mechanism or therapeutic target set.

**Keywords:** osteoarthritis; ovarian cancer; context-dependent convergence; cross-disease transcriptomics; single-cell RNA sequencing; pathway direction; evidence panel; genetic liability
"""


INTRODUCTION = """## Introduction

Osteoarthritis is a whole-joint disorder characterized by cartilage loss, extracellular-matrix remodeling, subchondral bone change, synovial inflammation, and altered chondrocyte states [1,3,17-19]. High-grade serous ovarian cancer is an aggressive malignancy marked by genomic instability, dissemination, and a heterogeneous stromal and immune microenvironment [2,14,20,21]. Their anatomy, comparator tissues, and clinical trajectories are fundamentally different. A cross-disease comparison must therefore ask a narrower question than whether both conditions contain inflammatory or remodeling-associated transcripts.

The premise of this study was context-dependent transcriptomic convergence: unrelated tissues may recruit partially overlapping stress, matrix, metabolic, or immune-response genes without engaging the same disease program. Four questions follow. First, do OA and OC share transcriptional alterations? Second, are the directions and pathway states concordant? Third, do the shared signals occupy comparable cellular contexts? Fourth, is the observed convergence consistent with a disease-to-disease inherited-liability effect? These questions require different evidence layers and cannot be answered by a DEG intersection alone.

Cross-disease transcriptomic studies are vulnerable to several interpretive shortcuts. Shared membership can conceal opposite fold-change signs; bulk overlap can reflect different cell populations; feature screening outside resampling can exaggerate stability; tumor-reference and chronic-degeneration contrasts do not measure one common classification task; and expression association does not establish genetic causality. Database-derived protein, regulatory, or communication networks add hypotheses but do not repair these limitations unless their inference boundaries are explicit.

We therefore analyzed OA and OC on separate disease-specific tracks and connected them only through prespecified gene-, pathway-, evidence-, and interpretation-level comparisons. The primary narrative follows direction-aware gene overlap, paired pathway direction, cellular localization, and the inherited-causality boundary. Co-expression, machine-learning stability, protein-association, CellChat/NicheNet, HPA, TCGA-OV, immune, and focused upstream analyses were retained as secondary or supplementary context. The goal was not to construct a universal diagnostic model, but to determine why molecular overlap can be shared yet nonidentical.
"""


STUDY_DESIGN = """### Study design and reproducibility

OA and OC were analyzed separately at every disease-specific stage and connected only by explicit gene-, pathway-, evidence-, and boundary-level comparisons. Figure 1 asks four questions: whether transcriptional alterations are shared, whether their directions and pathway states agree, whether they occupy similar cellular contexts, and whether disease-to-disease genetic liability is supported. The primary evidence chain therefore comprised bulk discovery and direction, paired Hallmark direction, a transparent evidence panel, cross-cohort molecular separability, five-atlas single-cell localization, and supplementary bidirectional MR. Disease-specific WGCNA, strict nested feature stability, STRING protein associations, sample-consensus CellChat, a bounded NicheNet prior overlay, HPA/TCGA audits, immune signatures, and focused upstream resources supplied secondary context. Fixed configurations, manifests, seeds, cached external responses, tests, result tables, exact figure source data, and a claim-evidence registry are stored in the one-command project; randomized analyses used seed 20260726.
"""


BULK_METHOD = """### Bulk transcriptomic cohorts, metadata, and unsupervised quality audit

Discovery cohorts were OA knee cartilage GSE114007 (20 OA, 18 non-OA) [3] and OC GSE18520 (53 advanced high-grade serous tumors, 10 normal ovarian surface-epithelium samples) [4]. External cohorts were GSE117999 (10 OA, 10 controls) and GSE82107 (10 OA, 7 controls) for OA, and GSE54388 (16 tumors, 6 normal samples) and GSE12470 (43 serous carcinomas, 10 normal peritoneal samples) for OC [5]. Expression and metadata were checked for identifier agreement, phenotype balance, duplicate features, and non-finite values; repository metadata, not expression clustering, defined groups. Uneven availability of age, sex, stage, grade, treatment, and technical variables prevented consistent patient-level covariate adjustment (Table S1).

Each discovery cohort underwent a separate unsupervised audit using the 1,000 genes with highest sample variance. Centered and scaled principal-component analysis was fitted without phenotype labels; labels were added only for display. Pearson correlations were calculated among samples using the same genes, and each sample's median and minimum correlation to all other samples were recorded. The three lowest median-correlation samples were labeled for audit, but no sample was removed solely because it occupied an extreme principal-component or correlation position. OA and OC were not combined, and this audit was not used as a batch-correction or outcome-informed filtering step (Figure S16; Table S31).
"""


PPI_METHOD = """### Supplementary protein-association landscape

The 286 shared DEGs were submitted without additional neighbors to STRING v12.0 for *Homo sapiens* identifier mapping and association retrieval [36]. The primary graph retained physical STRING associations with score >=0.700; a high-confidence functional-association graph was retained as sensitivity. Mapped isolates remained in all denominators. Degree, normalized betweenness, density, component count, and largest-component size were descriptive. Concordant and discordant induced subgraphs were compared with 10,000 fixed-size permutations of direction labels on the fixed mapped graph. Network topology was not used to choose the evidence panel. This supplementary analysis describes a database-derived protein-association landscape; it does not demonstrate physical binding, tissue-specific activity, or a shared OA-OC mechanism.
"""


CANDIDATE_METHOD = """### Evidence-panel construction, panel-size sensitivity, and strict nested stability

Candidate construction and model-performance estimation were treated as separate operations. Shared primary differentially expressed genes were ranked by adjusted significance and absolute effect size in both diseases and restricted to the union of disease-specific primary WGCNA modules when the prespecified minimum was met. Disease-specific LASSO and random forest supplied complementary model-vote evidence. SOX9 was the only cross-disease model consensus; the original ranked-vote rule completed the fixed ten-gene evidence panel (SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2). Table S16 reports direction, WGCNA membership, original LASSO/random-forest support, strict nested frequency, and single-cell context. The panel is an interpretable evidence summary, not an optimized predictive signature.

Sensitivity to the arbitrary panel size was evaluated without re-optimizing any endpoint. Deterministic top-5, top-10, and top-15 panels were taken from the same original four-model vote ranking; ranks 11-15 were KIT, MYZAP, NOD2, OGN, and RTN1. We compared concordant/discordant membership and absolute discovery effects. Descriptive Hallmark over-representation used the 286 shared DEGs as the background, Fisher/hypergeometric right-tail probabilities, and Benjamini-Hochberg correction within panel size. Pathway-profile similarity used Spearman correlation of fold-enrichment profiles and top-10-term Jaccard overlap. Detection-based cell-context sensitivity was calculated in the count-level GSE255460 OA and GSE154600 OC atlases: per-gene detection fractions were aggregated within exact source labels, labels with fewer than 100 cells or named Unassigned, Ambiguous, or Other were excluded, and panel means were ranked within atlas. This analysis tests dependence on panel size; it does not define pathway activity or select a new signature (Figure S15; Tables S29-S30).

Internal molecular separation was re-estimated independently of the fixed panel. In 50 repeated five-fold outer resamples, all measured genes were ranked by training-only two-group statistics and the top 100 entered modeling. Five-fold inner cross-validation selected LASSO lambda; random-forest mtry was selected from 2, 5, and 10 by training-fold out-of-bag AUC before fitting 300 trees [11,12]. Outer test samples were excluded from screening, tuning, scaling, and fitting. Balanced accuracy, Brier score, AUC, model-specific feature importance, and selection frequency were summarized across outer predictions [32]. AUC remained secondary.
"""


SEPARABILITY_METHOD = """### Disease-context-dependent separability of the evidence panel

External gene orientation was fixed from discovery log2 fold change. The signed molecular summary averaged standardized expression multiplied by discovery signs. AUC and DeLong 95% confidence intervals were calculated without choosing direction in validation data [13]. GSE54388 also underwent unsupervised PCA using the 2,000 genes with highest sample standard deviation, with phenotype labels used only for display. Robustness used 1,000 label permutations and sequential sample omission.

Because the same summary was evaluated in tasks with different tissues and comparator scales, each cohort was annotated by tissue, comparator, sample size, AUC interval, permutation result, and leave-one-out range. Hedges g summarized the standardized disease-reference score contrast. Descriptive random-effects summaries were calculated separately for the two OA and two OC cohorts; with two cohorts per disease, pooled effects and heterogeneity remained imprecise. The score was not a locked probability model. Leave-one-sample-out ridge recalibration was therefore retained only as a supplementary reliability sensitivity, not transported calibration or clinical utility. No nomogram or decision-curve analysis was performed because no clinical decision model or threshold was specified.
"""


SINGLE_CELL_METHOD = """### Single-cell localization, cell-context specificity, and UCell scoring

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17-21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. OA and OC were not integrated into a shared latent space, and labels remained dataset specific. Disease inference used pseudobulk only where biological replication and an interpretable contrast were available [23,24].

For each candidate and eligible exact source label, the Candidate Cell Context Specificity Score (CCSS) was calculated from within-dataset detection fractions. Labels named Unassigned, Ambiguous, or Other and strata below the cell threshold were excluded. Thresholds of 50, 100, and 200 cells were audited; 100 cells was primary. Medians were combined across datasets only for identical labels. Fibroblast was not relabeled as cancer-associated fibroblast, and Ovarian.cancer.cell was not relabeled as tumor epithelium.

UCell v2.14.0 scored the unsigned ten-gene evidence panel and a discovery-direction-compatible version from cell-level ranks with maxRank 1500 [37]. Scores were summarized first within biological sample and then within cell type and displayed as within-atlas ranks. The unsigned score represents joint rank abundance, not pathway activity, and neither score is a single-cell differential-expression effect. Exploratory marker-based functional joins and eligible pseudobulk contrasts remained descriptive and sample aware.
"""


COMMUNICATION_METHOD = """### Supplementary sample-consensus CellChat and bounded NicheNet context

Communication inference was restricted to datasets with sufficient sample-resolved expression and usable cell labels: five OC samples from GSE154600 and three control plus eight OA donors from GSE255460. CellChat v2 was fitted independently to each biological sample after deterministic cell-type-balanced subsampling and then aggregated using prespecified context-specific sample-consensus thresholds [38,39]. Cells were never treated as independent patient replicates. OA and OC probabilities were not compared because tissues, label systems, and comparator designs differ.

Consensus CellChat ligands were projected onto the official NicheNet v2 human ligand-target matrix [40,41]. Regulatory-potential weights linking those ligands to the fixed evidence panel were retained as an external prior-consistency overlay. Full ligand-activity inference was not performed because OC lacked a symmetric disease/reference receiver-cell contrast and the fixed target set was too small for an unbiased activity test. CellChat and NicheNet results are supplementary hypothesis-generating context; they do not establish spatial contact, signaling flux, regulation, mediation, or mechanism.
"""


RESULTS = """## Results

### Discovery-cohort audits showed phenotype-associated separation and substantial internal structure

The study included two bulk discovery cohorts, four external bulk cohorts, TCGA-OV, five single-cell datasets, HPA annotations, focused supplementary resources, and two GWAS datasets (Figure 1; Table S1). In separate unsupervised discovery-cohort analyses, OA PC1 and PC2 explained 36.0% and 13.8% of variance, whereas OC PC1 and PC2 explained 41.5% and 16.1% (Figure S16; Table S31). Both cohorts showed phenotype-associated separation together with marked within-group substructure. Median sample correlations ranged from 0.139 to 0.627 in OA and from 0.080 to 0.809 in OC. These audits motivated explicit cohort-composition limitations; no sample was removed on the basis of PCA position or correlation alone.

### The shared DEG set was directionally heterogeneous

At the primary threshold, 2,008 genes were differentially expressed in OA and 2,310 in OC, with 286 shared (Figure 2; Table S2). The four direction quadrants contained 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. Only 146/286 (51.0%) were concordant. Directional heterogeneity persisted across all six prespecified DEG thresholds (Table S3).

All ten evidence-panel genes met the primary rule in both diseases, but only AKAP12, JUNB, and DIRAS3 were concordant. SOX9 was lower in OA (log2 fold change -2.072) and higher in OC (+2.612), whereas EFEMP1 was higher in OA (+3.092) and lower in OC (-2.781). Shared membership was therefore not a uniform activation or suppression program.

### Pathway-level convergence was also directionally heterogeneous

Complete paired analysis retained all 50 Hallmark sets (Figure 3; Figure S7; Table S18). Twenty-two were significant in OA and 16 in OC; 10 reached FDR <0.05 in both diseases. Four had matching NES signs and six had opposite signs. Epithelial-mesenchymal transition was positive in OA (NES 2.347) and negative in OC (-1.833), whereas glycolysis was negative in OA (-1.410) and positive in OC (1.774). Paired membership therefore did not indicate the same pathway state or a conserved mechanism.

The supplementary STRING landscape mapped 275/286 shared DEGs and retained 11 as unmapped (Figure S12; Tables S25a-e). The high-confidence physical graph contained 62 edges among 46 connected products. Concordant labels were more densely organized than expected under fixed-size label permutation, whereas discordant edge count did not differ from its null. This was a database-association pattern, not evidence of tissue-specific protein interactions or a common OA-OC mechanism.

### The evidence panel summarized heterogeneous support and remained panel-size dependent

Disease-specific WGCNA module-trait correlations were -0.951 in OA and -0.879 in OC, with bootstrap sign stability of 1.000 across 2,000 resamples (Figure S8; Table S4). SOX9 was the only cross-disease original model consensus. The other genes entered by the prespecified ranked-vote completion, and strict nested selection frequencies were heterogeneous: maximum OA frequency among panel genes ranged from 0 to 0.940 and maximum OC frequency from 0 to 0.200 (Figure S8; Tables S5 and S16). The ten genes therefore did not form a uniformly stable predictive panel.

Panel-size sensitivity clarified what was and was not robust (Figure S15; Tables S29-S30). Discordant genes comprised 3/5 (60.0%), 7/10 (70.0%), and 10/15 (66.7%) of deterministic top-5, top-10, and top-15 panels, supporting the direction-heterogeneity conclusion. Descriptive Hallmark fold-enrichment profiles correlated from 0.540 to 0.883 across panel sizes, but no panel-level Hallmark term survived FDR correction and top-10-term Jaccard overlap ranged from 0.429 to 0.667. The pathway ordering of a small candidate panel was therefore not treated as a stable mechanistic result.

Candidate-centered residual associations for SOX9, DDIT3, BNC1, and AKAP12 provided disease-specific hypotheses (Figure S9; Table S20). They were group-adjusted co-expression contexts, not perturbation, mediation, or pathway-activation evidence.

### Disease context determined retrospective separability of the evidence panel

As a secondary reproducibility assessment, Figure 5 presents unsupervised structure, standardized score contrasts, permutation behavior, and ROC in that order. In GSE54388, PC1 and PC2 explained 33.0% and 15.1% of variance. Empirical permutation P values were 0.469 and 0.207 in the two OA cohorts and 0.001 in both OC cohorts. Leave-one-out AUC ranges were 0.456-0.589 and 0.583-0.698 in OA, remained 1.000 in GSE54388, and were 0.971-0.993 in GSE12470.

The direction-fixed summary yielded AUCs of 0.520 (95% CI 0.250-0.790) and 0.629 (0.350-0.907) in OA, compared with 1.000 and 0.979 (0.944-1.000) in OC (Figure 5D; Tables S6 and S21). Hedges g values were 0.015 and 0.306 in OA versus 5.243 and 2.692 in OC. These are not competing estimates of one universal model: OA contrasts chronic within-tissue degeneration, whereas OC contrasts malignant tissue with ovarian or peritoneal references. The result is disease-context-dependent retrospective molecular separability, not a diagnostic claim.

Cross-fitted calibration was stable only in GSE12470; the other three cohorts were numerically or directionally unstable (Figure S10; Table S22). Calibration did not rescue modest OA separation or convert the OC contrasts into clinical validation.

### Single-cell dissection revealed distinct and panel-size-sensitive cellular contexts

The five adapters audited 1,187,436 cells and retained 1,025,361 quality-control-pass cells (Figure 4; Figure S3; Table S9). CCSS and sample-aware UCell retained exact source labels and localized the evidence panel to HomC or preHTC contexts in OA atlases and to Fibroblast or Ovarian.cancer.cell contexts in OC atlases (Tables S24a-c). These labels were not asserted to be homologous, and scores were not compared numerically across diseases.

The count-level panel-size sensitivity showed high within-atlas profile correlations but a meaningful boundary. In GSE255460, top-5/top-10/top-15 localization-profile correlations ranged from 0.709 to 0.918 and HomC remained the highest mean-detection label. In GSE154600, correlations were 0.964-1.000, but the highest label shifted from Ovarian.cancer.cell for top 5 to Fibroblast for top 10 and top 15 (Figure S15; Table S30). Broad cellular structure was reproducible, whereas the single top OC label depended on panel composition.

The gene-cell-function matrix linked SOX9 to an OA hypertrophic-chondrocyte mineralization context and an OC Ovarian.cancer.cell translation context, and AKAP12 to OA preHTC extracellular-matrix organization and OC Fibroblast collagen-fibril organization (Table S19). These data support different cellular meanings, not cell-type-specific causal functions.

### Supplementary communication, regulatory, tissue, and clinical contexts remained bounded

CellChat completed independently for 16 biological samples and yielded 199 sample-consensus context-pathway records (Figure S13; Tables S26a-c). Although 526 consensus interactions overlapped a shared DEG, none contained a fixed panel gene as a direct ligand or receptor (Table S26d). NicheNet supplied prior paths from consensus ligands to all ten panel targets, but ligand activity was not estimable under the available design (Figure S14; Tables S27-S28). The communication layer therefore did not supply a direct signaling explanation for the evidence panel.

Focused KnockTF and miRTarBase matrices, HPA normal references, rank-based immune signatures, and TCGA-OV analyses remained supplementary context (Figures S4-S5 and S11; Tables S8, S14-S15, and S23). The age/stage-adjusted TCGA continuous-score association had a hazard ratio of 1.262 per standard deviation (95% CI 1.076-1.479), but optimism-corrected concordance was 0.586 and feature selection was unstable. This was not an independently validated prognostic model.

### Bidirectional MR defined the inherited-causality boundary

MR used ebi-a-GCST007092 and ieu-a-1120 with 21 and 11 instruments. Inverse-variance-weighted estimates were null for OA-to-OC (OR 1.015, 95% CI 0.900-1.144; P=0.811) and OC-to-OA (OR 1.040, 0.955-1.132; P=0.371); reverse MR also showed heterogeneity and MR-PRESSO outliers (Figure S2; Table S12). Bidirectional MR found no evidence that genetic liability to OA affected OC risk, or that genetic liability to OC affected OA risk, under the selected datasets, instruments, and assumptions. It does not test shared heritability, genetic correlation, or common susceptibility loci. Figure 6 integrates this boundary with the gene, pathway, and cellular evidence layers.
"""


DISCUSSION_AND_CONCLUSION = """## Discussion

### Shared genes are not shared programs

The central finding is the heterogeneity hidden within a 286-gene intersection. Only 51.0% of shared DEGs were directionally concordant, seven of ten evidence-panel genes changed oppositely, and six of ten jointly significant Hallmarks had opposite NES signs. The top-5/top-10/top-15 sensitivity preserved discordant predominance, so this conclusion was not unique to a ten-gene cutoff. At the same time, none of the small-panel Hallmark profiles survived FDR correction. Shared membership can therefore summarize parallel remodeling-associated responses without defining one conserved program.

The evidence panel was designed to make the candidate hierarchy inspectable, not to maximize an AUC. SOX9 was the only cross-disease model consensus, yet it changed in opposite bulk directions and occupied different candidate-centered pathway and cell contexts. DDIT3 and BNC1 had different disease-specific nested stability profiles, and AKAP12 was directionally concordant but mapped to different OA cartilage and OC Fibroblast contexts. These patterns explain why evidence can converge on the same gene without assigning the same biological role.

### Shared genes have different cellular meanings

Single-cell analysis moved the question from whether a gene was detected to which exact source-defined cellular context carried it. OA and OC were not forced into one latent space, and Fibroblast and Ovarian.cancer.cell labels were not relabeled into stronger biological identities. The panel-size audit reinforced this restraint: OA HomC localization was consistent, whereas the top OC label shifted with panel composition despite high overall profile correlations. Same gene therefore did not imply the same cellular meaning, and one top-scoring label should not be treated as a fixed disease mechanism.

CellChat and NicheNet were retained as supplementary boundary analyses. Their most informative result was negative: the broader shared-DEG set intersected consensus communication, but no fixed panel gene was a direct consensus ligand or receptor, and the available design did not support NicheNet ligand-activity inference. STRING and focused upstream resources likewise supplied database context rather than tissue-specific mechanism. These layers nominate experiments; they do not close a causal chain.

### Molecular separability depends on the validation task

The OA and OC AUCs should not be compared as estimates of one universal diagnostic model. The same signed summary was applied to chronic within-tissue OA contrasts and to malignant OC tumor-reference contrasts that differ in signal amplitude, comparator tissue, and cellular composition. Near-complete OC separation and modest OA separation are therefore properties of different retrospective tasks. Permutation, sample omission, Hedges g, and bounded calibration increase transparency but cannot make those tasks biologically equivalent.

Unsupervised discovery-cohort PCA and correlation audits also showed substantial internal structure. This finding argues against a simple gene-only explanation and strengthens the need for external-cohort and context-specific interpretation. It does not justify outcome-informed sample exclusion or post hoc batch correction without reliable technical covariates.

### Transcriptomic convergence does not imply genetic causality

Bidirectional MR found no evidence that genetic liability to OA affected OC risk, or that genetic liability to OC affected OA risk, under the selected datasets, instruments, and assumptions. MR was used to bound disease-to-disease inherited causality, not to test shared heritability or genetic correlation. The result does not damage the transcriptomic observations; it limits their interpretation. European ancestry predominance, hip/knee OA phenotypes, heterogeneity, and outlier sensitivity restrict generalizability, and MR cannot exclude environmentally mediated pathways.

### Parallel stress adaptation is a testable hypothesis

The observed convergence may reflect parallel adaptation to chronic stress rather than a shared disease program. Aging, inflammation, matrix remodeling, tissue injury, metabolic stress, and endocrine context could recruit overlapping transcripts in unrelated tissues while producing different directions, pathway states, and cellular consequences. This remains a hypothesis for spatial, longitudinal, protein-level, and perturbational studies, not a demonstrated mechanism.

### Strengths and limitations

Strengths include separate disease-specific discovery, two external cohorts per disease, explicit gene and pathway directions, prespecified DEG thresholds, 2,000 WGCNA bootstraps, strict outer-fold feature selection, transparent panel accounting, top-5/top-10/top-15 sensitivity, unsupervised discovery-cohort PCA and correlation audits, fixed-direction external scoring, five single-cell adapters, sample-aware CCSS/UCell/pseudobulk summaries, complete MR provenance, exact figure source data, and a claim-evidence registry. Null and unfavorable results were retained.

The study remains observational and retrospective. Tissues, platforms, processing, disease stages, comparator tissues, and clinical composition differed. Discovery-cohort PCA and correlation profiles revealed marked substructure, but complete technical covariates were unavailable for consistent blocking; residual technical and cell-composition confounding remains possible. DEG overlap was threshold dependent, while candidate-panel pathways and the exact top OC cell label were panel-size dependent. WGCNA lacked independent module preservation, only two external cohorts per disease limited meta-analysis, and the signed score was not a locked probability model. HPA/GTEx lacks articular cartilage; immune and TCGA context scores are transcriptomic proxies; TCGA survival analysis lacked independent validation. Single-cell datasets had unequal controls, annotations, and replication. CellChat was eligible in only two atlases and lacked spatial evidence; NicheNet ligand activity was not estimated; STRING and regulatory resources were database derived. MR was predominantly European and represented hip/knee OA. Protein, spatial, prospective, perturbational, and multi-cancer specificity validation were unavailable.

## Conclusions

Osteoarthritis and ovarian cancer show partial context-dependent transcriptomic convergence, but shared genes are not a shared program. The overlap is approximately half discordant, pathway states are frequently opposite, cellular localization is atlas and panel dependent, retrospective separability differs by comparator scale, and bidirectional MR supplies no evidence of a disease-to-disease inherited-liability effect under the selected assumptions. The ten-gene evidence panel is an auditable summary, not a universal biomarker, therapeutic target set, or established mechanism.
"""


SUPPLEMENTARY_INDEX = """## Supplementary table index

- **Table S1:** Datasets, cohorts, and clinical-metadata availability.
- **Table S2:** Primary shared OA-OC DEGs.
- **Table S3a-b:** DEG threshold sensitivity.
- **Table S4:** WGCNA stability.
- **Table S5:** Strict nested machine-learning resampling and feature stability.
- **Table S6a-c:** Direction-fixed external validation and unsupervised PCA source data.
- **Table S7:** TCGA-OV model sensitivity.
- **Table S8:** Immune-signature comparisons.
- **Table S9:** Single-cell QC and downstream status.
- **Table S10:** Single-cell candidate localization and eligible pseudobulk results.
- **Table S11a-c:** Hallmark, Gene Ontology, and KEGG enrichment.
- **Table S12a-b:** Bidirectional MR estimates, provenance, and diagnostics.
- **Table S13:** Supplementary interaction lookup data.
- **Table S14:** HPA normal-tissue and normal-cell context.
- **Table S15a-b:** TCGA-OV relative stromal/immune context scores and correlations.
- **Table S16:** Transparent ten-gene evidence-panel matrix.
- **Table S17:** Exploratory cell-type marker GO annotation.
- **Table S18:** Complete paired Hallmark direction matrix.
- **Table S19:** Candidate gene-cell-function context matrix.
- **Table S20:** Candidate-centered Hallmark association context.
- **Table S21:** Cross-cohort molecular-separability context.
- **Table S22a-b:** External signed-score effect sizes and cross-fitted calibration sensitivity.
- **Table S23a-b:** Focused KnockTF and miRTarBase context.
- **Table S24a-c:** Dataset/context CCSS, exact-label consensus, and sample-aware UCell summaries.
- **Table S25a-e:** STRING mapping, edges, node/subgraph topology, and label permutations.
- **Table S26a-d:** CellChat sample audit, interactions, consensus pathways, and shared-DEG anchoring.
- **Table S27:** NicheNet v2 prior-consistency overlay.
- **Table S28:** Communication-analysis feasibility decisions and inference boundaries.
- **Table S29a-d:** Evidence-panel composition, direction, Hallmark, and profile sensitivity.
- **Table S30a-c:** Extended gene-cell detection and panel-size localization sensitivity.
- **Table S31a-b:** Discovery-cohort unsupervised PCA and sample-correlation QC.
"""


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_inputs(panel_path: Path, localization_path: Path, bulk_path: Path) -> None:
    panel = read_rows(panel_path)
    if len(panel) != 15:
        raise ValueError(f"Expected 15 ranked genes, found {len(panel)}")
    expected = [
        "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1", "CFI", "DDIT3",
        "DIRAS3", "EFEMP1", "HK2", "KIT", "MYZAP", "NOD2", "OGN", "RTN1",
    ]
    if [row["gene"] for row in panel] != expected:
        raise ValueError("Panel ranking does not match the audited deterministic order")
    localization = read_rows(localization_path)
    top = {
        (row["dataset_id"], row["panel"]): row["cell_type"]
        for row in localization if row["top_cell_type"].upper() == "TRUE"
    }
    expected_top = {
        ("GSE255460", "Top 5"): "HomC",
        ("GSE255460", "Top 10"): "HomC",
        ("GSE255460", "Top 15"): "HomC",
        ("GSE154600", "Top 5"): "Ovarian.cancer.cell",
        ("GSE154600", "Top 10"): "Fibroblast",
        ("GSE154600", "Top 15"): "Fibroblast",
    }
    if top != expected_top:
        raise ValueError(f"Unexpected localization sensitivity: {top}")
    bulk = read_rows(bulk_path)
    if len(bulk) != 101:
        raise ValueError(f"Expected 101 bulk PCA samples, found {len(bulk)}")


def response_matrix() -> str:
    return """# OC-OA manuscript V3.1 revision decision matrix

The supplied Scientific Reports article was used as a visual and structural comparator, not as a template for nomograms, DCA, or hub-gene claims. Manuscript/review content was processed locally.

| Recommendation | V3.1 action | Location | Status |
|---|---|---|---|
| Center one scientific question | Reframed the study around context-dependent transcriptomic convergence and four explicit questions. | Title; Introduction; Figure 1 | Implemented |
| Reduce analysis accumulation | Main narrative now follows gene direction, pathway direction, cellular context, separability, and the genetic-liability boundary. | Figures 1-6; Results | Implemented |
| Explain why ten genes | Added deterministic top-5/top-10/top-15 sensitivity from the original ranking; no re-optimization. | Figure S15; Tables S29-S30 | Implemented |
| Report unfavorable sensitivity results | Retained zero FDR-significant small-panel Hallmarks and the OC top-label shift. | Abstract; Results; Limitations | Implemented |
| Add bulk PCA/QC | Added separate unsupervised PCA and sample-correlation audits for OA and OC discovery cohorts; no outcome-informed exclusion. | Figure S16; Table S31 | Implemented |
| Demote PPI | Renamed it protein-association landscape and moved it fully to the supplement. | Figure S12; Tables S25a-e | Implemented |
| Keep ML secondary | Retained strict nested stability and original LASSO/RF evidence in supplements; no new algorithm. | Figure S8; Tables S5/S16 | Implemented |
| Reframe AUC | Renamed the section disease-context-dependent separability; Figure 5 order is PCA, Hedges g, permutation, ROC. | Figure 5; Results/Discussion | Implemented |
| Prioritize single-cell context | CCSS/UCell/pseudobulk remain primary; exact labels retained; panel-size localization added. | Figure 4; Figure S15; Tables S24/S30 | Implemented |
| Keep CellChat/NicheNet bounded | Retained as supplementary, sample-aware, exploratory analyses; no cross-disease probability comparison or ligand-activity claim. | Figures S13-S14; Tables S26-S28 | Implemented |
| Preserve MR as a boundary | Used exact disease-to-disease genetic-liability wording; did not reinterpret MR as shared heritability. | Figure S2; Figure 6; Discussion | Implemented |
| Approach reference visual style | Adopted white backgrounds, bold panel letters, compact black hierarchy, and blue/coral disease cues while retaining colorblind-safe shape redundancy. | All rebuilt V3.1 figures | Implemented |
| Avoid low-value expansion | Added no DCA, nomogram, additional ML, drug prediction, or broad TF-miRNA network. | Full scope | Preserved |

## Author-controlled items still required

- Authors, affiliations, and corresponding-author details.
- Target journal and current journal-specific formatting.
- Funding, competing interests, and CRediT contributions.
- Public repository URL and archival DOI after deposition.
- Final human scientific, statistical, citation, language, and policy approval.
"""


def build(source: Path, legends_path: Path) -> str:
    text = source.read_text(encoding="utf-8")
    text = re.sub(r"^# .+$", TITLE, text, count=1, flags=re.MULTILINE)
    text = text.replace(
        "**Running title:** Context-dependent convergence and divergence",
        "**Running title:** Context-dependent transcriptomic convergence",
    )
    text = v22.replace_h2(text, "Abstract", ABSTRACT)
    text = v22.replace_h2(text, "Introduction", INTRODUCTION)
    text = v22.replace_h3(text, "Study design and reproducibility", STUDY_DESIGN)
    text = v22.replace_h3(text, "Bulk transcriptomic cohorts and clinical metadata", BULK_METHOD)
    text = v22.replace_h3(text, "Direction-aware protein-association context", PPI_METHOD)
    text = v22.replace_h3(text, "Candidate prioritization and strict nested machine learning", CANDIDATE_METHOD)
    text = v22.replace_h3(text, "Cross-cohort molecular separability", SEPARABILITY_METHOD)
    text = v22.replace_h3(text, "Single-cell localization, cell-context specificity, and UCell scoring", SINGLE_CELL_METHOD)
    text = v22.replace_h3(text, "Sample-consensus CellChat and bounded NicheNet prior context", COMMUNICATION_METHOD)
    text = v22.replace_h2(text, "Results", RESULTS)
    pattern = r"^## Discussion\n.*?(?=^## Data and code availability)"
    text, count = re.subn(
        pattern,
        DISCUSSION_AND_CONCLUSION.rstrip() + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one Discussion block, found {count}")
    text = text.replace("`run_submission_v30.ps1`", "`run_submission_v31.ps1`")
    text = text.replace("claim–evidence registry", "claim-evidence registry")
    text = text.replace("ten-gene set", "ten-gene evidence panel")
    text = text.replace("fixed ten-gene set", "fixed ten-gene evidence panel")
    legends = legends_path.read_text(encoding="utf-8").strip()
    legend_pattern = r"^## Figure legends\n.*?(?=^## Supplementary table index)"
    text, count = re.subn(
        legend_pattern,
        "## Figure legends\n\n" + legends + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one figure-legends block, found {count}")
    text = re.sub(r"^## Supplementary table index\n.*\Z", SUPPLEMENTARY_INDEX.rstrip() + "\n", text, flags=re.MULTILINE | re.DOTALL)
    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--figure-legends", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--localization", required=True, type=Path)
    parser.add_argument("--bulk-pca", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()
    validate_inputs(args.panel, args.localization, args.bulk_pca)
    manuscript = build(args.source, args.figure_legends)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(manuscript, encoding="utf-8")
    args.response_output.write_text(response_matrix(), encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.response_output}")


if __name__ == "__main__":
    main()
