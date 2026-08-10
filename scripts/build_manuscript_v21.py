from __future__ import annotations

import argparse
import re
from pathlib import Path


TITLE = (
    "# Context-dependent remodeling convergence between osteoarthritis and "
    "ovarian cancer without shared genetic causality: a multi-cohort "
    "transcriptomic and single-cell study"
)


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in fundamentally different tissues, yet both involve stromal, immune, metabolic, and stress-response remodeling. We tested the hypothesis that OA and OC partially converge on remodeling-associated molecular programs, but that these programs are context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality.

### Methods

We analyzed public human bulk and single-cell transcriptomic datasets with a reproducible workflow. Differential expression was estimated separately in OA cartilage (GSE114007; 20 OA and 18 non-OA samples) and ovarian tissue (GSE18520; 53 high-grade serous tumors and 10 normal samples). Direction-aware overlap was integrated with disease-specific weighted gene co-expression network analysis (WGCNA). Machine-learning robustness was tested with 50 repeated five-fold outer resamples in which screening across all measured genes was repeated inside each training fold. A discovery-oriented signed ten-gene score was examined in two external cohorts per disease using fixed ROC direction, unsupervised principal-component analysis (PCA), 1,000 label permutations, and leave-one-sample-out analyses. Human Protein Atlas (HPA), exploratory TCGA-OV survival and relative stromal/immune context scores, five single-cell datasets, and supplementary bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes (DEGs), with 286 shared. Only 146/286 (51.0%) were directionally concordant; the four OA/OC direction quadrants contained 112, 34, 86, and 54 genes. Ten genes were prioritized (SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2), but only AKAP12, JUNB, and DIRAS3 were concordant. The fixed score showed weak OA separation (AUC 0.520 and 0.629) and strong OC tumor–normal separation (AUC 1.000 and 0.979). In GSE54388, unsupervised PCA based on the 2,000 most variable genes visually separated the 6 normal from 16 tumor samples; permutation and sample-omission checks supported cohort-level molecular separation but not diagnostic utility. Relative TCGA-OV context scores showed modest candidate associations (maximum absolute Spearman ρ=0.393) and were not interpreted as histologic purity. Across five single-cell datasets, 1,025,361 QC-pass cells placed candidate expression in distinct cartilage, malignant, stromal, and immune contexts. MR using 21 OA-to-OC and 11 OC-to-OA instruments detected no causal effect in either direction.

### Conclusions

OA and OC partially converge on remodeling-associated molecular programs, but the convergence is directionally heterogeneous and cellularly context dependent. The evidence supports an auditable cross-disease molecular resource, not a clinically validated diagnostic panel, shared cellular state, therapeutic target set, or causal OA–OC relationship.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; single-cell RNA sequencing; WGCNA; nested cross-validation; directionality; reproducibility
"""


INTRODUCTION = """## Introduction

Osteoarthritis is a whole-joint disorder characterized by cartilage loss, extracellular-matrix remodeling, subchondral bone change, synovial inflammation, and altered chondrocyte states [1,3,17–19]. High-grade serous ovarian cancer is an aggressive malignancy marked by genomic instability, dissemination, and a heterogeneous stromal and immune microenvironment [2,14,20,21]. Their anatomy and clinical course are fundamentally different, so a molecular comparison requires a narrower premise than the observation that both involve inflammation.

We formulated the comparison around three intersecting axes. Aging-associated stress may connect chondrocyte senescence and matrix failure with malignant stress adaptation; mesenchymal remodeling may connect cartilage-matrix reorganization with tumor–stroma remodeling; and immune remodeling may connect inflammatory cartilage states with immune adaptation in the ovarian tumor microenvironment. Our unifying hypothesis was that OA and OC partially converge on remodeling-associated molecular programs, but that these programs are context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality.

Cross-disease transcriptomic studies commonly intersect threshold-defined DEGs and then add network, machine-learning, pathway, or ROC analyses. This can obscure three unresolved questions. First, does shared membership represent the same direction of change? Second, can bulk overlap be resolved to compatible cellular contexts rather than tissue-composition differences? Third, does transcriptomic overlap coincide with inherited causal liability? Feature screening before resampling and choosing ROC direction in each validation cohort can further inflate apparent reproducibility [6,10–13,32].

We therefore conducted an evidence-bounded multi-cohort analysis. The primary aims were to quantify direction-aware overlap, assess disease-specific co-expression evidence, repeat feature screening inside strict outer resampling, evaluate fixed-direction cross-cohort molecular separation, and localize candidates to dataset-specific cellular contexts. HPA, functional and immune annotation, exploratory TCGA-OV analyses, and supplementary bidirectional MR were used to define what the shared transcriptomic signal did—and did not—support.
"""


METHODS = """## Methods

### Study design and reproducibility

OA and OC were analyzed separately at every disease-specific stage and connected only by explicit gene-level comparisons. The workflow included bulk discovery, direction-aware overlap, disease-specific WGCNA, strict nested machine learning, external reproducibility, HPA context, functional and immune annotation, exploratory TCGA-OV analyses, single-cell localization, and supplementary bidirectional MR (Figure 1). Fixed configurations, manifests, seeds, tests, result tables, figure source data, and a claim–evidence registry are stored in the one-command project; randomized analyses used seed 20260726.

### Bulk transcriptomic cohorts and clinical metadata

Discovery cohorts were OA knee cartilage GSE114007 (20 OA, 18 non-OA) [3] and OC GSE18520 (53 advanced high-grade serous tumors, 10 normal ovarian surface-epithelium samples) [4]. External cohorts were GSE117999 (10 OA, 10 controls) and GSE82107 (10 OA, 7 controls) for OA, and GSE54388 (16 tumors, 6 normal samples) and GSE12470 (43 serous carcinomas, 10 normal peritoneal samples) for OC [5]. Expression and metadata were checked for identifier agreement, phenotype balance, duplicate features, and non-finite values; repository metadata, not expression clustering, defined groups. Uneven availability of age, sex, stage, grade, and treatment prevented consistent patient-level covariate adjustment (Table S1).

### Differential expression and direction-aware overlap

Differential expression was estimated separately with limma [6]. The primary rule was Benjamini–Hochberg FDR <0.05 and absolute log2 fold change ≥1 [7]. Shared DEGs met this rule in both diseases; concordance required equal fold-change signs. Prespecified sensitivity covered FDR 0.01/0.05 crossed with absolute log2 fold-change thresholds 0.5/1.0/1.5, recording DEG, overlap, concordance, Jaccard, and candidate-retention metrics without post hoc threshold selection.

### Functional and immune annotation

Ranked enrichment used Hallmark, canonical pathway, Gene Ontology, and oncogenic-signature collections [8,9]; Gene Ontology and KEGG over-representation used the primary shared set. Direction was interpreted separately by disease. Rank-based immune signatures were treated as expression associations, not measured cell counts.

### Co-expression analysis

Signed WGCNA networks were fitted separately [10]. Primary modules were OA green (269 genes; power 8) and OC brown (424 genes; power 16). Stability used power perturbation (−2/+2), 2,000 sample bootstraps of the fixed module, and leave-one-sample-out estimates. These checks quantify within-cohort association stability but do not replace independent preservation, particularly for the 38-sample OA cohort.

### Candidate prioritization and strict nested machine learning

The evidence synthesis prioritized SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2 and retained their signed score for transparent external assessment. Internal performance was re-estimated without restricting screening to this global set. In 50 repeated five-fold outer resamples, all measured genes were ranked by training-only two-group statistics and the top 100 entered modeling. Five-fold inner cross-validation selected LASSO lambda; random-forest mtry was chosen from 2, 5, and 10 by training-fold out-of-bag AUC before fitting 300 trees [11,12]. Outer test samples were excluded from screening, tuning, scaling, and fitting. AUC, balanced accuracy, Brier score, fit success, and selection frequency were summarized across outer predictions [32].

### Cross-cohort molecular reproducibility

External gene orientation was fixed from discovery log2 fold change. AUC and DeLong 95% confidence intervals were calculated without choosing direction in validation data [13]; the signed score averaged standardized expression multiplied by discovery signs. GSE54388 also underwent unsupervised PCA using the 2,000 genes with highest sample standard deviation, with gene-wise centering and unit-variance scaling; phenotype labels were used only to display the resulting sample coordinates. Robustness used 1,000 label permutations per cohort and sequential sample omission with score recomputation. These are retrospective molecular-separation analyses, not prospective diagnostic evaluations.

### Human Protein Atlas normal-reference context

HPA version 25.1/Ensembl 109 records captured normal RNA tissue specificity/distribution, enriched or enhanced tissues, single-cell-type specificity, and whether ovary was listed [31]. HPA integrates normal HPA/GTEx expression but lacks articular cartilage, so it cannot establish OA specificity or remove tissue-composition confounding.

### TCGA-OV exploratory analyses

TCGA-OV RNA-seq and clinical data provided exploratory context [14]. Penalized Cox regression selected AKAP12, BNC1, and DDIT3 [15,16]. Continuous risk was evaluated per standard deviation, including age/stage adjustment (n=303), 200 optimism bootstraps, 200 bootstrap LASSO fits, Schoenfeld diagnostics, and a time-varying-coefficient sensitivity model. Separately, published ESTIMATE rank-based stromal and immune signatures [33] generated relative within-cohort transcriptomic context scores. Candidate–score associations used Spearman correlation with FDR correction. Absolute tumor purity and histologic cell fractions were not inferred from RNA-seq.

### Single-cell datasets and inference

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17–21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. Candidate expression was summarized by annotated cell type using mean log-normalized expression and detection fraction. Disease inference used pseudobulk only with biological replication and an interpretable contrast [23]; other datasets contributed localization. OA and OC were not integrated into a shared latent space, and transferred labels remained conservative and dataset specific [24].

### Supplementary bidirectional MR

OpenGWAS datasets were hip or knee OA ebi-a-GCST007092 (European; n=417,596; 39,427 cases) and ovarian cancer ieu-a-1120 (European females; n=66,450; 25,509 cases), both GRCh37 [25,28,29]. Instruments used P<5×10−8, LD clumping r²<0.001 within 10,000 kb, and outcome harmonization, retaining 21 OA-to-OC and 11 OC-to-OA variants. Five estimators were reported with heterogeneity, Egger-intercept, leave-one-out, and MR-PRESSO diagnostics [26,27]. MR was independent of transcriptomic selection; a null estimate meant no detected inherited causal evidence under available instruments and assumptions.

### Statistical reporting

Tests were two sided unless prespecified otherwise. Benjamini–Hochberg correction was applied within analysis families, and 95% intervals are reported where available. Resampling intervals describe implemented resampling distributions. Null, heterogeneous, and direction-inconsistent findings were retained.
"""


RESULTS = """## Results

### The shared DEG set was directionally heterogeneous

The study included two bulk discovery cohorts, four external bulk cohorts, TCGA-OV, five single-cell datasets, HPA annotations, and two GWAS datasets (Figure 1; Table S1).

At the primary threshold, 2,008 genes were differentially expressed in OA and 2,310 in OC, with 286 shared (Figure 2; Table S2). The explicit quadrants contained 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. Thus, only 146/286 (51.0%) were concordant. Directional heterogeneity persisted across all six prespecified thresholds (Table S3), showing that shared membership was not a uniform activation program.

All ten prioritized genes met the primary rule in both diseases, but only AKAP12, JUNB, and DIRAS3 were concordant. SOX9 was lower in OA (log2 fold change −2.072) and higher in OC (+2.612), whereas EFEMP1 was higher in OA (+3.092) and lower in OC (−2.781). The panel therefore cannot be described as uniformly activated or suppressed.

### Functional annotations supported convergence without common direction

Hallmark enrichment identified inflammatory, hypoxic, metabolic, and stress-response contexts, but several directions differed between OA and OC (Figure 6; Table S11). Immune-signature associations also differed within diseases, yet could reflect abundance or within-cell transcription. Functional overlap therefore generated hypotheses rather than evidence of identical pathway activation.

### WGCNA associations were stable within the available cohorts

Primary module–trait correlations were −0.951 in OA and −0.879 in OC. Bootstrap sign stability was 1.000 in both across 2,000 resamples. Power perturbation retained absolute correlations of 0.904–0.948 in OA and 0.856–0.885 in OC; leave-one-sample-out estimates did not reverse signs (Figure 3; Figure S1; Table S4). These support stable within-cohort associations, with limited precision in the small OA network.

### Strict nested feature selection retained strong internal separation

Screening across all measured genes was repeated inside every outer training fold. Across 50 repeated five-fold resamples, median AUC was 1.000 for LASSO and random forest in both diseases; the OA LASSO 95% resampling range was 0.998–1.000 and remaining ranges were 1.000–1.000 (Table S5). This addresses leakage from global screening but does not convert a large retrospective disease–reference contrast into clinical accuracy or establish specificity of the ten-gene panel (Figure 3C).

### OC tumor–normal separation was robust but not diagnostic validation

The signed score yielded OA AUCs of 0.520 (95% CI 0.250–0.790) in GSE117999 and 0.629 (0.350–0.907) in GSE82107, with permutation P=0.469 and 0.207 and leave-one-out ranges 0.456–0.589 and 0.583–0.698. OA external reproducibility was weak and imprecise.

In OC, AUC was 1.000 in GSE54388 and 0.979 (0.944–1.000) in GSE12470; both empirical permutation P values were 0.001. Leave-one-out AUC remained 1.000 throughout GSE54388 and 0.971–0.993 in GSE12470 (Figure 4; Table S6). In an unsupervised GSE54388 PCA, PC1 and PC2 explained 33.0% and 15.1% of variance, and the 6 normal samples were visually separated from the 16 tumors without label-informed fitting. Together, these findings support a strong tumor–normal molecular contrast, while small normal groups, platform structure, comparator tissue, and retrospective preprocessing preclude a diagnostic claim.

### Normal-reference annotations did not support a uniformly tissue-specific panel

HPA classified AKAP12 and DDIT3 as low-specificity, six candidates as tissue enhanced, CFI as liver enriched, and DIRAS3 as group enriched in brain, ovary, and pituitary. Only DIRAS3 listed ovary among specific tissues (Figure S4; Table S14). Cartilage is absent from HPA, and bulk cohorts contain different cellular mixtures, so residual composition confounding remains unresolved.

### TCGA-OV provided exploratory prognostic and tissue-context information

The continuous TCGA risk score retained an age/stage-adjusted association (hazard ratio per standard deviation 1.262, 95% CI 1.076–1.479). Apparent concordance was 0.597 and optimism-corrected concordance 0.586; selection instability and proportional-hazards concerns remained (Figure 6; Figure S1; Table S7). The model is exploratory, not clinically ready.

Relative transcriptomic context scores provided an additional composition audit across 307 TCGA-OV samples. Nine of ten candidates were present in the matrix (EFEMP1 was unavailable). Associations were modest: the largest absolute correlation was BNC1 with stromal score (Spearman ρ=0.393, FDR=2.55×10−11), followed by BNC1 with the combined score (ρ=0.358) and ELF3 with immune score (ρ=0.297) (Figure S5; Table S15). These scores contextualize expression but do not measure absolute purity or histologic fractions.

### Single-cell data showed distinct cellular contexts

The five adapters audited 1,187,436 cells and retained 1,025,361 QC-pass cells (Figure 1; Table S9). The gene–disease context matrix localized candidates to cartilage chondrocyte/fibrochondrocyte states in OA and to malignant epithelial, fibroblast, endothelial, mast-cell, and other immune-associated contexts in OC (Figure 5; Figure S3; Table S10). Eligible OA pseudobulk contrasts supported cell-state-specific associations for CFI, DIRAS3, AKAP12, EFEMP1, and JUNB; OC datasets mainly contributed localization because controls or biological replication were inconsistent.

These findings do not define a conserved OA–OC cellular niche. The same bulk-level gene can carry different cellular meaning in cartilage and ovarian tumors, which helps explain disease-level directional discordance.

### Supplementary MR did not support inherited causal liability

MR used the prespecified ebi-a-GCST007092 and ieu-a-1120 datasets with 21 and 11 instruments. Inverse-variance-weighted estimates were null for OA-to-OC (OR 1.015, 95% CI 0.900–1.144; P=0.811) and OC-to-OA (OR 1.040, 0.955–1.132; P=0.371); reverse MR also showed heterogeneity and MR-PRESSO outliers (Figure S2; Table S12). Thus, transcriptomic overlap did not coincide with detected inherited genetic liability, although power and instrument assumptions limit the inference.
"""


DISCUSSION_AND_CONCLUSION = """## Discussion

Our unifying interpretation is that OA and OC partially converge on remodeling-associated molecular programs, but these programs are context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality. This framing narrows the study from a search for a universal cross-disease signature to an auditable comparison of molecular components and their boundaries.

First, shared genes are not identical changes. Only 146 of 286 shared DEGs were concordant, and seven of ten prioritized genes changed oppositely. A gene can participate in remodeling in both diseases while reflecting compensatory cartilage responses in one and malignant or stromal programs in the other. The explicit four-quadrant display therefore carries more biological information than the intersection alone.

Second, the same gene can have different cellular meaning. OA candidates localized to cartilage cell states, whereas OC expression appeared across malignant epithelial, fibroblast, endothelial, and immune-associated compartments. The cell-context matrix does not claim homologous cell states; it shows how bulk overlap can arise from distinct sources. A forced shared latent space would not establish conservation without matched controls, harmonized annotations, and prespecified state scores.

Third, transcriptomic overlap is not inherited genetic causality. Bidirectional MR was null under the available instruments, while transcriptomic and single-cell results remained descriptive associations. Environmental exposures, aging, tissue injury, or convergent remodeling can generate expression overlap without a causal disease-to-disease pathway.

The strong OC AUCs require the same restraint. Strict nested analysis addressed information leakage, and permutation, sample-omission, and unsupervised PCA checks support genuine retrospective cohort structure. They do not solve the malignant-versus-normal tissue contrast, small normal groups, platform effects, or prospective clinical utility. The correct term is reproducible molecular separation, not diagnostic performance or biomarker validation.

Tissue specificity remains a complementary constraint. HPA showed that most candidates are not uniformly ovary specific and cannot assess cartilage. Relative TCGA stromal/immune scores further showed that candidate expression is partly associated with transcriptome-derived context, but these scores are neither histologic purity nor measured cell fractions. Together with single-cell localization, this supports context-dependent interpretation rather than a disease-specific universal panel.

The practical contribution is a reproducible resource containing direction-aware shared DEGs, disease-specific network evidence, strict resampling outputs, fixed-score external checks, tissue and cell-context annotations, exact source data, and explicit negative evidence. The ten candidates are regulatory or contextual nodes for follow-up, not therapeutic targets. Spatial transcriptomics, protein assays, perturbation experiments, and prospective multi-cancer validation would be required before mechanism or clinical application could be inferred.

TCGA-OV and MR define useful outer boundaries. The survival model had modest optimism-corrected discrimination and imperfect assumptions; MR did not detect inherited causal liability. Neither weakens the direction-aware transcriptomic finding, but both prevent its promotion to a prognostic model or causal OA–OC relationship.

### Strengths

Strengths include disease-specific discovery, two external cohorts per disease, explicit directional quadrants, prespecified DEG sensitivity, 2,000 WGCNA bootstraps, strict outer-fold feature selection, unsupervised PCA, permutation and sample-influence checks, transparent HPA and TCGA context audits, five single-cell adapters, complete MR provenance, exact figure source data, and a claim–evidence registry. Null and unfavorable results were retained.

### Limitations

The study is observational and retrospective. Cohorts differed in platform, tissue, processing, disease stage, comparator tissue, and clinical composition, with incomplete patient-level metadata. DEG overlap remained threshold dependent. OA WGCNA used 38 samples and lacked independent module-preservation testing. Near-perfect internal and OC external AUCs may reflect strong tumor–normal or disease–reference contrasts and are not clinical diagnostic performance. HPA lacks cartilage; immune and ESTIMATE-derived scores are transcriptomic proxies rather than measured fractions. TCGA survival modeling was exploratory and lacked independent validation. Single-cell datasets had unequal controls, annotations, and replication, and OC results were mainly descriptive. MR power and validity depend on available instruments. Protein, spatial, perturbational, prospective, and multi-cancer specificity validation were unavailable.

## Conclusions

Osteoarthritis and ovarian cancer partially converge on remodeling-associated molecular programs, but the convergence is directionally heterogeneous, cellularly context dependent, and unsupported as shared genetic causality. The evidence provides a reproducible framework for studying shared molecular components without equating shared membership with identical direction, cell of origin, diagnostic utility, prognosis, therapeutic actionability, or causal disease linkage.
"""


FIGURE_AND_TABLE_TEXT = """## Figure legends

### Figure 1. Study design, biological hypothesis, and audited data resources

**A,** The study tests context-dependent convergence of aging-associated, mesenchymal-remodeling, and immune-remodeling programs while keeping OA and OC analyses separate. Bidirectional MR is supplementary null evidence. **B,** Bulk discovery and external-cohort sample counts. **C,** Total and QC-pass single-cell counts.

### Figure 2. Direction-aware cross-disease transcriptomic discovery

**A–B,** OA and OC differential-expression volcano plots. **C,** Explicit quadrants for disease-specific effects among commonly measured genes. The 286 primary shared DEGs are colored by direction: 112 were higher in both, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. **D,** Shared-gene counts under six prespecified thresholds.

### Figure 3. Robust identification of molecular candidates

**A,** WGCNA module–trait association under soft-power perturbation. **B,** Primary-module gene retention. **C,** Feature frequency under strict outer-fold screening and nested tuning. **D,** Disease-specific effects for the ten prioritized genes.

### Figure 4. Direction-fixed cross-cohort molecular reproducibility

**A,** Unsupervised PCA of GSE54388 using the 2,000 genes with highest sample standard deviation; phenotype labels were used only for display. **B,** ROC curves for the signed ten-gene score with direction fixed from discovery. **C,** Null AUC distributions from 1,000 label permutations. **D,** Leave-one-sample-out AUCs. These panels show retrospective molecular separation, not clinical diagnostic performance.

### Figure 5. Candidate-gene localization across distinct cellular contexts

**A,** QC-pass fractions. **B,** Cell-count-weighted candidate detection. **C,** Gene–disease matrix showing the annotated cell context with the highest weighted detection for each candidate within OA and OC atlases. Labels remain dataset specific and do not imply homologous cell states. **D,** Significant eligible pseudobulk effects. OA and OC datasets were not integrated into a shared latent space.

### Figure 6. Functional and exploratory prognostic context

**A,** Hallmark enrichment. **B,** Rank-based immune-signature differences. **C,** Exploratory TCGA-OV Cox estimates. **D,** Bootstrap concordance and optimism correction.

### Supplementary Figure 1. Prespecified non-MR sensitivity analyses

DEG threshold retention, WGCNA leave-one-out stability, strict nested feature frequency, and TCGA LASSO selection stability.

### Supplementary Figure 2. Negative bidirectional MR results

Five MR estimators and heterogeneity/pleiotropy diagnostics in both directions. Null estimates mean no causal effect was detected under the available instruments and assumptions, not proof of absence.

### Supplementary Figure 3. Dataset-specific single-cell embeddings

Released or recomputed embeddings with dataset-specific labels. Quantitative summaries used all QC-pass cells; OA and OC were not integrated into a shared latent space.

### Supplementary Figure 4. HPA normal-tissue context

Normal-tissue specificity and distribution categories. Cartilage is absent, so this audit cannot establish OA–OC specificity.

### Supplementary Figure 5. TCGA-OV relative stromal and immune context audit

**A,** Spearman correlations between candidate expression and published ESTIMATE rank-based stromal, immune, and combined scores. **B–C,** Representative strongest absolute stromal and immune associations. Scores are relative within-cohort transcriptomic proxies; absolute tumor purity and histologic cell fractions were not inferred from RNA-seq.

## Supplementary table index

- **Table S1:** Datasets, cohorts, and clinical-metadata availability.
- **Table S2:** Primary shared OA–OC DEGs.
- **Table S3a–b:** DEG threshold sensitivity.
- **Table S4:** WGCNA stability.
- **Table S5:** Strict nested machine-learning resampling and feature stability.
- **Table S6:** Fixed-direction external AUC, permutation, leave-one-out, sample scores, and PCA source data.
- **Table S7:** TCGA-OV model sensitivity.
- **Table S8:** Immune-signature comparisons.
- **Table S9:** Single-cell QC and downstream status.
- **Table S10:** Single-cell candidate localization and eligible pseudobulk results.
- **Table S11a–c:** Hallmark, Gene Ontology, and KEGG enrichment.
- **Table S12a–b:** Bidirectional MR estimates, dataset provenance, instrument parameters, and sensitivity diagnostics.
- **Table S13:** Regulatory and compound interaction lookups retained as supplementary lookup data only.
- **Table S14:** HPA normal-tissue and normal-cell context.
- **Table S15a–b:** TCGA-OV relative stromal/immune context scores and candidate correlations.
"""


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(
        rf"(?ms)^{re.escape(start)}\n.*?(?=^{re.escape(end)}\n)"
    )
    updated, count = pattern.subn(replacement.rstrip() + "\n\n", text)
    if count != 1:
        raise RuntimeError(
            f"Expected exactly one section {start!r} before {end!r}; found {count}."
        )
    return updated


def word_count(section: str) -> int:
    return len(re.findall(r"\b[\w×≥ρ²−–]+\b", section, flags=re.UNICODE))


def build(source: Path, destination: Path) -> None:
    text = source.read_text(encoding="utf-8")
    original_methods = re.search(
        r"(?ms)^## Methods\n.*?(?=^## Results\n)", text
    )
    if original_methods is None:
        raise RuntimeError("Could not locate original Methods section.")

    text = re.sub(r"(?m)^# .+$", TITLE, text, count=1)
    text = replace_section(text, "## Abstract", "## Introduction", ABSTRACT)
    text = replace_section(text, "## Introduction", "## Methods", INTRODUCTION)
    text = replace_section(text, "## Methods", "## Results", METHODS)
    text = replace_section(text, "## Results", "## Discussion", RESULTS)
    text = replace_section(
        text,
        "## Discussion",
        "## Data and code availability",
        DISCUSSION_AND_CONCLUSION,
    )
    text = text.replace(
        "Manuscript-facing outputs can be rebuilt with one command.",
        "Manuscript-facing outputs can be rebuilt with one command.",
    )
    text = text.replace(
        "The local reproducible project contains one-command runners, fixed "
        "configurations, environment and input manifests, processing and test "
        "code, logs, result tables, exact figure source data, and the "
        "claim–evidence registry.",
        "The local reproducible project contains one-command runners, including "
        "`run_submission_v21.ps1`, fixed configurations, environment and input "
        "manifests, processing and test code, logs, result tables, exact figure "
        "source data, and the claim–evidence registry.",
    )
    reference = (
        "33. Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour "
        "purity and stromal and immune cell admixture from expression data. "
        "*Nat Commun*. 2013;4:2612. doi:10.1038/ncomms3612."
    )
    text = text.replace(
        "32. Varma S, Simon R. Bias in error estimation when using "
        "cross-validation for model selection. *BMC Bioinformatics*. "
        "2006;7:91. doi:10.1186/1471-2105-7-91.\n",
        "32. Varma S, Simon R. Bias in error estimation when using "
        "cross-validation for model selection. *BMC Bioinformatics*. "
        "2006;7:91. doi:10.1186/1471-2105-7-91.\n"
        + reference
        + "\n",
    )
    text = re.sub(
        r"(?ms)^## Figure legends\n.*\Z",
        FIGURE_AND_TABLE_TEXT.rstrip() + "\n",
        text,
        count=1,
    )
    if "GCST90018888" in text or "GCST90038686" in text:
        raise RuntimeError("Incorrect simulated-review GWAS identifiers entered the manuscript.")
    for required in (
        "ebi-a-GCST007092",
        "ieu-a-1120",
        "P<5×10−8",
        "r²<0.001",
        "10,000 kb",
        "21 OA-to-OC",
        "11 OC-to-OA",
        "Figure S5",
        "Table S15",
    ):
        if required not in text:
            raise RuntimeError(f"Required V2.1 item is missing: {required}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8", newline="\n")

    original_words = word_count(original_methods.group(0))
    revised_words = word_count(METHODS)
    reduction = 100 * (original_words - revised_words) / original_words
    print(
        f"Methods words: {original_words} -> {revised_words} "
        f"({reduction:.1f}% reduction)"
    )
    print(f"Wrote {destination}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("manuscript/OC_OA_manuscript_revision_v2.md"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("manuscript/OC_OA_manuscript_revision_v21.md"),
    )
    args = parser.parse_args()
    build(args.source, args.output)


if __name__ == "__main__":
    main()
