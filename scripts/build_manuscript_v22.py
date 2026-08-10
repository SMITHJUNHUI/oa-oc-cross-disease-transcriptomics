from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TITLE = (
    "# Context-dependent remodeling-associated molecular convergence between "
    "osteoarthritis and ovarian cancer without shared genetic causality: a "
    "multi-cohort transcriptomic and single-cell study"
)


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in fundamentally different tissues, yet both involve stromal, immune, metabolic, and stress-response remodeling. We tested the hypothesis that OA and OC partially converge on remodeling-associated molecular programs, but that these programs are context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality.

### Methods

We analyzed public human bulk and single-cell transcriptomic datasets with a reproducible workflow. Differential expression was estimated separately in OA cartilage (GSE114007; 20 OA and 18 non-OA samples) and ovarian tissue (GSE18520; 53 high-grade serous tumors and 10 normal samples). Direction-aware overlap was integrated with disease-specific weighted gene co-expression network analysis (WGCNA). Candidate prioritization was enumerated gene by gene, while machine-learning robustness was evaluated separately with 50 repeated five-fold outer resamples in which screening across all measured genes was repeated inside each training fold. A discovery-oriented signed ten-gene score was examined in two external cohorts per disease using fixed ROC direction, unsupervised principal-component analysis (PCA), 1,000 label permutations, and leave-one-sample-out analyses. Human Protein Atlas (HPA), exploratory TCGA-OV survival and relative stromal/immune context scores, five single-cell datasets, descriptive cell-type functional annotation, and supplementary bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes (DEGs), with 286 shared. Only 146/286 (51.0%) were directionally concordant; the four OA/OC direction quadrants contained 112, 34, 86, and 54 genes. Ten genes were prioritized (SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2), but only AKAP12, JUNB, and DIRAS3 were concordant. The fixed score showed weak OA separation (AUC 0.520 and 0.629) and near-complete retrospective OC tumor-normal molecular separation (AUC 1.000 and 0.979). In GSE54388, unsupervised PCA visually separated the 6 normal from 16 tumor samples; permutation and sample-omission checks supported cohort-level molecular separation but not diagnostic utility. Across five single-cell datasets, 1,025,361 QC-pass cells placed candidate expression in distinct cartilage, malignant, stromal, and immune contexts; descriptive cluster-marker annotation reinforced cell-type-specific interpretation rather than a conserved mechanism. Relative TCGA-OV context associations were modest, and MR using 21 OA-to-OC and 11 OC-to-OA instruments detected no causal effect in either direction.

### Conclusions

OA and OC show partial remodeling-associated molecular convergence, but the overlap is directionally heterogeneous and cellularly context dependent. The evidence supports an auditable cross-disease molecular resource, not a clinically validated diagnostic panel, shared cellular state, therapeutic target set, biological remodeling mechanism, or causal OA-OC relationship.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; single-cell RNA sequencing; WGCNA; nested cross-validation; directionality; reproducibility
"""


INTRODUCTION = """## Introduction

Osteoarthritis is a whole-joint disorder characterized by cartilage loss, extracellular-matrix remodeling, subchondral bone change, synovial inflammation, and altered chondrocyte states [1,3,17-19]. High-grade serous ovarian cancer is an aggressive malignancy marked by genomic instability, dissemination, and a heterogeneous stromal and immune microenvironment [2,14,20,21]. Their anatomy and clinical course are fundamentally different, so a molecular comparison requires a narrower premise than the observation that both involve inflammation.

We formulated the comparison around three intersecting axes. Aging-associated stress may connect chondrocyte senescence and matrix failure with malignant stress adaptation; mesenchymal remodeling may connect cartilage-matrix reorganization with tumor-stroma remodeling; and immune remodeling may connect inflammatory cartilage states with immune adaptation in the ovarian tumor microenvironment. Our unifying hypothesis was that OA and OC partially converge on remodeling-associated molecular programs, but that these programs are context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality.

Past cross-disease transcriptomic studies have often treated DEG overlap as the endpoint. That approach leaves three gaps directly relevant to interpretation: shared membership can ignore opposite effect directions; bulk overlap does not identify the contributing cell type; and expression similarity does not establish inherited causal liability. Feature screening before resampling and choosing ROC direction in each validation cohort can further inflate apparent reproducibility [6,10-13,32]. Thus, a convincing cross-disease analysis must separate membership, direction, cellular origin, predictive behavior, and genetic causality instead of treating them as interchangeable evidence.

We therefore conducted an evidence-bounded multi-cohort analysis. The primary aims were to quantify direction-aware overlap, assess disease-specific co-expression evidence, make candidate prioritization transparent, repeat feature screening inside strict outer resampling, evaluate fixed-direction cross-cohort molecular separation, and localize candidates to dataset-specific cellular and functional contexts. HPA, functional and immune annotation, exploratory TCGA-OV analyses, and supplementary bidirectional MR were used to define what the shared transcriptomic signal did-and did not-support.
"""


CANDIDATE_METHOD = """### Candidate prioritization and strict nested machine learning

Candidate construction and model-performance estimation were treated as separate operations. Shared primary DEGs were ranked by a combined evidence score comprising adjusted significance and absolute effect size in both diseases. When at least five shared genes occurred in the disease-specific primary WGCNA modules, the candidate space was restricted to that module union. Disease-specific LASSO and random forest models then provided consensus and model-vote evidence. SOX9 was the cross-disease consensus; because this was below the prespecified minimum, the set was completed to ten genes by ranked multi-model votes. Table S16 reports every tier, disease-specific module membership, original model votes, strict nested selection frequencies, and single-cell support. The resulting signed score was retained only for transparent external assessment.

Internal molecular separation was re-estimated independently of that global candidate set. In 50 repeated five-fold outer resamples, all measured genes were ranked by training-only two-group statistics and the top 100 entered modeling. Five-fold inner cross-validation selected LASSO lambda; random-forest mtry was chosen from 2, 5, and 10 by training-fold out-of-bag AUC before fitting 300 trees [11,12]. Outer test samples were excluded from screening, tuning, scaling, and fitting. Balanced accuracy, Brier score, AUC, fit success, and feature-selection frequency were summarized across outer predictions [32]; AUC was retained in Table S5 rather than treated as the primary biological result.
"""


HPA_METHOD = """### Human Protein Atlas normal-reference context

HPA version 25.1/Ensembl 109 records captured normal RNA tissue specificity/distribution, enriched or enhanced tissues, single-cell-type specificity, and whether ovary was listed [31]. The HPA tissue resource integrates normal HPA and GTEx expression, but neither resource provides articular cartilage and normal ovary does not reproduce the ovarian-tumor microenvironment. These references therefore contextualize tissue specificity but cannot establish OA specificity, tumor specificity, or remove tissue-composition confounding.
"""


SINGLE_CELL_METHOD = """### Single-cell datasets and inference

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17-21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. Candidate expression was summarized by annotated cell type using mean log-normalized expression and detection fraction. Disease inference used pseudobulk only with biological replication and an interpretable contrast [23]; other datasets contributed localization. OA and OC were not integrated into a shared latent space, and transferred labels remained conservative and dataset specific [24].

For exploratory functional localization, the top 25 mean-expression contrast genes per analysis cluster in GSE104782 and GSE154600 were mapped to the majority dataset-specific cell label, pooled within cell type, and tested for Gene Ontology Biological Process over-representation against the annotated human-symbol universe. Benjamini-Hochberg correction was applied separately within each cell-type analysis. These cluster-marker ranks are descriptive rather than inferential differential-expression tests; their enrichment localizes functional themes but does not demonstrate conserved function or mechanism.
"""


MR_METHOD = """### Supplementary bidirectional MR

OpenGWAS datasets were hip or knee OA ebi-a-GCST007092 (European; n=417,596; 39,427 cases) and ovarian cancer ieu-a-1120 (European females; n=66,450; 25,509 cases), both GRCh37 [25,28,29]. A provenance audit matched these accessions across configuration, OpenGWAS metadata, result filenames, harmonized tables, and combined estimates. Instruments used P<5x10-8, LD clumping r2<0.001 within 10,000 kb, and outcome harmonization, retaining 21 OA-to-OC and 11 OC-to-OA variants. Five estimators were reported with heterogeneity, Egger-intercept, leave-one-out, and MR-PRESSO diagnostics [26,27]. MR was independent of transcriptomic selection; a null estimate meant no detected inherited causal evidence under available instruments and assumptions.
"""


RESULTS = """## Results

### The shared DEG set was directionally heterogeneous

The study included two bulk discovery cohorts, four external bulk cohorts, TCGA-OV, five single-cell datasets, HPA annotations, and two GWAS datasets (Figure 1; Table S1).

At the primary threshold, 2,008 genes were differentially expressed in OA and 2,310 in OC, with 286 shared (Figure 2; Table S2). The explicit quadrants contained 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. Thus, only 146/286 (51.0%) were concordant. Directional heterogeneity persisted across all six prespecified thresholds (Table S3), showing that shared membership was not a uniform activation program.

All ten prioritized genes met the primary rule in both diseases, but only AKAP12, JUNB, and DIRAS3 were concordant. SOX9 was lower in OA (log2 fold change -2.072) and higher in OC (+2.612), whereas EFEMP1 was higher in OA (+3.092) and lower in OC (-2.781). The candidate set therefore cannot be described as uniformly activated or suppressed.

### Functional annotations supported convergence without common direction

Hallmark enrichment identified inflammatory, hypoxic, metabolic, and stress-response contexts, but several directions differed between OA and OC (Figure 6; Table S11). Immune-signature associations also differed within diseases, yet could reflect abundance or within-cell transcription. Functional overlap therefore generated hypotheses rather than evidence of identical pathway activation.

### WGCNA associations were stable within the available cohorts

Primary module-trait correlations were -0.951 in OA and -0.879 in OC. Bootstrap sign stability was 1.000 in both across 2,000 resamples. Power perturbation retained absolute correlations of 0.904-0.948 in OA and 0.856-0.885 in OC; leave-one-sample-out estimates did not reverse signs (Figure 3; Figure S1; Table S4). These support stable within-cohort associations, with limited precision in the small OA network.

### Single-cell data resolved distinct cellular and functional contexts

The five adapters audited 1,187,436 cells and retained 1,025,361 QC-pass cells (Figure 1; Table S9). The gene-disease context matrix localized candidates to cartilage chondrocyte/fibrochondrocyte states in OA and to malignant epithelial, fibroblast, endothelial, mast-cell, and other immune-associated contexts in OC (Figure 5; Figure S3; Table S10). Eligible OA pseudobulk contrasts supported cell-state-specific associations for CFI, DIRAS3, AKAP12, EFEMP1, and JUNB; OC datasets mainly contributed localization because controls or biological replication were inconsistent.

Exploratory annotation of descriptive top cluster markers further separated functional contexts (Figure S6; Table S17). OA hypertrophic-chondrocyte markers were associated with bone mineralization and biomineral tissue development, and prehypertrophic-chondrocyte markers with extracellular-matrix organization. In OC, fibroblast markers were associated with collagen and extracellular-matrix organization, myeloid markers with antigen presentation, and malignant-cell markers with cytoplasmic translation and oxidative phosphorylation. These results do not test conserved OA-OC function; they show that superficially shared bulk genes occur within different cell-type programs.

### Candidate hierarchy and nested resampling emphasized feature stability

The prioritization matrix shows why each of the ten genes was retained (Table S16). SOX9 was the only cross-disease model consensus, while the remaining genes entered through the prespecified ranked-vote completion; five candidates belonged to the OA primary WGCNA module and five to the OC primary module. Strict outer-fold frequencies were heterogeneous rather than uniformly high: among the final genes, the maximum OA frequency ranged from 0 to 0.940 and the maximum OC frequency from 0 to 0.200 (Figure 3C; Table S5).

Screening across all measured genes was repeated inside every outer training fold. Median balanced accuracy was 1.000 for three disease-model combinations and 0.991 for the OC random forest; median Brier scores ranged from 0.0005 to 0.0198. These near-ceiling within-cohort contrasts address leakage from global screening but do not establish specificity of the ten-gene score, diagnostic accuracy, or prospective transportability. AUC distributions are reported in Table S5 as secondary model metrics.

### OC showed near-complete retrospective molecular separation, not diagnostic validation

The signed score yielded OA AUCs of 0.520 (95% CI 0.250-0.790) in GSE117999 and 0.629 (0.350-0.907) in GSE82107, with permutation P=0.469 and 0.207 and leave-one-out ranges 0.456-0.589 and 0.583-0.698. OA external reproducibility was weak and imprecise.

In OC, AUC was 1.000 in GSE54388 and 0.979 (0.944-1.000) in GSE12470; both empirical permutation P values were 0.001. Leave-one-out AUC remained 1.000 throughout GSE54388 and 0.971-0.993 in GSE12470 (Figure 4; Table S6). In an unsupervised GSE54388 PCA, PC1 and PC2 explained 33.0% and 15.1% of variance, and the 6 normal samples were visually separated from the 16 tumors without label-informed fitting. Together, these findings support near-complete retrospective tumor-normal molecular separation, while small normal groups, platform structure, comparator tissue, and retrospective preprocessing preclude a diagnostic claim.

### Normal-reference annotations did not support a uniformly tissue-specific panel

HPA classified AKAP12 and DDIT3 as low-specificity, six candidates as tissue enhanced, CFI as liver enriched, and DIRAS3 as group enriched in brain, ovary, and pituitary. Only DIRAS3 listed ovary among specific tissues (Figure S4; Table S14). GTEx-integrated HPA lacks articular cartilage, and normal ovary cannot reproduce tumor composition; residual tissue-composition confounding therefore remains unresolved.

### TCGA-OV showed an exploratory survival association and modest context correlations

The continuous TCGA risk score retained an age/stage-adjusted association (hazard ratio per standard deviation 1.262, 95% CI 1.076-1.479), but apparent concordance was 0.597 and optimism-corrected concordance 0.586; selection instability and proportional-hazards concerns remained (Figure 6; Figure S1; Table S7). This was an exploratory association, not a prognostic model.

Relative transcriptomic context scores provided a composition audit across 307 TCGA-OV samples. Nine of ten candidates were present in the matrix (EFEMP1 was unavailable). Associations were modest: the largest absolute correlation was BNC1 with stromal score (Spearman rho=0.393, FDR=2.55x10-11), followed by BNC1 with the combined score (rho=0.358) and ELF3 with immune score (rho=0.297) (Figure S5; Table S15). These scores contextualize expression but do not measure absolute purity or histologic fractions.

### Supplementary MR did not support inherited causal liability

MR used the verified ebi-a-GCST007092 and ieu-a-1120 datasets with 21 and 11 instruments. Inverse-variance-weighted estimates were null for OA-to-OC (OR 1.015, 95% CI 0.900-1.144; P=0.811) and OC-to-OA (OR 1.040, 0.955-1.132; P=0.371); reverse MR also showed heterogeneity and MR-PRESSO outliers (Figure S2; Table S12). Thus, transcriptomic overlap did not coincide with detected inherited genetic liability, although power and instrument assumptions limit the inference.
"""


DISCUSSION_TO_CONCLUSION = """## Discussion

Our unifying interpretation is that OA and OC show partial remodeling-associated molecular convergence, but the overlap is context dependent, directionally heterogeneous, and shaped by distinct cellular environments rather than shared genetic causality. This framing deliberately distinguishes transcriptomic association from biological remodeling mechanism. The observed convergence may reflect parallel adaptation to chronic stress rather than a shared disease program.

First, shared genes are not identical changes. Only 146 of 286 shared DEGs were concordant, and seven of ten prioritized genes changed oppositely. A gene can participate in remodeling-associated expression in both diseases while reflecting compensatory cartilage responses in one and malignant or stromal programs in the other. The explicit four-quadrant display therefore carries more biological information than the intersection alone.

Second, the same gene can have different cellular meaning. OA candidates localized to cartilage cell states, whereas OC expression appeared across malignant epithelial, fibroblast, endothelial, and immune-associated compartments. Exploratory cluster-marker annotation added functional context-ECM and mineralization themes in OA versus stromal, antigen-presentation, translational, and metabolic themes in OC-without demonstrating conserved function. A forced shared latent space or descriptive enrichment cannot establish mechanistic conservation.

Third, transcriptomic overlap is not inherited genetic causality. Bidirectional MR was null under the available instruments, while transcriptomic and single-cell results remained descriptive associations. Environmental exposures, aging, tissue injury, or parallel stress adaptation can generate expression overlap without a causal disease-to-disease pathway.

The machine-learning results require similar restraint. Transparent candidate accounting showed a single cross-disease model consensus and heterogeneous strict nested feature frequencies. Near-ceiling internal resampling metrics and near-complete external OC tumor-normal separation can arise from large retrospective disease-reference contrasts, small normal groups, comparator-tissue differences, and platform structure. They support cohort-level molecular separation, not a diagnostic panel.

Tissue and TCGA analyses define additional boundaries. HPA/GTEx context does not cover articular cartilage, and most candidates were not uniformly ovary specific. The TCGA survival association had modest optimism-corrected discrimination, while relative stromal/immune scores were transcriptomic proxies rather than histologic fractions. Neither supports prognostic deployment or absolute composition inference.

The practical contribution is a reproducible resource containing direction-aware shared DEGs, disease-specific network evidence, a gene-level candidate hierarchy, strict resampling outputs, fixed-score external checks, tissue and cell-context annotations, exact source data, and explicit negative evidence. The ten candidates are regulatory or contextual nodes for follow-up, not therapeutic targets. Protein assays, spatial transcriptomics, perturbation experiments, prospective multi-cancer validation, and matched cell-state comparisons would be required before mechanism or clinical application could be inferred.

### Strengths

Strengths include disease-specific discovery, two external cohorts per disease, explicit directional quadrants, prespecified DEG sensitivity, 2,000 WGCNA bootstraps, strict outer-fold feature selection, a transparent candidate-prioritization matrix, unsupervised PCA, permutation and sample-influence checks, HPA and TCGA context audits, five single-cell adapters, exploratory cell-type functional annotation with explicit boundaries, complete MR provenance, exact figure source data, and a claim-evidence registry. Null and unfavorable results were retained.

### Limitations

The study is observational and retrospective. Cohorts differed in platform, tissue, processing, disease stage, comparator tissue, and clinical composition, with incomplete patient-level metadata. DEG overlap remained threshold dependent. OA WGCNA used 38 samples and lacked independent module-preservation testing. Candidate completion by ranked model votes and heterogeneous strict nested selection frequencies limit claims of a stable panel. Near-perfect internal and OC external AUCs may reflect strong tumor-normal or disease-reference contrasts and are not clinical diagnostic performance. HPA/GTEx lacks cartilage; immune and ESTIMATE-derived scores are transcriptomic proxies rather than measured fractions. TCGA survival analysis lacked independent validation. Single-cell datasets had unequal controls, annotations, and replication; cluster-marker enrichment was descriptive, and OC results were mainly localization. MR power and validity depend on available instruments. Protein, spatial, perturbational, prospective, and multi-cancer specificity validation were unavailable.

## Conclusions

Osteoarthritis and ovarian cancer show partial remodeling-associated molecular convergence, but the overlap is directionally heterogeneous, cellularly context dependent, and unsupported as shared genetic causality. The evidence provides a reproducible framework for studying shared molecular components without equating transcriptomic association with identical direction, conserved cell function, biological mechanism, diagnostic utility, prognosis, therapeutic actionability, or causal disease linkage.
"""


FIGURE_S6 = """### Supplementary Figure 6. Exploratory cell-type functional annotation

**A,** Gene Ontology Biological Process over-representation among top cluster markers from the GSE104782 OA cartilage atlas. **B,** Corresponding annotation in GSE154600 ovarian tumors. Clusters were mapped to their majority dataset-specific cell label; the top two terms per label are shown. The dashed line marks BH FDR 0.05. Cluster-marker ranks are descriptive, so these panels localize functional themes but do not demonstrate conserved function or mechanism.
"""


def replace_h2(text: str, heading: str, replacement: str) -> str:
    pattern = rf"^## {re.escape(heading)}\n.*?(?=^## |\Z)"
    updated, count = re.subn(
        pattern,
        replacement.rstrip() + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one H2 section {heading!r}, found {count}")
    return updated


def replace_h3(text: str, heading: str, replacement: str) -> str:
    pattern = rf"^### {re.escape(heading)}\n.*?(?=^### |^## |\Z)"
    updated, count = re.subn(
        pattern,
        replacement.rstrip() + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one H3 section {heading!r}, found {count}")
    return updated


def validate_candidate_table(path: Path) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 10:
        raise ValueError(f"Expected ten candidates, found {len(rows)}")
    if sum(
        row["selection_stage"] == "cross-disease model consensus"
        for row in rows
    ) != 1:
        raise ValueError("Expected one cross-disease model consensus candidate")


def build_manuscript(source: Path, candidate_table: Path) -> str:
    validate_candidate_table(candidate_table)
    text = source.read_text(encoding="utf-8")
    text = re.sub(r"^# .+$", TITLE, text, count=1, flags=re.MULTILINE)
    text = replace_h2(text, "Abstract", ABSTRACT)
    text = replace_h2(text, "Introduction", INTRODUCTION)
    text = replace_h3(
        text,
        "Candidate prioritization and strict nested machine learning",
        CANDIDATE_METHOD,
    )
    text = replace_h3(
        text,
        "Human Protein Atlas normal-reference context",
        HPA_METHOD,
    )
    text = replace_h3(
        text,
        "Single-cell datasets and inference",
        SINGLE_CELL_METHOD,
    )
    text = replace_h3(text, "Supplementary bidirectional MR", MR_METHOD)
    text = replace_h2(text, "Results", RESULTS)

    discussion_pattern = r"^## Discussion\n.*?(?=^## Data and code availability)"
    text, count = re.subn(
        discussion_pattern,
        DISCUSSION_TO_CONCLUSION.rstrip() + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one Discussion-to-Conclusion block, found {count}")

    text = text.replace(
        "`run_submission_v21.ps1`",
        "`run_submission_v22.ps1`",
    )
    text = text.replace(
        "### Figure 6. Functional and exploratory prognostic context",
        "### Figure 6. Functional context and exploratory survival association",
    )
    text = text.replace(
        "**C,** Exploratory TCGA-OV Cox estimates for the three selected genes and the age/stage-adjusted continuous risk score. **D,** Bootstrap concordance and optimism correction.",
        "**C,** Exploratory TCGA-OV survival associations for the three selected genes and the age/stage-adjusted continuous score. **D,** Bootstrap concordance and optimism correction; these panels do not constitute a prognostic model.",
    )

    table_index = "## Supplementary table index"
    if FIGURE_S6.splitlines()[0] not in text:
        text = text.replace(
            table_index,
            FIGURE_S6.rstrip() + "\n\n" + table_index,
        )
    text = text.replace(
        "- **Table S15a–b:** TCGA-OV relative stromal/immune context scores and candidate correlations.",
        "\n".join(
            [
                "- **Table S15a–b:** TCGA-OV relative stromal/immune context scores and candidate correlations.",
                "- **Table S16:** Candidate prioritization matrix spanning shared-DEG, WGCNA, model-vote, strict nested frequency, and single-cell evidence.",
                "- **Table S17:** Exploratory GO Biological Process annotation of descriptive top cluster markers by dataset-specific cell type.",
            ]
        ),
    )
    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def response_matrix() -> str:
    return """# OC-OA manuscript V2.2: response and change matrix

This document is a working revision record for author verification. The simulated external-review text was processed locally. No unpublished manuscript or review content was sent to an external service.

| Reviewer request | V2.2 response | Manuscript/output location | Resolution |
|---|---|---|---|
| Narrow the meaning of “remodeling convergence” | Added “associated” to the title and consistently distinguished transcriptomic association from biological remodeling mechanism. | Title; Abstract; Discussion; Conclusions | Implemented |
| Explain why the ten genes were chosen | Added a gene-level hierarchy covering shared DEGs, disease-specific primary WGCNA modules, original model consensus/votes, strict nested feature frequency, single-cell context, and eligible pseudobulk evidence. The text now states that SOX9 was the sole cross-disease consensus and that the remainder followed the prespecified ranked-vote completion. | Methods; Results; Table S16 | Implemented |
| Reduce the prominence of AUC=1 in nested ML | Removed “median AUC was 1.000” from the main Results narrative. Main text now emphasizes heterogeneous selection frequency, balanced accuracy, and Brier score; AUC remains a secondary metric in Table S5. | Abstract; ML Methods/Results; Table S5; Figure 3C | Implemented |
| Avoid implying a prognostic model from TCGA | Renamed the Figure 6 context to “exploratory survival association,” shortened the Discussion, and explicitly stated that the analysis is not a prognostic model. | TCGA Results; Discussion; Figure 6 legend | Implemented |
| Add cell-type functional context without claiming mechanism | Added exploratory GO Biological Process annotation of descriptive top cluster markers from GSE104782 and GSE154600. The method and legend state that this localizes functional themes but does not demonstrate conserved function or mechanism. | Single-cell Methods/Results; Figure S6; Table S17 | Implemented with bounded inference |
| Revise Abstract wording for OC separation | Replaced “strong OC tumor-normal separation” with “near-complete retrospective OC tumor-normal molecular separation” and retained the non-diagnostic boundary. | Abstract; external-reproducibility Results | Implemented |
| Strengthen the Introduction research gap | Made the three gaps explicit: direction ignored, cell origin unresolved, and genetic causality unknown; also separated candidate selection from predictive behavior. | Introduction, paragraph 3 | Implemented |
| Verify MR accession provenance | Re-audited configuration, local OpenGWAS metadata, filenames, harmonized tables, and combined estimates. All consistently use OA `ebi-a-GCST007092` and OC `ieu-a-1120`, with 21 and 11 instruments. Unsupported alternative accessions were not introduced. | MR Methods/Results; Table S12a-b; audit report | Verified; reviewer’s hypothetical replacement not adopted |
| Reorder Results so biology precedes ML | Moved single-cell localization and functional context before the nested-resampling and external-score sections. | Results section order | Implemented |
| Add the chronic-stress alternative explanation | Added: “The observed convergence may reflect parallel adaptation to chronic stress rather than a shared disease program.” | Discussion, paragraph 1 | Implemented |
| Clarify GTEx/HPA limitations | Stated that the HPA tissue resource integrates HPA/GTEx normal expression, lacks articular cartilage, and cannot reproduce the ovarian-tumor microenvironment. | HPA Methods/Results; Limitations | Implemented |
| Avoid unnecessary virtual knockout, PPI, or drug-prediction expansion | No such analyses were added; existing interaction lookups remain supplementary and hypothesis-generating only. | Table S13 boundary; Discussion | Preserved scope |

## Required author completion before submission

- Add author names, affiliations, and corresponding-author details.
- Select the target journal and apply its current author instructions.
- Add funding, competing interests, and the final CRediT statement.
- Add the public GitHub URL and archival DOI after deposition.
- Confirm institutional requirements for secondary analysis of public data.
- Perform final human scientific, statistical, citation, and journal-policy approval.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()

    manuscript = build_manuscript(
        args.source.resolve(),
        args.candidate_table.resolve(),
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(manuscript, encoding="utf-8")

    response_output = args.response_output.resolve()
    response_output.parent.mkdir(parents=True, exist_ok=True)
    response_output.write_text(response_matrix(), encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Wrote {response_output}")


if __name__ == "__main__":
    main()
