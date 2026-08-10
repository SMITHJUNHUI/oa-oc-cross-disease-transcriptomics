from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import build_manuscript_v22 as v22


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in fundamentally different tissues, yet both involve stromal, immune, metabolic, and stress-response remodeling. We tested whether they show context-dependent molecular convergence characterized by directional and cellular heterogeneity rather than a shared disease program or shared genetic causality.

### Methods

We analyzed public human bulk and single-cell transcriptomic datasets with a reproducible workflow. Differential expression was estimated separately in OA cartilage (GSE114007; 20 OA and 18 non-OA samples) and ovarian tissue (GSE18520; 53 high-grade serous tumors and 10 normal samples). Direction-aware overlap was integrated with disease-specific weighted gene co-expression network analysis (WGCNA). The ten-gene set was treated as an interpretable evidence summary, while machine-learning robustness was evaluated separately with 50 repeated five-fold outer resamples in which screening was repeated inside each training fold. Secondary external-cohort analyses used unsupervised principal-component analysis (PCA), label permutation, leave-one-sample-out assessment, and a direction-fixed score. Complete paired Hallmark enrichment, Human Protein Atlas (HPA), exploratory TCGA-OV associations, five single-cell datasets, descriptive gene-cell-function mapping, and supplementary bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes (DEGs), with 286 shared. Only 146/286 (51.0%) were directionally concordant. Among 50 paired Hallmark sets, 10 were significant in both diseases; 6 were directionally discordant and 4 concordant. Ten genes were prioritized (SOX9, ELF3, JUNB, AKAP12, BNC1, CFI, DDIT3, DIRAS3, EFEMP1, and HK2), but only AKAP12, JUNB, and DIRAS3 were concordant. Secondary assessment showed weak OA separation and near-complete retrospective OC tumor-normal molecular separation, without supporting diagnostic utility. Across five single-cell datasets, 1,025,361 QC-pass cells placed candidate expression in distinct cartilage, malignant, stromal, and immune contexts; gene-cell-function mapping reinforced that the same gene did not imply the same biological meaning. Relative TCGA-OV context associations were modest, and MR using 21 OA-to-OC and 11 OC-to-OA instruments detected no causal effect in either direction.

### Conclusions

OA and OC show partial context-dependent molecular convergence, but the overlap is directionally heterogeneous, cellularly contingent, and unsupported as shared genetic causality. The evidence supports an auditable cross-disease molecular resource, not a shared disease mechanism, diagnostic signature, prognostic model, therapeutic target set, or causal OA-OC relationship.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; single-cell RNA sequencing; WGCNA; nested cross-validation; pathway direction; reproducibility
"""


FUNCTIONAL_METHOD = """### Functional, pathway-direction, and immune annotation

Ranked enrichment used Hallmark, canonical pathway, Gene Ontology, and oncogenic-signature collections [8,9]; Gene Ontology and KEGG over-representation used the primary shared set. Direction was interpreted separately by disease. Rank-based immune signatures were treated as expression associations, not measured cell counts.

For the V2.3 secondary pathway-direction analysis, all 50 Hallmark sets were retained by re-estimating enrichment from the original complete OA and OC log2-fold-change rankings with the same gene-set definitions, size limits, GSEA implementation, and Benjamini-Hochberg correction. OA and OC normalized enrichment scores (NES) were paired by Hallmark identifier. Same-sign NES values were classified as concordant and opposite signs as discordant. Their product was reported as a descriptive paired direction index: positive values denote concordance and negative values discordance. This post hoc index compares independently estimated directions and is not a test of a shared pathway mechanism.
"""


CANDIDATE_METHOD = """### Candidate prioritization and strict nested machine learning

Candidate construction and model-performance estimation were treated as separate operations. Shared primary DEGs were ranked by a combined evidence score comprising adjusted significance and absolute effect size in both diseases. When at least five shared genes occurred in the disease-specific primary WGCNA modules, the candidate space was restricted to that module union. Disease-specific LASSO and random forest models then provided consensus and model-vote evidence. SOX9 was the cross-disease consensus; because this was below the prespecified minimum, the set was completed to ten genes by ranked multi-model votes. The ten-gene set was selected as an interpretable evidence summary rather than an optimized predictive signature. Table S16 reports shared-DEG status, direction, WGCNA support, separate LASSO and random-forest support, strict nested selection frequencies, and single-cell context for every gene. The resulting signed score was retained only for transparent secondary external assessment.

Internal molecular separation was re-estimated independently of that global candidate set. In 50 repeated five-fold outer resamples, all measured genes were ranked by training-only two-group statistics and the top 100 entered modeling. Five-fold inner cross-validation selected LASSO lambda; random-forest mtry was chosen from 2, 5, and 10 by training-fold out-of-bag AUC before fitting 300 trees [11,12]. Outer test samples were excluded from screening, tuning, scaling, and fitting. Balanced accuracy, Brier score, AUC, fit success, and feature-selection frequency were summarized across outer predictions [32]; AUC was retained in Table S5 as a secondary model metric.
"""


SINGLE_CELL_METHOD = """### Single-cell datasets and inference

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17-21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. Candidate expression was summarized by annotated cell type using mean log-normalized expression and detection fraction. Disease inference used pseudobulk only with biological replication and an interpretable contrast [23]; other datasets contributed localization. OA and OC were not integrated into a shared latent space, and transferred labels remained conservative and dataset specific [24].

For exploratory functional localization, the top 25 mean-expression contrast genes per analysis cluster in GSE104782 and GSE154600 were mapped to the majority dataset-specific cell label, pooled within cell type, and tested for Gene Ontology Biological Process over-representation against the annotated human-symbol universe. Benjamini-Hochberg correction was applied separately within each cell-type analysis. The gene-cell-function matrix then joined each candidate's highest-detection disease-specific cell context to the lowest-FDR available term for that label. Labels not represented in the two marker-annotation datasets were retained as unavailable rather than imputed. These cluster-marker ranks and joins are descriptive; they localize functional themes but do not demonstrate conserved function or mechanism.
"""


RESULTS = """## Results

### The shared DEG set was directionally heterogeneous

The study included two bulk discovery cohorts, four external bulk cohorts, TCGA-OV, five single-cell datasets, HPA annotations, and two GWAS datasets (Figure 1; Table S1).

At the primary threshold, 2,008 genes were differentially expressed in OA and 2,310 in OC, with 286 shared (Figure 2; Table S2). The explicit quadrants contained 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. Thus, only 146/286 (51.0%) were concordant. Directional heterogeneity persisted across all six prespecified thresholds (Table S3), showing that shared membership was not a uniform activation program.

All ten prioritized genes met the primary rule in both diseases, but only AKAP12, JUNB, and DIRAS3 were concordant. SOX9 was lower in OA (log2 fold change -2.072) and higher in OC (+2.612), whereas EFEMP1 was higher in OA (+3.092) and lower in OC (-2.781). The candidate set therefore cannot be described as uniformly activated or suppressed.

### Pathway-level convergence was also directionally heterogeneous

Complete paired analysis retained all 50 Hallmark sets (Figure S7; Table S18). Twenty-two were significant in OA and 16 in OC; 10 reached FDR <0.05 in both diseases. Of these, four were concordant (E2F targets, G2M checkpoint, apoptosis, and mitotic spindle), whereas six were discordant (epithelial-mesenchymal transition, coagulation, KRAS signaling up, mTORC1 signaling, complement, and glycolysis). For example, epithelial-mesenchymal transition was positive in OA (NES 2.347) and negative in OC (-1.833), while glycolysis was negative in OA (-1.410) and positive in OC (1.774). These paired directions reinforce that pathway membership does not imply the same pathway state or a shared mechanism.

### WGCNA associations were stable within the available cohorts

Primary module-trait correlations were -0.951 in OA and -0.879 in OC. Bootstrap sign stability was 1.000 in both across 2,000 resamples. Power perturbation retained absolute correlations of 0.904-0.948 in OA and 0.856-0.885 in OC; leave-one-sample-out estimates did not reverse signs (Figure 3; Figure S1; Table S4). These support stable within-cohort associations, with limited precision in the small OA network.

### Single-cell data resolved distinct cellular and functional contexts

The five adapters audited 1,187,436 cells and retained 1,025,361 QC-pass cells (Figure 1; Table S9). The gene-disease context matrix localized candidates to cartilage chondrocyte/fibrochondrocyte states in OA and to malignant epithelial, fibroblast, endothelial, mast-cell, and other immune-associated contexts in OC (Figure 5; Figure S3; Table S10). Eligible OA pseudobulk contrasts supported cell-state-specific associations for CFI, DIRAS3, AKAP12, EFEMP1, and JUNB; OC datasets mainly contributed localization because controls or biological replication were inconsistent.

Exploratory annotation of descriptive top cluster markers further separated functional contexts (Figure S6; Table S17). The candidate gene-cell-function matrix made this difference explicit (Table S19). SOX9 was most detected in OA hypertrophic chondrocytes linked descriptively to bone mineralization, but in OC its top context was ovarian cancer cells linked to cytoplasmic translation. AKAP12 mapped to OA prehypertrophic chondrocytes with extracellular-matrix organization and to OC fibroblasts with collagen-fibril organization. Several top contexts were not represented in the two marker-annotation datasets and were reported as unavailable rather than inferred. These mappings do not test conserved OA-OC function; they show why the same gene need not have the same biological meaning.

### Candidate hierarchy and nested resampling emphasized feature stability

The prioritization matrix shows why each of the ten genes was retained (Table S16). It separately displays shared-DEG status, direction, WGCNA membership, LASSO support, random-forest support, strict nested frequency, and single-cell context. SOX9 was the only cross-disease model consensus, while the remaining genes entered through the prespecified ranked-vote completion; five candidates belonged to the OA primary WGCNA module and five to the OC primary module. The set is therefore an interpretable evidence summary, not an optimized predictive signature.

Strict outer-fold frequencies were heterogeneous rather than uniformly high: among the final genes, the maximum OA frequency ranged from 0 to 0.940 and the maximum OC frequency from 0 to 0.200 (Figure 3C; Figure S8; Table S5). Screening across all measured genes was repeated inside every outer training fold. Median balanced accuracy was 1.000 for three disease-model combinations and 0.991 for the OC random forest; median Brier scores ranged from 0.0005 to 0.0198. These near-ceiling within-cohort contrasts address leakage from global screening but do not establish specificity of the ten-gene summary, diagnostic accuracy, or prospective transportability.

### OC showed near-complete retrospective molecular separation, not diagnostic validation

As a secondary reproducibility assessment, Figure 4 first presents unsupervised structure, permutation behavior, and sample influence before the ROC display. In GSE54388, PC1 and PC2 explained 33.0% and 15.1% of variance, and the 6 normal samples were visually separated from the 16 tumors without label-informed fitting. Empirical permutation P values were 0.469 and 0.207 in the two OA cohorts and 0.001 in both OC cohorts. Leave-one-out ranges were 0.456-0.589 and 0.583-0.698 in OA, remained 1.000 throughout GSE54388, and were 0.971-0.993 in GSE12470.

The direction-fixed score yielded OA AUCs of 0.520 (95% CI 0.250-0.790) and 0.629 (0.350-0.907), compared with 1.000 in GSE54388 and 0.979 (0.944-1.000) in GSE12470 (Figure 4D; Table S6). These findings support near-complete retrospective OC tumor-normal molecular separation, while small normal groups, platform structure, comparator tissue, and retrospective preprocessing preclude a diagnostic claim.

### Normal-reference annotations did not support a uniformly tissue-specific panel

HPA classified AKAP12 and DDIT3 as low-specificity, six candidates as tissue enhanced, CFI as liver enriched, and DIRAS3 as group enriched in brain, ovary, and pituitary. Only DIRAS3 listed ovary among specific tissues (Figure S4; Table S14). GTEx-integrated HPA lacks articular cartilage, and normal ovary cannot reproduce tumor composition; residual tissue-composition confounding therefore remains unresolved.

### TCGA-OV showed an exploratory survival association and modest context correlations

The continuous TCGA score retained an age/stage-adjusted association (hazard ratio per standard deviation 1.262, 95% CI 1.076-1.479), but apparent concordance was 0.597 and optimism-corrected concordance 0.586; selection instability and proportional-hazards concerns remained (Figure 6; Figure S1; Table S7). This was an exploratory association, not a prognostic model.

Relative transcriptomic context scores provided a composition audit across 307 TCGA-OV samples. Nine of ten candidates were present in the matrix (EFEMP1 was unavailable). Associations were modest: the largest absolute correlation was BNC1 with stromal score (Spearman rho=0.393, FDR=2.55x10-11), followed by BNC1 with the combined score (rho=0.358) and ELF3 with immune score (rho=0.297) (Figure S5; Table S15). These scores contextualize expression but do not measure absolute purity or histologic fractions.

### Supplementary MR did not support inherited causal liability

MR used the verified ebi-a-GCST007092 and ieu-a-1120 datasets with 21 and 11 instruments. Inverse-variance-weighted estimates were null for OA-to-OC (OR 1.015, 95% CI 0.900-1.144; P=0.811) and OC-to-OA (OR 1.040, 0.955-1.132; P=0.371); reverse MR also showed heterogeneity and MR-PRESSO outliers (Figure S2; Table S12). Thus, transcriptomic convergence did not coincide with detected inherited genetic liability under the available instruments and assumptions.
"""


DISCUSSION_TO_CONCLUSION = """## Discussion

### Shared genes are not shared programs

The central result is not the existence of 286 overlapping DEGs, but the heterogeneity hidden within that overlap. Only 51.0% were directionally concordant, and seven of ten prioritized genes changed oppositely. The pathway-direction analysis extended the same observation above the gene level: among 10 Hallmark sets significant in both diseases, 6 had opposite NES signs. Shared membership can therefore mark participation in remodeling-associated biology without indicating an equivalent disease program.

This distinction is especially important for broad terms such as epithelial-mesenchymal transition, glycolysis, coagulation, complement, or mTORC1 signaling. Their opposing OA and OC enrichment directions are independently estimated cohort-level states. They do not establish pathway inhibition, activation in a specific cell type, or a shared pathway mechanism.

### Shared genes have different cellular meanings

Single-cell localization showed that candidate expression arose in cartilage cell states in OA but in malignant, stromal, endothelial, mast-cell, and other immune-associated contexts in OC. The gene-cell-function matrix provided a transparent descriptive bridge: SOX9 linked to an OA hypertrophic-chondrocyte mineralization context but an OC malignant-cell translation context, whereas AKAP12 linked to OA extracellular-matrix organization and OC fibroblast collagen organization. Same gene therefore did not equal same biological meaning.

These joins remain constrained by dataset-specific annotations and by the two datasets with marker-based functional annotation. Unmatched labels were reported rather than imputed, OA and OC were not forced into a shared latent space, and cell-type enrichment was not treated as perturbational evidence.

### Transcriptomic convergence does not imply genetic causality

Bidirectional MR did not detect inherited causal liability in either direction, while transcriptomic, pathway, and single-cell results remained associations. First, the available GWAS datasets were predominantly derived from European populations, limiting generalizability. Second, OA GWAS mainly represented hip/knee phenotypes and may not capture all OA subtypes. Third, MR cannot exclude environmentally mediated pathways, such as aging-related inflammation or tissue injury responses.

The HPA/GTEx and TCGA analyses reinforce this boundary. HPA/GTEx does not cover articular cartilage, normal ovary does not reproduce a tumor microenvironment, and relative stromal/immune scores are not histologic fractions. The TCGA result is an exploratory survival association with modest optimism-corrected discrimination, not a validated prognostic model.

### Parallel stress adaptation is a hypothesis, not a mechanism

The observed convergence may reflect parallel adaptation to chronic stress rather than a shared disease program. Aging, persistent inflammation, matrix remodeling, tissue injury, metabolic stress, and endocrine context could produce partially overlapping expression responses in anatomically unrelated diseases. This is a hypothesis for future spatial, longitudinal, protein-level, and perturbational work; the current design does not establish a common mechanism.

The machine-learning results require the same restraint. The ten-gene set was an interpretable evidence summary rather than an optimized predictive signature. Heterogeneous nested selection frequencies, weak OA external separation, strong retrospective tumor-normal contrasts, small normal groups, and platform or comparator-tissue structure limit transportability. Figure 4 therefore places PCA, permutation, and sample-influence checks before ROC curves, and the detailed OA/OC feature heatmaps are supplementary.

### Strengths

Strengths include disease-specific discovery, two external cohorts per disease, explicit gene and pathway directions, prespecified DEG sensitivity, 2,000 WGCNA bootstraps, strict outer-fold feature selection, transparent candidate accounting, unsupervised PCA, permutation and sample-influence checks, HPA and TCGA context audits, five single-cell adapters, a bounded gene-cell-function matrix, complete MR provenance, exact figure source data, and a claim-evidence registry. Null and unfavorable results were retained.

### Limitations

The study is observational and retrospective. Cohorts differed in platform, tissue, processing, disease stage, comparator tissue, and clinical composition, with incomplete patient-level metadata. DEG overlap remained threshold dependent. The complete Hallmark direction comparison was a secondary post hoc analysis, and the paired direction index is descriptive rather than inferential. OA WGCNA used 38 samples and lacked independent module-preservation testing. Candidate completion by ranked model votes and heterogeneous strict nested selection frequencies limit claims of a stable panel. Near-perfect internal and OC external AUCs may reflect strong tumor-normal or disease-reference contrasts and are not clinical diagnostic performance. HPA/GTEx lacks cartilage; immune and ESTIMATE-derived scores are transcriptomic proxies rather than measured fractions. TCGA survival analysis lacked independent validation. Single-cell datasets had unequal controls, annotations, and replication; some candidate contexts lacked a matched marker-enrichment label, and OC results were mainly localization. MR was predominantly European, represented hip/knee OA, and cannot exclude environmentally mediated pathways. Protein, spatial, perturbational, prospective, and multi-cancer specificity validation were unavailable.

## Conclusions

Osteoarthritis and ovarian cancer show partial context-dependent molecular convergence, but the overlap is shared without being identical, associated without being causal, and cellularly contingent rather than mechanistically conserved. The evidence provides a reproducible framework for studying cross-disease molecular overlap without presenting the ten genes as an optimized signature, therapeutic targets, diagnostic markers, or a shared disease mechanism.
"""


FIGURE_2 = """### Figure 2. Direction-aware cross-disease transcriptomic discovery

**A-B,** OA and OC differential-expression volcano plots. **C,** All commonly measured genes with the 286 primary shared DEGs colored by direction; the inset isolates the ten-gene interpretable evidence summary. **D,** Shared-gene counts under six prespecified thresholds.
"""


FIGURE_3 = """### Figure 3. Robust identification of molecular candidates

**A,** WGCNA module-trait association under soft-power perturbation. **B,** Primary-module gene retention. **C,** Summary of strict nested feature-selection frequency across the ten candidates; points are genes and diamonds are within-model maxima. **D,** Disease-specific discovery effects. Detailed OA and OC model heatmaps are in Figure S8.
"""


FIGURE_4 = """### Figure 4. Direction-fixed cross-cohort molecular reproducibility

**A,** Unsupervised PCA of GSE54388 using the 2,000 genes with highest sample standard deviation; group labels were used only for display. **B,** Null AUC distributions from 1,000 label permutations. **C,** Leave-one-sample-out AUCs. **D,** ROC curves for the fixed signed ten-gene score, deliberately shown last as a secondary display. These panels show retrospective molecular separation, not clinical diagnostic performance.
"""


FIGURE_S7 = """### Supplementary Figure 7. Hallmark pathway direction across diseases

**A,** OA and OC normalized enrichment scores (NES) for all 50 Hallmark sets. Same-sign estimates indicate directional concordance and opposite signs indicate discordance. **B,** Paired NES for the strongest pathways significant in both diseases. The paired direction index is descriptive and does not establish a shared mechanism.
"""


FIGURE_S8 = """### Supplementary Figure 8. Detailed strict nested feature stability

**A,** OA LASSO and random-forest selection frequencies for the ten candidates. **B,** Corresponding OC frequencies. These detailed heatmaps are separated from the main feature-stability summary.
"""


def validate_table(path: Path, expected_rows: int, required: set[str]) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != expected_rows:
        raise ValueError(f"{path.name}: expected {expected_rows} rows, found {len(rows)}")
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")


def build_manuscript(
    source: Path,
    candidate_table: Path,
    pathway_table: Path,
    gene_cell_table: Path,
) -> str:
    validate_table(
        candidate_table,
        10,
        {
            "shared_DEG",
            "direction",
            "WGCNA_support",
            "LASSO_support",
            "random_forest_support",
            "strict_nested_frequency",
            "single_cell_context",
            "ten_gene_set_role",
        },
    )
    validate_table(
        pathway_table,
        50,
        {"OA_NES", "OC_NES", "direction_class", "both_significant"},
    )
    validate_table(
        gene_cell_table,
        10,
        {
            "OA_cell_context",
            "OA_functional_theme",
            "OC_cell_context",
            "OC_functional_theme",
        },
    )

    text = source.read_text(encoding="utf-8")
    text = v22.replace_h2(text, "Abstract", ABSTRACT)
    text = v22.replace_h3(
        text,
        "Functional and immune annotation",
        FUNCTIONAL_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Candidate prioritization and strict nested machine learning",
        CANDIDATE_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Single-cell datasets and inference",
        SINGLE_CELL_METHOD,
    )
    text = v22.replace_h2(text, "Results", RESULTS)
    discussion_pattern = r"^## Discussion\n.*?(?=^## Data and code availability)"
    text, count = re.subn(
        discussion_pattern,
        DISCUSSION_TO_CONCLUSION.rstrip() + "\n\n",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"Expected one Discussion-to-Conclusion block, found {count}")

    text = v22.replace_h3(
        text,
        "Figure 2. Direction-aware cross-disease transcriptomic discovery",
        FIGURE_2,
    )
    text = v22.replace_h3(
        text,
        "Figure 3. Robust identification of molecular candidates",
        FIGURE_3,
    )
    text = v22.replace_h3(
        text,
        "Figure 4. Direction-fixed cross-cohort molecular reproducibility",
        FIGURE_4,
    )
    table_index = "## Supplementary table index"
    if FIGURE_S7.splitlines()[0] not in text:
        text = text.replace(
            table_index,
            FIGURE_S7.rstrip()
            + "\n\n"
            + FIGURE_S8.rstrip()
            + "\n\n"
            + table_index,
        )
    text = text.replace(
        "`run_submission_v22.ps1`",
        "`run_submission_v23.ps1`",
    )
    text = text.replace(
        "- **Table S16:** Candidate prioritization matrix spanning shared-DEG, WGCNA, model-vote, strict nested frequency, and single-cell evidence.",
        "- **Table S16:** Candidate evidence-summary matrix spanning shared-DEG status, direction, WGCNA, LASSO, random forest, strict nested frequency, and single-cell context.",
    )
    text = text.replace(
        "- **Table S17:** Exploratory GO Biological Process annotation of descriptive top cluster markers by dataset-specific cell type.",
        "\n".join(
            [
                "- **Table S17:** Exploratory GO Biological Process annotation of descriptive top cluster markers by dataset-specific cell type.",
                "- **Table S18:** Complete OA-OC Hallmark pathway-direction matrix with paired NES, FDR, and descriptive direction index.",
                "- **Table S19:** Candidate gene-cell-function context matrix with explicit unavailable states for unmatched labels.",
            ]
        ),
    )
    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def response_matrix() -> str:
    return """# OC-OA manuscript V2.3: pre-submission refinement matrix

This working record maps the third-round recommendations to implemented changes. Processing remained local. The accountable authors must verify the final scientific wording and journal policy before submission.

| Recommendation | V2.3 action | Location | Status |
|---|---|---|---|
| Fix the article around context-dependent molecular convergence | Retained the V2.2 title and rewrote the Abstract, Discussion, and Conclusion around “shared but not identical; associated but not causal.” | Title; Abstract; Discussion; Conclusion | Implemented |
| Explain why there are ten genes | Upgraded Table S16 with explicit shared-DEG, direction, WGCNA, separate LASSO/RF, strict nested frequency, and single-cell columns. Added the exact boundary that the set is an interpretable evidence summary rather than an optimized predictive signature. | Methods; Results; Table S16 | Implemented |
| Reduce AUC centrality | Removed AUC values from the Abstract, introduced external validation as secondary, reordered Figure 4 to PCA, permutation, leave-one-out, then ROC, and kept exact AUC values only in the detailed Results/Table S6. | Abstract; Results; Figure 4; Table S6 | Implemented |
| Add three MR limitations without rerunning MR | Added population ancestry, hip/knee OA phenotype, and environmentally mediated pathway limitations. MR estimates were not rerun. | Discussion; Limitations | Implemented |
| Add pathway-direction analysis | Re-estimated the complete 50-set Hallmark output from the original full rankings, paired OA/OC NES, classified concordance/discordance, and added a descriptive paired direction index. | Methods; Results; Figure S7; Table S18 | Implemented |
| Add a gene-cell-function matrix | Joined each candidate’s top disease-specific cell context to the lowest-FDR available cell-type marker term; unmatched labels remain explicitly unavailable. | Methods; Results; Table S19 | Implemented with missing-state protection |
| Improve Figure 2C | Main panel now shows all genes without candidate labels; a separate inset isolates the ten candidates. | Figure 2C | Implemented |
| Simplify Figure 3C | Replaced the main heatmap with a feature-stability summary and moved separate OA/OC detailed heatmaps to Figure S8. | Figure 3C; Figure S8 | Implemented |
| Reorganize Discussion | Rebuilt the interpretation around four themes: shared genes are not shared programs; different cellular meanings; no genetic causality; parallel stress adaptation as a hypothesis. | Discussion | Implemented |
| Standardize cautious terminology | Positive claims use association, context, convergence, localization, or descriptive direction. Mechanism, diagnostic, prognostic, and target terms are retained only in explicit negations or limitations. | Full manuscript; claim-evidence registry | Implemented |
| Do not add PPI, drug prediction, virtual knockout, or more ML models | Added none of these analyses. Existing interaction tables remain supplementary and hypothesis-generating. | Project scope; Table S13 boundary | Preserved |
| Complete submission information | No author, affiliation, funding, conflict, CRediT, journal, or DOI information was invented. These remain author-controlled pending fields. | Title page; availability/declaration sections | Awaiting author input |

## Author-controlled items still required

- Authors, affiliations, and corresponding-author details.
- Target journal and current journal-specific formatting.
- Funding, competing interests, and CRediT contributions.
- Public repository URL and archival DOI after deposition.
- Final citation, language, statistical, and scientific approval.
- Journal-required disclosure of AI-assisted editing, if applicable.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--candidate-table", required=True, type=Path)
    parser.add_argument("--pathway-table", required=True, type=Path)
    parser.add_argument("--gene-cell-table", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()

    manuscript = build_manuscript(
        args.source.resolve(),
        args.candidate_table.resolve(),
        args.pathway_table.resolve(),
        args.gene_cell_table.resolve(),
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
