from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import build_manuscript_v22 as v22


ABSTRACT = """## Abstract

### Background

Osteoarthritis (OA) and ovarian cancer (OC) arise in fundamentally different tissues. We tested whether their transcriptomes show context-dependent molecular convergence characterized by directional, cellular, and regulatory heterogeneity rather than a shared disease program or inherited causal relationship.

### Methods

We analyzed public human bulk and single-cell transcriptomic datasets with a reproducible workflow. Differential expression was estimated separately in OA cartilage (GSE114007; 20 OA and 18 non-OA samples) and ovarian tissue (GSE18520; 53 high-grade serous tumors and 10 normal samples). Direction-aware overlap was integrated with disease-specific weighted gene co-expression network analysis. The ten-gene set was treated as an interpretable evidence summary; LASSO, random forest, and strict nested resampling provided complementary feature evidence. External cohorts were evaluated as biologically different molecular-separability tasks using a direction-fixed score, unsupervised principal-component analysis, label permutation, sample omission, standardized effect sizes, and cross-fitted calibration sensitivity. Paired Hallmark direction, exploratory candidate-centered enrichment, focused KnockTF/miRTarBase context, Human Protein Atlas, TCGA-OV, five single-cell datasets, and supplementary bidirectional Mendelian randomization (MR) constrained interpretation.

### Results

We identified 2,008 OA and 2,310 OC differentially expressed genes, with 286 shared; only 146/286 (51.0%) were directionally concordant. Among 50 paired Hallmark sets, 10 were significant in both diseases, including 6 with discordant directions. Ten genes were prioritized, but only AKAP12, JUNB, and DIRAS3 were concordant, and strict nested selection frequencies were heterogeneous. The same molecular summary showed modest and imprecise separation in OA cohorts but near-complete retrospective tumor-reference separation in OC cohorts. Candidate-centered analyses of SOX9, DDIT3, BNC1, and AKAP12 showed disease-specific pathway associations without establishing single-gene mechanisms. Across five single-cell datasets, 1,025,361 quality-control-pass cells placed candidates in distinct cartilage, malignant, stromal, and immune contexts. Focused regulator matrices remained hypothesis-generating, TCGA-OV context associations were modest, and MR detected no causal effect in either direction.

### Conclusions

OA and OC show partial context-dependent molecular convergence, but the overlap is shared without being identical, differently separable across validation tasks, cellularly contingent, and unsupported as shared genetic causality. The evidence supports an auditable cross-disease molecular resource, not a shared mechanism, diagnostic signature, prognostic model, therapeutic target set, or causal OA-OC relationship.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; single-cell RNA sequencing; WGCNA; nested cross-validation; pathway direction; molecular separability
"""


STUDY_DESIGN_METHOD = """### Study design and reproducibility

OA and OC were analyzed separately at every disease-specific stage and connected only by explicit gene- and pathway-level comparisons. The linear evidence chain comprised bulk discovery, direction-aware overlap, paired pathway direction, disease-specific WGCNA, transparent candidate evidence, strict nested feature stability, cross-cohort molecular separability, single-cell localization, HPA/TCGA context audits, focused upstream regulatory context, and supplementary bidirectional MR (Figure 1). Each layer answered a distinct question and was assigned an explicit inference boundary. Fixed configurations, manifests, seeds, tests, result tables, figure source data, and a claim-evidence registry are stored in the one-command project; randomized analyses used seed 20260726.
"""


FUNCTIONAL_METHOD = """### Functional, pathway-direction, and immune annotation

Ranked enrichment used Hallmark, canonical pathway, Gene Ontology, and oncogenic-signature collections [8,9]; Gene Ontology and KEGG over-representation used the primary shared set. Direction was interpreted separately by disease. Rank-based immune signatures were treated as expression associations, not measured cell counts.

All 50 Hallmark sets were retained by re-estimating enrichment from the original complete OA and OC log2-fold-change rankings with the same gene-set definitions, size limits, GSEA implementation, and Benjamini-Hochberg correction. OA and OC normalized enrichment scores (NES) were paired by identifier. Same-sign NES values were classified as concordant and opposite signs as discordant. Their product was reported as a descriptive paired direction index; it is not a test of a shared pathway mechanism.

Four transparent post hoc exemplars were chosen for candidate-centered context: SOX9 was the only cross-disease model consensus, DDIT3 and BNC1 had the highest strict nested selection frequencies in OA and OC, respectively, and AKAP12 combined concordant bulk direction with non-zero OA stability and explicit OA/OC cell-function contexts. Within each discovery cohort, every transcript and the candidate were residualized on disease/reference status. Residual partial-association t statistics ranked genes for complete Hallmark GSEA. Missing estimates caused by insufficiently balanced residual rankings were retained explicitly. This analysis tests candidate-associated pathway context after group adjustment; it is not single-gene perturbation, mediation, pathway activation, or mechanism.
"""


CANDIDATE_METHOD = """### Candidate prioritization and strict nested machine learning

Candidate construction and model-performance estimation were treated as separate operations. Shared primary differentially expressed genes were ranked by adjusted significance and absolute effect size in both diseases, restricted to the union of disease-specific primary WGCNA modules when the prespecified minimum was met. Disease-specific LASSO and random forest then supplied complementary model-vote evidence. SOX9 was the only cross-disease consensus; the set was completed to ten genes by the prespecified ranked multi-model vote rule. Table S16 reports shared-DEG status, direction, WGCNA membership, original LASSO support, original random-forest support, strict nested frequency, and single-cell context. Thus, LASSO and random forest support feature prioritization, whereas nested resampling estimates stability; their intersection does not replace the nested design. The ten-gene set is an interpretable evidence summary rather than an optimized predictive signature.

Internal molecular separation was re-estimated independently of the global candidate set. In 50 repeated five-fold outer resamples, all measured genes were ranked by training-only two-group statistics and the top 100 entered modeling. Five-fold inner cross-validation selected LASSO lambda; random-forest mtry was selected from 2, 5, and 10 by training-fold out-of-bag AUC before fitting 300 trees [11,12]. Outer test samples were excluded from screening, tuning, scaling, and fitting. Balanced accuracy, Brier score, AUC, fit success, model-specific feature importance, and selection frequency were summarized across outer predictions [32]. AUC remained a secondary model metric.
"""


CROSS_COHORT_METHOD = """### Cross-cohort molecular separability

External gene orientation was fixed from discovery log2 fold change. AUC and DeLong 95% confidence intervals were calculated without choosing direction in validation data [13]; the signed molecular summary averaged standardized expression multiplied by discovery signs. GSE54388 also underwent unsupervised PCA using the 2,000 genes with highest sample standard deviation, with phenotype labels used only for display. Robustness used 1,000 label permutations and sequential sample omission.

Because the same summary was evaluated in validation tasks with different tissues and comparator scales, each cohort was annotated by tissue, comparator, sample size, AUC interval, permutation result, and leave-one-out range. Hedges g summarized the standardized disease-reference score contrast. Descriptive random-effects summaries were calculated separately for the two OA and two OC cohorts; with only two cohorts per disease, pooled effects and heterogeneity remained imprecise.

The signed score was not a locked probability model. Calibration was therefore restricted to a supplementary sensitivity analysis: each sample received a probability from leave-one-sample-out ridge logistic recalibration of the fixed score (slope L2 penalty 0.1). Brier score and calibration intercept/slope were reported only when numerically and directionally stable. This is cross-fitted cohort-specific reliability assessment, not transported external probability calibration or clinical utility analysis. Decision-curve analysis and a nomogram were not performed because no clinical decision model was specified.
"""


SINGLE_CELL_AND_REGULATORY_METHOD = """### Single-cell datasets and inference

OA GSE104782, GSE169454, and GSE255460 and OC GSE154600 and GSE180661 were processed with dataset-specific adapters [17-21]. Released count layers and metadata were preserved; count, feature, and mitochondrial-fraction outliers were assessed within reliable sample partitions, and scDblFinder was required where count-level data and suitable partitions allowed it [22]. Candidate expression was summarized by annotated cell type. Disease inference used pseudobulk only with biological replication and an interpretable contrast [23]; other datasets contributed localization. OA and OC were not integrated into a shared latent space, and transferred labels remained conservative and dataset specific [24].

For exploratory functional localization, the top 25 mean-expression contrast genes per analysis cluster in GSE104782 and GSE154600 were mapped to the majority dataset-specific cell label, pooled within cell type, and tested for Gene Ontology Biological Process over-representation. The gene-cell-function matrix joined each candidate's highest-detection disease-specific cell context to the lowest-FDR available term for that label. Unmatched labels remained unavailable rather than imputed. These joins localize functional themes but do not demonstrate conserved function or mechanism.

### Focused upstream regulatory context

The regulatory analysis was designed to ask whether candidate genes occurred within different upstream evidence contexts, not to construct a generic network. KnockTF 2.0 perturbational top-differential target sets were tested against the complete OA and OC discovery rankings [34]. Only transcription factors whose target sets contained SOX9, DDIT3, BNC1, or AKAP12 were retained and paired across diseases. Target-set NES does not estimate transcription-factor activity or direct regulation in OA or OC.

Curated miRTarBase interactions were filtered to the same four candidates [35]. Multi-candidate miRNAs were summarized as an interaction-coverage matrix. miRNA expression and activity were not measured in the analyzed cohorts; these records are testable upstream hypotheses, not disease-specific regulatory inference. No TF-miRNA causal network was asserted.
"""


RESULTS = """## Results

### The shared DEG set was directionally heterogeneous

The study included two bulk discovery cohorts, four external bulk cohorts, TCGA-OV, five single-cell datasets, HPA annotations, focused regulatory resources, and two GWAS datasets (Figure 1; Table S1).

At the primary threshold, 2,008 genes were differentially expressed in OA and 2,310 in OC, with 286 shared (Figure 2; Table S2). The four explicit direction quadrants contained 112 genes higher in both diseases, 34 lower in both, 86 higher in OA/lower in OC, and 54 lower in OA/higher in OC. Only 146/286 (51.0%) were concordant. Directional heterogeneity persisted across all six prespecified thresholds (Table S3).

All ten candidates met the primary rule in both diseases, but only AKAP12, JUNB, and DIRAS3 were concordant. SOX9 was lower in OA (log2 fold change -2.072) and higher in OC (+2.612), whereas EFEMP1 was higher in OA (+3.092) and lower in OC (-2.781). Shared membership was therefore not a uniform activation or suppression program.

### Pathway-level convergence was also directionally heterogeneous

Complete paired analysis retained all 50 Hallmark sets (Figure S7; Table S18). Twenty-two were significant in OA and 16 in OC; 10 reached FDR <0.05 in both diseases. Four were concordant and six discordant. Epithelial-mesenchymal transition was positive in OA (NES 2.347) and negative in OC (-1.833), whereas glycolysis was negative in OA (-1.410) and positive in OC (1.774). Paired pathway membership therefore did not indicate the same pathway state or a shared mechanism.

### WGCNA and multi-algorithm evidence supported a transparent, heterogeneous candidate hierarchy

Primary module-trait correlations were -0.951 in OA and -0.879 in OC. Bootstrap sign stability was 1.000 in both across 2,000 resamples; soft-power perturbation and leave-one-sample-out checks did not reverse signs (Figure 3; Figure S1; Table S4).

Table S16 separates WGCNA, original LASSO, original random-forest, strict nested, direction, and single-cell evidence. SOX9 was the only cross-disease LASSO/random-forest consensus; the remaining genes entered by the prespecified ranked-vote completion. Strict outer-fold frequencies were heterogeneous: among the ten candidates, maximum OA frequency ranged from 0 to 0.940 and maximum OC frequency from 0 to 0.200 (Figure 3C; Figure S8; Table S5). LASSO and random forest therefore supplied complementary prioritization evidence, while strict nested frequencies showed that the ten genes were not a uniformly stable predictive panel.

### Candidate-centered pathway contexts illustrated disease-specific associations

The exploratory residual-association analysis retained 400 candidate-disease-Hallmark rows, including explicit unavailable estimates for 7 candidate-disease pathways (Figure S9; Table S20). SOX9 in OC was positively associated with MYC targets, oxidative phosphorylation, mTORC1 signaling, protein secretion, and DNA repair. DDIT3 in OA was negatively associated with G2M checkpoint, mitotic spindle, epithelial-mesenchymal transition, and protein secretion, while positively associated with oxidative phosphorylation. BNC1 showed strong positive OC associations with epithelial-mesenchymal transition and interferon responses but negative OA associations with several proliferative/stress-related sets. AKAP12 was inversely associated with oxidative phosphorylation in OA but positively associated with oxidative phosphorylation, MYC targets, and E2F targets in OC. These are group-adjusted co-expression contexts, not candidate-gene perturbation effects.

### The same molecular summary had different cross-cohort separability

As a secondary reproducibility assessment, Figure 4 presents unsupervised structure, permutation behavior, and sample influence before ROC. In GSE54388, PC1 and PC2 explained 33.0% and 15.1% of variance, and the 6 reference samples were visually separated from the 16 tumors without label-informed fitting. Empirical permutation P values were 0.469 and 0.207 in the two OA cohorts and 0.001 in both OC cohorts. Leave-one-out AUC ranges were 0.456-0.589 and 0.583-0.698 in OA, remained 1.000 in GSE54388, and were 0.971-0.993 in GSE12470.

The direction-fixed summary yielded AUCs of 0.520 (95% CI 0.250-0.790) and 0.629 (0.350-0.907) in OA, compared with 1.000 and 0.979 (0.944-1.000) in OC (Figure 4D; Tables S6 and S21). The molecular summary demonstrated strong retrospective separation in ovarian cancer cohorts, whereas discrimination was modest in osteoarthritis cohorts, reflecting differences in tissue context, disease biology, comparator scale, and cohort composition.

Standardized score contrasts showed the same pattern: cohort Hedges g values were 0.015 and 0.306 in OA versus 5.243 and 2.692 in OC. Descriptive two-cohort random-effects summaries were 0.146 (95% CI -0.505 to 0.797) in OA and 3.826 (1.341 to 6.311) in OC (Figure S10A; Table S22a). These summaries are imprecise because each disease contributed only two heterogeneous cohorts.

Cross-fitted Brier scores were 0.316 and 0.305 in OA versus 0.010 and 0.057 in OC (Figure S10B; Table S22b). Calibration intercept/slope was estimable only in GSE12470 (0.008 and 0.919); the other three cohorts were numerically or directionally unstable. Calibration therefore did not rescue the weak OA task or convert the OC contrast into clinical validation.

### Single-cell data resolved distinct cellular and functional contexts

The five adapters audited 1,187,436 cells and retained 1,025,361 quality-control-pass cells (Figure 1; Table S9). Candidates localized to cartilage chondrocyte/fibrochondrocyte states in OA and to malignant epithelial, fibroblast, endothelial, mast-cell, and other immune-associated contexts in OC (Figure 5; Figure S3; Table S10). Eligible OA pseudobulk contrasts supported cell-state-specific associations for CFI, DIRAS3, AKAP12, EFEMP1, and JUNB; OC datasets mainly contributed localization because controls or biological replication were inconsistent.

The gene-cell-function matrix made contextual differences explicit (Table S19). SOX9 mapped to an OA hypertrophic-chondrocyte mineralization context but an OC malignant-cell translation context. AKAP12 mapped to OA prehypertrophic chondrocytes with extracellular-matrix organization and OC fibroblasts with collagen-fibril organization. Same gene therefore did not imply the same cellular meaning.

### Focused upstream evidence defined hypotheses rather than a regulatory mechanism

Among 157 KnockTF target-set rows linked to the four exemplars, several target sets showed concordant enrichment, whereas MAP3K7, MEIS2, DLX1, RNF2, TLE2, and MSX1 were positive in OA and negative in OC (Figure S11A; Table S23a). These directions describe enrichment of perturbational target sets and do not estimate transcription-factor activity. miRTarBase contributed 363 interaction-coverage rows; hsa-miR-335-5p was the only miRNA linked to all four exemplars in the local release (Figure S11B; Table S23b). Because cohort miRNA abundance was unavailable, no disease-specific miRNA direction was inferred.

### Normal-reference, TCGA, and immune analyses constrained interpretation

HPA classified AKAP12 and DDIT3 as low-specificity, six candidates as tissue enhanced, CFI as liver enriched, and DIRAS3 as group enriched in brain, ovary, and pituitary. Only DIRAS3 listed ovary among specific tissues (Figure S4; Table S14). HPA/GTEx lacks articular cartilage, and normal ovary cannot reproduce tumor composition.

The continuous TCGA score retained an age/stage-adjusted association (hazard ratio per standard deviation 1.262, 95% CI 1.076-1.479), but optimism-corrected concordance was 0.586 and selection instability remained (Figure 6; Figure S1; Table S7). Relative context-score associations were modest; the largest absolute correlation was BNC1 with stromal score (Spearman rho=0.393, FDR=2.55x10^-11) (Figure S5; Table S15). Rank-based immune associations remained expression proxies rather than measured cell fractions (Table S8).

### Supplementary MR defined an inherited-causality boundary

MR used ebi-a-GCST007092 and ieu-a-1120 with 21 and 11 instruments. Inverse-variance-weighted estimates were null for OA-to-OC (OR 1.015, 95% CI 0.900-1.144; P=0.811) and OC-to-OA (OR 1.040, 0.955-1.132; P=0.371); reverse MR also showed heterogeneity and MR-PRESSO outliers (Figure S2; Table S12). MR was performed not to make OA or OC the presumed cause of the other disease, but to test whether transcriptomic convergence was supported by inherited genetic liability. It was not detected under the available instruments and assumptions.
"""


DISCUSSION_TO_CONCLUSION = """## Discussion

### Shared genes are not shared programs

The central result is not the existence of 286 overlapping differentially expressed genes, but the heterogeneity hidden within that overlap. Only 51.0% were directionally concordant, seven of ten candidates changed oppositely, and six of ten Hallmark sets significant in both diseases had opposite NES signs. Shared membership can therefore mark participation in remodeling-associated biology without indicating an equivalent disease program.

The four candidate exemplars sharpen this distinction without creating a new signature. SOX9 represented the only cross-disease model consensus but had opposite bulk directions and different candidate-centered pathway and cell contexts. DDIT3 was selected most frequently within OA strict resampling and was associated with an OA stress/proliferation pattern that differed from its OC context. BNC1 carried the strongest OC nested frequency and the largest TCGA stromal correlation among available candidates, but its OA and OC pathway associations diverged. AKAP12 was bulk-concordant yet mapped to different OA cartilage and OC fibroblast functions and opposite oxidative-phosphorylation associations. These profiles explain why evidence can converge on the same gene without implying the same biological role.

### Molecular separability is a property of the validation task

The OA and OC AUCs should not be interpreted as competing estimates of one universal model. The same signed molecular summary was applied to within-tissue chronic degenerative or inflammatory contrasts in OA and to malignant tumor-reference contrasts in OC. These tasks differ in signal amplitude, tissue composition, comparator tissue, and disease biology. Near-complete OC separation and modest OA separation therefore support context dependence rather than a clinical prediction claim.

Permutation, sample omission, standardized score contrasts, and cross-fitted calibration add dimensions of reproducibility but do not erase the task difference. The unstable OA calibration and extremely large OC standardized effects make a nomogram, decision-curve analysis, or clinical-utility claim inappropriate. The result is cross-cohort molecular separability, not diagnostic performance.

### Shared genes have different cellular and regulatory meanings

Single-cell localization showed candidate expression in cartilage cell states in OA but malignant, stromal, endothelial, mast-cell, and other immune-associated contexts in OC. The gene-cell-function matrix provided a transparent descriptive bridge, while the candidate-centered Hallmark analysis showed different residual association structures after disease-status adjustment. Same gene did not equal same cellular or pathway meaning.

The focused upstream analysis adds hypotheses rather than a completed mechanism. KnockTF target-set direction cannot be read as transcription-factor activity, and miRTarBase interaction presence cannot establish regulation in an unmeasured tissue context. The value of these matrices is to nominate experimentally testable axes while displaying exactly which evidence is absent. They complement rather than replace the stronger transcriptomic, cellular, and genetic boundary analyses.

### Transcriptomic convergence does not imply genetic causality

Bidirectional MR did not detect inherited causal liability in either direction, while transcriptomic, pathway, single-cell, and regulatory results remained associations. MR was retained as an explicit causality boundary, not as the core mechanism analysis. First, the available GWAS datasets were predominantly derived from European populations, limiting generalizability. Second, OA GWAS mainly represented hip/knee phenotypes and may not capture all OA subtypes. Third, MR cannot exclude environmentally mediated pathways, including aging-related inflammation or tissue injury responses.

HPA/GTEx and TCGA reinforce the same boundary. HPA/GTEx does not cover articular cartilage, normal ovary does not reproduce a tumor microenvironment, relative stromal/immune scores are not histologic fractions, and the TCGA association is not a validated prognostic model.

### Parallel stress adaptation is a hypothesis, not a mechanism

The observed convergence may reflect parallel adaptation to chronic stress rather than a shared disease program. Aging, persistent inflammation, matrix remodeling, tissue injury, metabolic stress, and endocrine context could produce partially overlapping expression responses in unrelated tissues. This remains a hypothesis for spatial, longitudinal, protein-level, and perturbational work.

### Strengths

Strengths include disease-specific discovery, two external cohorts per disease, explicit gene and pathway directions, prespecified DEG sensitivity, 2,000 WGCNA bootstraps, strict outer-fold feature selection, separate LASSO/random-forest evidence, transparent candidate accounting, unsupervised PCA, permutation and sample-influence checks, cross-cohort effect sizes, bounded cross-fitted calibration, HPA and TCGA context audits, five single-cell adapters, a gene-cell-function matrix, complete MR provenance, focused perturbational/curated regulatory matrices, exact figure source data, and a claim-evidence registry. Null and unfavorable results were retained.

### Limitations

The study is observational and retrospective. Cohorts differed in platform, tissue, processing, disease stage, comparator tissue, and clinical composition. DEG overlap remained threshold dependent. The paired pathway-direction, candidate-centered enrichment, external effect-size summaries, and upstream regulatory analyses were secondary or post hoc. Only two external cohorts per disease precluded reliable meta-regression and made pooled heterogeneity imprecise. The fixed signed score was not a locked probability model; cross-fitted calibration was cohort-specific and unstable in three cohorts. Candidate completion by ranked model votes and heterogeneous strict nested frequencies limit claims of a stable panel. HPA/GTEx lacks cartilage; immune and TCGA context scores are transcriptomic proxies. TCGA survival analysis lacked independent validation. Single-cell datasets had unequal controls, annotations, and replication. KnockTF target-set enrichment does not measure transcription-factor activity, and miRNA abundance was unavailable. MR was predominantly European, represented hip/knee OA, and cannot exclude environmentally mediated pathways. Protein, spatial, prospective, perturbational, and multi-cancer specificity validation were unavailable.

## Conclusions

Osteoarthritis and ovarian cancer show partial context-dependent molecular convergence, but the overlap is shared without being identical, differently separable across validation tasks, associated without being causal, and cellularly contingent rather than mechanistically conserved. The ten genes are an interpretable evidence summary, and the focused pathway/regulatory analyses define testable contexts rather than diagnostic markers, therapeutic targets, or an established shared mechanism.
"""


FIGURE_1 = """### Figure 1. Linear study design and audited evidence boundaries

**A,** Linear evidence chain from disease-specific discovery through direction-aware overlap, pathway comparison, WGCNA and nested-model stability, external molecular separability, single-cell localization, context audits, and bidirectional MR. **B,** Audited resource scale, intended analytic role, and inference boundary. Each layer answers a distinct question; no layer is treated as proof of a shared mechanism.
"""


FIGURE_4 = """### Figure 4. Direction-fixed cross-cohort molecular separability

**A,** Unsupervised PCA of GSE54388 using the 2,000 most variable genes; group labels were used only for display. **B,** Null AUC distributions from 1,000 label permutations. **C,** Leave-one-sample-out AUCs. **D,** ROC curves for the fixed signed ten-gene molecular summary, shown last as a secondary display. The same summary was evaluated in tasks with different tissues and comparator scales (Table S21); panels show retrospective molecular separability, not clinical performance.
"""


FIGURE_S9 = """### Supplementary Figure 9. Candidate-centered Hallmark contexts

Disease-status-adjusted residual association rankings for SOX9, DDIT3, BNC1, and AKAP12 were analyzed against all 50 Hallmark sets. Six strongest pathways per candidate are shown; filled points denote FDR <0.05. This transparent post hoc analysis is exploratory and does not represent single-gene perturbation, mediation, pathway activation, or mechanism.
"""


FIGURE_S10 = """### Supplementary Figure 10. External evaluation context

**A,** Direction-fixed signed-score standardized mean differences by cohort with descriptive two-cohort random-effects summaries within disease. **B,** Three-bin calibration sensitivity from leave-one-sample-out ridge recalibration of the fixed score. A locked probability model was not transported from discovery; this is a cross-fitted sensitivity analysis rather than external clinical calibration.
"""


FIGURE_S11 = """### Supplementary Figure 11. Focused upstream regulatory context

**A,** OA and OC enrichment of KnockTF perturbational target sets for transcription factors connected to four candidate exemplars; NES does not estimate transcription-factor activity. **B,** Curated miRTarBase interaction coverage for multi-candidate miRNAs. miRNA abundance and activity were not measured, so these matrices define testable regulatory context rather than an inferred disease mechanism.
"""


def validate_table(path: Path, expected_rows: int, required: set[str]) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != expected_rows:
        raise ValueError(
            f"{path.name}: expected {expected_rows} rows, found {len(rows)}"
        )
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")


def build_manuscript(
    source: Path,
    candidate_pathways: Path,
    validation_context: Path,
    calibration: Path,
    tf_context: Path,
    mirna_context: Path,
) -> str:
    validate_table(
        candidate_pathways,
        400,
        {
            "candidate",
            "disease",
            "pathway_id",
            "NES",
            "FDR",
            "calculation_status",
            "inference_boundary",
        },
    )
    validate_table(
        validation_context,
        4,
        {
            "tissue_and_comparator",
            "validation_task_scale",
            "fixed_direction_AUC",
            "inference_boundary",
        },
    )
    validate_table(
        calibration,
        4,
        {
            "Brier_score",
            "calibration_status",
            "inference_boundary",
        },
    )
    with tf_context.open(encoding="utf-8-sig", newline="") as handle:
        tf_rows = list(csv.DictReader(handle))
    with mirna_context.open(encoding="utf-8-sig", newline="") as handle:
        mirna_rows = list(csv.DictReader(handle))
    if len(tf_rows) != 157 or len(mirna_rows) != 363:
        raise ValueError(
            f"Unexpected regulatory row counts: TF={len(tf_rows)}, "
            f"miRNA={len(mirna_rows)}"
        )

    text = source.read_text(encoding="utf-8")
    text = v22.replace_h2(text, "Abstract", ABSTRACT)
    text = v22.replace_h3(
        text,
        "Study design and reproducibility",
        STUDY_DESIGN_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Functional, pathway-direction, and immune annotation",
        FUNCTIONAL_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Candidate prioritization and strict nested machine learning",
        CANDIDATE_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Cross-cohort molecular reproducibility",
        CROSS_COHORT_METHOD,
    )
    text = v22.replace_h3(
        text,
        "Single-cell datasets and inference",
        SINGLE_CELL_AND_REGULATORY_METHOD,
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
        raise ValueError(
            f"Expected one Discussion-to-Conclusion block, found {count}"
        )

    text = v22.replace_h3(
        text,
        "Figure 1. Study design, biological hypothesis, and audited data resources",
        FIGURE_1,
    )
    text = v22.replace_h3(
        text,
        "Figure 4. Direction-fixed cross-cohort molecular reproducibility",
        FIGURE_4,
    )
    table_index = "## Supplementary table index"
    if FIGURE_S9.splitlines()[0] not in text:
        text = text.replace(
            table_index,
            FIGURE_S9.rstrip()
            + "\n\n"
            + FIGURE_S10.rstrip()
            + "\n\n"
            + FIGURE_S11.rstrip()
            + "\n\n"
            + table_index,
        )

    text = text.replace("`run_submission_v23.ps1`", "`run_submission_v24.ps1`")
    text = text.replace(
        "- **Table S19:** Candidate gene-cell-function context matrix with explicit unavailable states for unmatched labels.",
        "\n".join(
            [
                "- **Table S19:** Candidate gene-cell-function context matrix with explicit unavailable states for unmatched labels.",
                "- **Table S20:** Disease-status-adjusted candidate-centered Hallmark association context for four transparent exemplars.",
                "- **Table S21:** External-cohort tissue, comparator, task scale, AUC interval, permutation, sample-omission, and interpretation matrix.",
                "- **Table S22a-b:** Direction-fixed signed-score effect sizes and cross-fitted calibration sensitivity.",
                "- **Table S23a-b:** Focused KnockTF target-set and miRTarBase interaction context with explicit activity boundaries.",
            ]
        ),
    )

    reference_anchor = (
        "33. Yoshihara K, Shahmoradgoli M, Mart\u00ednez E, et al. "
        "Inferring tumour purity and stromal and immune cell admixture from "
        "expression data. *Nat Commun*. 2013;4:2612. "
        "doi:10.1038/ncomms3612."
    )
    references = "\n".join(
        [
            reference_anchor,
            "34. Feng C, Song C, Song S, et al. KnockTF 2.0: a comprehensive gene expression profile database with knockdown/knockout of transcription (co-)factors in multiple species. *Nucleic Acids Res*. 2024;52:D183-D193. doi:10.1093/nar/gkad1016.",
            "35. Huang HY, Lin YCD, Cui S, et al. miRTarBase update 2022: an informative resource for experimentally validated miRNA-target interactions. *Nucleic Acids Res*. 2022;50:D222-D230. doi:10.1093/nar/gkab1079.",
        ]
    )
    if "doi:10.1093/nar/gkad1016" not in text:
        text = text.replace(reference_anchor, references)

    return re.sub(r"\n{3,}", "\n\n", text).rstrip() + "\n"


def response_matrix() -> str:
    return """# OC-OA manuscript V2.4: Scientific Reports-informed strengthening matrix

This record maps the supplied RA-RF article and the subsequent reviewer-style recommendations to V2.4. The reference article was used as a structural comparator, not as a template for copying claims or low-evidence modules.

| Recommendation | V2.4 action | Location | Status |
|---|---|---|---|
| Learn from the RA-RF article's complete discovery-to-validation narrative | Rebuilt Figure 1 and Results as a linear evidence chain with an explicit role and inference boundary for every data layer. | Figure 1; Methods; Results | Implemented |
| Explain why RA-RF and OA-OC AUC patterns are not directly comparable | Added a four-cohort validation-task matrix covering tissue, comparator, sample size, AUC interval, permutation, leave-one-out range, and biological task scale. | Results; Discussion; Table S21 | Implemented |
| Reframe model performance as cross-cohort molecular separability | Renamed the section and Figure 4; added the exact distinction between modest OA separation and strong retrospective OC separation. | Results; Figure 4; Discussion | Implemented |
| Retain nested ML while adding LASSO/RF as auxiliary evidence | Kept strict nested screening/tuning; explicitly separated LASSO/RF prioritization evidence from nested stability estimation. No non-nested intersection replaced the primary design. | Methods; Figure 3; Figure S8; Tables S5 and S16 | Implemented |
| Add feature importance and selection frequency | Retained original model support, detailed disease/model selection frequencies, and the main stability summary without duplicating them in Figure 4. | Figure 3C; Figure S8; Tables S5 and S16 | Already present and clarified |
| Add calibration cautiously | Added leave-one-sample-out ridge recalibration of the fixed score, Brier scores, and calibration status. Unstable intercept/slopes are marked non-interpretable. | Figure S10B; Table S22b | Implemented as bounded sensitivity |
| Add external validation meta-analysis | Added cohort Hedges g and descriptive two-cohort random-effects summaries separately for OA and OC. | Figure S10A; Table S22a | Implemented with k=2 limitation |
| Do not let calibration imply a clinical model | Stated that no locked probability model was transported; did not add DCA, clinical utility claims, or a nomogram. | Methods; Discussion; Figure S10 legend | Preserved |
| Add single-gene GSEA that directly tests context dependence | Added disease-status-adjusted candidate-centered Hallmark analysis for SOX9, DDIT3, BNC1, and AKAP12, selected by transparent evidence axes. | Results; Figure S9; Table S20 | Implemented as exploratory |
| Add gene-by-gene interpretation | Added evidence-based profiles for SOX9, DDIT3, BNC1, and AKAP12 using direction, nested frequency, pathway, cell, and TCGA evidence. RELA, TNF, and IL6 were not promoted because they are outside the audited ten-gene set. | Discussion | Implemented |
| Retain MR but reduce its centrality | Kept MR as the final inherited-causality boundary and explicitly explained why it was performed. Estimates were not rerun. | Results last subsection; Discussion; Figure S2 | Implemented |
| Add TF/miRNA only if it explains context rather than making a decorative network | Added paired KnockTF perturbational target-set matrices and miRTarBase interaction coverage for four candidates. No TF activity, miRNA direction, or causal network was inferred. | Methods; Results; Figure S11; Tables S23a-b | Implemented with explicit missing-evidence boundary |
| Keep immune analysis | Retained rank-based immune associations and their proxy limitation; no additional immune estimator was added. | Results; Table S8 | Preserved |
| Avoid traditional hub-gene expansion | Added no nomogram, DCA, PPI, drug prediction, virtual knockout, SVM, XGBoost, LightGBM, or broad TF-miRNA network. | Full project scope | Preserved |
| Complete submission information | No author, affiliation, funding, conflict, CRediT, target-journal, repository URL, or DOI information was invented. | Title page and declarations | Awaiting author input |

## Author-controlled items still required

- Authors, affiliations, and corresponding-author details.
- Target journal and journal-specific formatting.
- Funding, competing interests, and CRediT contributions.
- Public repository URL and archival DOI after deposition.
- Final citation, language, statistical, and scientific approval.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--candidate-pathways", required=True, type=Path)
    parser.add_argument("--validation-context", required=True, type=Path)
    parser.add_argument("--calibration", required=True, type=Path)
    parser.add_argument("--tf-context", required=True, type=Path)
    parser.add_argument("--mirna-context", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()

    manuscript = build_manuscript(
        args.source.resolve(),
        args.candidate_pathways.resolve(),
        args.validation_context.resolve(),
        args.calibration.resolve(),
        args.tf_context.resolve(),
        args.mirna_context.resolve(),
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
