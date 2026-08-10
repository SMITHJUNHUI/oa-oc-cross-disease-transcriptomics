from __future__ import annotations

import argparse
import re
from pathlib import Path


MANUSCRIPT = r"""# Shared molecular features between osteoarthritis and ovarian cancer across tissue, cellular and systemic contexts

**Article type:** Original Research  
**Running title:** Shared molecular features across OA and OC  
**Authors:** [Author names to be added]  
**Affiliations:** [Affiliations to be added]  
**Corresponding author:** [Name, postal address and email to be added]  
**Target journal:** [To be selected before submission]

## Abstract

**Background.** Osteoarthritis (OA) and ovarian cancer (OC) arise in different tissues and follow distinct clinical courses. Both diseases nevertheless involve extracellular-matrix remodelling, inflammatory signalling and cellular stress. Whether their transcriptomic overlap is reproducible across tissue cohorts, cellular compartments and peripheral blood remains unclear.

**Methods.** OA and OC datasets were analysed on separate disease tracks. Discovery tissue transcriptomes were used to identify shared differentially expressed genes and classify their direction of change. Independent tissue cohorts evaluated replication. Representative genes were localized in separate single-cell atlases, and a prespecified peripheral-blood screen tested whether tissue-concordant signals persisted in both diseases.

**Results.** The discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes. Of these, 146 changed in the same direction and 140 changed in opposite directions. The shared set was enriched for extracellular-matrix, immune, stress and cell-cycle functions, while six of ten Hallmark pathways significant in both diseases had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 localized to different source-defined cell populations in OA and OC atlases. The blood screen retained one gene after independent false-discovery-rate control in both cohorts. G0S2 was lower in OA and OC tissues, OA peripheral-blood mononuclear cells and the OC blood-cell fraction.

**Conclusions.** OA and OC share a finite transcriptional background, but the direction, reproducibility and cellular location of these signals depend on disease context. G0S2 represents a candidate systemic molecular signal that requires prospective and protein-level confirmation.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; external validation; single-cell RNA sequencing; peripheral blood

## Introduction

Osteoarthritis is a whole-joint disorder characterized by cartilage loss, subchondral-bone remodelling, synovial changes and low-grade inflammatory signalling [1-7]. Transcriptomic and single-cell studies have further shown that diseased cartilage contains heterogeneous chondrocyte states and interacting stromal compartments [8-14]. Peripheral-blood expression studies suggest that selected OA-associated changes can also be detected outside the joint, although circulating profiles remain sensitive to cell composition and cohort design [15].

Ovarian cancer is a heterogeneous malignancy in which genomic instability, molecular subtype, stromal remodelling and immune context shape disease behaviour [16-25]. Bulk and single-cell studies of high-grade serous ovarian cancer have resolved diverse malignant, fibroblast, endothelial and immune populations across primary and metastatic sites [26-30]. Blood-based studies have also detected tumour-associated transcriptional variation, but these measurements integrate tumour effects with systemic and haematological responses [31].

Cross-disease transcriptomics can identify biological processes that recur across anatomically distant conditions [32-39]. Gene-list overlap alone, however, does not establish that two diseases use the same transcriptional state. A shared gene may change in opposite directions, fail to reproduce in another cohort or arise from different cellular sources. Single-cell analysis provides one route to separate bulk-tissue recurrence from cell-specific context [40-55]. Peripheral blood offers a further test of whether a tissue-derived signal has a detectable systemic component, but normal inter-individual variation and blood-cell composition can strongly influence such measurements [56-60].

This study addressed three questions. First, which transcriptional alterations were shared between OA cartilage and OC tissue, and were their directions concordant? Second, how consistently did these features recur in independent tissue cohorts and source-defined cellular populations? Third, did any tissue-concordant feature remain detectable in peripheral blood from both diseases? OA and OC were analysed separately at every stage, and evidence was integrated only after disease-specific estimation.

## Methods

### Study design

This exploratory observational study used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets. OA and OC were processed as separate disease tracks (Figure 1). The evidence sequence comprised tissue discovery, shared-gene classification, external tissue replication, functional analysis, cellular localization and conditional blood evaluation. Five genes were retained as representative molecular features after the full evidence chain was assembled. No clinical prediction model was developed.

### Tissue transcriptomic cohorts and preprocessing

Public datasets were obtained from the NCBI Gene Expression Omnibus (GEO) [78,79]. The OA discovery cohort was knee cartilage GSE114007, with 20 OA and 18 non-OA samples [8]. The OC discovery cohort was GSE18520, with 53 advanced high-grade serous tumours and 10 normal ovarian surface-epithelium samples [23]. External OA cohorts were GSE117999 (10 OA and 10 controls) and GSE82107 (10 OA and 7 controls). External OC cohorts were GSE54388 (16 tumours and 6 normal samples) and GSE12470 (43 serous carcinomas and 10 normal peritoneal samples) [24]. Dataset roles, platforms and group definitions are listed in Table S1.

Each expression matrix was checked for identifier consistency, duplicated features, missing or non-finite values and group balance. Repository metadata defined disease and reference groups. Multiple probes for one gene were collapsed according to the dataset-specific processing record. Normalization and modelling were performed within each cohort to avoid pooling measurements across platforms. Principal-component and sample-correlation displays used the 1,000 genes with the highest variance and served as quality-control summaries (Supplementary Figure 4). Bioconductor resources supported the analysis environment [80], and potential batch effects were considered during cohort-specific processing [73].

### Differential expression analysis

Differential expression was estimated separately in each cohort with limma [65,66]. The primary threshold was Benjamini-Hochberg false-discovery rate (FDR) <0.05 and absolute log2 fold change >=1 [67]. Sensitivity analysis combined FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Supplementary Figure 1; Tables S3a-b).

### Shared-gene direction classification

A gene entered the shared set only when it met the primary threshold in both discovery cohorts. Shared genes were divided into four sign classes: higher in both diseases, lower in both, higher in OA and lower in OC, or lower in OA and higher in OC. The first two classes were designated concordant and the latter two discordant. This classification was based on within-disease effects and did not require cross-platform expression-scale alignment.

### External tissue replication

Each external cohort was modelled independently with the same disease-versus-reference orientation as its corresponding discovery cohort. For every measurable shared gene, replication was summarized by sign agreement with the relevant discovery effect. An exact binomial test compared the agreement proportion with 0.5, and Spearman correlation described gene-wise concordance between discovery and external log2 fold changes. External FDR support and the intersection of FDR significance with sign agreement were also recorded. Expression values were not pooled across cohorts.

### Functional analysis

Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation analyses were applied to the 286 shared genes using clusterProfiler [74,75] and the corresponding knowledge resources [76,77]. The main display contains ten significant, non-redundant GO terms representing the major biological categories; complete results are provided in Table S12. All FDR-significant KEGG results are listed in Table S13.

Ranked Hallmark gene-set enrichment analysis was performed independently on the complete OA and OC discovery statistics [68,69]. A pathway was shared at this level when FDR <0.05 in both diseases. Paired normalized enrichment scores described concordant or discordant pathway states. The same analysis was applied to external tissue cohorts when a suitable ranked statistic was available. Enrichment results were interpreted as transcriptional associations.

### Single-cell analysis

Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29,30]. Count-level datasets underwent cell- and feature-level quality control, mitochondrial-content assessment, doublet handling when compatible sample-level counts were available, normalization and dimensional reduction. TPM-only input was not treated as raw counts. The workflow followed established single-cell analysis principles [40-49], with scDblFinder used for compatible doublet assessment and UCell retained for sample-aware score auditing [52,53]. Compatibility and quality-control outcomes are reported in Table S8a.

The main cellular analysis used the count-level OA atlas GSE255460 and OC atlas GSE154600. The datasets remained separate because they differed in tissue, platform, cellular composition and source annotation. Exact source labels were retained. UMAPs used deterministic display subsampling. For each representative gene, the fraction of cells with detected counts and the mean UMI count per cell were calculated within source labels. Detection fractions were interpreted within each atlas. All five embeddings are shown in Supplementary Figure 3, and detection summaries are reported in Table S8b.

### Peripheral blood evaluation of systemic molecular signals

Official GEO metadata were used to identify eligible disease-versus-control blood cohorts. The OA cohort was GSE48556, an Illumina HumanHT-12 V3.0 peripheral-blood mononuclear-cell dataset containing 106 OA cases and 33 healthy controls [15]. The OC cohort was GSE31682, an ABI Human Genome Survey Microarray V2 blood-cell-fraction dataset containing 48 epithelial ovarian cancer cases and 20 healthy controls [31]. The cohorts differed in platform and blood fraction and were therefore modelled separately (Table S9).

Official platform annotations were used for probe-to-gene mapping. Ambiguous mappings were discarded. When several unambiguous probes mapped to one gene, the probe with the highest interquartile range was retained. Limma estimated separate disease-control contrasts. Hedges g and 95% confidence intervals summarized standardized effects within each blood cohort; effects were not pooled across diseases.

The screen applied five sequential criteria: membership in the shared tissue set, concordant OA and OC tissue direction, measurement in both blood cohorts, the same sign in all four tissue and blood contrasts, and FDR <0.05 independently in each blood cohort. Denominators were retained for every stage. This sequence evaluated whether a tissue-derived signal also had a reproducible systemic component.

### Secondary network analysis

Signed weighted gene co-expression networks were fitted separately in the two discovery cohorts [70]. The primary disease-associated module was assessed by soft-power perturbation, 2,000 fixed-module bootstrap resamples and leave-one-sample-out estimates. This analysis described within-cohort co-expression stability.

The 286 shared genes were mapped to STRING v12.0 for Homo sapiens without adding neighbouring proteins [71]. The supplementary network retained high-confidence physical associations (score >=0.700). Mapped isolates remained in denominator-based summaries. Network topology was not used to select the representative genes.

### Statistical reporting and reproducibility

Tests were two-sided unless otherwise stated. Multiple-testing procedures and direction definitions were fixed before biological interpretation. Sample and cell units are reported for each evidence layer, and denominators accompany all filtering steps. Randomized procedures used seed 20260726. The project retains source manifests, parameters, figure source data, checksums and audit reports in accordance with FAIR principles [72]. The V4.0 manuscript and figures can be rebuilt with `run_submission_v40.ps1`.

## Results

### Shared transcriptomic alterations were identified in OA and OC tissues

OA and OC discovery tissues shared 286 differentially expressed genes at the primary threshold (Figure 2A-C; Table S2). The OA analysis identified 2,008 genes and the OC analysis identified 2,310. Shared-gene membership persisted across all six prespecified threshold combinations, although the total varied with the cut-off (Supplementary Figure 1; Tables S3a-b). Quality-control displays showed phenotype-associated structure together with within-group heterogeneity (Supplementary Figure 4).

### Shared genes converged on matrix, immune, stress and cell-cycle functions

The shared set was functionally concentrated in extracellular-matrix, immune, stress and cell-cycle processes (Figure 3; Tables S12-S14). GO terms included chromosome segregation, matrix organization, cytokine responses, hypoxia and integrated stress signalling. KEGG identified Cell cycle as the only FDR-significant pathway. Ten Hallmark pathways were significant in both discovery cohorts. Four had matching normalized enrichment-score signs, whereas six had opposite signs.

### Directional heterogeneity divided the shared set almost evenly

Directional heterogeneity was a defining property of the 286-gene overlap (Figure 2D). There were 112 genes higher in both diseases and 34 lower in both. A further 86 genes were higher in OA and lower in OC, while 54 were lower in OA and higher in OC. The resulting split was 146 concordant genes (51.0%) and 140 discordant genes (49.0%).

### External tissue replication was disease dependent

External replication was stronger in OC than in OA (Figure 4A; Tables S4-S6). GSE117999 reproduced the discovery sign for 143 of 280 measurable genes (51.1%; binomial P=0.765; Spearman rho=0.061). GSE82107 reproduced 178 of 286 signs (62.2%; P=4.15 x 10^-5; rho=0.219), although no individual shared gene reached external FDR <0.05. The two OC cohorts showed 260 of 286 (90.9%; P=1.01 x 10^-49; rho=0.847) and 179 of 226 (79.2%; P=2.47 x 10^-19; rho=0.628) matching signs. The corresponding numbers meeting both external FDR and sign criteria were 179 and 134.

G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 summarized complementary tissue-replication patterns (Figure 4B-C; Table S7). EFEMP1 and AKAP12 matched the discovery direction in all four external cohorts. G0S2, SOX9 and DDIT3 matched in three of four. These genes were carried forward as interpretable molecular features spanning concordant and discordant disease effects.

### Single-cell localization placed representative genes in distinct compartments

Representative genes localized to different source-defined populations in OA and OC atlases (Figure 5; Tables S8a-b). Across five datasets, 1,187,436 cells were audited and 1,025,361 passed dataset-specific quality-control criteria. The main display preserved six OA and six OC source labels from GSE255460 and GSE154600, respectively.

G0S2 had its highest detection fraction in ProC in the OA atlas (0.340) and Myeloid.cell in the OC atlas (0.409). EFEMP1 was most frequent in preHTC in OA (0.455) and Fibroblast in OC (0.179). AKAP12 was highest in preInfC in OA (0.165) and Fibroblast in OC (0.203). SOX9 was frequent in HTC in OA (0.783) and detected in Ovarian.cancer.cell in OC (0.072). DDIT3 was most frequent in HomC in OA (0.705) and Ovarian.cancer.cell in OC (0.187).

### Peripheral blood evaluation retained G0S2 as a candidate systemic signal

The prespecified blood screen reduced the tissue-derived set to one candidate systemic signal (Figure 6; Tables S9-S11). The sequence retained 146 tissue-concordant genes, 127 genes measured in both blood cohorts, 38 genes with the same sign in all four contrasts and three genes with nominal significance in both blood cohorts. Only G0S2 met FDR <0.05 independently in each blood cohort.

G0S2 was lower in OA tissue (log2 fold change -2.027; FDR=0.0028), OC tissue (-1.219; FDR=1.70 x 10^-5), OA peripheral-blood mononuclear cells (-0.115; FDR=0.0048) and the OC blood-cell fraction (-1.236; FDR=0.0102). Standardized blood effects were Hedges g=-0.78 (95% CI, -1.18 to -0.38) in OA and -0.80 (95% CI, -1.33 to -0.26) in OC.

### Secondary network analyses supported within-cohort structure

Secondary network analyses provided supporting context without altering the main evidence chain. The primary disease-associated WGCNA module-trait correlations were -0.951 in OA and -0.879 in OC, with stable signs across bootstrap and leave-one-out analyses (Supplementary Figure 1; Table S15). STRING mapped 275 of 286 shared genes and identified 62 high-confidence physical edges among 46 connected products (Supplementary Figure 5; Tables S16a-c).

## Discussion

The main finding is that OA and OC share a reproducible but heterogeneous transcriptional background. The 286-gene intersection persisted across threshold choices, yet its almost even concordant-discordant split shows that overlap did not correspond to one common expression state. External replication strengthened this interpretation. OC cohorts reproduced most discovery directions, whereas OA replication was weaker and more cohort dependent. Differences in tissue structure, disease severity and reference composition offer plausible explanations for this contrast [1-8,16-25,32-39].

The functional and cellular results indicate that shared genes can participate in recurring biological themes without carrying identical biological meaning. Matrix remodelling, immune responses, stress and cell-cycle regulation were represented in both diseases, but six jointly significant Hallmark pathways changed in opposite directions. The single-cell atlases added a cellular boundary: G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 were concentrated in different source-defined populations across OA and OC. This pattern is consistent with the broader observation that bulk transcriptional states reflect both cell-intrinsic regulation and cellular composition [9-14,26-30,40-55].

Most tissue-derived signals did not propagate through the blood screen. This attrition is biologically plausible because local tissue expression is diluted by blood-cell composition, systemic physiology and inter-individual variation [15,31,56-60]. G0S2 was the single feature that remained lower across both tissues and both blood cohorts after independent FDR control. Its known roles in lipid mobilization, quiescence and growth regulation provide a basis for follow-up, but they do not specify the source of the circulating transcriptional change [61-64]. The present evidence therefore positions G0S2 as a candidate systemic molecular signal for prospective evaluation.

Several limitations define the scope of these findings. The analysis used retrospective public datasets with incomplete harmonization of age, sex, disease stage, treatment, joint site and tumour histology. OA and OC tissue contrasts were anatomically distinct, and external OA replication was modest. Single-cell datasets differed in protocol and annotation, so detection fractions were interpreted only within each atlas. The blood cohorts used different platforms and cell fractions, and a second compatible OA-OC blood pair was unavailable. WGCNA lacked external module-preservation testing, and STRING associations do not establish disease-specific protein interactions. Clinical performance, tissue origin of the blood signal and protein-level regulation require independent prospective or experimental study.

## Conclusions

OA and OC share 286 tissue transcriptional alterations spanning matrix, immune, stress and cell-cycle functions. Their direction, external reproducibility and cellular location were strongly context dependent. G0S2 was the only feature retained by the independent dual-cohort blood screen. Together, these observations support partial molecular convergence across tissue, cellular and systemic contexts and define a focused candidate for further validation.

## Declarations

### Ethics approval and consent to participate

This study used de-identified public datasets. Ethical approval and informed consent were obtained by the original studies as described in their repository records and publications. No new participants were recruited.

### Consent for publication

Not applicable.

### Data availability

All datasets are available from the NCBI Gene Expression Omnibus under accessions GSE114007, GSE117999, GSE82107, GSE18520, GSE54388, GSE12470, GSE104782, GSE169454, GSE255460, GSE154600, GSE180661, GSE48556 and GSE31682. Derived tables, figure source data, parameter manifests, checksums and audit reports are included in the reproducible project. A public repository DOI or archived release URL should be inserted before submission.

### Code availability

The complete V4.0 workflow is available in the accompanying project and can be executed from `run_submission_v40.ps1`. The entry point records software versions and regenerates the figures, manuscript files and submission audit. A public repository URL and immutable release identifier should be added before submission.

### Competing interests

The authors declare no competing interests.

### Funding

[Funding information to be added.]

### Author contributions

[CRediT author-contribution statement to be added after the author list is finalized.]

### Acknowledgements

The authors thank the investigators and participants who generated and shared the public datasets used in this study.
"""


SUPPLEMENTARY_INDEX = r"""## Supplementary table index

- **Table S1:** Tissue cohorts, platforms, sample groups and analysis roles.
- **Table S2:** The 286 shared tissue differentially expressed genes with OA and OC effects and direction classes.
- **Tables S3a-b:** Differential-expression threshold sensitivity and membership.
- **Table S4:** External tissue gene-level effects.
- **Table S5:** External tissue direction-agreement summary.
- **Table S6:** External tissue Hallmark enrichment.
- **Table S7:** Five-gene descriptive evidence summary.
- **Table S8a:** Single-cell compatibility, quality control and permitted analysis layers.
- **Table S8b:** Representative-gene detection by exact source label.
- **Table S9:** Peripheral-blood cohort audit.
- **Table S10:** FDR-supported blood-replicated systemic signal.
- **Table S11:** Prespecified blood-screen attrition.
- **Table S12:** GO enrichment of shared tissue genes.
- **Table S13:** KEGG enrichment of shared tissue genes.
- **Table S14:** Discovery-cohort Hallmark direction matrix.
- **Table S15:** WGCNA soft-power, bootstrap and leave-one-out stability.
- **Tables S16a-c:** STRING mapping, physical associations and node topology.
"""


RESPONSE = r"""# V4.0 final-polish revision matrix

| V4.0 objective | Action | Evidence boundary |
|---|---|---|
| Fix the research position | Reframed the manuscript around reproducible shared molecular features across tissue, cellular and systemic contexts | Shared expression does not imply identical regulation |
| Simplify the abstract | Used a structured background-methods-results-conclusions format and retained only decisive numbers | No new analysis or unsupported quantitative claim |
| Strengthen the Introduction | Added a direct gap and three study questions, with an 80-reference verified bibliography | Prior work is used for context, not as study data |
| Reorder Methods | Aligned Methods with the evidence sequence from tissue discovery to blood evaluation and secondary networks | Each cohort remains independently modelled |
| Rebuild Results | Opened each subsection with its conclusion and reduced numerical repetition | All retained numbers trace to existing V3.4 tables |
| Focus the Discussion | Organized interpretation around partial overlap, functional/cellular context, G0S2 and study limits | Causal and clinical inferences remain outside scope |
| Reduce defensive language | Removed repeated classifier, signature and shared-mechanism disclaimers from the main narrative | Necessary boundaries are concentrated in Study design and Limitations |
| De-emphasize networks | Reduced WGCNA and STRING to one short Results subsection plus supplementary evidence | Networks do not select the five representative genes |
| Harmonize terminology | Locked shared molecular features, candidate systemic molecular signal and source-defined cell populations | Avoided biomarker panel, therapeutic target and common pathogenic mechanism wording |
| Standardize figures | Retained the validated 183-mm R figure system and revised Figure 4 language | No visual change alters data or uncertainty |
"""


def normalize(text: str) -> str:
    return re.sub(r"\n{3,}", "\n\n", text.strip()) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--references", required=True, type=Path)
    parser.add_argument("--figure-legends", required=True, type=Path)
    parser.add_argument("--template-output", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()

    references = args.references.read_text(encoding="utf-8").strip()
    reference_count = len(re.findall(r"^\d+\.\s", references, flags=re.MULTILINE))
    if reference_count != 80:
        raise RuntimeError(f"Expected 80 references, found {reference_count}")
    legends = args.figure_legends.read_text(encoding="utf-8").strip()
    legend_count = len(re.findall(r"^## (?:Supplementary )?Figure \d+\.", legends, flags=re.MULTILINE))
    if legend_count != 12:
        raise RuntimeError(f"Expected 12 figure legends, found {legend_count}")

    template = normalize(MANUSCRIPT + "\n\n## References\n\n" + references)
    manuscript = normalize(template + "\n\n## Figure legends\n\n" + legends + "\n\n" + SUPPLEMENTARY_INDEX)

    for path in (args.template_output, args.output, args.response_output):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.template_output.write_text(template, encoding="utf-8")
    args.output.write_text(manuscript, encoding="utf-8")
    args.response_output.write_text(normalize(RESPONSE), encoding="utf-8")
    print(f"Wrote {args.template_output}")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.response_output}")


if __name__ == "__main__":
    main()
