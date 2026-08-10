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

**Background.** Osteoarthritis (OA) and ovarian cancer (OC) arise in different tissues but both involve matrix remodelling, inflammatory signalling and cellular stress. It remains unclear whether their transcriptomic overlap persists across independent tissues, source-defined cell populations and peripheral blood.

**Methods.** OA and OC datasets were analysed on separate disease tracks. Discovery tissue transcriptomes defined shared differentially expressed genes and their direction of change. Independent tissue cohorts evaluated replication. Five illustrative genes were localized in separate single-cell atlases, and a prespecified blood screen tested whether tissue-concordant signals persisted in both diseases.

**Results.** The discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes. Of these, 146 changed in the same direction and 140 changed in opposite directions. Six of ten Hallmark pathways significant in both diseases also had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. The five illustrative genes showed different source-defined localization patterns in OA and OC atlases. Only G0S2 passed independent false-discovery-rate control in both blood cohorts and was lower in all four tissue and blood contrasts.

**Conclusions.** OA and OC share a limited transcriptional background whose direction, reproducibility and cellular localization depend on disease context. G0S2 is a candidate systemic molecular signal that requires prospective and protein-level confirmation.

**Keywords:** osteoarthritis; ovarian cancer; cross-disease transcriptomics; external validation; single-cell RNA sequencing; peripheral blood

## Introduction

Osteoarthritis is a whole-joint disorder characterized by cartilage loss, subchondral-bone remodelling, synovial change and low-grade inflammatory signalling [1-7]. Transcriptomic and single-cell studies have identified heterogeneous chondrocyte and stromal states in diseased cartilage [8-14]. Selected OA-associated expression changes are also detectable in blood, although circulating profiles depend on cohort design and cell composition [15].

Ovarian cancer is a heterogeneous malignancy shaped by genomic instability, molecular subtype, stromal remodelling and immune context [16-25]. Bulk and single-cell studies of high-grade serous ovarian cancer have resolved malignant, fibroblast, endothelial and immune populations across disease sites [26-30,44,45]. Blood-based measurements capture tumour-associated variation together with systemic and haematological responses [31].

Cross-disease transcriptomics can identify molecular features that recur across anatomically distant conditions [32-34]. Gene-list overlap alone does not show that two diseases share one transcriptional state. Shared genes may change in opposite directions, fail to replicate in another cohort or arise from different cellular sources. Single-cell analysis can distinguish bulk-tissue recurrence from cell-associated localization [35-38]. Peripheral blood provides a further test of systemic persistence, but normal inter-individual variation and blood-cell composition can obscure local tissue signals [46-50].

We therefore asked three questions. Which tissue transcriptional alterations were shared between OA and OC, and were their directions concordant? How consistently did these features recur in independent tissue cohorts and source-defined cell populations? Did any tissue-concordant feature remain detectable in peripheral blood from both diseases? Each disease was estimated separately, and evidence was integrated only after disease-specific analysis.

## Methods

### Study design

This exploratory observational study used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets. OA and OC were processed as separate disease tracks (Figure 1). The evidence sequence comprised tissue discovery, direction classification, functional analysis, external tissue replication, cellular localization and blood evaluation. Five genes were chosen after these analyses to illustrate distinct evidence roles. They were not used to define the shared set or the blood-screen criteria.

### Tissue transcriptomic cohorts and preprocessing

Datasets were obtained from the NCBI Gene Expression Omnibus [68,69]. The OA discovery cohort was knee cartilage GSE114007, with 20 OA and 18 non-OA samples [8]. The OC discovery cohort was GSE18520, with 53 advanced high-grade serous tumours and 10 normal ovarian surface-epithelium samples [23]. External OA cohorts were GSE117999 (10 OA and 10 controls) and GSE82107 (10 OA and 7 controls). External OC cohorts were GSE54388 (16 tumours and 6 normal samples) and GSE12470 (43 serous carcinomas and 10 normal peritoneal samples) [24]. Cohort roles, platforms and group definitions are listed in Table S1.

Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Multiple probes for one gene were collapsed according to the dataset-specific record. Normalization and modelling were performed within each cohort, with batch structure considered during processing [63,70]. Principal-component and sample-correlation summaries are provided in Supplementary Figure 4.

### Differential expression and direction classification

Differential expression was estimated separately with limma [55,56]. The primary threshold was Benjamini-Hochberg false-discovery rate (FDR) <0.05 and absolute log2 fold change >=1 [57]. Sensitivity analysis combined FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Supplementary Figure 1; Tables S3a-b).

A gene entered the shared set only when it met the primary threshold in both discovery cohorts. Shared genes were classified as higher in both diseases, lower in both, higher in OA and lower in OC, or lower in OA and higher in OC. The first two classes were concordant and the latter two discordant. Classification used within-disease effects and did not align expression scales across platforms.

### External tissue replication

Each external cohort was modelled independently with the disease-versus-reference orientation of its discovery cohort. Replication was summarized by sign agreement for every measurable shared gene. An exact binomial test compared agreement with 0.5, and Spearman correlation described gene-wise concordance between discovery and external log2 fold changes. External FDR support was recorded without pooling expression values across cohorts.

### Functional analysis

Gene Ontology and Kyoto Encyclopedia of Genes and Genomes over-representation analyses were applied to the 286 shared genes with clusterProfiler [64-67]. Figure 3 shows ten significant, non-redundant terms representing the principal categories; complete results are provided in Tables S12-S13.

Ranked Hallmark gene-set enrichment analysis was performed independently on the complete OA and OC discovery statistics [58,59]. A pathway was shared when FDR <0.05 in both diseases. Paired normalized enrichment scores described matching or opposite pathway directions. External tissue cohorts were analysed when a suitable ranked statistic was available (Tables S6 and S14).

### Single-cell localization

Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29,30]. Count-level datasets underwent cell- and feature-level quality control, mitochondrial-content assessment, compatible doublet handling, normalization and dimensional reduction. TPM-only input was not treated as raw counts. The workflow followed established single-cell analysis guidance [35-39]. Doublet methods were reviewed for compatibility, and scDblFinder was used when sample-level counts permitted [40-42]. UCell supported sample-aware score auditing [43]. Dataset eligibility and quality-control outcomes are reported in Table S8a.

The main analysis used count-level OA GSE255460 and OC GSE154600. The atlases remained separate because their tissues, platforms, cellular composition and annotations differed. Exact source labels were retained. For each illustrative gene, detection fraction and mean UMI count were calculated within source labels. Fractions were interpreted within each atlas, not as direct OA-OC comparisons. All five embeddings and gene-level summaries are provided in Supplementary Figure 3 and Table S8b.

### Peripheral blood evaluation

Eligible disease-versus-control blood cohorts were identified from official repository metadata. OA GSE48556 contained peripheral-blood mononuclear cells from 106 cases and 33 healthy controls [15]. OC GSE31682 contained a blood-cell fraction from 48 cases and 20 healthy controls [31]. The platforms and blood fractions differed, so cohorts were modelled separately (Table S9).

Official platform annotations were used for probe-to-gene mapping. Ambiguous mappings were discarded, and the probe with the highest interquartile range represented genes with multiple unambiguous probes. Limma estimated separate contrasts. Hedges g and 95% confidence intervals summarized within-cohort effects; effects were not pooled across diseases.

The prespecified screen required membership in the shared tissue set, concordant tissue direction, measurement in both blood cohorts, the same sign in all four contrasts and FDR <0.05 independently in each blood cohort. Denominators were retained at every step (Tables S9-S11).

### Supportive analyses and statistical reporting

Co-expression stability and high-confidence STRING associations were evaluated as supportive analyses and are reported only in Supplementary Figure 1, Supplementary Figure 5 and Tables S15-S16 [60,61]. These analyses did not select the illustrative genes.

Tests were two-sided unless otherwise stated. Multiple-testing procedures and direction definitions were fixed before biological interpretation. Sample and cell units are reported for each evidence layer, and denominators accompany filtering steps. Randomized procedures used seed 20260726. Full parameters, source manifests and executable scripts are included in the reproducible project.

## Results

### Tissue discovery identified 286 shared alterations with balanced direction classes

The OA discovery analysis identified 2,008 differentially expressed genes and the OC analysis identified 2,310. Their intersection contained 286 genes (Figure 2A-C; Table S2). Shared-gene membership persisted across all six prespecified threshold combinations, although the total changed with the cut-off (Supplementary Figure 1; Tables S3a-b).

Direction divided the overlap almost evenly (Figure 2D). There were 112 genes higher in both diseases and 34 lower in both. Another 86 genes were higher in OA and lower in OC, while 54 were lower in OA and higher in OC. The shared set therefore contained 146 concordant genes (51.0%) and 140 discordant genes (49.0%).

### Shared genes mapped to recurring functions but not uniform pathway states

The shared set was enriched for extracellular-matrix, immune, stress and cell-cycle processes (Figure 3; Tables S12-S14). Gene Ontology terms included chromosome segregation, matrix organization, cytokine responses, hypoxia and integrated stress signalling. Cell cycle was the only FDR-significant Kyoto Encyclopedia of Genes and Genomes pathway. Ten Hallmark pathways were significant in both discovery cohorts. Four had matching normalized enrichment-score signs, whereas six had opposite signs.

### External direction replicated more strongly in OC cohorts

External direction agreement was higher in OC than in OA (Figure 4A; Tables S4-S6). GSE117999 reproduced 143 of 280 measurable OA signs (51.1%; binomial P=0.765; Spearman rho=0.061). GSE82107 reproduced 178 of 286 signs (62.2%; P=4.15 x 10^-5; rho=0.219), although no shared gene reached external FDR <0.05. The OC cohorts reproduced 260 of 286 signs (90.9%; P=1.01 x 10^-49; rho=0.847) and 179 of 226 signs (79.2%; P=2.47 x 10^-19; rho=0.628). The corresponding numbers meeting both external FDR and sign criteria were 179 and 134.

G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 were used as illustrative genes because they captured different evidence patterns (Figure 4B-C; Table S7). EFEMP1 and AKAP12 matched the discovery direction in all four external cohorts. G0S2, SOX9 and DDIT3 matched in three of four. The set also spanned concordant and discordant disease effects, source-defined cellular localization and the blood-screen outcome.

### Illustrative genes showed source-defined cellular localization patterns

Across the five audited single-cell datasets, 1,187,436 cells were assessed and 1,025,361 passed dataset-specific quality-control criteria. Figure 5 shows six exact OA and six exact OC source labels from GSE255460 and GSE154600; all dataset embeddings are provided in Supplementary Figure 3.

The highest detection fraction for each illustrative gene occurred in different source labels across the two atlases (Figure 5; Table S8b). G0S2 was highest in OA ProC (0.340) and OC Myeloid.cell (0.409). EFEMP1 was highest in OA preHTC (0.455) and OC Fibroblast (0.179), while AKAP12 was highest in OA preInfC (0.165) and OC Fibroblast (0.203). SOX9 was frequent in OA HTC (0.783) and detected in OC Ovarian.cancer.cell (0.072). DDIT3 was highest in OA HomC (0.705) and OC Ovarian.cancer.cell (0.187). These values describe localization within each source-defined atlas and do not support functional inference.

### Peripheral blood evaluation retained only G0S2

The blood screen reduced the tissue-derived set to one candidate systemic molecular signal (Figure 6; Tables S9-S11). The sequence retained 146 tissue-concordant genes, 127 genes measured in both blood cohorts, 38 genes with the same sign in all four contrasts and three genes with nominal significance in both blood cohorts. Only G0S2 met FDR <0.05 independently in each cohort.

G0S2 was lower in OA tissue (log2 fold change -2.027; FDR=0.0028), OC tissue (-1.219; FDR=1.70 x 10^-5), OA peripheral-blood mononuclear cells (-0.115; FDR=0.0048) and the OC blood-cell fraction (-1.236; FDR=0.0102). Standardized blood effects were Hedges g=-0.78 (95% CI, -1.18 to -0.38) in OA and -0.80 (95% CI, -1.33 to -0.26) in OC. The result was consistent across direction and FDR criteria despite the different blood fractions.

## Discussion

OA and OC shared a limited but heterogeneous transcriptional background. The 286-gene overlap persisted across threshold choices, yet its nearly equal concordant and discordant components did not define one common expression state. External replication reinforced this boundary. OC cohorts reproduced most discovery directions, whereas OA replication was weaker and more cohort dependent. Tissue composition, disease severity and reference-group differences may contribute to this contrast [2,3,8,16,18,23,24,32-34].

Recurring functions also retained disease and cellular context. Matrix, immune, stress and cell-cycle categories were represented in both diseases, but six jointly significant Hallmark pathways had opposite directions. The five illustrative genes localized to different source-defined populations across separate OA and OC atlases. This observation supports cell-associated localization differences, not different gene functions or homologous cell-state correspondence [9-14,26-30,35,36,38,44,45].

Most tissue-concordant genes did not persist through the blood screen, consistent with dilution by blood-cell composition, systemic physiology and inter-individual variation [15,31,46-50]. G0S2 alone remained lower in both tissues and both blood cohorts after independent FDR control. Its reported roles in lipolysis, quiescence and growth regulation motivate follow-up but do not identify the source or function of the circulating change [51-54]. G0S2 should therefore be treated as a candidate systemic molecular signal rather than a clinical marker.

The study has several limitations. Public cohorts differed in age, sex, disease stage, treatment, joint site, tumour histology and reference tissue, and these variables could not be harmonized completely. External OA replication was modest. Single-cell protocols and annotations differed, so detection fractions were interpreted only within each atlas. The blood cohorts used different platforms and cell fractions, and no second compatible OA-OC blood pair was available. The analysis was retrospective and transcriptomic; tissue origin, protein-level regulation and prospective reproducibility remain unresolved.

## Conclusions

OA and OC shared 286 tissue transcriptional alterations, but their direction, external reproducibility and cellular localization were context dependent. Only G0S2 persisted through the independent dual-cohort blood screen. These findings define partial molecular convergence across tissue, cellular and systemic contexts while setting clear limits on mechanistic and translational interpretation.

## Declarations

### Ethics approval and consent to participate

This study used de-identified public datasets. Ethical approval and informed consent were obtained by the original studies as described in their repository records and publications. No new participants were recruited.

### Consent for publication

Not applicable.

### Data availability

All datasets are available from the NCBI Gene Expression Omnibus under accessions GSE114007, GSE117999, GSE82107, GSE18520, GSE54388, GSE12470, GSE104782, GSE169454, GSE255460, GSE154600, GSE180661, GSE48556 and GSE31682. Derived tables, figure source data, parameters and audit reports are included in the reproducible project to support FAIR reuse [62]. A public repository DOI or archived release URL should be inserted before submission.

### Code availability

The complete V4.1 workflow is included in the accompanying project. A public repository URL and immutable release identifier should be added before submission.

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
- **Table S7:** Five-gene illustrative evidence summary.
- **Table S8a:** Single-cell compatibility, quality control and permitted analysis layers.
- **Table S8b:** Illustrative-gene detection by exact source label.
- **Table S9:** Peripheral-blood cohort audit.
- **Table S10:** FDR-supported blood-replicated systemic signal.
- **Table S11:** Prespecified blood-screen attrition.
- **Table S12:** Gene Ontology enrichment of shared tissue genes.
- **Table S13:** Kyoto Encyclopedia of Genes and Genomes enrichment of shared tissue genes.
- **Table S14:** Discovery-cohort Hallmark direction matrix.
- **Table S15:** WGCNA soft-power, bootstrap and leave-one-out stability.
- **Tables S16a-c:** STRING mapping, physical associations and node topology.
"""


REVISION_MATRIX = r"""# V4.1 strategic-compression revision matrix

| Review concern | V4.1 action | Resolution test |
|---|---|---|
| Analytical breadth diluted the main story | Removed WGCNA and STRING results from the main Results and reduced Methods to one supportive sentence | Main evidence chain reads tissue overlap to direction, external replication, cellular localization and blood persistence |
| Five-gene set resembled a selected panel | Renamed the set illustrative genes and stated its post-analysis descriptive role | No optimized, predictive or clinical claim is attached to the set |
| Introduction and single-cell citations were dense | Reduced the verified bibliography from 80 to 70 and redistributed citations by claim | Every retained reference is cited and no sentence carries a broad citation block |
| Methods mixed scientific procedure with project-management detail | Removed checksum, FAIR and executable-name detail from Methods | Reproducibility remains documented in Data and Code availability |
| Single-cell interpretation exceeded localization evidence | Replaced functional language with source-defined localization language | Results explicitly state that localization does not establish different functions |
| External validation needed greater prominence | Preserved all cohort denominators, agreement rates, correlations and FDR-supported counts | Figure 4 and the corresponding Results subsection remain central |
| G0S2 needed a clearer climax and boundary | Retained four-context effects and two Hedges g estimates, while preserving candidate wording | Blood result is prominent without becoming a clinical claim |
| Integrated model was not sequential | Rebuilt Figure 7 as a linear evidence sequence | Diagram ends with one G0S2 candidate and a bounded interpretation |
"""


def normalize(text: str) -> str:
    return re.sub(r"\n{3,}", "\n\n", text.strip()) + "\n"


def cited_numbers(text: str) -> set[int]:
    result: set[int] = set()
    for match in re.finditer(r"\[([0-9,\-]+)\]", text):
        for part in match.group(1).split(","):
            if "-" in part:
                start, end = (int(value) for value in part.split("-", 1))
                result.update(range(start, end + 1))
            else:
                result.add(int(part))
    return result


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
    if reference_count != 70:
        raise RuntimeError(f"Expected 70 references, found {reference_count}")
    cited = cited_numbers(MANUSCRIPT)
    if cited != set(range(1, 71)):
        raise RuntimeError(f"Citation coverage mismatch: missing={sorted(set(range(1,71))-cited)} extra={sorted(cited-set(range(1,71)))}")

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
    args.response_output.write_text(normalize(REVISION_MATRIX), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
