from __future__ import annotations

import argparse
import re
from pathlib import Path


def expand_citation(content: str) -> list[int]:
    values: list[int] = []
    for part in content.split(","):
        if "-" in part:
            start, end = (int(value) for value in part.split("-", 1))
            values.extend(range(start, end + 1))
        else:
            values.append(int(part))
    return values


def compress_numbers(numbers: list[int]) -> str:
    ordered = sorted(set(numbers))
    if not ordered:
        return ""
    ranges: list[str] = []
    start = previous = ordered[0]
    for number in ordered[1:]:
        if number == previous + 1:
            previous = number
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = number
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def remap_citations(text: str) -> str:
    def replace(match: re.Match[str]) -> str:
        old = [value for value in expand_citation(match.group(1)) if value != 60]
        mapped = [value if value < 60 else value - 1 for value in old]
        if not mapped:
            return ""
        return f"[{compress_numbers(mapped)}]"

    return re.sub(r"\[([0-9,\-]+)\]", replace, text)


def transform_references(text: str) -> str:
    corrections = {
        1: "Glyn-Jones S, Palmer AJR, Agricola R, Price AJ, Vincent TL, Weinans H, et al. Osteoarthritis. The Lancet. 2015;386(9991):376-387. doi:10.1016/s0140-6736(14)60802-3.",
        10: "Fu W, Hettinghouse A, Chen Y, Hu W, Ding X, Chen M, et al. 14-3-3 epsilon is an intracellular component of the TNFR2 receptor complex and its activation protects against osteoarthritis. Annals of the Rheumatic Diseases. 2021;80(12):1615-1627. doi:10.1136/annrheumdis-2021-220000.",
        19: "Cancer Genome Atlas Research Network. Integrated genomic analyses of ovarian carcinoma. Nature. 2011;474(7353):609-615. doi:10.1038/nature10166.",
        22: "Konecny GE, Wang C, Hamidi H, Winterhoff B, Kalli KR, Dering J, et al. Prognostic and therapeutic relevance of molecular subtypes in high-grade serous ovarian cancer. JNCI: Journal of the National Cancer Institute. 2014;106(10):dju249. doi:10.1093/jnci/dju249.",
        32: "Dudley JT, Tibshirani R, Deshpande T, Butte AJ. Disease signatures are robust across tissues and experiments. Molecular Systems Biology. 2009;5:307. doi:10.1038/msb.2009.66.",
        36: "Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular Systems Biology. 2019;15(6):e8746. doi:10.15252/msb.20188746.",
        48: "Pepe MS, Etzioni R, Feng Z, Potter JD, Thompson ML, Thornquist M, et al. Phases of biomarker development for early detection of cancer. JNCI: Journal of the National Cancer Institute. 2001;93(14):1054-1061. doi:10.1093/jnci/93.14.1054.",
        55: "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research. 2015;43(7):e47. doi:10.1093/nar/gkv007.",
        57: "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B (Methodological). 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.",
        67: "Barrett T, Wilhite SE, Ledoux P, Evangelista C, Kim IF, Tomashevsky M, et al. NCBI GEO: archive for functional genomics data sets-update. Nucleic Acids Research. 2013;41(D1):D991-D995. doi:10.1093/nar/gks1193.",
    }
    records: list[str] = []
    for line in text.splitlines():
        match = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if not match:
            continue
        number = int(match.group(1))
        if number == 60:
            continue
        new_number = number if number < 60 else number - 1
        record = match.group(2)
        record = record.replace("B枚hringer", "Böhringer").replace("B枚hm", "Böhm")
        record = corrections.get(new_number, record)
        record = re.sub(
            r"doi:(10\.\S+?)(?=\.$)",
            lambda doi_match: "doi:" + doi_match.group(1).lower(),
            record,
        )
        records.append(f"{new_number}. {record}")
    if len(records) != 69:
        raise RuntimeError(f"Expected 69 BMC references, found {len(records)}")
    return "\n".join(records)


def ensure_cross_references(text: str) -> str:
    text = text.replace(
        "Eligible external cohorts were analysed similarly (Additional file 1: Tables S6 and S14).",
        "Eligible external cohorts were analysed similarly (Additional file 2: Figure S2; Additional file 1: Tables S6 and S14).",
    )
    text = text.replace(
        "The similar directions and standardized effects across the two blood datasets made G0S2 the only blood-persistent feature that met the prespecified rule.",
        "The similar directions and standardized effects across the two blood datasets made G0S2 the only blood-persistent feature that met the prespecified rule. The multi-layer evidence sequence is summarized in Figure 7.",
    )
    return text


def renumber_by_first_appearance(body: str, references: str) -> tuple[str, str]:
    records: dict[int, str] = {}
    for line in references.splitlines():
        match = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if match:
            records[int(match.group(1))] = match.group(2)

    order: list[int] = []
    for match in re.finditer(r"\[([0-9,\-]+)\]", body):
        for number in expand_citation(match.group(1)):
            if number not in order:
                order.append(number)
    expected = set(records)
    if set(order) != expected or len(order) != len(records):
        raise RuntimeError(
            f"Cannot renumber references by first appearance: order={len(order)} records={len(records)}"
        )
    mapping = {old: new for new, old in enumerate(order, start=1)}

    def replace(match: re.Match[str]) -> str:
        mapped = [mapping[number] for number in expand_citation(match.group(1))]
        return f"[{compress_numbers(mapped)}]"

    body = re.sub(r"\[([0-9,\-]+)\]", replace, body)
    references = "\n".join(
        f"{new}. {records[old]}" for new, old in enumerate(order, start=1)
    )
    return body, references


def replace_required(text: str, old: str, new: str) -> str:
    if old not in text:
        if new in text:
            return text
        raise RuntimeError(f"Required polishing source text was not found: {old[:80]}")
    return text.replace(old, new)


def polish_bmc_narrative(text: str) -> str:
    replacements = [
        (
            "# Directionally heterogeneous transcriptomic overlap between osteoarthritis and ovarian cancer across tissue, cellular and blood contexts",
            "# Shared and context-dependent transcriptomic features between osteoarthritis and ovarian cancer across tissue, single-cell and blood datasets",
        ),
        (
            "**Running title:** Directionally heterogeneous OA-OC overlap",
            "**Running title:** Context-dependent OA-OC transcriptomic overlap",
        ),
        (
            "**Background.** Osteoarthritis (OA) and ovarian cancer (OC) arise in different tissues but both involve matrix remodelling, inflammatory signalling and cellular stress. It remains unclear whether their transcriptomic overlap defines a reproducible cross-disease state across independent tissues, source-defined cell populations and peripheral blood.",
            "**Background.** Osteoarthritis (OA) and ovarian cancer (OC) arise in distinct tissues, yet both engage matrix remodelling, inflammatory signalling and cellular stress. The extent to which their transcriptomic overlap reflects conserved or context-dependent signals across tissues, source-defined cell populations and peripheral blood remains unresolved.",
        ),
        (
            "**Methods.** OA and OC datasets were analysed on separate disease tracks. Discovery tissue transcriptomes defined shared differentially expressed genes and their direction of change. Independent tissue cohorts evaluated replication. Five illustrative genes were localized in separate single-cell atlases, and a prespecified blood screen tested whether tissue-concordant signals persisted in both diseases.",
            "**Methods.** OA and OC datasets were analysed on separate disease tracks. Discovery tissue transcriptomes defined shared differentially expressed genes and their direction of change. Independent tissue cohorts evaluated replication, and separate single-cell atlases localized five illustrative genes. A prespecified blood screen then tested whether tissue-concordant signals persisted in both diseases.",
        ),
        (
            "**Results.** The discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes. Of these, 146 changed in the same direction and 140 changed in opposite directions. Six of ten Hallmark pathways significant in both diseases also had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. The five illustrative genes showed different source-defined localization patterns in OA and OC atlases. Only G0S2 passed independent false-discovery-rate control in both blood cohorts and was lower in all four tissue and blood contrasts.",
            "**Results.** The discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes. Of these, 146 changed in the same direction and 140 changed in opposite directions. Six of ten Hallmark pathways significant in both diseases also had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts, revealing disease-specific stability. The five illustrative genes occupied different source-defined contexts in the OA and OC atlases. Under independent false-discovery-rate control, G0S2 was the sole blood-persistent feature and was lower in all four tissue and blood contrasts.",
        ),
        (
            "**Conclusions.** OA and OC share a directionally heterogeneous transcriptomic overlap rather than a uniform disease state. Replication and cellular localization remain context dependent. G0S2 is a blood-persistent transcriptomic association that requires prospective and protein-level confirmation.",
            "**Conclusions.** Multi-context integration distinguished a restricted conserved component from broader context-dependent transcriptomic overlap between OA and OC. Nearly balanced gene and pathway directions, together with disease-specific replication and cellular localization, defined the boundary of this convergence. G0S2 emerged as a blood-persistent association for prospective cell-resolved and protein-level validation.",
        ),
        (
            "Cross-disease transcriptomics can identify molecular features that recur across anatomically distant conditions [32-34]. Gene-list overlap alone does not show that two diseases share one transcriptional state. Shared genes may change in opposite directions, fail to replicate in another cohort or arise from different cellular sources. Single-cell analysis can distinguish bulk-tissue recurrence from cell-associated localization [35-38]. Peripheral blood provides a further test of systemic persistence, but normal inter-individual variation and blood-cell composition can obscure local tissue signals [46-50].",
            "Cross-disease transcriptomics can quantify how molecular features recur across anatomically distant conditions [32-34]. Three dimensions determine whether an overlap is conserved: direction of change, reproducibility across cohorts and cellular context. Single-cell analysis localizes bulk-tissue signals to source-defined populations [35-38], while peripheral blood tests whether a tissue-concordant component also appears in circulation [46-50]. Together, these layers distinguish broadly recurrent signals from tissue- and context-dependent responses.",
        ),
        (
            "We therefore asked three questions. Which tissue transcriptional alterations were shared between OA and OC, and were their directions concordant? How consistently did these features recur in independent tissue cohorts and source-defined cell populations? Did any tissue-concordant feature remain detectable in peripheral blood from both diseases? Each disease was estimated separately, and evidence was integrated only after disease-specific analysis. This design treated overlap as a hypothesis to be tested rather than evidence of common pathogenesis.",
            "We characterized tissue-level overlap between OA and OC, classified each feature by direction, evaluated external replication and localized illustrative genes in separate single-cell atlases. We then applied a prespecified dual-cohort blood screen. Each disease was analysed independently before evidence integration. This design distinguished reproducibly conserved signals from overlap that depended on tissue, cohort or cellular context.",
        ),
        (
            "Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Public processed matrices were log2-transformed only when indicated by their scale, and multiple probes for one gene were collapsed according to the dataset-specific record. Reliable technical batch covariates were not consistently available. No explicit batch-correction method was applied, and cohorts were not merged or adjusted with ComBat. Each cohort was modelled separately with disease group as the only design factor. Principal-component and sample-correlation summaries assessed cohort structure but did not trigger outcome-informed sample removal or adjustment (Additional file 2: Figure S4) [62].",
            "Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Public processed matrices were log2-transformed only when indicated by their scale, and multiple probes for one gene were collapsed according to the dataset-specific record. Because reliable technical batch covariates were unavailable and cohorts differed in platform and tissue source, each cohort was modelled separately. Cross-cohort ComBat adjustment and expression-level pooling were deliberately avoided. Disease group was the sole design factor within each cohort. Principal-component and sample-correlation summaries assessed cohort structure without outcome-informed sample removal or adjustment (Additional file 2: Figure S4) [62].",
        ),
        (
            "Classification used within-disease effects and did not align expression scales across platforms.",
            "Classification preserved within-disease effect scales rather than aligning expression values across platforms.",
        ),
        (
            "Fractions were interpreted within each atlas, not as direct OA-OC comparisons.",
            "Fractions were summarized within each atlas to preserve its source-defined cellular composition.",
        ),
        (
            "Thus, gene-list overlap did not define a uniform cross-disease expression state.",
            "This near-even partition showed that tissue-level recurrence contained both conserved and disease-divergent components.",
        ),
        (
            "### Shared genes mapped to recurring functions but not uniform pathway states",
            "### Shared genes converged on recurring functions while pathway direction remained heterogeneous",
        ),
        (
            "### External direction replicated more strongly in OC cohorts",
            "### External replication revealed disease-specific transcriptomic stability",
        ),
        (
            "The highest detection fraction for each illustrative gene occurred in different source labels across the two atlases (Figure 5; Additional file 1: Table S8b). G0S2 was highest in OA ProC (0.340) and OC Myeloid.cell (0.409). EFEMP1 was highest in OA preHTC (0.455) and OC Fibroblast (0.179), while AKAP12 was highest in OA preInfC (0.165) and OC Fibroblast (0.203). SOX9 was frequent in OA HTC (0.783) and detected in OC Ovarian.cancer.cell (0.072). DDIT3 was highest in OA HomC (0.705) and OC Ovarian.cancer.cell (0.187). These values describe localization within each source-defined atlas and do not support functional inference.",
            "The highest detection fraction for each illustrative gene occurred in different source labels across the two atlases (Figure 5; Additional file 1: Table S8b). G0S2 was highest in OA ProC (0.340) and OC Myeloid.cell (0.409). EFEMP1 was highest in OA preHTC (0.455) and OC Fibroblast (0.179), while AKAP12 was highest in OA preInfC (0.165) and OC Fibroblast (0.203). SOX9 was frequent in OA HTC (0.783) and detected in OC Ovarian.cancer.cell (0.072). DDIT3 was highest in OA HomC (0.705) and OC Ovarian.cancer.cell (0.187). These atlas-specific distributions placed each feature within a distinct source-defined context across OA and OC.",
        ),
        (
            "### Peripheral-blood evaluation retained only G0S2 under the dual-cohort FDR rule",
            "### Stringent blood validation identified G0S2 as a restricted systemic component",
        ),
        (
            "The blood screen reduced the tissue-derived set to one blood-persistent transcriptomic association (Figure 6; Additional file 1: Tables S9-S11). The sequence retained 146 tissue-concordant genes, 127 genes measured in both blood cohorts, 38 genes with the same sign in all four contrasts and three genes with nominal significance in both blood cohorts. Only G0S2 met FDR <0.05 independently in each cohort.",
            "The prespecified blood screen progressively resolved the tissue-derived set to one blood-persistent transcriptomic association (Figure 6; Additional file 1: Tables S9-S11). The sequence retained 146 tissue-concordant genes, 127 genes measured in both blood cohorts, 38 genes with the same sign in all four contrasts and three genes with nominal significance in both blood cohorts. G0S2 was the sole gene meeting FDR <0.05 independently in each cohort.",
        ),
        (
            "G0S2 was lower in OA tissue (log2 fold change -2.027; FDR=0.0028), OC tissue (-1.219; FDR=1.70 x 10^-5), OA peripheral-blood mononuclear cells (-0.115; FDR=0.0048) and the OC blood-cell fraction (-1.236; FDR=0.0102).",
            "G0S2 was lower in OA tissue (log2 fold change -2.027; FDR=0.0028) and OC tissue (-1.219; FDR=1.70 x 10^-5). It was also lower in OA peripheral-blood mononuclear cells (-0.115; FDR=0.0048) and the OC blood-cell fraction (-1.236; FDR=0.0102).",
        ),
        (
            "The result met the prespecified direction and FDR criteria despite the different blood fractions, but it did not establish a shared cellular source.",
            "The concordant direction and comparable standardized effects therefore identified G0S2 as the blood-persistent component of the tissue-derived set across the two blood fractions.",
        ),
        (
            "OA and OC showed a limited transcriptomic overlap rather than a uniform shared disease state. The 286-gene intersection was stable to prespecified thresholds, but concordant and discordant directions were nearly balanced. External replication further constrained the result. OC cohorts reproduced most discovery directions, whereas OA replication was weaker and more cohort dependent. Differences in cartilage source, joint site, disease severity, platform and reference groups may contribute to this contrast [2-3,8,16,18,23-24,32-34].",
            "This study defines the extent and boundary of transcriptomic overlap between OA and OC. The 286-gene intersection was robust to prespecified thresholds, yet its near-even directional partition showed that recurrence included both conserved and disease-divergent components. External cohorts further separated a strongly conserved OC component from a more heterogeneous OA component. Variation in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
        ),
        (
            "Recurring matrix, immune, stress and cell-cycle categories therefore should not be interpreted as common pathogenesis. Six jointly significant Hallmark pathways had opposite enrichment directions, and the illustrative genes localized to different source-defined populations in separate OA and OC atlases. These observations identify contextual heterogeneity but do not establish different gene functions or homologous cell states [9-14,26-30,35-36,38,44-45].",
            "Functional enrichment and single-cell localization clarified why shared membership was not equivalent to uniform biology. Matrix, immune, stress and cell-cycle categories recurred across diseases, but six jointly significant Hallmark pathways had opposite enrichment directions. Illustrative genes also localized to different source-defined populations in the separate OA and OC atlases. Direction and cellular context are therefore integral properties of the observed convergence [9-14,26-30,35-36,38,44-45].",
        ),
        (
            "The blood screen retained G0S2 because it was lower in both tissues, had the same sign in both blood cohorts and passed FDR control independently in each cohort. The other same-sign genes did not meet this complete rule. This attrition may reflect biological heterogeneity, blood-cell composition, measurement sensitivity or cohort differences. G0S2 was most frequently detected in OA ProC and OC Myeloid.cell populations, but neither atlas identifies the cells responsible for the blood association. Reported roles in lipolysis, quiescence and growth regulation provide biological context, not a mechanism for the four observed contrasts [51-54]. G0S2 is therefore a blood-persistent transcriptomic association that requires cell-resolved, protein-level and prospective confirmation.",
            "Stringent systemic validation reduced the broad tissue overlap to a restricted circulating component. G0S2 was lower in both tissues, showed the same direction in both blood cohorts and passed FDR control independently in each cohort. The attrition from 286 shared tissue genes to one dual-FDR result indicates that circulating conservation was substantially more selective than tissue recurrence. G0S2 was most frequently detected in OA ProC and OC Myeloid.cell populations, while reported roles in lipolysis, quiescence and growth regulation provide a biological basis for targeted follow-up [51-54]. Cell-resolved, protein-level and prospective studies are now needed to define its source and functional relevance.",
        ),
        (
            "The study is limited by retrospective public data and incomplete covariate annotation. Reliable technical batch covariates were not consistently available. No explicit batch-correction method was applied, and cohorts were not merged or adjusted with ComBat. Age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be harmonized completely. Single-cell protocols and annotations differed, and the two blood cohorts used different platforms and cell fractions. No second compatible OA-OC blood pair was available. These constraints limit mechanistic and translational interpretation.",
            "Several design features define the scope of these findings. The study used retrospective public data with incomplete covariate annotation, and age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be fully harmonized. Cohorts were therefore analysed separately without cross-cohort ComBat adjustment. Single-cell protocols and annotations differed, while the two blood cohorts used different platforms and cell fractions. These factors may contribute to the observed heterogeneity. Prospective tissue-blood pairing, independent blood replication and protein measurements represent the next validation steps.",
        ),
        (
            "OA and OC shared 286 tissue transcriptional alterations, but the near-balanced direction classes did not define a uniform shared disease state. External reproducibility and cellular localization remained context dependent. Only G0S2 passed the independent dual-cohort blood rule. This result is a blood-persistent association, not evidence of common pathogenesis or a validated circulating marker.",
            "Multi-context analysis distinguished a reproducible but restricted conserved component from broader context-dependent overlap between OA and OC. The 286 shared tissue genes were nearly evenly divided by direction, and external replication and single-cell localization varied by disease context. Stringent dual-cohort blood validation identified G0S2 as the sole blood-persistent association, providing a focused candidate for prospective cell-resolved and protein-level validation.",
        ),
    ]
    for old, new in replacements:
        text = replace_required(text, old, new)
    return text


def apply_discovery_framing(text: str) -> str:
    replacements = [
        (
            "# Shared and context-dependent transcriptomic features between osteoarthritis and ovarian cancer across tissue, single-cell and blood datasets",
            "# Shared molecular features between osteoarthritis and ovarian cancer revealed by multi-layer transcriptomic analyses",
        ),
        (
            "**Running title:** Context-dependent OA-OC transcriptomic overlap",
            "**Running title:** Shared molecular features in OA and OC",
        ),
        (
            "**Background.** Osteoarthritis (OA) and ovarian cancer (OC) arise in distinct tissues, yet both engage matrix remodelling, inflammatory signalling and cellular stress. The extent to which their transcriptomic overlap reflects conserved or context-dependent signals across tissues, source-defined cell populations and peripheral blood remains unresolved.",
            "**Background.** Osteoarthritis (OA) and ovarian cancer (OC) are distinct age-associated diseases that both involve matrix remodelling, inflammatory signalling and cellular stress. Whether these common biological pressures generate detectable molecular features across tissues, source-defined cell populations and peripheral blood remains unclear.",
        ),
        (
            "**Results.** The discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes. Of these, 146 changed in the same direction and 140 changed in opposite directions. Six of ten Hallmark pathways significant in both diseases also had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts, revealing disease-specific stability. The five illustrative genes occupied different source-defined contexts in the OA and OC atlases. Under independent false-discovery-rate control, G0S2 was the sole blood-persistent feature and was lower in all four tissue and blood contrasts.",
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. The shared set comprised 146 concordant and 140 discordant genes, and six of ten jointly significant Hallmark pathways had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
        ),
        (
            "**Conclusions.** Multi-context integration distinguished a restricted conserved component from broader context-dependent transcriptomic overlap between OA and OC. Nearly balanced gene and pathway directions, together with disease-specific replication and cellular localization, defined the boundary of this convergence. G0S2 emerged as a blood-persistent association for prospective cell-resolved and protein-level validation.",
            "**Conclusions.** OA and OC share a measurable molecular landscape involving matrix remodelling, immune regulation and cellular stress. Multi-layer analysis mapped these features across independent cohorts and disease-specific cellular contexts, while revealing context-dependent transcriptional regulation. G0S2 emerged as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
        ),
        (
            "Cross-disease transcriptomics can quantify how molecular features recur across anatomically distant conditions [32-34]. Three dimensions determine whether an overlap is conserved: direction of change, reproducibility across cohorts and cellular context. Single-cell analysis localizes bulk-tissue signals to source-defined populations [35-38], while peripheral blood tests whether a tissue-concordant component also appears in circulation [46-50]. Together, these layers distinguish broadly recurrent signals from tissue- and context-dependent responses.",
            "Despite their distinct anatomy and pathology, OA and OC share links to ageing, matrix remodelling, inflammatory signalling and cellular stress [1-7,16-25]. Cross-disease transcriptomics can test whether these recurring biological pressures leave a detectable molecular imprint across anatomically distant conditions [32-34]. Direction, reproducibility and cellular context then show how shared features are deployed within each disease. Single-cell analysis maps bulk-tissue signals to source-defined populations [35-38], while peripheral blood evaluates their potential systemic persistence [46-50].",
        ),
        (
            "We characterized tissue-level overlap between OA and OC, classified each feature by direction, evaluated external replication and localized illustrative genes in separate single-cell atlases. We then applied a prespecified dual-cohort blood screen. Each disease was analysed independently before evidence integration. This design distinguished reproducibly conserved signals from overlap that depended on tissue, cohort or cellular context.",
            "We therefore asked whether OA and OC share detectable molecular features and how these signals are represented across disease contexts. We defined tissue-level overlap, characterized shared biological programs, evaluated external replication and mapped illustrative genes in separate single-cell atlases. A prespecified dual-cohort blood screen then tested for systemic persistence. Each disease was analysed independently before the evidence layers were integrated.",
        ),
        (
            "### Tissue discovery identified 286 shared alterations with balanced direction classes",
            "### Tissue discovery revealed a substantial shared transcriptomic set",
        ),
        (
            "This near-even partition showed that tissue-level recurrence contained both conserved and disease-divergent components.",
            "Thus, a substantial shared transcriptomic set coexisted with disease-context-dependent regulation.",
        ),
        (
            "### Shared genes converged on recurring functions while pathway direction remained heterogeneous",
            "### Shared genes converged on matrix, immune, stress and cell-cycle programs",
        ),
        (
            "### External replication revealed disease-specific transcriptomic stability",
            "### External cohorts validated disease-specific components of the shared landscape",
        ),
        (
            "### Illustrative genes showed source-defined cellular localization patterns",
            "### Single-cell analysis mapped shared transcripts to disease-specific cellular contexts",
        ),
        (
            "### Stringent blood validation identified G0S2 as a restricted systemic component",
            "### Stringent blood validation identified G0S2 as a candidate systemic feature",
        ),
        (
            "This study defines the extent and boundary of transcriptomic overlap between OA and OC. The 286-gene intersection was robust to prespecified thresholds, yet its near-even directional partition showed that recurrence included both conserved and disease-divergent components. External cohorts further separated a strongly conserved OC component from a more heterogeneous OA component. Variation in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
            "The principal finding of this study is a reproducible shared molecular landscape between OA and OC. The 286-gene intersection was robust to prespecified thresholds and converged on matrix, immune, stress and cell-cycle programs. Directional analysis added an important layer of resolution: concordant and discordant genes were nearly balanced, indicating that shared programs were regulated differently across disease contexts. External cohorts supported this interpretation, with strong conservation in OC and greater cohort dependence in OA. Variation in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
        ),
        (
            "Functional enrichment and single-cell localization clarified why shared membership was not equivalent to uniform biology. Matrix, immune, stress and cell-cycle categories recurred across diseases, but six jointly significant Hallmark pathways had opposite enrichment directions. Illustrative genes also localized to different source-defined populations in the separate OA and OC atlases. Direction and cellular context are therefore integral properties of the observed convergence [9-14,26-30,35-36,38,44-45].",
            "Single-cell mapping extended the tissue-level findings by locating shared transcripts within disease-specific cellular environments. Illustrative genes were distributed across chondrocyte-related populations in OA and malignant, stromal or myeloid populations in OC. Together with the mixed Hallmark directions, these patterns support a model in which common biological pressures engage overlapping molecular features through different tissue microenvironments [9-14,26-30,35-36,38,44-45].",
        ),
        (
            "Stringent systemic validation reduced the broad tissue overlap to a restricted circulating component. G0S2 was lower in both tissues, showed the same direction in both blood cohorts and passed FDR control independently in each cohort. The attrition from 286 shared tissue genes to one dual-FDR result indicates that circulating conservation was substantially more selective than tissue recurrence. G0S2 was most frequently detected in OA ProC and OC Myeloid.cell populations, while reported roles in lipolysis, quiescence and growth regulation provide a biological basis for targeted follow-up [51-54]. Cell-resolved, protein-level and prospective studies are now needed to define its source and functional relevance.",
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Its persistence across tissue and blood compartments identifies G0S2 as a focused candidate systemic feature. Detection in OA ProC and OC Myeloid.cell populations, together with reported roles in lipolysis, quiescence and growth regulation, provides a biological basis for targeted follow-up [51-54]. Prospective cell-resolved and protein-level studies can now test its source and functional relevance.",
        ),
        (
            "Several design features define the scope of these findings. The study used retrospective public data with incomplete covariate annotation, and age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be fully harmonized. Cohorts were therefore analysed separately without cross-cohort ComBat adjustment. Single-cell protocols and annotations differed, while the two blood cohorts used different platforms and cell fractions. These factors may contribute to the observed heterogeneity. Prospective tissue-blood pairing, independent blood replication and protein measurements represent the next validation steps.",
            "The study also has defined boundaries. It used retrospective public data with incomplete covariate annotation, and age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be fully harmonized. Cohorts were therefore analysed separately without cross-cohort ComBat adjustment. Single-cell protocols and annotations differed, while the two blood cohorts used different platforms and cell fractions. These factors may contribute to between-cohort variability. Prospective tissue-blood pairing, independent blood replication and protein measurements represent the next validation steps.",
        ),
        (
            "Multi-context analysis distinguished a reproducible but restricted conserved component from broader context-dependent overlap between OA and OC. The 286 shared tissue genes were nearly evenly divided by direction, and external replication and single-cell localization varied by disease context. Stringent dual-cohort blood validation identified G0S2 as the sole blood-persistent association, providing a focused candidate for prospective cell-resolved and protein-level validation.",
            "Integrated multi-layer transcriptomic analysis identified a shared molecular landscape between OA and OC. The 286 shared tissue genes converged on matrix, immune, stress and cell-cycle programs, while external cohorts and single-cell mapping resolved their disease-specific regulation and cellular contexts. Stringent dual-cohort blood validation prioritized G0S2 as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
        ),
    ]
    for old, new in replacements:
        text = replace_required(text, old, new)
    return text


def apply_final_iteration(text: str) -> str:
    """Apply the final discovery-led, reviewer-bounded narrative pass.

    This pass runs after citation remapping so its source strings match the
    submission-ready manuscript. It preserves directional stratification as an
    analytical contribution while keeping shared molecular features as the
    principal finding.
    """
    replacements = [
        (
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. The shared set comprised 146 concordant and 140 discordant genes, and six of ten jointly significant Hallmark pathways had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. Directional stratification showed that the shared set contained both conserved and disease-context-dependent expression patterns; six of ten jointly significant Hallmark pathways had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
        ),
        (
            "**Conclusions.** OA and OC share a measurable molecular landscape involving matrix remodelling, immune regulation and cellular stress. Multi-layer analysis mapped these features across independent cohorts and disease-specific cellular contexts, while revealing context-dependent transcriptional regulation. G0S2 emerged as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
            "**Conclusions.** We identified a shared molecular landscape between OA and OC involving extracellular-matrix remodelling, immune regulation and cellular stress. Directional stratification and multi-layer mapping resolved how these signals were represented across tissues, source-defined cellular populations and peripheral blood. G0S2 emerged as a candidate systemic feature requiring prospective cell-resolved and protein-level validation.",
        ),
        (
            "This exploratory observational study used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets. OA and OC were processed as separate disease tracks (Figure 1). The evidence sequence comprised tissue discovery, direction classification, functional analysis, external tissue replication, cellular localization and blood evaluation. Five genes were chosen after these analyses to illustrate distinct evidence roles. They were not used to define the shared set or the blood-screen criteria.",
            "This exploratory observational study used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets. OA and OC were analysed as separate disease tracks (Figure 1). The sequence comprised tissue discovery, directional stratification, functional analysis, external replication, cellular localization and blood evaluation. Five genes were selected after the primary analyses only to illustrate distinct evidence patterns; they did not define the shared set or blood-screen criteria.",
        ),
        (
            "Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Public processed matrices were log2-transformed only when indicated by their scale, and multiple probes for one gene were collapsed according to the dataset-specific record. Because reliable technical batch covariates were unavailable and cohorts differed in platform and tissue source, each cohort was modelled separately. Cross-cohort ComBat adjustment and expression-level pooling were deliberately avoided. Disease group was the sole design factor within each cohort. Principal-component and sample-correlation summaries assessed cohort structure without outcome-informed sample removal or adjustment (Additional file 2: Figure S4) [61].",
            "Processed matrices were checked for identifier integrity, duplicated features, non-finite values and group balance. Log2 transformation and probe collapse followed dataset-specific scales and annotations. Each cohort was modelled separately because reliable technical batch covariates were unavailable and platforms or tissue sources differed. Cross-cohort ComBat adjustment and expression-level pooling were not performed. Principal-component and sample-correlation summaries assessed cohort structure without outcome-informed exclusion (Additional file 2: Figure S4) [61].",
        ),
        (
            "Differential expression was estimated separately with limma [55-56]. The primary threshold was Benjamini-Hochberg false-discovery rate (FDR) <0.05 and absolute log2 fold change >=1 [57]. Sensitivity analysis combined FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Additional file 2: Figure S1; Additional file 1: Tables S3a and S3b).\n\nA gene entered the shared set only when it met the primary threshold in both discovery cohorts. Shared genes were classified as higher in both diseases, lower in both, higher in OA and lower in OC, or lower in OA and higher in OC. The first two classes were concordant and the latter two discordant. Classification preserved within-disease effect scales rather than aligning expression values across platforms.",
            "Limma estimated each disease-versus-reference contrast separately [55-56]. The primary threshold was Benjamini-Hochberg FDR <0.05 and absolute log2 fold change >=1 [57]. Sensitivity analysis combined FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Additional file 2: Figure S1; Additional file 1: Tables S3a and S3b).\n\nGenes meeting the primary threshold in both discovery cohorts entered the shared set. Matching OA and OC signs defined concordant genes; opposing signs defined disease-context-dependent patterns. Classification retained within-disease effect scales without cross-platform expression alignment.",
        ),
        (
            "Gene Ontology and Kyoto Encyclopedia of Genes and Genomes over-representation analyses were applied to the 286 shared genes with clusterProfiler [62-65]. Figure 3 shows ten significant, non-redundant terms representing the principal categories; complete results are provided in Additional file 1: Tables S12 and S13.\n\nRanked Hallmark gene-set enrichment analysis was performed independently on the complete OA and OC discovery statistics [58-59]. A pathway was shared when FDR <0.05 in both diseases. Paired normalized enrichment scores described matching or opposite pathway directions. External tissue cohorts were analysed when a suitable ranked statistic was available (Additional file 1: Tables S6 and S14).",
            "Gene Ontology and Kyoto Encyclopedia of Genes and Genomes over-representation analyses were applied to the 286 shared genes with clusterProfiler [62-65]. Figure 3 shows ten significant, non-redundant terms; complete results are provided in Additional file 1: Tables S12 and S13.\n\nHallmark gene-set enrichment analysis used the complete ranked OA and OC discovery statistics [58-59]. Pathways with FDR <0.05 in both diseases were compared by normalized enrichment-score direction. External cohorts were analysed when a suitable ranked statistic was available (Additional file 1: Tables S6 and S14).",
        ),
        (
            "Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29-30]. Count-level datasets underwent cell- and feature-level quality control, mitochondrial-content assessment, compatible doublet handling, normalization and dimensional reduction. TPM-only input was not treated as raw counts. The workflow followed established single-cell analysis guidance [35-39]. Doublet methods were reviewed for compatibility, and scDblFinder was used when sample-level counts permitted [40-42]. UCell supported sample-aware score auditing [43]. Dataset eligibility and quality-control outcomes are reported in Additional file 1: Table S8a.\n\nThe main analysis used count-level OA GSE255460 and OC GSE154600. The atlases remained separate because their tissues, platforms, cellular composition and annotations differed. Exact source labels were retained. For each illustrative gene, detection fraction and mean UMI count were calculated within source labels. Fractions were summarized within each atlas to preserve its source-defined cellular composition. All five embeddings and gene-level summaries are provided in Additional file 2: Figure S3 and Additional file 1: Table S8b.",
            "Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29-30]. Eligible count-level datasets underwent cell- and feature-level quality control, mitochondrial-content assessment, compatible doublet handling, normalization and dimensional reduction; TPM-only input was not treated as raw counts [35-42]. Dataset-specific thresholds, scDblFinder use and sample-aware UCell audits are documented in Additional file 1: Table S8a [43].\n\nThe main localization analysis used count-level OA GSE255460 and OC GSE154600, which remained separate because their tissues, platforms and source annotations differed. Detection fractions and mean UMI counts for the illustrative genes were summarized within exact source labels. All embeddings and gene-level summaries are provided in Additional file 2: Figure S3 and Additional file 1: Table S8b.",
        ),
        (
            "Official platform annotations were used for probe-to-gene mapping. Ambiguous mappings were discarded, and the probe with the highest interquartile range represented genes with multiple unambiguous probes. Limma estimated separate contrasts. Hedges g and 95% confidence intervals summarized within-cohort effects; effects were not pooled across diseases.\n\nThe prespecified screen required membership in the shared tissue set, concordant tissue direction, measurement in both blood cohorts, the same sign in all four contrasts and FDR <0.05 independently in each blood cohort. Denominators were retained at every step (Additional file 1: Tables S9-S11).",
            "Official platform annotations supported probe-to-gene mapping; ambiguous probes were removed, and the highest-interquartile-range probe represented multiply mapped genes. Limma estimated separate contrasts. Hedges g and 95% confidence intervals summarized within-cohort effects without cross-disease pooling.\n\nThe prespecified screen required shared-set membership, concordant tissue direction, measurement in both blood cohorts, the same sign in all four contrasts and FDR <0.05 in each blood cohort. Denominators were retained throughout (Additional file 1: Tables S9-S11).",
        ),
        (
            "Direction divided the overlap almost evenly (Figure 2D).",
            "### Directional analysis refined the shared molecular landscape\n\nDirectional stratification divided the overlap almost evenly (Figure 2D).",
        ),
        (
            "### External cohorts validated disease-specific components of the shared landscape",
            "### External cohorts evaluated the reproducibility of shared transcriptional features",
        ),
        (
            "External direction agreement was higher in OC than in OA (Figure 4A; Additional file 1: Tables S4-S6). GSE117999 reproduced 143 of 280 measurable OA signs (51.1%; binomial P=0.765; Spearman rho=0.061). GSE82107 reproduced 178 of 286 signs (62.2%; P=4.15 x 10^-5; rho=0.219), although no shared gene reached external FDR <0.05. The OC cohorts reproduced 260 of 286 signs (90.9%; P=1.01 x 10^-49; rho=0.847) and 179 of 226 signs (79.2%; P=2.47 x 10^-19; rho=0.628). The corresponding numbers meeting both external FDR and sign criteria were 179 and 134.",
            "External direction agreement varied by disease and cohort (Figure 4A; Additional file 1: Tables S4-S6). The OA cohorts reproduced 143 of 280 (51.1%; Spearman rho=0.061) and 178 of 286 (62.2%; rho=0.219) measurable signs. The OC cohorts reproduced 260 of 286 (90.9%; rho=0.847) and 179 of 226 (79.2%; rho=0.628) signs. The corresponding numbers meeting both external FDR and sign criteria were 179 and 134 in OC. Full binomial and FDR results are reported in Additional file 1: Tables S4-S6.",
        ),
        (
            "G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 were used as illustrative genes because they captured different evidence patterns (Figure 4B-C; Additional file 1: Table S7).",
            "G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 were used as illustrative genes because they captured different evidence patterns (Figure 4B; Additional file 2: Figure S5; Additional file 1: Table S7).",
        ),
        (
            "The principal finding of this study is a reproducible shared molecular landscape between OA and OC. The 286-gene intersection was robust to prespecified thresholds and converged on matrix, immune, stress and cell-cycle programs. Directional analysis added an important layer of resolution: concordant and discordant genes were nearly balanced, indicating that shared programs were regulated differently across disease contexts. External cohorts supported this interpretation, with strong conservation in OC and greater cohort dependence in OA. Variation in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
            "The principal finding of this study is a measurable shared molecular landscape between OA and OC. The 286-gene intersection was robust to prespecified thresholds and converged on matrix, immune, stress and cell-cycle programs. Directional stratification provided additional resolution by separating conserved from disease-context-dependent expression patterns. External cohorts evaluated reproducibility, revealing high agreement in OC and greater variability in OA. Differences in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
        ),
        (
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Its persistence across tissue and blood compartments identifies G0S2 as a focused candidate systemic feature. Detection in OA ProC and OC Myeloid.cell populations, together with reported roles in lipolysis, quiescence and growth regulation, provides a biological basis for targeted follow-up [51-54]. Prospective cell-resolved and protein-level studies can now test its source and functional relevance.",
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Its persistence across tissue and blood compartments identifies G0S2 as a focused candidate systemic feature. Reported roles in lipolysis, quiescence and growth regulation offer a hypothesis for targeted follow-up, potentially reflecting age-associated systemic transcriptional regulation involving cell quiescence and metabolic pathways [51-54]. This interpretation remains provisional pending prospective cell-resolved and protein-level validation.",
        ),
        (
            "Integrated multi-layer transcriptomic analysis identified a shared molecular landscape between OA and OC. The 286 shared tissue genes converged on matrix, immune, stress and cell-cycle programs, while external cohorts and single-cell mapping resolved their disease-specific regulation and cellular contexts. Stringent dual-cohort blood validation prioritized G0S2 as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
            "We identified a shared molecular landscape between OA and OC involving extracellular-matrix remodelling, immune regulation and cellular stress. Directional stratification, external cohorts and single-cell mapping resolved how these signals were represented across disease contexts. Stringent dual-cohort blood validation prioritized G0S2 as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
        ),
    ]
    for old, new in replacements:
        text = replace_required(text, old, new)
    return text


def apply_acceptance_polish(text: str) -> str:
    """Apply the final BMC-focused acceptance polish.

    The pass keeps all quantitative findings while reducing negative framing,
    shortening the main-text Methods and expanding the evidence-bounded G0S2
    interpretation. Detailed single-cell audit parameters remain in the
    supplementary tables.
    """
    replacements = [
        (
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. Directional stratification showed that the shared set contained both conserved and disease-context-dependent expression patterns; six of ten jointly significant Hallmark pathways had opposite enrichment directions. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. Directional stratification further characterized disease-context-dependent regulation among the shared programs. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
        ),
        (
            "Despite their distinct anatomy and pathology, OA and OC share links to ageing, matrix remodelling, inflammatory signalling and cellular stress [1-7,16-25]. Cross-disease transcriptomics can test whether these recurring biological pressures leave a detectable molecular imprint across anatomically distant conditions [32-34]. Direction, reproducibility and cellular context then show how shared features are deployed within each disease. Single-cell analysis maps bulk-tissue signals to source-defined populations [35-38], while peripheral blood evaluates their potential systemic persistence [46-50].",
            "Despite their distinct anatomy and pathology, OA and OC share links to ageing, matrix remodelling, inflammatory signalling and cellular stress [1-7,16-25]. Cross-disease transcriptomics can test whether these recurring biological pressures leave a detectable molecular imprint across anatomically distant conditions [32-34]. Such comparative analyses may reveal conserved molecular processes that are not apparent when diseases are studied independently. Direction, reproducibility and cellular context then show how shared features are represented within each disease. Single-cell analysis maps bulk-tissue signals to source-defined populations [35-38], while peripheral blood evaluates their potential systemic persistence [46-50].",
        ),
        (
            "This exploratory observational study used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets. OA and OC were analysed as separate disease tracks (Figure 1). The sequence comprised tissue discovery, directional stratification, functional analysis, external replication, cellular localization and blood evaluation. Five genes were selected after the primary analyses only to illustrate distinct evidence patterns; they did not define the shared set or blood-screen criteria.",
            "This exploratory observational study analysed public bulk-tissue, single-cell and peripheral-blood transcriptomes on separate OA and OC tracks (Figure 1). Evidence was integrated across tissue discovery, directional stratification, functional analysis, external replication, cellular localization and blood evaluation. Five genes illustrated evidence patterns after the primary analyses and did not define the shared set or blood screen.",
        ),
        (
            "Datasets were obtained from the NCBI Gene Expression Omnibus [66-67]. The OA discovery cohort was knee cartilage GSE114007, with 20 OA and 18 non-OA samples [8]. The OC discovery cohort was GSE18520, with 53 advanced high-grade serous tumours and 10 normal ovarian surface-epithelium samples [23]. External OA cohorts were GSE117999 (10 OA and 10 controls) and GSE82107 (10 OA and 7 controls). External OC cohorts were GSE54388 (16 tumours and 6 normal samples) and GSE12470 (43 serous carcinomas and 10 normal peritoneal samples) [24]. Cohort roles, platforms and group definitions are listed in Additional file 1: Table S1.",
            "Datasets were obtained from the NCBI Gene Expression Omnibus [66-67]. Discovery cohorts were knee cartilage GSE114007 (20 OA, 18 non-OA) [8] and GSE18520 (53 advanced high-grade serous tumours, 10 normal ovarian surface-epithelium samples) [23]. External cohorts were OA GSE117999 (10 OA, 10 controls), OA GSE82107 (10 OA, 7 controls), OC GSE54388 (16 tumours, 6 normal samples) and OC GSE12470 (43 serous carcinomas, 10 normal peritoneal samples) [24]. Roles, platforms and group definitions are provided in Additional file 1: Table S1.",
        ),
        (
            "Processed matrices were checked for identifier integrity, duplicated features, non-finite values and group balance. Log2 transformation and probe collapse followed dataset-specific scales and annotations. Each cohort was modelled separately because reliable technical batch covariates were unavailable and platforms or tissue sources differed. Cross-cohort ComBat adjustment and expression-level pooling were not performed. Principal-component and sample-correlation summaries assessed cohort structure without outcome-informed exclusion (Additional file 2: Figure S4) [61].",
            "Processed matrices were checked for identifiers, duplicate features, non-finite values and group balance. Transformation and probe collapse followed dataset-specific scales and annotations. Cohorts were modelled separately because reliable batch covariates were unavailable and platforms or tissues differed. No cross-cohort ComBat adjustment or expression pooling was performed. Principal-component and sample-correlation summaries assessed structure without outcome-informed exclusion (Additional file 2: Figure S4) [61].",
        ),
        (
            "Limma estimated each disease-versus-reference contrast separately [55-56]. The primary threshold was Benjamini-Hochberg FDR <0.05 and absolute log2 fold change >=1 [57]. Sensitivity analysis combined FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Additional file 2: Figure S1; Additional file 1: Tables S3a and S3b).",
            "Limma estimated each disease-versus-reference contrast [55-56]. The primary threshold was Benjamini-Hochberg FDR <0.05 and absolute log2 fold change >=1 [57]. Sensitivity analysis crossed FDR thresholds of 0.01 and 0.05 with absolute log2-fold-change thresholds of 0.5, 1.0 and 1.5 (Additional file 2: Figure S1; Additional file 1: Tables S3a and S3b).",
        ),
        (
            "Genes meeting the primary threshold in both discovery cohorts entered the shared set. Matching OA and OC signs defined concordant genes; opposing signs defined disease-context-dependent patterns. Classification retained within-disease effect scales without cross-platform expression alignment.",
            "Genes meeting the primary threshold in both discovery cohorts entered the shared set. Matching signs defined conserved direction; differing signs defined disease-context-dependent direction. Classification retained within-disease effect scales without cross-platform expression alignment.",
        ),
        (
            "Each external cohort was modelled independently with the disease-versus-reference orientation of its discovery cohort. Replication was summarized by sign agreement for every measurable shared gene. An exact binomial test compared agreement with 0.5, and Spearman correlation described gene-wise concordance between discovery and external log2 fold changes. External FDR support was recorded without pooling expression values across cohorts.",
            "External cohorts used the disease-versus-reference orientation of their discovery cohort. Replication was summarized by sign agreement for each measurable shared gene. Exact binomial tests compared agreement with 0.5, and Spearman correlations related discovery and external log2 fold changes. External FDR support was recorded without pooling expression values.",
        ),
        (
            "Gene Ontology and Kyoto Encyclopedia of Genes and Genomes over-representation analyses were applied to the 286 shared genes with clusterProfiler [62-65]. Figure 3 shows ten significant, non-redundant terms; complete results are provided in Additional file 1: Tables S12 and S13.",
            "ClusterProfiler tested Gene Ontology and Kyoto Encyclopedia of Genes and Genomes over-representation among the 286 shared genes [62-65]. Figure 3 shows ten significant, non-redundant terms; complete results are provided in Additional file 1: Tables S12 and S13.",
        ),
        (
            "Hallmark gene-set enrichment analysis used the complete ranked OA and OC discovery statistics [58-59]. Pathways with FDR <0.05 in both diseases were compared by normalized enrichment-score direction. External cohorts were analysed when a suitable ranked statistic was available (Additional file 1: Tables S6 and S14).",
            "Hallmark gene-set enrichment analysis used complete ranked discovery statistics [58-59]. Pathways with FDR <0.05 in both diseases were compared by normalized enrichment-score direction. Eligible external cohorts were analysed similarly (Additional file 1: Tables S6 and S14).",
        ),
        (
            "Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29-30]. Eligible count-level datasets underwent cell- and feature-level quality control, mitochondrial-content assessment, compatible doublet handling, normalization and dimensional reduction; TPM-only input was not treated as raw counts [35-42]. Dataset-specific thresholds, scDblFinder use and sample-aware UCell audits are documented in Additional file 1: Table S8a [43].",
            "Five public single-cell datasets were audited with dataset-specific adapters: OA GSE104782, GSE169454 and GSE255460, and OC GSE154600 and GSE180661 [9-12,29-30]. Eligible count-level data underwent quality control, normalization and dimensional reduction following established guidance [35-43]. Dataset eligibility, input type, thresholds, doublet handling and score audits are reported in Additional file 1: Table S8a.",
        ),
        (
            "The main localization analysis used count-level OA GSE255460 and OC GSE154600, which remained separate because their tissues, platforms and source annotations differed. Detection fractions and mean UMI counts for the illustrative genes were summarized within exact source labels. All embeddings and gene-level summaries are provided in Additional file 2: Figure S3 and Additional file 1: Table S8b.",
            "The main localization analysis used OA GSE255460 and OC GSE154600 separately because their tissues, platforms and annotations differed. Detection fractions and mean UMI counts for illustrative genes were summarized within source labels. Embeddings and gene-level summaries are provided in Additional file 2: Figure S3 and Additional file 1: Table S8b.",
        ),
        (
            "Eligible disease-versus-control blood cohorts were identified from official repository metadata. OA GSE48556 contained peripheral-blood mononuclear cells from 106 cases and 33 healthy controls [15]. OC GSE31682 contained a blood-cell fraction from 48 cases and 20 healthy controls [31]. The platforms and blood fractions differed, so cohorts were modelled separately (Additional file 1: Table S9).",
            "Official metadata identified eligible disease-versus-control blood cohorts. OA GSE48556 contained peripheral-blood mononuclear cells from 106 cases and 33 healthy controls [15]. OC GSE31682 contained a blood-cell fraction from 48 cases and 20 controls [31]. Different platforms and blood fractions required separate models (Additional file 1: Table S9).",
        ),
        (
            "Official platform annotations supported probe-to-gene mapping; ambiguous probes were removed, and the highest-interquartile-range probe represented multiply mapped genes. Limma estimated separate contrasts. Hedges g and 95% confidence intervals summarized within-cohort effects without cross-disease pooling.",
            "Official annotations supported probe-to-gene mapping. Ambiguous probes were removed, and the highest-interquartile-range probe represented genes with multiple probes. Limma estimated separate contrasts; Hedges g and 95% confidence intervals summarized effects without cross-disease pooling.",
        ),
        (
            "The prespecified screen required shared-set membership, concordant tissue direction, measurement in both blood cohorts, the same sign in all four contrasts and FDR <0.05 in each blood cohort. Denominators were retained throughout (Additional file 1: Tables S9-S11).",
            "The prespecified screen required shared-set membership, conserved tissue direction, measurement in both blood cohorts, a consistent sign across all four contrasts, and FDR <0.05 in each blood cohort. Denominators were retained throughout (Additional file 1: Tables S9-S11).",
        ),
        (
            "Tests were two-sided unless otherwise stated. Multiple-testing procedures and direction definitions were fixed before biological interpretation. Sample and cell units are reported for each evidence layer, and denominators accompany filtering steps. Randomized procedures used seed 20260726. Bioconductor resources supported the analysis environment [68]. Full parameters, source manifests and executable scripts are included in the reproducible project.",
            "Tests were two-sided unless stated otherwise. Multiple-testing procedures and direction definitions preceded biological interpretation. Sample and cell units and filtering denominators are reported for each evidence layer. Randomized procedures used seed 20260726. Bioconductor supported the analysis environment [68]. Full parameters, manifests and scripts are included in the reproducible project.",
        ),
        (
            "Directional stratification divided the overlap almost evenly (Figure 2D). There were 112 genes higher in both diseases and 34 lower in both. Another 86 genes were higher in OA and lower in OC, while 54 were lower in OA and higher in OC. The shared set therefore contained 146 concordant genes (51.0%) and 140 discordant genes (49.0%). Thus, a substantial shared transcriptomic set coexisted with disease-context-dependent regulation.",
            "Directional stratification further characterized the shared molecular landscape (Figure 2D). There were 112 genes higher in both diseases and 34 lower in both. Another 86 genes were higher in OA and lower in OC, while 54 were lower in OA and higher in OC. The shared set therefore contained 146 genes with conserved direction (51.0%) and 140 with disease-context-dependent direction (49.0%). Thus, a substantial shared transcriptomic set coexisted with context-dependent regulation.",
        ),
        (
            "The shared set was enriched for extracellular-matrix, immune, stress and cell-cycle processes (Figure 3; Additional file 1: Tables S12-S14). Gene Ontology terms included chromosome segregation, matrix organization, cytokine responses, hypoxia and integrated stress signalling. Cell cycle was the only FDR-significant Kyoto Encyclopedia of Genes and Genomes pathway. Ten Hallmark pathways were significant in both discovery cohorts. Four had matching normalized enrichment-score signs, whereas six had opposite signs.",
            "The shared set was enriched for extracellular-matrix, immune, stress and cell-cycle processes (Figure 3; Additional file 1: Tables S12-S14). Gene Ontology terms included chromosome segregation, matrix organization, cytokine responses, hypoxia and integrated stress signalling. Cell cycle was the only FDR-significant Kyoto Encyclopedia of Genes and Genomes pathway. Ten Hallmark pathways were significant in both discovery cohorts. Four had matching normalized enrichment-score signs, whereas six differed between diseases.",
        ),
        (
            "External direction agreement varied by disease and cohort (Figure 4A; Additional file 1: Tables S4-S6). The OA cohorts reproduced 143 of 280 (51.1%; Spearman rho=0.061) and 178 of 286 (62.2%; rho=0.219) measurable signs. The OC cohorts reproduced 260 of 286 (90.9%; rho=0.847) and 179 of 226 (79.2%; rho=0.628) signs. The corresponding numbers meeting both external FDR and sign criteria were 179 and 134 in OC. Full binomial and FDR results are reported in Additional file 1: Tables S4-S6.",
            "External cohorts demonstrated variable reproducibility across disease contexts (Figure 4A; Additional file 1: Tables S4-S6). The OA cohorts reproduced 143 of 280 (51.1%; Spearman rho=0.061) and 178 of 286 (62.2%; rho=0.219) measurable signs. The OC cohorts reproduced 260 of 286 (90.9%; rho=0.847) and 179 of 226 (79.2%; rho=0.628) signs. The corresponding numbers meeting both external FDR and sign criteria were 179 and 134 in OC. Full binomial and FDR results are reported in Additional file 1: Tables S4-S6.",
        ),
        (
            "The set also spanned concordant and discordant disease effects, source-defined cellular localization and the blood-screen outcome.",
            "The set also spanned conserved and disease-context-dependent effects, source-defined cellular localization and the blood-screen outcome.",
        ),
        (
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Its persistence across tissue and blood compartments identifies G0S2 as a focused candidate systemic feature. Reported roles in lipolysis, quiescence and growth regulation offer a hypothesis for targeted follow-up, potentially reflecting age-associated systemic transcriptional regulation involving cell quiescence and metabolic pathways [51-54]. This interpretation remains provisional pending prospective cell-resolved and protein-level validation.",
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Recurrence across cartilage, tumour tissue and two blood fractions is notable because these compartments differ markedly in cellular composition. Prior studies link G0S2 to lipolysis, cellular quiescence and growth regulation [51-54]. The observed pattern is compatible with a systemic response to age-associated metabolic and inflammatory pressures, but it does not identify a common pathological driver. This interpretation remains hypothesis-generating. Prospective tissue-blood pairing, cell-resolved analysis and protein measurement are needed to determine the cellular source and biological relevance of this signal.",
        ),
    ]
    for old, new in replacements:
        text = replace_required(text, old, new)
    return text


def apply_human_style_polish(text: str) -> str:
    """Reduce formulaic phrasing without changing claims, data or terminology.

    This deliberately limited pass varies sentence structure in the abstract,
    background and discussion, replaces noun-heavy workflow summaries with
    direct prose, and keeps the Methods and Results quantitatively unchanged.
    """
    replacements = [
        (
            "**Background.** Osteoarthritis (OA) and ovarian cancer (OC) are distinct age-associated diseases that both involve matrix remodelling, inflammatory signalling and cellular stress. Whether these common biological pressures generate detectable molecular features across tissues, source-defined cell populations and peripheral blood remains unclear.",
            "**Background.** Although osteoarthritis (OA) and ovarian cancer (OC) arise through different pathological processes, both involve extracellular-matrix alteration, inflammatory regulation and cellular stress. It remains unclear whether these shared biological pressures leave reproducible molecular signals in tissue, source-defined cell populations and peripheral blood.",
        ),
        (
            "**Methods.** OA and OC datasets were analysed on separate disease tracks. Discovery tissue transcriptomes defined shared differentially expressed genes and their direction of change. Independent tissue cohorts evaluated replication, and separate single-cell atlases localized five illustrative genes. A prespecified blood screen then tested whether tissue-concordant signals persisted in both diseases.",
            "**Methods.** We analysed OA and OC separately before comparing their transcriptional changes. Discovery cohorts defined the shared genes and their direction of change, while independent cohorts tested reproducibility. We then examined five illustrative genes in separate single-cell atlases and screened tissue-concordant signals in two blood datasets.",
        ),
        (
            "**Results.** Discovery analyses identified 2,008 OA and 2,310 OC differentially expressed genes, including 286 shared genes enriched for matrix, immune, stress and cell-cycle processes. Directional stratification further characterized disease-context-dependent regulation among the shared programs. External directional agreement was 51.1-62.2% in OA cohorts and 79.2-90.9% in OC cohorts. Single-cell analysis mapped five illustrative genes to disease-specific cellular contexts. Stringent blood validation identified G0S2 as a candidate systemic feature that was lower in all four tissue and blood contrasts.",
            "**Results.** The OA and OC discovery analyses identified 2,008 and 2,310 differentially expressed genes, respectively, with 286 genes in common. These genes were enriched for matrix, immune, stress and cell-cycle processes, but their direction of change varied by disease context. External direction agreement ranged from 51.1% to 62.2% in OA and from 79.2% to 90.9% in OC. Single-cell data placed five illustrative genes in disease-specific cellular contexts. G0S2 alone passed false-discovery-rate control in both blood cohorts and was lower in all four tissue and blood contrasts.",
        ),
        (
            "**Conclusions.** We identified a shared molecular landscape between OA and OC involving extracellular-matrix remodelling, immune regulation and cellular stress. Directional stratification and multi-layer mapping resolved how these signals were represented across tissues, source-defined cellular populations and peripheral blood. G0S2 emerged as a candidate systemic feature requiring prospective cell-resolved and protein-level validation.",
            "**Conclusions.** OA and OC share a measurable, but context-dependent, transcriptional component involving extracellular-matrix remodelling, immune regulation and cellular stress. G0S2 warrants further study as a possible systemic feature, although its cellular source and protein-level relevance remain unknown.",
        ),
        (
            "Osteoarthritis is a whole-joint disorder characterized by cartilage loss, subchondral-bone remodelling, synovial change and low-grade inflammatory signalling [1-7]. Transcriptomic and single-cell studies have identified heterogeneous chondrocyte and stromal states in diseased cartilage [8-14]. Selected OA-associated expression changes are also detectable in blood, although circulating profiles depend on cohort design and cell composition [15].",
            "Osteoarthritis affects the entire joint. Its pathology includes cartilage loss, remodelling of subchondral bone, synovial changes and low-grade inflammatory signalling [1-7]. Transcriptomic and single-cell studies have described diverse chondrocyte and stromal states in diseased cartilage [8-14]. Some OA-related expression changes are detectable in blood, although results vary with cohort design and cellular composition [15].",
        ),
        (
            "Ovarian cancer is a heterogeneous malignancy shaped by genomic instability, molecular subtype, stromal remodelling and immune context [16-25]. Bulk and single-cell studies of high-grade serous ovarian cancer have resolved malignant, fibroblast, endothelial and immune populations across disease sites [26-30,44-45]. Blood-based measurements capture tumour-associated variation together with systemic and haematological responses [31].",
            "Ovarian cancer is also heterogeneous, although its transcriptional variation is shaped by genomic instability, molecular subtype, stromal remodelling and immune context [16-25]. Studies of high-grade serous ovarian cancer at bulk and single-cell resolution have resolved malignant, fibroblast, endothelial and immune populations across disease sites [26-30,44-45]. Blood measurements combine tumour-associated variation with systemic and haematological responses [31].",
        ),
        (
            "Despite their distinct anatomy and pathology, OA and OC share links to ageing, matrix remodelling, inflammatory signalling and cellular stress [1-7,16-25]. Cross-disease transcriptomics can test whether these recurring biological pressures leave a detectable molecular imprint across anatomically distant conditions [32-34]. Such comparative analyses may reveal conserved molecular processes that are not apparent when diseases are studied independently. Direction, reproducibility and cellular context then show how shared features are represented within each disease. Single-cell analysis maps bulk-tissue signals to source-defined populations [35-38], while peripheral blood evaluates their potential systemic persistence [46-50].",
            "OA and OC differ in anatomy and pathogenesis, yet both are linked to ageing, matrix remodelling, inflammation and cellular stress [1-7,16-25]. Comparing them may expose molecular processes that are difficult to recognize when each disease is examined alone [32-34]. Direction of change, reproducibility and cellular location matter because a shared gene need not represent the same biological state in both tissues. Single-cell data can localize bulk signals [35-38], while blood data provide a separate test of systemic persistence [46-50].",
        ),
        (
            "We therefore asked whether OA and OC share detectable molecular features and how these signals are represented across disease contexts. We defined tissue-level overlap, characterized shared biological programs, evaluated external replication and mapped illustrative genes in separate single-cell atlases. A prespecified dual-cohort blood screen then tested for systemic persistence. Each disease was analysed independently before the evidence layers were integrated.",
            "We asked which transcriptional changes were shared by OA cartilage and OC tissue, whether they recurred in independent cohorts, and where selected genes were detected at single-cell resolution. A prespecified screen then tested whether any tissue-concordant signal was also present in blood. OA and OC were analysed separately throughout; comparisons were made only after within-disease estimation.",
        ),
        (
            "This exploratory observational study analysed public bulk-tissue, single-cell and peripheral-blood transcriptomes on separate OA and OC tracks (Figure 1). Evidence was integrated across tissue discovery, directional stratification, functional analysis, external replication, cellular localization and blood evaluation. Five genes illustrated evidence patterns after the primary analyses and did not define the shared set or blood screen.",
            "We used public bulk-tissue, single-cell and peripheral-blood transcriptomic datasets in an exploratory observational design (Figure 1). OA and OC were analysed separately. Tissue-level discovery was followed by functional analysis and external replication; selected genes were then examined in single-cell and blood datasets. The five illustrative genes were chosen after the primary analyses and did not determine the shared set or blood-screen criteria.",
        ),
        (
            "Directional stratification further characterized the shared molecular landscape (Figure 2D). There were 112 genes higher in both diseases and 34 lower in both. Another 86 genes were higher in OA and lower in OC, while 54 were lower in OA and higher in OC. The shared set therefore contained 146 genes with conserved direction (51.0%) and 140 with disease-context-dependent direction (49.0%). Thus, a substantial shared transcriptomic set coexisted with context-dependent regulation.",
            "Direction of change added detail to the shared set (Figure 2D). Of the 286 genes, 112 were higher in both diseases and 34 were lower in both. Eighty-six were higher in OA and lower in OC, whereas 54 showed the reverse pattern. Thus, 146 genes (51.0%) had a conserved direction and 140 (49.0%) were disease-context-dependent.",
        ),
        (
            "The concordant direction and comparable standardized effects therefore identified G0S2 as the blood-persistent component of the tissue-derived set across the two blood fractions.",
            "The similar directions and standardized effects across the two blood datasets made G0S2 the only blood-persistent feature that met the prespecified rule.",
        ),
        (
            "The principal finding of this study is a measurable shared molecular landscape between OA and OC. The 286-gene intersection was robust to prespecified thresholds and converged on matrix, immune, stress and cell-cycle programs. Directional stratification provided additional resolution by separating conserved from disease-context-dependent expression patterns. External cohorts evaluated reproducibility, revealing high agreement in OC and greater variability in OA. Differences in cartilage source, joint site, disease severity, platform and reference group may contribute to the OA pattern [2-3,8,16,18,23-24,32-34].",
            "The comparison identified 286 genes shared by OA and OC, with enrichment for matrix, immune, stress and cell-cycle processes. This overlap was not uniform: approximately half the genes changed in the same direction, and external agreement was high in OC but more variable in OA. Differences in cartilage source, joint site, disease severity, platform and reference group are plausible contributors to the OA variation [2-3,8,16,18,23-24,32-34].",
        ),
        (
            "Single-cell mapping extended the tissue-level findings by locating shared transcripts within disease-specific cellular environments. Illustrative genes were distributed across chondrocyte-related populations in OA and malignant, stromal or myeloid populations in OC. Together with the mixed Hallmark directions, these patterns support a model in which common biological pressures engage overlapping molecular features through different tissue microenvironments [9-14,26-30,35-36,38,44-45].",
            "The single-cell data provide a second level of context. The illustrative genes occurred in chondrocyte-related populations in OA, but in malignant, stromal or myeloid populations in OC. Hallmark enrichment also differed by disease. These observations are consistent with similar biological pressures acting through different tissue environments; they do not show that a shared gene has the same function in OA and OC [9-14,26-30,35-36,38,44-45].",
        ),
        (
            "Among the shared transcripts, G0S2 showed the strongest systemic consistency. It was lower in both discovery tissues, changed in the same direction in both blood cohorts and passed independent FDR control in each disease. Recurrence across cartilage, tumour tissue and two blood fractions is notable because these compartments differ markedly in cellular composition. Prior studies link G0S2 to lipolysis, cellular quiescence and growth regulation [51-54]. The observed pattern is compatible with a systemic response to age-associated metabolic and inflammatory pressures, but it does not identify a common pathological driver. This interpretation remains hypothesis-generating. Prospective tissue-blood pairing, cell-resolved analysis and protein measurement are needed to determine the cellular source and biological relevance of this signal.",
            "G0S2 was the most consistent cross-compartment result. Its expression was lower in both discovery tissues and both blood datasets, and it passed FDR control in each blood cohort. This recurrence is notable because cartilage, tumour tissue and the two blood fractions have markedly different cellular compositions. Previous work has linked G0S2 to lipolysis, cellular quiescence and growth regulation [51-54]. One possibility is that its lower expression reflects a systemic response to metabolic or inflammatory stress rather than a common disease driver. Paired tissue and blood samples will be needed to distinguish blood-cell changes from signals related to the diseased tissues. Protein measurements would provide an additional test of this interpretation.",
        ),
        (
            "The study also has defined boundaries. It used retrospective public data with incomplete covariate annotation, and age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be fully harmonized. Cohorts were therefore analysed separately without cross-cohort ComBat adjustment. Single-cell protocols and annotations differed, while the two blood cohorts used different platforms and cell fractions. These factors may contribute to between-cohort variability. Prospective tissue-blood pairing, independent blood replication and protein measurements represent the next validation steps.",
            "Several limitations affect the interpretation. The public datasets were retrospective and lacked a common set of covariates. Age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be fully harmonized. We therefore analysed cohorts separately and did not apply cross-cohort ComBat correction. Single-cell protocols and annotations differed, as did the platforms and cell fractions used for blood profiling. Independent blood cohorts, ideally with paired tissue samples, are needed before G0S2 can be considered reproducible beyond the datasets analysed here.",
        ),
        (
            "We identified a shared molecular landscape between OA and OC involving extracellular-matrix remodelling, immune regulation and cellular stress. Directional stratification, external cohorts and single-cell mapping resolved how these signals were represented across disease contexts. Stringent dual-cohort blood validation prioritized G0S2 as a candidate systemic feature for prospective cell-resolved and protein-level validation.",
            "OA and OC shared 286 transcriptional changes related mainly to matrix, immune, stress and cell-cycle processes. Their direction, external reproducibility and cellular localization varied with disease context. G0S2 was the only feature supported in both blood cohorts, but its source and biological significance require prospective validation.",
        ),
        (
            "Full parameters, manifests and scripts are included in the reproducible project.",
            "Full parameters, source manifests and executable scripts are supplied in Additional file 3.",
        ),
    ]
    for old, new in replacements:
        text = replace_required(text, old, new)
    return text


def transform_main_body(text: str) -> str:
    text = text.replace("**Article type:** Original Research", "**Article type:** Research Article")
    text = text.replace("**Target journal:** [To be selected before submission]", "**Target journal:** BMC Medical Genomics")
    text = re.sub(
        r"\*\*Authors:\*\*.*",
        "**Authors:** Junhui Shi¹†, Mengxiang Liu¹†, Repkat Inayatilla¹, Ke Li¹, Lei Chen¹*  ",
        text,
        count=1,
    )
    text = text.replace("**Affiliations:** [Affiliations to be added]", "**Affiliation:** ¹Department of Orthopedics Center, The First Affiliated Hospital of Shihezi University, Shihezi University, Shihezi 832008, China")
    text = text.replace("**Corresponding author:** [Name, postal address and email to be added]", "**Corresponding author:** *Lei Chen, Department of Orthopedics Center, The First Affiliated Hospital of Shihezi University, Shihezi University, Shihezi 832008, China; email: 564386249@qq.com")
    if "**Equal contribution:**" not in text:
        text = text.replace("**Target journal:** BMC Medical Genomics", "**Equal contribution:** †Junhui Shi and Mengxiang Liu contributed equally and share first authorship.\n**Target journal:** BMC Medical Genomics")
    text = text.replace("## Introduction", "## Background")

    text = re.sub(
        r"### Supportive analyses and statistical reporting\n\n"
        r"High-confidence STRING associations.*?define the evidence sequence\.\n\n",
        "### Statistical analysis and reproducibility\n\n",
        text,
        flags=re.S,
    )

    # Cite the BMC additional files in the order in which they first appear.
    table_replacements = {
        "Table S1": "Additional file 1: Table S1",
        "Table S2": "Additional file 1: Table S2",
        "Tables S3a-b": "Additional file 1: Tables S3a and S3b",
        "Tables S12-S13": "Additional file 1: Tables S12 and S13",
        "Tables S6 and S14": "Additional file 1: Tables S6 and S14",
        "Table S8a": "Additional file 1: Table S8a",
        "Table S8b": "Additional file 1: Table S8b",
        "Table S9": "Additional file 1: Table S9",
        "Tables S9-S11": "Additional file 1: Tables S9-S11",
        "Tables S4-S6": "Additional file 1: Tables S4-S6",
        "Table S7": "Additional file 1: Table S7",
        "Tables S12-S14": "Additional file 1: Tables S12-S14",
    }
    # Longest keys first prevents partial substitutions.
    for old in sorted(table_replacements, key=len, reverse=True):
        text = text.replace(old, table_replacements[old])
    for number in range(1, 15):
        text = re.sub(
            rf"(?<!Additional file 1: )Table S{number}(?![a-zA-Z0-9])",
            f"Additional file 1: Table S{number}",
            text,
        )
    for number in range(1, 5):
        text = text.replace(
            f"Supplementary Figure {number}",
            f"Additional file 2: Figure S{number}",
        )

    text = apply_discovery_framing(polish_bmc_narrative(text))

    text = re.sub(r"\n## Declarations\n.*$", "", text, flags=re.S)
    abbreviations = """## Abbreviations

CI: confidence interval; FDR: false discovery rate; OA: osteoarthritis; OC: ovarian cancer; PCA: principal-component analysis; UMAP: Uniform Manifold Approximation and Projection; UMI: unique molecular identifier.

## Declarations

### Ethics approval and consent to participate

This secondary analysis used de-identified datasets available in public repositories. No new participants were recruited and no new ethics approval or consent was required for this study. Ethics approval and informed-consent procedures for the original studies are described in the corresponding repository records and publications.

### Consent for publication

Not applicable.

### Availability of data and materials

All datasets analysed in this study are publicly available from the NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/) under accessions GSE114007, GSE117999, GSE82107, GSE18520, GSE54388, GSE12470, GSE104782, GSE169454, GSE255460, GSE154600, GSE180661, GSE48556 and GSE31682. The derived data supporting the conclusions are supplied in Additional file 1, and the supplementary figures are supplied in Additional file 2. Versioned analysis scripts, configuration templates, source manifests and execution instructions are supplied in Additional file 3 and are publicly available at https://github.com/SMITHJUNHUI/oa-oc-cross-disease-transcriptomics. The immutable v1.0.0 release is archived in Zenodo (doi:10.5281/zenodo.21876012) [70]. No new restricted-access dataset was generated for this study. The reporting structure supports findable, accessible, interoperable and reusable data stewardship [61].

### Competing interests

The authors declare that they have no competing interests.

### Funding

This research was supported by the Science and Technology Program of XPCC (2023CB008-34), the Scientific Research Project of Shihezi University (ZZZC2023050), the Tianshan Talents Training Program for High-Level Medical and Health Talents (CZ001222), and the Clinical Basic Research Project of the First Affiliated Hospital of Shihezi University (LC2024009). The funders had no role in the study design, data analysis, interpretation, manuscript preparation or decision to submit.

### Authors' contributions

JHS and MXL contributed equally to this work and share first authorship. JHS contributed to conceptualization, methodology, software, formal analysis, data curation, visualization and writing: original draft. MXL contributed to methodology, validation, formal analysis, investigation, visualization and writing: review and editing. RI contributed to investigation, data curation, validation, literature review and writing: review and editing. KL contributed to methodology, validation, data curation and writing: review and editing. LC contributed to resources, supervision, project administration, funding acquisition and writing: review and editing. All authors read and approved the final manuscript.

### Acknowledgements

The authors thank the investigators and participants who generated and shared the public datasets used in this study.
"""
    text = text.rstrip() + "\n\n" + abbreviations.strip() + "\n"
    text = apply_human_style_polish(
        apply_acceptance_polish(apply_final_iteration(remap_citations(text)))
    )
    return ensure_cross_references(text)


def transform_legends(text: str) -> str:
    text = re.sub(r"(?<!Additional file 1: )Table S5", "Additional file 1: Table S5", text)
    text = re.sub(r"(?<!Additional file 1: )Table S8b", "Additional file 1: Table S8b", text)
    text = re.sub(r"(?<!Additional file 1: )Tables S9-S11", "Additional file 1: Tables S9-S11", text)
    text = replace_required(
        text,
        "**A-B,** OA GSE255460 and OC GSE154600 UMAPs. Six exact source labels are highlighted in each atlas; remaining source labels are grey. **C,** Fraction of cells with detected G0S2, EFEMP1, AKAP12, SOX9 or DDIT3 expression within selected source labels. Fractions are interpreted within each atlas and do not support functional inference. G0S2 was most frequently detected in ProC in the OA atlas and Myeloid.cell in the OC atlas. These atlas-specific observations do not identify the cell population responsible for the blood association.",
        "**A-B,** OA GSE255460 and OC GSE154600 UMAPs. Six exact source labels are highlighted in each atlas; remaining source labels are grey. **C,** Fraction of cells with detected G0S2, EFEMP1, AKAP12, SOX9 or DDIT3 expression within selected source labels. Fractions summarize localization within each atlas. G0S2 was most frequently detected in ProC in the OA atlas and Myeloid.cell in the OC atlas. These atlas-specific distributions provide cellular context for targeted follow-up.",
    )
    text = replace_required(
        text,
        "## Figure 6. Peripheral-blood screen retains G0S2 under independent FDR control",
        "## Figure 6. Stringent peripheral-blood validation identifies G0S2 as a restricted systemic signal",
    )
    text = replace_required(
        text,
        "**A,** Prespecified attrition from 286 shared tissue genes to one gene meeting the independent dual-cohort blood FDR rule. **B,** G0S2 log2 fold changes in OA tissue, OC tissue, OA peripheral-blood mononuclear cells and the OC blood-cell fraction. **C,** Within-cohort standardized blood effects shown as Hedges g with 95% confidence intervals. Effects were not pooled across diseases, and the analysis does not establish a shared circulating-cell source.",
        "**A,** Prespecified progression from 286 shared tissue genes to one gene meeting the independent dual-cohort blood FDR rule. **B,** G0S2 log2 fold changes in OA tissue, OC tissue, OA peripheral-blood mononuclear cells and the OC blood-cell fraction. **C,** Within-cohort standardized blood effects shown as Hedges g with 95% confidence intervals. Effects were analysed within each disease because the blood fractions differed.",
    )
    text = replace_required(
        text,
        "## Figure 7. Evidence sequence constrains cross-disease interpretation",
        "## Figure 7. Multi-context evidence distinguishes conserved from context-dependent transcriptomic features",
    )
    text = replace_required(
        text,
        "The sequential interpretation moves from tissue overlap to directional heterogeneity, external replication, source-defined cellular localization and blood persistence. G0S2 was the only gene retained after independent FDR control in both blood cohorts. The diagram treats this result as a blood-persistent association requiring confirmation, not as evidence of common pathogenesis or a validated circulating marker.",
        "The evidence sequence moves from tissue overlap to directional heterogeneity, external replication, source-defined cellular localization and blood persistence. Under independent FDR control, G0S2 was the sole feature retained in both blood cohorts. This restricted persistence defines the systemic component supported by the current data and prioritizes targeted prospective validation.",
    )
    text = text.replace(
        "## Figure 6. Stringent peripheral-blood validation identifies G0S2 as a restricted systemic signal",
        "## Figure 6. Stringent peripheral-blood validation identifies G0S2 as a candidate systemic feature",
    )
    text = text.replace(
        "## Figure 7. Multi-context evidence distinguishes conserved from context-dependent transcriptomic features",
        "## Figure 7. Multi-layer evidence maps shared molecular features across tissue, cellular and blood contexts",
    )
    text = text.replace(
        "The evidence sequence moves from tissue overlap to directional heterogeneity, external replication, source-defined cellular localization and blood persistence. Under independent FDR control, G0S2 was the sole feature retained in both blood cohorts. This restricted persistence defines the systemic component supported by the current data and prioritizes targeted prospective validation.",
        "The evidence sequence moves from shared tissue programs to external replication, disease-specific cellular mapping and systemic persistence. Directional analysis resolves how the shared features are regulated within each disease. Under independent FDR control, G0S2 emerged as the candidate systemic feature supported in both blood cohorts and prioritized for prospective validation.",
    )
    legacy_figure4 = "**A,** Proportion of measurable shared genes with the same direction as the corresponding OA or OC discovery contrast. Labels above bars report percentages; cohort-specific denominators are provided in Additional file 1: Table S5. **B,** Discovery and external-cohort log2 fold changes for G0S2, EFEMP1, AKAP12, SOX9 and DDIT3. **C,** Descriptive evidence summary showing external direction agreement, source-defined single-cell localization, dual-blood FDR status and interpretive role. The five genes illustrate distinct evidence patterns and were not used to define the primary analyses."
    compact_figure4 = "**A,** Proportion of measurable shared genes with the same direction as the corresponding OA or OC discovery contrast. Labels above bars report percentages; cohort-specific denominators are provided in Additional file 1: Table S5. **B,** Discovery and external-cohort log2 fold changes for G0S2, EFEMP1, AKAP12, SOX9 and DDIT3. The cross-layer descriptive evidence summary is provided in Additional file 2: Figure S5."
    text = text.replace(legacy_figure4, compact_figure4)
    text = text.replace(
        "The cross-layer descriptive evidence summary is provided in Supplementary Figure S5.",
        "The cross-layer descriptive evidence summary is provided in Additional file 2: Figure S5.",
    )
    if "The cross-layer descriptive evidence summary is provided in Additional file 2: Figure S5." not in text:
        raise RuntimeError("Figure 4 legend did not resolve to the final two-panel description")
    text = text.replace(
        "Directional analysis resolves how the shared features are regulated within each disease.",
        "Directional stratification resolves how the shared features are regulated within each disease.",
    )
    text = text.replace(
        "connector colour indicates matching or opposite signs.",
        "connector colour indicates matching or differing signs.",
    )
    text = text.replace(
        "Directional stratification resolves how the shared features are regulated within each disease.",
        "Directional stratification shows how the shared features are represented within each disease context.",
    )
    return text.strip()


def cited_numbers(text: str) -> set[int]:
    result: set[int] = set()
    for match in re.finditer(r"\[([0-9,\-]+)\]", text):
        result.update(expand_citation(match.group(1)))
    return result


def normalize(text: str) -> str:
    return re.sub(r"\n{3,}", "\n\n", text.strip()) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cover-letter", required=True, type=Path)
    parser.add_argument("--checklist", required=True, type=Path)
    args = parser.parse_args()

    source = args.source.read_text(encoding="utf-8")
    body, tail = source.split("\n## References\n", 1)
    references, legends_tail = tail.split("\n## Figure legends\n", 1)
    main_legends = legends_tail.split("\n## Supplementary Figure 1.", 1)[0]

    body = transform_main_body(body)
    references = transform_references(references)
    body, references = renumber_by_first_appearance(body, references)
    legends = transform_legends(main_legends)
    additional_files = """## Additional files

**Additional file 1.** XLSX. *Supplementary tables.* Tables S1-S14 provide cohort metadata, differential-expression results, threshold sensitivity, external replication, functional enrichment, illustrative-gene evidence, single-cell quality control and localization, and the peripheral-blood screening record.

**Additional file 2.** PDF. *Supplementary figures.* Figure S1 shows differential-expression threshold sensitivity; Figure S2 shows Hallmark states in discovery and external tissue cohorts; Figure S3 shows dataset-specific single-cell embeddings; Figure S4 shows discovery-cohort PCA and sample-correlation quality control; and Figure S5 summarizes the cross-layer evidence for the five illustrative genes.

**Additional file 3.** ZIP. *Reproducible analysis code.* Versioned R and Python scripts, configuration templates, source manifests, dependency records, tests and execution instructions required to regenerate the reported analyses from the public GEO inputs.
"""
    manuscript = normalize(
        body + "\n\n## References\n\n" + references + "\n\n## Figure legends\n\n" + legends + "\n\n" + additional_files
    )

    citations = cited_numbers(body)
    expected = set(range(1, 70))
    if citations != expected:
        raise RuntimeError(
            f"Citation coverage mismatch: missing={sorted(expected-citations)} extra={sorted(citations-expected)}"
        )
    if "STRING" in manuscript or "Supplementary table index" in manuscript:
        raise RuntimeError("Excluded supplementary-directory or STRING content remains")
    if len(re.findall(r"^## Figure \d+\.", manuscript, flags=re.M)) != 7:
        raise RuntimeError("Expected seven main figure legends")
    if len(re.findall(r"^\*\*(?:Background|Methods|Results|Conclusions)\.\*\*", manuscript, flags=re.M)) != 4:
        raise RuntimeError("Structured abstract labels are incomplete")

    cover_letter = """# Cover letter

11 August 2026

Editorial Office  
BMC Medical Genomics

Dear Editors,

Please consider our Research Article, “Shared molecular features between osteoarthritis and ovarian cancer revealed by multi-layer transcriptomic analyses,” for publication in BMC Medical Genomics.

We analysed osteoarthritis and ovarian cancer separately before comparing their transcriptional changes. The analysis identified 286 genes shared by the two discovery tissues, with enrichment for matrix, immune, stress and cell-cycle processes. Independent cohorts tested reproducibility, and single-cell data placed selected genes in disease-specific cellular environments. G0S2 was the only tissue-concordant gene that passed false-discovery-rate control in both blood cohorts.

The manuscript is relevant to BMC Medical Genomics because it examines a cross-disease molecular signal across tissue, single-cell and blood datasets while keeping the two diseases analytically separate. Source-data tables, supplementary figures and dataset accessions accompany the submission. The reproducible analysis code is public on GitHub and the immutable v1.0.0 release is archived in Zenodo (doi:10.5281/zenodo.21876012).

The manuscript has not been published previously and is not under consideration by another journal. All authors have read and approved the submitted version and agree to its submission to BMC Medical Genomics. The authors declare no competing interests. The work uses de-identified public datasets and did not recruit new participants. We are not aware of any issues relating to journal policies.

Thank you for your consideration.

Sincerely,

Lei Chen, Professor  
Department of Orthopedics Center, The First Affiliated Hospital of Shihezi University  
Shihezi University, Shihezi 832008, China  
Email: 564386249@qq.com  
Telephone: 13579758836
"""

    checklist = """# BMC Medical Genomics pre-submission checklist

Checked against the official BMC Medical Genomics Research Article and submission-guidelines pages on 11 August 2026.

- Research Article requirements: https://link.springer.com/journal/12920/submission-guidelines/research-article
- Journal submission guidelines: https://link.springer.com/journal/12920/submission-guidelines

## Completed in this package

- Research Article structure: Background, Methods, Results, Discussion and Conclusions.
- Structured abstract with Background, Methods, Results and Conclusions; below the 350-word limit.
- Six keywords (within the required range of 3-10).
- Exact BMC declaration headings, including “Not applicable” where appropriate.
- Abbreviations list.
- Seven main figures supplied as separate composite PDF files; legends remain in the main manuscript.
- Supplementary materials consolidated into three sequentially cited additional files.
- Supplementary figure directory and unprovided/low-priority STRING network materials removed.
- Main manuscript formatted as a double-spaced, page-numbered and continuously line-numbered DOCX.
- References follow the BMC Vancouver sequence, are numbered in order of first appearance, and all 69 DOI records have been verified against Crossref or DataCite.
- Main Figures 1-7, Supplementary Figures S1-S5 and Additional files 1-3 are all cited in the manuscript.
- GEO accessions are stated, and a peer-review-ready reproducible code archive is supplied as Additional file 3.
- Each additional file is below the 20 MB limit; each main figure is a one-page composite PDF below 10 MB.
- Main-figure titles and legends are retained in the manuscript; title and legend lengths meet the journal limits.
- The current article-processing charge and any institutional agreement, waiver or country-tier eligibility must be checked in the submission portal because the applicable amount is set by the acceptance date.

## Author sign-off required before upload

- Public repository and archive verified: https://github.com/SMITHJUNHUI/oa-oc-cross-disease-transcriptomics; immutable v1.0.0 DOI: https://doi.org/10.5281/zenodo.21876012; all-version concept DOI: https://doi.org/10.5281/zenodo.21876011.
- [ ] Funding: Lei Chen confirms that 2023CB008-34, ZZZC2023050, CZ001222 and LC2024009 supported this OA-OC study and that the no-role statement is accurate. Signature/date: __________
- [x] Junhui Shi (2223727941@qq.com): final manuscript approval, authorship contribution and no competing interests confirmed by the authors.
- [x] Mengxiang Liu (2877992646@qq.com): final manuscript approval, equal first authorship, authorship contribution and no competing interests confirmed by the authors.
- [x] Repkat Inayatilla (2047733903@qq.com): final manuscript approval, authorship contribution and no competing interests confirmed by the authors.
- [x] Ke Li (1413458714@qq.com): final manuscript approval, authorship contribution and no competing interests confirmed by the authors.
- [x] Lei Chen (564386249@qq.com): final manuscript approval, corresponding-authorship details and no competing interests confirmed by the authors.
- [ ] Confirm whether the target journal requests suggested/excluded reviewers in the submission portal.

## Upload map

1. Main manuscript: `OC_OA_BMC_Medical_Genomics_main_manuscript.docx`
2. Main figures: `Figure_1.pdf` to `Figure_7.pdf`
3. Additional file 1: `Additional_file_1_supplementary_tables.xlsx`
4. Additional file 2: `Additional_file_2_supplementary_figures.pdf`
5. Additional file 3: `Additional_file_3_reproducible_code.zip`
6. Cover letter: `OC_OA_BMC_Medical_Genomics_cover_letter.docx`
7. Author details reference (not uploaded unless requested): `OC_OA_BMC_Author_Information.docx`
"""

    for path in (args.output, args.cover_letter, args.checklist):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(manuscript, encoding="utf-8")
    args.cover_letter.write_text(normalize(cover_letter), encoding="utf-8")
    args.checklist.write_text(normalize(checklist), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
