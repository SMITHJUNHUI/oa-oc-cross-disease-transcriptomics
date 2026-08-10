from __future__ import annotations

import argparse
import re
from pathlib import Path

import build_manuscript_v41 as v41


DISCUSSION_V42 = r"""## Discussion

OA and OC showed a limited transcriptomic overlap rather than a uniform shared disease state. The 286-gene intersection was stable to prespecified thresholds, but concordant and discordant directions were nearly balanced. External replication further constrained the result. OC cohorts reproduced most discovery directions, whereas OA replication was weaker and more cohort dependent. Differences in cartilage source, joint site, disease severity, platform and reference groups may contribute to this contrast [2,3,8,16,18,23,24,32-34].

Recurring matrix, immune, stress and cell-cycle categories therefore should not be interpreted as common pathogenesis. Six jointly significant Hallmark pathways had opposite enrichment directions, and the illustrative genes localized to different source-defined populations in separate OA and OC atlases. These observations identify contextual heterogeneity but do not establish different gene functions or homologous cell states [9-14,26-30,35,36,38,44,45].

The blood screen retained G0S2 because it was lower in both tissues, had the same sign in both blood cohorts and passed FDR control independently in each cohort. The other same-sign genes did not meet this complete rule. This attrition may reflect biological heterogeneity, blood-cell composition, measurement sensitivity or cohort differences. G0S2 was most frequently detected in OA ProC and OC Myeloid.cell populations, but neither atlas identifies the cells responsible for the blood association. Reported roles in lipolysis, quiescence and growth regulation provide biological context, not a mechanism for the four observed contrasts [51-54]. G0S2 is therefore a blood-persistent transcriptomic association that requires cell-resolved, protein-level and prospective confirmation.

The study is limited by retrospective public data and incomplete covariate annotation. Reliable technical batch covariates were not consistently available. No explicit batch-correction method was applied, and cohorts were not merged or adjusted with ComBat. Age, sex, disease stage, treatment, joint site, tumour histology and reference tissue could not be harmonized completely. Single-cell protocols and annotations differed, and the two blood cohorts used different platforms and cell fractions. No second compatible OA-OC blood pair was available. These constraints limit mechanistic and translational interpretation.
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
- **Table S10:** FDR-supported blood-persistent association.
- **Table S11:** Prespecified blood-screen attrition.
- **Table S12:** Gene Ontology enrichment of shared tissue genes.
- **Table S13:** Kyoto Encyclopedia of Genes and Genomes enrichment of shared tissue genes.
- **Table S14:** Discovery-cohort Hallmark direction matrix.
- **Tables S15a-c:** STRING mapping, physical associations and node topology.
"""


REVISION_MATRIX = r"""# V4.2 reviewer-iteration revision matrix

| Review concern | V4.2 action | Resolution test |
|---|---|---|
| Shared overlap could still be read as common pathogenesis | Added an explicit distinction between gene-list overlap and a uniform shared disease state in the Abstract, Results, Discussion and Conclusions | Every headline interpretation treats directional heterogeneity as the boundary of the overlap |
| Batch handling was ambiguous | Stated that cohorts were processed and modelled separately, no ComBat or cross-cohort batch correction was applied, and the disease group was the only design factor | Methods now distinguishes scale transformation, unsupervised QC and statistical adjustment |
| G0S2 was over-emphasized | Reframed G0S2 as a blood-persistent transcriptomic association; added its atlas-specific localization, alternative explanations and unresolved circulating-cell source | No biomarker, causal, shared-cell-origin or mechanism claim remains |
| WGCNA and STRING increased analytical weight | Removed WGCNA from the manuscript, supplementary figure package and table index; retained STRING only as auxiliary database context | The submission contains one threshold-sensitivity panel and no WGCNA table or figure |
| Discussion exceeded the evidence | Rewrote the Discussion as four compact paragraphs covering the central boundary, context, G0S2 and limitations | Interpretation is shorter and each paragraph has one evidentiary job |
"""


def expand_citation(content: str) -> list[int]:
    numbers: list[int] = []
    for part in content.split(","):
        if "-" in part:
            start, end = (int(value) for value in part.split("-", 1))
            numbers.extend(range(start, end + 1))
        else:
            numbers.append(int(part))
    return numbers


def compress_numbers(numbers: list[int]) -> str:
    ordered = sorted(set(numbers))
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
        old_numbers = [number for number in expand_citation(match.group(1)) if number != 60]
        mapped = [number if number < 60 else number - 1 for number in old_numbers]
        if not mapped:
            raise RuntimeError("A citation contained only the removed WGCNA reference.")
        return f"[{compress_numbers(mapped)}]"

    return re.sub(r"\[([0-9,\-]+)\]", replace, text)


def build_manuscript_text() -> str:
    text = v41.MANUSCRIPT
    text = text.replace(
        "# Shared molecular features between osteoarthritis and ovarian cancer across tissue, cellular and systemic contexts",
        "# Directionally heterogeneous transcriptomic overlap between osteoarthritis and ovarian cancer across tissue, cellular and blood contexts",
    )
    text = text.replace(
        "**Running title:** Shared molecular features across OA and OC",
        "**Running title:** Directionally heterogeneous OA-OC overlap",
    )
    text = text.replace(
        "It remains unclear whether their transcriptomic overlap persists across independent tissues, source-defined cell populations and peripheral blood.",
        "It remains unclear whether their transcriptomic overlap defines a reproducible cross-disease state across independent tissues, source-defined cell populations and peripheral blood.",
    )
    text = text.replace(
        "OA and OC share a limited transcriptional background whose direction, reproducibility and cellular localization depend on disease context. G0S2 is a candidate systemic molecular signal that requires prospective and protein-level confirmation.",
        "OA and OC share a directionally heterogeneous transcriptomic overlap rather than a uniform disease state. Replication and cellular localization remain context dependent. G0S2 is a blood-persistent transcriptomic association that requires prospective and protein-level confirmation.",
    )
    text = text.replace(
        "Each disease was estimated separately, and evidence was integrated only after disease-specific analysis.",
        "Each disease was estimated separately, and evidence was integrated only after disease-specific analysis. This design treated overlap as a hypothesis to be tested rather than evidence of common pathogenesis.",
    )
    text = text.replace(
        "Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Multiple probes for one gene were collapsed according to the dataset-specific record. Normalization and modelling were performed within each cohort, with batch structure considered during processing [63,70]. Principal-component and sample-correlation summaries are provided in Supplementary Figure 4.",
        "Expression matrices were checked for identifier consistency, duplicated features, non-finite values and group balance. Public processed matrices were log2-transformed only when indicated by their scale, and multiple probes for one gene were collapsed according to the dataset-specific record. Reliable technical batch covariates were not consistently available. No explicit batch-correction method was applied, and cohorts were not merged or adjusted with ComBat. Each cohort was modelled separately with disease group as the only design factor. Principal-component and sample-correlation summaries assessed cohort structure but did not trigger outcome-informed sample removal or adjustment (Supplementary Figure 4) [63].",
    )
    text = text.replace(
        "Co-expression stability and high-confidence STRING associations were evaluated as supportive analyses and are reported only in Supplementary Figure 1, Supplementary Figure 5 and Tables S15-S16 [60,61]. These analyses did not select the illustrative genes.",
        "High-confidence STRING associations were evaluated only as auxiliary database context and are reported in Supplementary Figure 5 and Tables S15a-c [61]. This analysis did not select the illustrative genes or define the evidence sequence.",
    )
    text = text.replace(
        "Tests were two-sided unless otherwise stated. Multiple-testing procedures and direction definitions were fixed before biological interpretation. Sample and cell units are reported for each evidence layer, and denominators accompany filtering steps. Randomized procedures used seed 20260726. Full parameters, source manifests and executable scripts are included in the reproducible project.",
        "Tests were two-sided unless otherwise stated. Multiple-testing procedures and direction definitions were fixed before biological interpretation. Sample and cell units are reported for each evidence layer, and denominators accompany filtering steps. Randomized procedures used seed 20260726. Bioconductor resources supported the analysis environment [70]. Full parameters, source manifests and executable scripts are included in the reproducible project.",
    )
    text = text.replace(
        "The shared set therefore contained 146 concordant genes (51.0%) and 140 discordant genes (49.0%).",
        "The shared set therefore contained 146 concordant genes (51.0%) and 140 discordant genes (49.0%). Thus, gene-list overlap did not define a uniform cross-disease expression state.",
    )
    text = text.replace(
        "### Peripheral blood evaluation retained only G0S2",
        "### Peripheral-blood evaluation retained only G0S2 under the dual-cohort FDR rule",
    )
    text = text.replace(
        "The blood screen reduced the tissue-derived set to one candidate systemic molecular signal (Figure 6; Tables S9-S11).",
        "The blood screen reduced the tissue-derived set to one blood-persistent transcriptomic association (Figure 6; Tables S9-S11).",
    )
    text = text.replace(
        "The result was consistent across direction and FDR criteria despite the different blood fractions.",
        "The result met the prespecified direction and FDR criteria despite the different blood fractions, but it did not establish a shared cellular source.",
    )
    text = re.sub(r"## Discussion\n.*?\n## Conclusions\n", DISCUSSION_V42 + "\n## Conclusions\n", text, flags=re.S)
    text = text.replace(
        "OA and OC shared 286 tissue transcriptional alterations, but their direction, external reproducibility and cellular localization were context dependent. Only G0S2 persisted through the independent dual-cohort blood screen. These findings define partial molecular convergence across tissue, cellular and systemic contexts while setting clear limits on mechanistic and translational interpretation.",
        "OA and OC shared 286 tissue transcriptional alterations, but the near-balanced direction classes did not define a uniform shared disease state. External reproducibility and cellular localization remained context dependent. Only G0S2 passed the independent dual-cohort blood rule. This result is a blood-persistent association, not evidence of common pathogenesis or a validated circulating marker.",
    )
    text = text.replace("The complete V4.1 workflow", "The complete V4.2 workflow")
    return remap_citations(text)


def normalize(text: str) -> str:
    return re.sub(r"\n{3,}", "\n\n", text.strip()) + "\n"


def cited_numbers(text: str) -> set[int]:
    result: set[int] = set()
    for match in re.finditer(r"\[([0-9,\-]+)\]", text):
        result.update(expand_citation(match.group(1)))
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
    if reference_count != 69:
        raise RuntimeError(f"Expected 69 references, found {reference_count}")
    manuscript_body = build_manuscript_text()
    cited = cited_numbers(manuscript_body)
    if cited != set(range(1, 70)):
        raise RuntimeError(f"Citation coverage mismatch: missing={sorted(set(range(1,70))-cited)} extra={sorted(cited-set(range(1,70)))}")

    legends = args.figure_legends.read_text(encoding="utf-8").strip()
    legend_count = len(re.findall(r"^## (?:Supplementary )?Figure \d+\.", legends, flags=re.MULTILINE))
    if legend_count != 12:
        raise RuntimeError(f"Expected 12 figure legends, found {legend_count}")

    template = normalize(manuscript_body + "\n\n## References\n\n" + references)
    manuscript = normalize(template + "\n\n## Figure legends\n\n" + legends + "\n\n" + SUPPLEMENTARY_INDEX)
    for path in (args.template_output, args.output, args.response_output):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.template_output.write_text(template, encoding="utf-8")
    args.output.write_text(manuscript, encoding="utf-8")
    args.response_output.write_text(normalize(REVISION_MATRIX), encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
