from __future__ import annotations

import argparse
import re
from pathlib import Path


SUPPLEMENTARY_INDEX = """
## Supplementary table index

- **Table S1:** Tissue cohorts, platforms, sample groups and analysis roles.
- **Table S2:** The 286 shared tissue DEGs with OA and OC effects and direction classes.
- **Tables S3a–b:** DEG-threshold sensitivity and membership.
- **Table S4:** External tissue gene-level effects.
- **Table S5:** External tissue direction-agreement summary.
- **Table S6:** External tissue Hallmark enrichment.
- **Table S7:** Five-gene descriptive evidence summary.
- **Table S8a:** Single-cell compatibility, quality control and permitted analysis layers.
- **Table S8b:** Representative-gene detection by exact source label.
- **Table S9:** Peripheral-blood cohort audit.
- **Table S10:** FDR-supported blood-replicated systemic component.
- **Table S11:** Prespecified blood-screen attrition.
- **Table S12:** GO enrichment of shared tissue DEGs.
- **Table S13:** KEGG enrichment of shared tissue DEGs.
- **Table S14:** Discovery-cohort Hallmark direction matrix.
- **Table S15:** WGCNA soft-power, bootstrap and leave-one-out stability.
- **Tables S16a–c:** STRING mapping, physical associations and node topology.
"""


RESPONSE = """
# V3.4 converged revision matrix

V3.4 restructures the study around a single evidence chain: shared tissue alterations, recurring biological themes, external tissue replication, cellular localization and conditional blood validation. Modules that did not strengthen this chain were removed from the submission narrative.

| Recommendation | V3.4 action | Location | Boundary retained |
|---|---|---|---|
| Stop adding loosely connected analyses | Reduced the manuscript to seven main figures and five supplementary figures; excluded prediction and speculative regulatory modules from the submitted evidence chain | Entire manuscript; scope ledger | Analysis volume is no longer used as a proxy for evidence strength |
| Reframe the scientific question | Replaced “shared mechanism” language with partial molecular convergence across tissue, cellular and blood contexts | Title, Abstract, Introduction, Discussion, Figure 7 | Observational transcriptomic overlap is not causal equivalence |
| Add external tissue validation | Recomputed all four external cohorts using gene-wise direction, external FDR and discovery–validation correlation | Figure 4; Tables S4–S6 | OA replication is explicitly reported as weaker and more cohort dependent than OC |
| Reduce the ten-gene panel | Replaced it with G0S2, EFEMP1, AKAP12, SOX9 and DDIT3 as a transparent descriptive evidence summary | Figure 4; Table S7 | The five genes are not an optimized signature and have no attached AUC |
| Strengthen biological interpretation | Consolidated GO and paired Hallmark results into matrix, immune, stress and cell-cycle themes | Figure 3; Tables S12–S14 | Enrichment is a transcriptional association, not protein-level activity |
| Make single-cell analysis answer a focused question | Used separate OA and OC count atlases, preserved exact source labels and displayed within-atlas detection fractions | Figure 5; Supplementary Figure 3; Tables S8a–b | No forced OA–OC integration or relabelling as homologous cell types |
| Keep blood as systemic-component validation | Applied a prespecified tissue-to-blood filter to GSE48556 and GSE31682; retained only the dual-FDR result | Figure 6; Tables S9–S11 | G0S2 is a blood-replicated component, not a biomarker or diagnostic panel |
| Search for independent blood confirmation without forcing it | Audited GSE163552 but did not add it: the record contains only three OA and three healthy PBMC RNA-seq samples focused on lncRNA and no compatible independent OC blood counterpart was identified | Audit trail; limitations | An incompatible or unpaired cohort is not treated as validation |
| De-emphasize WGCNA and STRING | Moved both analyses to supplementary evidence and bounded their claims | Supplementary Figures 1 and 5; Tables S15–S16 | Co-expression and database associations do not establish shared mechanism |
| Standardize figures | Rebuilt main figures at 183-mm width with consistent typography, restrained blue/orange encoding, fixed panel spacing and 600-dpi TIFF plus PDF/SVG/PNG exports | Figures 1–7 | Styling follows common biomedical-journal conventions rather than imitating one paper |
| Reduce formulaic or AI-like prose | Rewrote transitions, removed repetitive defensive lists, varied paragraph structure and kept claims close to data | Entire manuscript | Language remains precise without overstating certainty |

## Repository records consulted for the blood audit

- GSE48556: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48556
- GSE31682: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31682
- GSE163552: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163552
"""


def build(template: Path, legends_path: Path) -> str:
    body = template.read_text(encoding="utf-8").strip()
    legends = legends_path.read_text(encoding="utf-8").strip()
    body += "\n\n## Figure legends\n\n" + legends
    body += "\n\n" + SUPPLEMENTARY_INDEX.strip() + "\n"
    return re.sub(r"\n{3,}", "\n\n", body)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--figure-legends", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--response-output", required=True, type=Path)
    args = parser.parse_args()
    manuscript = build(args.template, args.figure_legends)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(manuscript, encoding="utf-8")
    args.response_output.write_text(RESPONSE.strip() + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.response_output}")


if __name__ == "__main__":
    main()
