from __future__ import annotations

import argparse
import re
from pathlib import Path


SUPPLEMENTARY_INDEX = """
## Supplementary table index

- **Table S1:** Tissue data sources, cohorts, platforms, and analysis roles.
- **Table S2:** The 286 shared tissue DEGs with OA and OC effects and direction classes.
- **Tables S3a-b:** DEG-threshold sensitivity and gene membership.
- **Table S4:** WGCNA soft-power, bootstrap, and leave-one-out stability.
- **Table S9:** Single-cell QC, compatibility, and permitted analysis layers.
- **Table S10:** Representative-gene single-cell evidence and eligible pseudobulk summaries.
- **Tables S11a-c:** Hallmark, GO, and KEGG enrichment results.
- **Table S12:** Representative shared-gene direction and cell-context summary.
- **Table S13:** Officially audited blood cohorts, platforms, groups, and mapping totals.
- **Table S14:** The FDR-supported systemic-component result (G0S2).
- **Table S15:** Prespecified blood-screen attrition counts.
- **Table S18:** Complete Hallmark pathway-direction matrix.
- **Table S19:** Gene-cell-function context matrix using exact source labels.
- **Tables S24a-c:** Cell-context specificity and sample-aware UCell summaries.
- **Tables S25a-e:** STRING mapping, edges, topology, and direction-label permutation audit.
- **Table S30a:** Extended representative-gene cell-detection results.
- **Tables S31a-b:** Bulk discovery PCA and sample-correlation QC.
"""


RESPONSE = """
# V3.3 targeted revision response matrix

V3.3 implements the final rule that peripheral blood is a conditional systemic-component validation layer, not a biomarker-discovery, diagnostic-model, or independent research direction. Official GEO metadata were treated as authoritative when the proposed accession descriptions conflicted with repository records.

| Recommendation or issue | V3.3 action | Location | Scientific boundary |
|---|---|---|---|
| Add blood only as systemic-component validation | Added a prespecified tissue-to-blood filter after the direction analysis; no blood feature discovery, classifier, AUC, nomogram, or decision curve was performed | Methods; Results; Figure 5; Tables S13-S15 | Blood does not become a separate disease-biomarker track |
| Use OA peripheral blood | Retained GSE48556 after official audit: PBMC, GPL6947, 106 OA cases and 33 healthy controls | Methods; Figure 5A; Table S13 | Corrected the proposed description: this is PBMC, not whole blood or GPL10558 |
| Use OC peripheral blood | Replaced ineligible GSE69428 with GSE31682: blood-cell fraction, GPL2986, 48 epithelial OC cases and 20 healthy controls | Methods; Figure 5A; Table S13 | Platforms and blood fractions differ, so cohorts were modeled separately and not pooled |
| Audit proposed alternatives | Rejected GSE51588 (OA subchondral bone), GSE69428 (fallopian-tube/ovarian tissue), and GSE112790 (hepatocellular-carcinoma tissue) from the blood module | Methods; project audit trail | No result was forced from a biologically ineligible accession |
| Compare tissue and blood direction | Screened only shared tissue DEGs and required concordant OA/OC tissue signs, measurement in both blood cohorts, and matching signs across all four contexts | Methods; Figure 5B-C; Table S15 | Cross-context sign agreement is descriptive, not evidence of a common cause |
| Retain only positive results | Named only the result passing FDR <0.05 separately in both blood cohorts; nominal-only genes remain internal | Results; Figure 5; Table S14 | Stage-wise denominators are still reported to prevent selective-reporting bias |
| Report the final positive result | G0S2 was lower in both tissues and both blood cohorts and showed moderate standardized blood effects in each disease | Results; Figure 5C-D; Table S14 | G0S2 is a replication candidate, not a biomarker, mechanism, or target |
| Preserve the tissue-divergence story | Kept 286 shared DEGs, four direction quadrants, GO/KEGG, and paired Hallmark directions as the primary evidence | Figures 2-4; Tables S2, S11, and S18 | 140/286 genes and 6/10 joint Hallmarks were discordant; broad conservation is not claimed |
| Keep single-cell interpretation exact | Retained separate OA and OC atlases and exact source labels; no forced cross-disease integration or biological relabeling | Figure 6; Supplementary Figure 3; Tables S9, S10, S19, S24, and S30 | Fibroblast is not relabeled as CAF; Ovarian.cancer.cell is not relabeled as tumor epithelium |
| Stop analysis accumulation | Excluded classification/AUC, MR, CellChat/NicheNet, TF-miRNA, DCA/nomogram, and other nonessential exploratory modules from V3.3 submission | Scope decisions; manuscript | Earlier outputs remain preserved for provenance but do not enter the submitted evidence chain |
| Align figures with common journal style | Used seven 185-mm main figures and five supplementary figures with white backgrounds, restrained typography, consistent panel lettering, and color-blind-aware blue/orange encoding | Figures 1-7 and Supplementary Figures 1-5 | Styling is general journal format; target-journal production requirements remain to be applied after journal selection |
| Keep claims proportional | Rewrote the title, abstract, discussion, limitations, and conclusion around tissue divergence plus one limited blood-replicated component | Entire manuscript; Figure 7 | No shared mechanism, causal relationship, diagnostic panel, or therapeutic target claim is made |

## Accession audit links

- GSE48556: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48556
- GSE31682: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31682
- GSE51588: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51588
- GSE69428: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69428
- GSE112790: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112790
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
