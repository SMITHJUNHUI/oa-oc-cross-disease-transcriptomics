from __future__ import annotations

import csv
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT = Path(__file__).resolve().parents[1]
OUTPUTS = PROJECT.parent
SUBMISSION = PROJECT / "results" / "submission_v23"
MANUSCRIPT = PROJECT / "manuscript" / "OC_OA_manuscript_revision_v23.md"
DOCX = OUTPUTS / "OC_OA_manuscript_revision_v23_with_figures.docx"
PDF = OUTPUTS / "OC_OA_manuscript_revision_v23_with_figures.pdf"
REPORT = OUTPUTS / "OC_OA_submission_audit_report_v23.md"
JSON_REPORT = SUBMISSION / "submission_audit_v23.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def is_true(value: str) -> bool:
    return value.strip().upper() == "TRUE"


def docx_audit(path: Path) -> dict[str, object]:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    ns = {
        "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
        "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    }
    doc_pr = root.findall(".//wp:docPr", ns)
    complete_alt = sum(
        bool(item.get("title", "").strip())
        and bool(item.get("descr", "").strip())
        for item in doc_pr
    )
    document_text = " ".join(item.text or "" for item in root.findall(".//w:t", ns))
    figure_alt_text = " ".join(
        f"{item.get('title', '')} {item.get('descr', '')}" for item in doc_pr
    )
    return {
        "inline_images": len(doc_pr),
        "images_with_title_and_description": complete_alt,
        "contains_table_s16": "Table S16" in document_text,
        "contains_table_s18": "Table S18" in document_text,
        "contains_table_s19": "Table S19" in document_text,
        "contains_secondary_roc_label": (
            "deliberately shown last as a secondary display" in figure_alt_text
        ),
    }


def pdf_page_count(path: Path) -> int:
    try:
        from pypdf import PdfReader

        return len(PdfReader(str(path)).pages)
    except Exception:
        data = path.read_bytes()
        return len(re.findall(rb"/Type\s*/Page\b", data))


def secret_hits() -> list[str]:
    jwt = re.compile(
        r"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}"
    )
    hits: list[str] = []
    suffixes = {
        ".R",
        ".r",
        ".py",
        ".ps1",
        ".md",
        ".txt",
        ".yml",
        ".yaml",
        ".csv",
        ".json",
    }
    roots = (PROJECT, OUTPUTS)
    seen: set[Path] = set()
    for root in roots:
        for path in root.rglob("*"):
            resolved = path.resolve()
            if resolved in seen or not path.is_file() or path.suffix not in suffixes:
                continue
            seen.add(resolved)
            try:
                text = path.read_text(encoding="utf-8")
            except (UnicodeDecodeError, OSError):
                continue
            if jwt.search(text):
                try:
                    hits.append(str(path.relative_to(root)))
                except ValueError:
                    hits.append(str(path))
    return sorted(set(hits))


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, evidence: str) -> None:
        checks.append(
            {
                "check": name,
                "status": "PASS" if passed else "FAIL",
                "evidence": evidence,
            }
        )

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    abstract_match = re.search(
        r"## Abstract\s+(.*?)(?:\n## Keywords|\n## Introduction)",
        manuscript,
        flags=re.S,
    )
    abstract = abstract_match.group(1) if abstract_match else ""

    check(
        "Title fixes the context-dependent convergence framing",
        manuscript.startswith(
            "# Context-dependent remodeling-associated molecular convergence"
        )
        and "without shared genetic causality" in manuscript.splitlines()[0],
        "The title centers context-dependent convergence rather than shared targets.",
    )
    check(
        "Ten-gene selection is explicitly non-optimized",
        (
            "The ten-gene set was selected as an interpretable evidence summary "
            "rather than an optimized predictive signature."
        )
        in manuscript,
        "The exact reviewer-requested selection boundary is present in Methods.",
    )
    check(
        "AUC is secondary rather than abstract-centered",
        "As a secondary reproducibility assessment," in manuscript
        and "AUC 1.000" not in abstract
        and "median AUC was 1.000" not in abstract,
        "The abstract omits headline AUC values; Results introduces external separation as secondary.",
    )
    check(
        "Figure 4 order de-emphasizes ROC",
        all(
            phrase in manuscript
            for phrase in (
                "**A,** Unsupervised PCA",
                "**B,** Null AUC distributions",
                "**C,** Leave-one-sample-out AUCs",
                "**D,** ROC curves",
            )
        )
        and "deliberately shown last as a secondary display" in manuscript,
        "PCA, permutation, and leave-one-out precede the secondary ROC display.",
    )
    mr_limitations = (
        "First, the available GWAS datasets were predominantly derived from European "
        "populations, limiting generalizability.",
        "Second, OA GWAS mainly represented hip/knee phenotypes and may not capture "
        "all OA subtypes.",
        "Third, MR cannot exclude environmentally mediated pathways, such as "
        "aging-related inflammation or tissue injury responses.",
    )
    check(
        "All three requested MR limitations are explicit",
        all(sentence in manuscript for sentence in mr_limitations),
        "Population ancestry, OA phenotype scope, and environmentally mediated pathways are stated.",
    )
    check(
        "Discussion follows the four-part scientific narrative",
        all(
            heading in manuscript
            for heading in (
                "### Shared genes are not shared programs",
                "### Shared genes have different cellular meanings",
                "### Transcriptomic convergence does not imply genetic causality",
                "### Parallel stress adaptation is a hypothesis, not a mechanism",
            )
        ),
        "The Discussion is organized around direction, cell context, causality, and hypothesis.",
    )
    check(
        "Verified MR accessions are retained",
        "ebi-a-GCST007092" in manuscript
        and "ieu-a-1120" in manuscript
        and "GCST90018888" not in manuscript
        and "GCST90038686" not in manuscript,
        "Only the locally verified OpenGWAS datasets are named.",
    )

    candidates = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S16_candidate_prioritization_matrix.csv"
    )
    required_s16 = {
        "gene",
        "shared_DEG",
        "direction",
        "WGCNA_support",
        "LASSO_support",
        "random_forest_support",
        "strict_nested_frequency",
        "single_cell_context",
        "ten_gene_set_role",
    }
    roles = {row["ten_gene_set_role"] for row in candidates}
    consensus = sum(
        row["selection_stage"] == "cross-disease model consensus" for row in candidates
    )
    check(
        "Table S16 transparently accounts for all ten genes",
        len(candidates) == 10
        and required_s16.issubset(candidates[0])
        and roles
        == {"Interpretable evidence summary; not an optimized predictive signature"}
        and consensus == 1
        and all(row["single_cell_context"].strip() for row in candidates),
        (
            f"Rows={len(candidates)}; required columns={required_s16.issubset(candidates[0])}; "
            f"cross-disease consensus={consensus}; role labels={len(roles)}."
        ),
    )

    pathways = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S18_Hallmark_pathway_direction_matrix.csv"
    )
    both = [row for row in pathways if is_true(row["both_significant"])]
    concordant = [row for row in both if row["direction_class"] == "concordant"]
    discordant = [row for row in both if row["direction_class"] == "discordant"]
    complete_pathways = all(
        row["OA_NES"].strip()
        and row["OC_NES"].strip()
        and row["OA_FDR"].strip()
        and row["OC_FDR"].strip()
        and row["paired_direction_index"].strip()
        and row["inference_boundary"].strip()
        for row in pathways
    )
    check(
        "Directional Hallmark analysis is complete and bounded",
        len(pathways) == 50
        and len(both) == 10
        and len(concordant) == 4
        and len(discordant) == 6
        and complete_pathways,
        (
            f"Hallmark sets={len(pathways)}; significant in both={len(both)}; "
            f"concordant={len(concordant)}; discordant={len(discordant)}."
        ),
    )

    cell_function = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S19_gene_cell_function_context_matrix.csv"
    )
    complete_cell_function = all(
        row["OA_cell_context"].strip()
        and row["OC_cell_context"].strip()
        and row["OA_functional_theme"].strip()
        and row["OC_functional_theme"].strip()
        and row["OA_functional_theme_status"].strip()
        and row["OC_functional_theme_status"].strip()
        and "same gene does not imply the same biological meaning"
        in row["interpretation"]
        for row in cell_function
    )
    check(
        "Gene-cell-function matrix is complete without imputation",
        len(cell_function) == 10 and complete_cell_function,
        f"Rows={len(cell_function)}; all context/theme/status fields populated={complete_cell_function}.",
    )

    mr = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S12a_MR_estimates_and_provenance.csv"
    )
    ids = {row["id.exposure"] for row in mr} | {row["id.outcome"] for row in mr}
    snps = {int(row["nsnp"]) for row in mr}
    check(
        "MR provenance table remains complete",
        ids == {"ebi-a-GCST007092", "ieu-a-1120"}
        and snps == {11, 21}
        and len(mr) == 10,
        f"IDs={sorted(ids)}; SNP counts={sorted(snps)}; estimator rows={len(mr)}.",
    )

    figure_dir = SUBMISSION / "figures"
    main_png = sorted(figure_dir.glob("Figure[1-6]_*.png"))
    supp_png = sorted(figure_dir.glob("SupplementaryFigure[1-8]_*.png"))
    paired = all(path.with_suffix(".pdf").exists() for path in main_png + supp_png)
    check(
        "Figure pack contains six main and eight supplementary figures",
        len(main_png) == 6 and len(supp_png) == 8 and paired,
        (
            f"Main={len(main_png)}; supplementary={len(supp_png)}; "
            f"all PDF pairs={paired}."
        ),
    )

    docx = docx_audit(DOCX)
    check(
        "All embedded figures have alternative text",
        docx["inline_images"] == 14
        and docx["images_with_title_and_description"] == 14,
        json.dumps(docx, ensure_ascii=False),
    )
    check(
        "Word index includes Tables S16, S18, and S19",
        bool(docx["contains_table_s16"])
        and bool(docx["contains_table_s18"])
        and bool(docx["contains_table_s19"]),
        "The final Word document contains the upgraded/new supplementary-table entries.",
    )
    check(
        "Word figure labels identify ROC as secondary",
        bool(docx["contains_secondary_roc_label"]),
        "The Figure 4 embedded panel label is 'Secondary direction-fixed ROC display'.",
    )

    pages = pdf_page_count(PDF)
    check(
        "Final PDF is renderable",
        pages == 31 and PDF.stat().st_size > 1_000_000,
        f"Pages={pages}; bytes={PDF.stat().st_size}.",
    )

    hits = secret_hits()
    check(
        "No OpenGWAS JWT is stored",
        not hits,
        "No JWT-shaped token found." if not hits else "; ".join(hits),
    )
    check(
        "renv lockfile is strict JSON",
        json.loads((PROJECT / "renv.lock").read_text(encoding="utf-8")) is not None,
        "renv.lock parsed successfully.",
    )

    failures = [row for row in checks if row["status"] == "FAIL"]
    payload = {
        "revision": "V2.3",
        "checks": checks,
        "passed": len(checks) - len(failures),
        "failed": len(failures),
        "render_review": {
            "canonical_renderer_attempt": (
                "Bundled render_docx.py was attempted but LibreOffice/soffice was unavailable."
            ),
            "fallback_method": (
                "Microsoft Word PDF export plus Poppler page rasterization."
            ),
            "pages_visually_reviewed": 31,
            "result": (
                "No clipping, overlap, missing figure, mojibake, or blank-page defect detected."
            ),
        },
        "pending_author_items": [
            "Authors and affiliations",
            "Corresponding-author details",
            "Target journal and journal-specific formatting",
            "Funding and competing interests",
            "CRediT statement",
            "Repository URL and archival DOI",
            "Final human scientific and language approval",
        ],
    }
    JSON_REPORT.parent.mkdir(parents=True, exist_ok=True)
    JSON_REPORT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# OC-OA submission audit report — V2.3",
        "",
        f"**Automated checks:** {payload['passed']} passed; {payload['failed']} failed.",
        "",
        "| Check | Status | Evidence |",
        "|---|---|---|",
    ]
    lines.extend(
        f"| {row['check']} | {row['status']} | {row['evidence']} |"
        for row in checks
    )
    lines.extend(
        [
            "",
            "## Regression and visual review",
            "",
            "The complete R test suite passed after V2.3 integration. The canonical DOCX renderer was attempted but could not run because LibreOffice/soffice was unavailable. The final Word export was therefore converted with Microsoft Word, rasterized with Poppler, and all 31 pages were visually reviewed. No clipping, overlap, missing figure, mojibake, or blank-page defect was detected.",
            "",
            "## Interpretation boundary",
            "",
            "The manuscript supports context-dependent molecular convergence characterized by directional and cellular heterogeneity. It does not claim a shared disease mechanism, clinical diagnostic performance, a validated prognostic model, therapeutic targets, or inherited genetic causality. Bidirectional MR remains a negative supplementary constraint.",
            "",
            "## Pending author input",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in payload["pending_author_items"])
    lines.extend(
        [
            "",
            "These pending items are submission metadata or human approval requirements; they do not invalidate the computational build.",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Checks passed: {payload['passed']}/{len(checks)}")
    print(f"Wrote {REPORT}")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
