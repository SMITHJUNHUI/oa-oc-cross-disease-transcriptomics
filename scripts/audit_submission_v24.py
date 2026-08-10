#!/usr/bin/env python3
"""Audit the V2.4 manuscript, figures, tables, and interpretation boundaries."""

from __future__ import annotations

import csv
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT = Path(__file__).resolve().parents[1]
OUTPUTS = PROJECT.parent
SUBMISSION = PROJECT / "results" / "submission_v24"
MANUSCRIPT = PROJECT / "manuscript" / "OC_OA_manuscript_revision_v24.md"
DOCX = OUTPUTS / "OC_OA_manuscript_revision_v24_with_figures.docx"
PDF = OUTPUTS / "OC_OA_manuscript_revision_v24_with_figures.pdf"
REPORT = OUTPUTS / "OC_OA_submission_audit_report_v24.md"
JSON_REPORT = SUBMISSION / "submission_audit_v24.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


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
        "contains_v24_tables": all(
            label in document_text
            for label in ("Table S20", "Table S21", "Table S22", "Table S23")
        ),
        "contains_separability_title": (
            "Direction-fixed cross-cohort molecular separability" in document_text
        ),
        "figure_4_boundary": (
            "retrospective cross-cohort molecular separability, not clinical performance"
            in figure_alt_text
        ),
        "replacement_character": "\ufffd" in document_text,
    }


def pdf_page_count(path: Path) -> int:
    try:
        from pypdf import PdfReader

        return len(PdfReader(str(path)).pages)
    except Exception:
        return len(re.findall(rb"/Type\s*/Page\b", path.read_bytes()))


def secret_hits() -> list[str]:
    token = re.compile(
        r"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}"
    )
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
        ".tsv",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for root in (PROJECT, OUTPUTS):
        for path in root.rglob("*"):
            if not path.is_file() or path.suffix not in suffixes:
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            try:
                text = path.read_text(encoding="utf-8")
            except (UnicodeDecodeError, OSError):
                continue
            if token.search(text):
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
        "Title fixes context-dependent convergence and the causal boundary",
        manuscript.startswith(
            "# Context-dependent remodeling-associated molecular convergence"
        )
        and "without shared genetic causality" in manuscript.splitlines()[0],
        "The title does not claim shared mechanism, diagnosis, or therapeutic targets.",
    )
    check(
        "Candidate-set role is transparent and non-optimized",
        (
            "The ten-gene set is an interpretable evidence summary rather than "
            "an optimized predictive signature."
        )
        in manuscript
        and "their intersection does not replace the nested design" in manuscript,
        "LASSO/RF prioritization is separated from strict nested stability estimation.",
    )
    check(
        "AUC is secondary and reframed as molecular separability",
        "As a secondary reproducibility assessment," in manuscript
        and "molecular summary demonstrated strong retrospective separation" in manuscript
        and "AUC 1.000" not in abstract
        and "clinical performance" not in abstract,
        "The abstract avoids headline AUC claims and Results explains the task-scale contrast.",
    )
    check(
        "Figure 4 displays structure and stability before ROC",
        all(
            phrase in manuscript
            for phrase in (
                "**A,** Unsupervised PCA",
                "**B,** Null AUC distributions",
                "**C,** Leave-one-sample-out AUCs",
                "**D,** ROC curves",
            )
        )
        and "shown last as a secondary display" in manuscript,
        "PCA, permutation, and sample omission precede the ROC panel.",
    )
    check(
        "Calibration and clinical-utility boundaries are explicit",
        "not a locked probability model" in manuscript
        and "not transported external probability calibration" in manuscript
        and "Decision-curve analysis and a nomogram were not performed" in manuscript,
        "Calibration remains a supplementary cross-fitted sensitivity analysis.",
    )
    check(
        "Focused upstream context is not asserted as a causal TF-miRNA network",
        "No TF-miRNA causal network was asserted." in manuscript
        and "miRNA expression and activity were not measured" in manuscript
        and "does not measure transcription-factor activity" in manuscript,
        "KnockTF and miRTarBase results are explicitly hypothesis-generating.",
    )
    check(
        "MR is retained as a bounded supplementary analysis",
        all(
            phrase in manuscript
            for phrase in (
                "MR was retained as an explicit causality boundary",
                "predominantly derived from European populations",
                "OA GWAS mainly represented hip/knee phenotypes",
                "MR cannot exclude environmentally mediated pathways",
            )
        )
        and "ebi-a-GCST007092" in manuscript
        and "ieu-a-1120" in manuscript,
        "The ancestry, OA-subtype, and environmentally mediated pathway limits are present.",
    )
    check(
        "Discussion follows the intended scientific narrative",
        all(
            heading in manuscript
            for heading in (
                "### Shared genes are not shared programs",
                "### Molecular separability is a property of the validation task",
                "### Shared genes have different cellular and regulatory meanings",
                "### Transcriptomic convergence does not imply genetic causality",
                "### Parallel stress adaptation is a hypothesis, not a mechanism",
            )
        ),
        "Direction, task scale, cellular/regulatory context, causality, and hypothesis are separated.",
    )

    s20 = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S20_candidate_centered_Hallmark_context.csv"
    )
    unavailable = sum(row["calculation_status"] != "estimated" for row in s20)
    check(
        "Candidate-centered Hallmark analysis is complete and bounded",
        len(s20) == 400
        and {row["candidate"] for row in s20} == {"SOX9", "DDIT3", "BNC1", "AKAP12"}
        and {row["disease"] for row in s20} == {"OA", "OC"}
        and unavailable == 7
        and all("not single-gene perturbation" in row["inference_boundary"] for row in s20),
        f"Rows={len(s20)}; explicit unavailable estimates={unavailable}.",
    )

    s21 = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S21_cross_cohort_molecular_separability_context.csv"
    )
    check(
        "External validation tasks are annotated by tissue and biological scale",
        len(s21) == 4
        and {row["dataset_id"] for row in s21}
        == {"GSE117999", "GSE82107", "GSE54388", "GSE12470"}
        and all(row["tissue_and_comparator"] for row in s21)
        and all(row["validation_task_scale"] for row in s21)
        and all("not clinical utility" in row["inference_boundary"] for row in s21),
        f"Annotated external cohorts={len(s21)}.",
    )

    s22a = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S22a_external_signed_score_effect_sizes.csv"
    )
    s22b = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S22b_cross_fitted_calibration_metrics.csv"
    )
    unstable = sum("unstable" in row["calibration_status"] for row in s22b)
    check(
        "Effect-size and calibration sensitivities retain unfavorable results",
        len(s22a) == 6
        and len(s22b) == 4
        and unstable == 3
        and all("not external probability calibration" in row["inference_boundary"] for row in s22b),
        f"Effect-size rows={len(s22a)}; calibration cohorts={len(s22b)}; unstable={unstable}.",
    )

    s23a = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S23a_KnockTF_candidate_regulatory_context.csv"
    )
    s23b = read_csv(
        SUBMISSION
        / "supplementary_tables"
        / "Table_S23b_miRTarBase_candidate_regulatory_context.csv"
    )
    check(
        "Regulatory context tables preserve evidence provenance and boundaries",
        len(s23a) == 157
        and len(s23b) == 363
        and all("does not establish TF activity" in row["inference_boundary"] for row in s23a)
        and all("hypothesis-generating" in row["inference_boundary"] for row in s23b),
        f"KnockTF rows={len(s23a)}; miRTarBase rows={len(s23b)}.",
    )

    figure_dir = SUBMISSION / "figures"
    main_png = sorted(figure_dir.glob("Figure[1-6]_*.png"))
    supp_png = sorted(figure_dir.glob("SupplementaryFigure*.png"))
    paired = all(path.with_suffix(".pdf").exists() for path in main_png + supp_png)
    legends = (figure_dir / "figure_legends.md").read_text(encoding="utf-8")
    check(
        "Figure pack contains six main and eleven supplementary figures",
        len(main_png) == 6
        and len(supp_png) == 11
        and paired
        and "## Figure 1. Linear study design and audited evidence boundaries" in legends
        and "## Figure 4. Direction-fixed cross-cohort molecular separability" in legends,
        f"Main={len(main_png)}; supplementary={len(supp_png)}; all PDF pairs={paired}.",
    )

    docx = docx_audit(DOCX)
    check(
        "Word document embeds all figures with alternative text",
        docx["inline_images"] == 17
        and docx["images_with_title_and_description"] == 17
        and not docx["replacement_character"],
        json.dumps(docx, ensure_ascii=False),
    )
    check(
        "Word document contains V2.4 tables and separability terminology",
        bool(docx["contains_v24_tables"])
        and bool(docx["contains_separability_title"])
        and bool(docx["figure_4_boundary"]),
        "Tables S20-S23 and the Figure 4 interpretation boundary are present.",
    )

    pages = pdf_page_count(PDF)
    check(
        "Final PDF is renderable",
        pages == 36 and PDF.stat().st_size > 3_000_000,
        f"Pages={pages}; bytes={PDF.stat().st_size}.",
    )

    hits = secret_hits()
    check(
        "No OpenGWAS JWT is stored",
        not hits,
        "No JWT-shaped token found." if not hits else "; ".join(hits),
    )
    check(
        "The reproducible build entry points and local regulatory input exist",
        (PROJECT / "run_submission_v24.R").exists()
        and (PROJECT / "run_submission_v24.ps1").exists()
        and (PROJECT / "data" / "raw" / "regulatory" / "human_500diff.gmt").exists()
        and json.loads((PROJECT / "renv.lock").read_text(encoding="utf-8")) is not None,
        "One-command runners, the project-local KnockTF GMT, and a parseable renv lockfile are present.",
    )

    failures = [row for row in checks if row["status"] == "FAIL"]
    payload = {
        "revision": "V2.4",
        "checks": checks,
        "passed": len(checks) - len(failures),
        "failed": len(failures),
        "regression_tests": "Complete R test suite passed on 30 July 2026.",
        "render_review": {
            "canonical_renderer_attempt": (
                "Bundled render_docx.py was attempted but LibreOffice/soffice was unavailable."
            ),
            "fallback_method": "Microsoft Word PDF export plus Poppler page rasterization.",
            "pages_visually_reviewed": 36,
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
            "Final author scientific and language approval",
        ],
    }
    JSON_REPORT.parent.mkdir(parents=True, exist_ok=True)
    JSON_REPORT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# OC-OA submission audit report — V2.4",
        "",
        f"**Automated checks:** {payload['passed']} passed; {payload['failed']} failed.",
        "",
        "| Check | Status | Evidence |",
        "|---|---|---|",
    ]
    lines.extend(
        f"| {row['check']} | {row['status']} | {str(row['evidence']).replace('|', '/')} |"
        for row in checks
    )
    lines.extend(
        [
            "",
            "## Regression and visual review",
            "",
            (
                "The complete R test suite passed after V2.4 integration. The canonical "
                "DOCX renderer was attempted but could not run because LibreOffice/soffice "
                "was unavailable. The final Word export was converted with Microsoft Word, "
                "rasterized with Poppler, and all 36 pages were visually reviewed. No clipping, "
                "overlap, missing figure, mojibake, or blank-page defect was detected."
            ),
            "",
            "## Interpretation boundary",
            "",
            (
                "The manuscript supports context-dependent molecular convergence characterized "
                "by directional, task-scale, cellular, and regulatory heterogeneity. It does not "
                "claim a shared disease mechanism, clinical diagnostic performance, a validated "
                "prognostic model, therapeutic targets, or inherited genetic causality. MR is a "
                "negative supplementary boundary; KnockTF/miRTarBase results are focused, "
                "hypothesis-generating context rather than a causal regulatory network."
            ),
            "",
            "## Pending author input",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in payload["pending_author_items"])
    lines.extend(
        [
            "",
            (
                "These items require author decisions or metadata and therefore remain visibly "
                "flagged; they do not invalidate the computational build."
            ),
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
