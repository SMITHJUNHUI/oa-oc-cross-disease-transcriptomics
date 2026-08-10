from __future__ import annotations

import argparse
import csv
import json
import re
import zipfile
from pathlib import Path

from pypdf import PdfReader


MAIN_FIGURES = [
    "Figure1_study_design",
    "Figure2_shared_tissue_transcriptomics",
    "Figure3_shared_biological_themes",
    "Figure4_external_tissue_and_illustrative_genes",
    "Figure5_single_cell_localization",
    "Figure6_peripheral_blood_validation",
    "Figure7_integrated_interpretation",
]

SUPPLEMENTARY_FIGURES = [
    "SupplementaryFigure1_DEG_threshold_sensitivity",
    "SupplementaryFigure2_external_tissue_Hallmark",
    "SupplementaryFigure3_all_single_cell_UMAPs",
    "SupplementaryFigure4_bulk_PCA_QC",
    "SupplementaryFigure5_direction_aware_STRING_network",
]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def pdf_text(path: Path) -> str:
    return "\n".join((page.extract_text() or "") for page in PdfReader(path).pages)


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


def audit(project_root: Path, docx: Path, pdf: Path) -> dict:
    out = project_root / "results" / "submission_v42"
    figures = out / "figures"
    tables = out / "supplementary_tables"
    manuscript_path = project_root / "manuscript" / "OC_OA_manuscript_revision_v42.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    main_text = manuscript.split("\n## References\n", 1)[0]
    discussion = main_text.split("\n## Discussion\n", 1)[1].split("\n## Conclusions\n", 1)[0]
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    expected_figures = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing = [
        f"{stem}.{ext}"
        for stem in expected_figures
        for ext in ("png", "pdf")
        if not (figures / f"{stem}.{ext}").exists()
    ]
    check("figure package", not missing, "12 PNG and 12 PDF figures" if not missing else str(missing))
    check(
        "WGCNA figure removed",
        not list(figures.glob("SupplementaryFigure1_core_sensitivity.*")),
        "old WGCNA-containing composite absent",
    )

    production_missing = [
        f"{stem}.{ext}"
        for stem in MAIN_FIGURES
        for ext in ("svg", "tiff")
        if not (figures / f"{stem}.{ext}").exists()
    ]
    check("main figure exports", not production_missing, "seven main figures include SVG and TIFF" if not production_missing else str(production_missing))

    legend_text = (figures / "figure_legends.md").read_text(encoding="utf-8")
    legend_count = len(re.findall(r"^## (?:Supplementary )?Figure \d+\.", legend_text, flags=re.MULTILINE))
    check("figure legends", legend_count == 12, f"legend headings={legend_count}")
    check("single-cell boundary", "do not support functional inference" in legend_text and "do not identify the cell population responsible for the blood association" in legend_text, "localization and blood-source limits are explicit")

    figure7 = pdf_text(figures / "Figure7_integrated_interpretation.pdf")
    sequence = ["Tissue overlap", "Directional heterogeneity", "External replication", "Source-defined cell localization", "Blood persistence", "G0S2 association"]
    positions = [figure7.find(item) for item in sequence]
    check("Figure 7 evidence sequence", all(value >= 0 for value in positions) and positions == sorted(positions), "linear evidence sequence present")
    figure7_flat = re.sub(r"\s+", " ", figure7)
    check("Figure 7 G0S2 boundary", "source and function unresolved" in figure7_flat, "G0S2 endpoint remains bounded")

    actual_tables = sorted(path.name for path in tables.glob("*.csv"))
    check("supplementary tables", len(actual_tables) == 19, f"CSV tables={len(actual_tables)}")
    check("WGCNA table removed", not (tables / "Table_S15_WGCNA_stability.csv").exists(), "WGCNA table absent")
    string_tables = [tables / f"Table_S15{suffix}" for suffix in ("a_STRING_mapping_audit.csv", "b_STRING_edges.csv", "c_STRING_node_topology.csv")]
    check("STRING tables renumbered", all(path.exists() for path in string_tables), "STRING tables are S15a-c")

    shared = rows(tables / "Table_S2_shared_tissue_DEGs.csv")
    concordant = sum(float(row["logFC_OA"]) * float(row["logFC_OC"]) > 0 for row in shared)
    check("shared gene counts", len(shared) == 286 and concordant == 146, f"shared={len(shared)} concordant={concordant} discordant={len(shared)-concordant}")
    blood = rows(tables / "Table_S10_FDR_supported_blood_component.csv")
    check("blood result", len(blood) == 1 and blood[0].get("gene") == "G0S2", f"rows={len(blood)} genes={[row.get('gene') for row in blood]}")

    reference_count = len(re.findall(r"^\d+\.\s", manuscript, flags=re.MULTILINE))
    reference_audit = json.loads((out / "reference_audit" / "reference_audit_v42.json").read_text(encoding="utf-8"))
    status_ok = all(record["status"] in {"crossref_verified", "manual_verified"} for record in reference_audit)
    citations = cited_numbers(main_text)
    check("references", reference_count == 69 and len(reference_audit) == 69 and status_ok and citations == set(range(1, 70)), f"manuscript={reference_count} audit={len(reference_audit)} cited={len(citations)} verified={status_ok}")

    lower = main_text.lower()
    check("batch transparency", "no explicit batch-correction method was applied" in lower and "not merged or adjusted with combat" in lower and "disease group as the only design factor" in lower, "batch strategy and model design are explicit")
    check("WGCNA scope", "wgcna" not in lower and "co-expression stability" not in lower, "WGCNA absent from main manuscript")
    prohibited = ["systemic signature", "systemic biomarker", "diagnostic signature", "biomarker panel", "therapeutic target", "shared disease mechanism", "different biological function"]
    found = [term for term in prohibited if term in lower]
    check("terminology", not found, "bounded terminology present" if not found else str(found))
    check("G0S2 boundary", "blood-persistent transcriptomic association" in lower and "did not establish a shared cellular source" in lower, "G0S2 is statistical and source-bounded")
    discussion_words = len(re.findall(r"\b[\w-]+\b", discussion))
    check("Discussion compression", discussion_words <= 390, f"Discussion words={discussion_words}")
    check("formulaic transitions", not any(term in lower for term in ("interestingly", "notably", "importantly")), "formulaic emphasis transitions absent")
    check("ASCII punctuation", not re.search(r"[\u2010-\u2015]", manuscript), "no non-ASCII dash characters")

    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
    check("embedded figures", len(media) == 12, f"embedded media={len(media)}")
    pages = len(PdfReader(pdf).pages)
    check("PDF render", pages >= 24 and pdf.stat().st_size > 1_000_000, f"pages={pages} bytes={pdf.stat().st_size}")
    return {"passed": all(item["passed"] for item in checks), "checks": checks, "pdf_pages": pages}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    args = parser.parse_args()
    report = audit(args.project_root.resolve(), args.docx.resolve(), args.pdf.resolve())
    output_dir = args.project_root.resolve() / "results" / "submission_v42"
    (output_dir / "submission_audit_v42.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# V4.2 submission audit", "", f"Overall: {'PASS' if report['passed'] else 'FAIL'}", ""]
    lines.extend(f"- {'PASS' if item['passed'] else 'FAIL'} - {item['name']}: {item['detail']}" for item in report["checks"])
    (output_dir / "submission_audit_v42.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"V4.2 audit: {'PASS' if report['passed'] else 'FAIL'}")
    for item in report["checks"]:
        print(f"{'PASS' if item['passed'] else 'FAIL'} {item['name']}: {item['detail']}")
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
