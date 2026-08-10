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
    "Figure4_external_tissue_and_representative_genes",
    "Figure5_single_cell_localization",
    "Figure6_peripheral_blood_validation",
    "Figure7_integrated_interpretation",
]

SUPPLEMENTARY_FIGURES = [
    "SupplementaryFigure1_core_sensitivity",
    "SupplementaryFigure2_external_tissue_Hallmark",
    "SupplementaryFigure3_all_single_cell_UMAPs",
    "SupplementaryFigure4_bulk_PCA_QC",
    "SupplementaryFigure5_direction_aware_STRING_network",
]

EXPECTED_TABLES = [
    "Table_S1_tissue_data_sources_and_cohorts.csv",
    "Table_S2_shared_tissue_DEGs.csv",
    "Table_S3a_DEG_threshold_summary.csv",
    "Table_S3b_DEG_threshold_membership.csv",
    "Table_S4_external_tissue_gene_effects.csv",
    "Table_S5_external_tissue_direction_summary.csv",
    "Table_S6_external_tissue_Hallmark_GSEA.csv",
    "Table_S7_candidate_evidence_summary.csv",
    "Table_S8a_single_cell_QC_and_status.csv",
    "Table_S8b_representative_gene_single_cell_detection.csv",
    "Table_S9_blood_dataset_audit.csv",
    "Table_S10_FDR_supported_blood_component.csv",
    "Table_S11_blood_screen_attrition.csv",
    "Table_S12_GO_shared_genes.csv",
    "Table_S13_KEGG_shared_genes.csv",
    "Table_S14_discovery_Hallmark_direction.csv",
    "Table_S15_WGCNA_stability.csv",
    "Table_S16a_STRING_mapping_audit.csv",
    "Table_S16b_STRING_edges.csv",
    "Table_S16c_STRING_node_topology.csv",
]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def pdf_text(path: Path) -> str:
    return "\n".join((page.extract_text() or "") for page in PdfReader(path).pages)


def audit(project_root: Path, docx: Path, pdf: Path) -> dict:
    out = project_root / "results" / "submission_v40"
    figures = out / "figures"
    tables = out / "supplementary_tables"
    manuscript_path = project_root / "manuscript" / "OC_OA_manuscript_revision_v40.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    main_text = manuscript.split("\n## References\n", 1)[0]
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    expected_figures = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing_figures = [
        f"{stem}.{extension}"
        for stem in expected_figures
        for extension in ("png", "pdf")
        if not (figures / f"{stem}.{extension}").exists()
    ]
    check("figure package", not missing_figures, "12 PNG and 12 PDF figures" if not missing_figures else str(missing_figures))

    production_missing = [
        f"{stem}.{extension}"
        for stem in MAIN_FIGURES
        for extension in ("svg", "tiff")
        if not (figures / f"{stem}.{extension}").exists()
    ]
    check("main figure exports", not production_missing, "seven main figures include SVG and TIFF" if not production_missing else str(production_missing))

    legend_text = (figures / "figure_legends.md").read_text(encoding="utf-8")
    legend_count = len(re.findall(r"^## (?:Supplementary )?Figure \d+\.", legend_text, flags=re.MULTILINE))
    check("figure legends", legend_count == 12, f"legend headings={legend_count}")

    figure4_text = pdf_text(figures / "Figure4_external_tissue_and_representative_genes.pdf")
    check(
        "Figure 4 terminology",
        "Candidate molecular features across evidence layers" in figure4_text and "predictive signature" not in figure4_text.lower(),
        "candidate-feature wording present; predictive-signature wording absent",
    )

    actual_tables = sorted(path.name for path in tables.glob("*.csv"))
    check(
        "supplementary tables",
        sorted(EXPECTED_TABLES) == actual_tables,
        "20 uniquely numbered CSV tables" if sorted(EXPECTED_TABLES) == actual_tables else f"actual={actual_tables}",
    )

    shared = rows(tables / "Table_S2_shared_tissue_DEGs.csv")
    concordant = sum(float(row["logFC_OA"]) * float(row["logFC_OC"]) > 0 for row in shared)
    check("shared gene counts", len(shared) == 286 and concordant == 146, f"shared={len(shared)} concordant={concordant} discordant={len(shared)-concordant}")

    external = {row["dataset_id"]: row for row in rows(tables / "Table_S5_external_tissue_direction_summary.csv")}
    expected_external = {"GSE117999": (143, 280), "GSE82107": (178, 286), "GSE54388": (260, 286), "GSE12470": (179, 226)}
    external_ok = set(external) == set(expected_external) and all(
        (int(external[key]["direction_matches"]), int(external[key]["measured_shared_genes"])) == value
        for key, value in expected_external.items()
    )
    check("external tissue replication", external_ok, json.dumps(expected_external, ensure_ascii=False))

    blood = rows(tables / "Table_S10_FDR_supported_blood_component.csv")
    check("blood result", len(blood) == 1 and blood[0].get("gene") == "G0S2", f"rows={len(blood)} genes={[row.get('gene') for row in blood]}")

    reference_count = len(re.findall(r"^\d+\.\s", manuscript, flags=re.MULTILINE))
    reference_audit = json.loads((out / "reference_audit" / "reference_audit_v40.json").read_text(encoding="utf-8"))
    reference_status = all(record["status"] in {"crossref_verified", "manual_verified"} for record in reference_audit)
    check("references", reference_count == 80 and len(reference_audit) == 80 and reference_status, f"manuscript={reference_count} audit={len(reference_audit)} verified={reference_status}")

    prohibited = [
        "predictive signature", "diagnostic signature", "biomarker panel",
        "therapeutic target", "shared disease mechanism", "common pathogenic mechanism",
    ]
    found = [term for term in prohibited if term in main_text.lower()]
    check("terminology", not found, "prohibited main-text terms absent" if not found else str(found))

    formulaic = [term for term in ("Interestingly", "Notably", "Importantly") if term.lower() in main_text.lower()]
    check("formulaic transitions", not formulaic, "formulaic emphasis transitions absent" if not formulaic else str(formulaic))
    check("ASCII punctuation", not re.search(r"[‐‑‒–—−]", manuscript), "no non-ASCII dash characters")

    required_headings = [
        "### Study design",
        "### External tissue replication",
        "### Peripheral blood evaluation of systemic molecular signals",
        "### Directional heterogeneity divided the shared set almost evenly",
        "### Peripheral blood evaluation retained G0S2 as a candidate systemic signal",
    ]
    missing_headings = [heading for heading in required_headings if heading not in manuscript]
    check("section architecture", not missing_headings, "V4.0 evidence sequence present" if not missing_headings else str(missing_headings))

    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
    check("embedded figures", len(media) == 12, f"embedded media={len(media)}")
    pages = len(PdfReader(pdf).pages)
    check("PDF render", pages >= 25 and pdf.stat().st_size > 1_000_000, f"pages={pages} bytes={pdf.stat().st_size}")

    return {"passed": all(item["passed"] for item in checks), "checks": checks, "pdf_pages": pages}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    args = parser.parse_args()

    report = audit(args.project_root.resolve(), args.docx.resolve(), args.pdf.resolve())
    output_dir = args.project_root.resolve() / "results" / "submission_v40"
    json_path = output_dir / "submission_audit_v40.json"
    md_path = output_dir / "submission_audit_v40.md"
    json_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# V4.0 submission audit", "", f"Overall: {'PASS' if report['passed'] else 'FAIL'}", ""]
    for item in report["checks"]:
        lines.append(f"- {'PASS' if item['passed'] else 'FAIL'} - {item['name']}: {item['detail']}")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"V4.0 audit: {'PASS' if report['passed'] else 'FAIL'}")
    for item in report["checks"]:
        print(f"{'PASS' if item['passed'] else 'FAIL'} {item['name']}: {item['detail']}")
    if not report["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
