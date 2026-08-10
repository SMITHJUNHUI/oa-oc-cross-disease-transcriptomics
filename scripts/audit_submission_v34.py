from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import zipfile
from pathlib import Path


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


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def page_count(pdfinfo: Path, pdf: Path) -> int:
    result = subprocess.run([str(pdfinfo), str(pdf)], capture_output=True, text=True, check=True)
    match = re.search(r"^Pages:\s+(\d+)", result.stdout, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("pdfinfo did not report a page count")
    return int(match.group(1))


def audit(project_root: Path, docx: Path, pdf: Path, pdfinfo: Path) -> dict:
    out = project_root / "results" / "submission_v34"
    figures = out / "figures"
    tables = out / "supplementary_tables"
    manuscript_path = project_root / "manuscript" / "OC_OA_manuscript_revision_v34.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    expected_figures = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing_figures = [
        f"{stem}.{ext}"
        for stem in expected_figures
        for ext in ("png", "pdf")
        if not (figures / f"{stem}.{ext}").exists()
    ]
    check("figure package", not missing_figures, "12 PNG and 12 PDF figures present" if not missing_figures else str(missing_figures))

    new_main_missing = [f"{stem}.{ext}" for stem in MAIN_FIGURES for ext in ("svg", "tiff") if not (figures / f"{stem}.{ext}").exists()]
    check("main-figure production exports", not new_main_missing, "seven main figures include SVG and 600-dpi TIFF" if not new_main_missing else str(new_main_missing))

    legends = (figures / "figure_legends.md").read_text(encoding="utf-8")
    legend_count = len(re.findall(r"^## (?:Supplementary )?Figure \d+\.", legends, flags=re.MULTILINE))
    check("figure legends", legend_count == 12, f"legend headings={legend_count}")

    expected_tables = [
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
    missing_tables = [name for name in expected_tables if not (tables / name).exists()]
    extra_tables = sorted(path.name for path in tables.glob("*.csv") if path.name not in expected_tables)
    check("supplementary table numbering", not missing_tables and not extra_tables, "20 uniquely numbered CSV tables" if not missing_tables and not extra_tables else f"missing={missing_tables}; extra={extra_tables}")

    shared = rows(tables / "Table_S2_shared_tissue_DEGs.csv")
    concordant = sum(float(row["logFC_OA"]) * float(row["logFC_OC"]) > 0 for row in shared)
    quadrants = {
        "up_up": sum(float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) > 0 for r in shared),
        "down_down": sum(float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) < 0 for r in shared),
        "up_down": sum(float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) < 0 for r in shared),
        "down_up": sum(float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) > 0 for r in shared),
    }
    check("shared tissue counts", len(shared) == 286 and concordant == 146 and quadrants == {"up_up": 112, "down_down": 34, "up_down": 86, "down_up": 54}, f"n={len(shared)}; concordant={concordant}; quadrants={quadrants}")

    hallmark = rows(tables / "Table_S14_discovery_Hallmark_direction.csv")
    joint = [r for r in hallmark if truth(r["both_significant"])]
    discordant = [r for r in joint if r["direction_class"].lower() == "discordant"]
    check("discovery Hallmark direction", len(joint) == 10 and len(discordant) == 6, f"joint={len(joint)}; discordant={len(discordant)}")

    external = {r["dataset_id"]: r for r in rows(tables / "Table_S5_external_tissue_direction_summary.csv")}
    external_expected = {"GSE117999": (143, 280), "GSE82107": (178, 286), "GSE54388": (260, 286), "GSE12470": (179, 226)}
    external_ok = set(external) == set(external_expected) and all((int(external[k]["direction_matches"]), int(external[k]["measured_shared_genes"])) == v for k, v in external_expected.items())
    check("external tissue direction", external_ok, json.dumps({k: [external[k]["direction_matches"], external[k]["measured_shared_genes"]] for k in external}, ensure_ascii=False))

    candidate = rows(tables / "Table_S7_candidate_evidence_summary.csv")
    genes = [row["gene"] for row in candidate]
    check("representative evidence set", genes == ["G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3"] and sum(truth(row["dual_blood_FDR"]) for row in candidate) == 1, f"genes={genes}")

    sc_all = rows(tables / "Table_S8a_single_cell_QC_and_status.csv")
    sc_qc = [row for row in sc_all if row["source_table"] == "qc_status"]
    audited = sum(int(float(row["cells"])) for row in sc_qc)
    passed = sum(int(float(row["qc_pass"])) for row in sc_qc)
    sc_detection = rows(tables / "Table_S8b_representative_gene_single_cell_detection.csv")
    check("single-cell audit", len(sc_qc) == 5 and audited == 1187436 and passed == 1025361, f"datasets={len(sc_qc)}; audited={audited}; passed={passed}")
    check("single-cell representative genes", sorted({row["gene"] for row in sc_detection}) == sorted(genes), f"genes={sorted({row['gene'] for row in sc_detection})}")

    attrition = {r["stage"]: int(r["genes"]) for r in rows(tables / "Table_S11_blood_screen_attrition.csv")}
    expected_attrition = {
        "Shared tissue DEGs": 286,
        "Tissue-concordant DEGs": 146,
        "Measured in both blood cohorts": 127,
        "All four effects same direction": 38,
        "Nominal P<0.05 in both blood cohorts": 3,
        "FDR<0.05 in both blood cohorts": 1,
    }
    check("blood attrition", attrition == expected_attrition, json.dumps(attrition))

    positive = rows(tables / "Table_S10_FDR_supported_blood_component.csv")
    positive_ok = len(positive) == 1 and positive[0]["gene"] == "G0S2" and truth(positive[0]["both_blood_fdr"]) and all(float(positive[0][column]) < 0 for column in ("logFC_OA", "logFC_OC", "OA_blood_logFC", "OC_blood_logFC"))
    check("blood FDR result", positive_ok, f"rows={len(positive)}; genes={[r['gene'] for r in positive]}")

    required_phrases = [
        "286 genes shared between contrasts",
        "146 had the same direction and 140 had opposite directions",
        "79.2–90.9%",
        "51.1–62.2%",
        "one gene meeting fdr <0.05 independently in both blood cohorts",
        "descriptive evidence summary rather than a predictive signature",
        "partial molecular convergence",
    ]
    missing_claims = [phrase for phrase in required_phrases if phrase not in manuscript.lower()]
    check("bounded manuscript claims", not missing_claims, "required claims and boundaries present" if not missing_claims else f"missing={missing_claims}")

    title = manuscript.splitlines()[0].lstrip("# ")
    check("title scope", "across tissue, cellular and blood contexts" in title.lower() and "mechanism" not in title.lower(), title)
    check("nominal-only genes not promoted", "KPNA2" not in manuscript and "PRKX" not in manuscript, "nominal-only blood genes absent")

    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
    check("DOCX embedded figures", len(media) == 12, f"embedded media={len(media)}")

    try:
        pages = page_count(pdfinfo, pdf)
        check("PDF rendering", pages >= 24, f"pages={pages}")
    except Exception as exc:
        check("PDF rendering", False, str(exc))

    jwt = re.compile(rb"eyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}")
    secret_hits: list[str] = []
    for root in (project_root / "manuscript", out, project_root / "scripts", project_root / "R"):
        for path in root.rglob("*"):
            if path.is_file() and path.stat().st_size <= 10_000_000:
                try:
                    if jwt.search(path.read_bytes()):
                        secret_hits.append(str(path.relative_to(project_root)))
                except OSError:
                    pass
    check("secret scan", not secret_hits, "no JWT-like token found" if not secret_hits else f"hits={secret_hits}")
    check("one-command entry point", (project_root / "run_submission_v34.ps1").exists(), "run_submission_v34.ps1 present")

    status = "PASS" if all(item["passed"] for item in checks) else "FAIL"
    return {"status": status, "checks": checks}


def write_report(result: dict, markdown: Path, json_path: Path) -> None:
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# V3.4 submission audit", "", f"**Overall status: {result['status']}**", "", "| Check | Status | Detail |", "|---|---:|---|"]
    for item in result["checks"]:
        detail = str(item["detail"]).replace("|", "\\|")
        lines.append(f"| {item['name']} | {'PASS' if item['passed'] else 'FAIL'} | {detail} |")
    markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--pdfinfo", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = audit(root, args.docx.resolve(), args.pdf.resolve(), args.pdfinfo.resolve())
    out = root / "results" / "submission_v34"
    write_report(result, out / "submission_audit_v34.md", out / "submission_audit_v34.json")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
