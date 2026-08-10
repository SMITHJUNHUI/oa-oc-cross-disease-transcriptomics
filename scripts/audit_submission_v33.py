from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import zipfile
from pathlib import Path


MAIN_FIGURES = [
    "Figure1_study_workflow_with_blood",
    "Figure2_shared_transcriptomic_alterations",
    "Figure3_functional_characterization",
    "Figure4_transcriptional_divergence",
    "Figure5_peripheral_blood_systemic_component",
    "Figure6_single_cell_localization",
    "Figure7_integrated_summary_model",
]

SUPPLEMENTARY_FIGURES = [
    "SupplementaryFigure1_core_sensitivity",
    "SupplementaryFigure2_complete_Hallmark_direction",
    "SupplementaryFigure3_all_single_cell_UMAPs",
    "SupplementaryFigure4_bulk_PCA_QC",
    "SupplementaryFigure5_direction_aware_STRING_network",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def pdf_pages(pdfinfo: Path, pdf: Path) -> int:
    process = subprocess.run(
        [str(pdfinfo), str(pdf)], capture_output=True, text=True, check=True
    )
    match = re.search(r"^Pages:\s+(\d+)", process.stdout, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("pdfinfo did not return a Pages field")
    return int(match.group(1))


def audit(project_root: Path, docx: Path, pdf: Path, pdfinfo: Path) -> dict:
    out = project_root / "results" / "submission_v33"
    figures = out / "figures"
    tables = out / "supplementary_tables"
    manuscript_path = project_root / "manuscript" / "OC_OA_manuscript_revision_v33.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    expected = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing = [
        f"{stem}.{ext}"
        for stem in expected
        for ext in ("png", "pdf")
        if not (figures / f"{stem}.{ext}").exists()
    ]
    check(
        "figure files",
        not missing,
        "12 PNG and 12 PDF figures present" if not missing else f"missing: {missing}",
    )

    legends = (figures / "figure_legends.md").read_text(encoding="utf-8")
    headings = re.findall(
        r"^## (?:Supplementary )?Figure \d+\.", legends, flags=re.MULTILINE
    )
    check("figure legends", len(headings) == 12, f"legend headings={len(headings)}")

    shared = read_csv(tables / "Table_S2_shared_differentially_expressed_genes.csv")
    concordant = sum(
        float(row["logFC_OA"]) * float(row["logFC_OC"]) > 0 for row in shared
    )
    quadrants = {
        "OA+/OC+": sum(
            float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) > 0 for r in shared
        ),
        "OA-/OC-": sum(
            float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) < 0 for r in shared
        ),
        "OA+/OC-": sum(
            float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) < 0 for r in shared
        ),
        "OA-/OC+": sum(
            float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) > 0 for r in shared
        ),
    }
    check(
        "shared-gene counts",
        len(shared) == 286 and concordant == 146,
        f"n={len(shared)}; concordant={concordant}; discordant={len(shared)-concordant}",
    )
    check(
        "direction quadrants",
        quadrants
        == {"OA+/OC+": 112, "OA-/OC-": 34, "OA+/OC-": 86, "OA-/OC+": 54},
        json.dumps(quadrants),
    )

    hallmark = read_csv(tables / "Table_S18_Hallmark_pathway_direction_matrix.csv")
    joint = [r for r in hallmark if truth(r["both_significant"])]
    opposite = [r for r in joint if r["direction_class"].lower() == "discordant"]
    check(
        "Hallmark direction",
        len(joint) == 10 and len(opposite) == 6,
        f"joint={len(joint)}; opposite={len(opposite)}",
    )

    blood_audit = read_csv(tables / "Table_S13_blood_dataset_audit.csv")
    blood = {r["dataset_id"]: r for r in blood_audit}
    blood_ok = (
        set(blood) == {"GSE48556", "GSE31682"}
        and blood["GSE48556"]["biospecimen"] == "PBMC"
        and int(blood["GSE48556"]["disease_n"]) == 106
        and int(blood["GSE48556"]["control_n"]) == 33
        and blood["GSE31682"]["biospecimen"] == "blood cell fraction"
        and int(blood["GSE31682"]["disease_n"]) == 48
        and int(blood["GSE31682"]["control_n"]) == 20
    )
    check("blood cohort audit", blood_ok, json.dumps(blood, ensure_ascii=False))

    attrition_rows = read_csv(tables / "Table_S15_blood_screen_attrition.csv")
    attrition = {r["stage"]: int(r["genes"]) for r in attrition_rows}
    expected_attrition = {
        "Shared tissue DEGs": 286,
        "Tissue-concordant DEGs": 146,
        "Measured in both blood cohorts": 127,
        "All four effects same direction": 38,
        "Nominal P<0.05 in both blood cohorts": 3,
        "FDR<0.05 in both blood cohorts": 1,
    }
    check("blood attrition", attrition == expected_attrition, json.dumps(attrition))

    positive = read_csv(tables / "Table_S14_FDR_supported_systemic_component.csv")
    positive_ok = (
        len(positive) == 1
        and positive[0]["gene"] == "G0S2"
        and truth(positive[0]["both_blood_fdr"])
        and all(
            float(positive[0][column]) < 0
            for column in ("logFC_OA", "logFC_OC", "OA_blood_logFC", "OC_blood_logFC")
        )
        and abs(float(positive[0]["OA_blood_hedges_g"]) + 0.7776584) < 1e-5
        and abs(float(positive[0]["OC_blood_hedges_g"]) + 0.7950483) < 1e-5
    )
    check("blood FDR result", positive_ok, f"rows={len(positive)}; genes={[r['gene'] for r in positive]}")

    mapping = read_csv(tables / "Table_S25a_STRING_mapping_audit.csv")
    mapped = sum(bool(r.get("stringId", "").strip()) for r in mapping)
    topology = read_csv(tables / "Table_S25c_STRING_node_topology.csv")
    physical = [r for r in topology if r["network_type"] == "high-confidence physical"]
    connected = [r for r in physical if float(r["degree"]) > 0]
    connected_candidates = sorted(r["gene"] for r in connected if truth(r["candidate"]))
    check("STRING audit", mapped == 275 and len(connected) == 46, f"mapped={mapped}; connected={len(connected)}")
    check("STRING boundary", connected_candidates == ["JUNB"], f"connected representatives={connected_candidates}")

    title = manuscript.splitlines()[0].lower()
    required = [
        "286 shared tissue degs",
        "146 (51.0%)",
        "140 (49.0%)",
        "one with fdr <0.05 in each blood cohort",
        "the sole fdr-supported gene was g0s2",
        "fibroblast was not relabeled",
        "ovarian.cancer.cell was not relabeled",
    ]
    missing_claims = [phrase for phrase in required if phrase not in manuscript.lower()]
    check(
        "manuscript claims",
        not missing_claims,
        "bounded claims present" if not missing_claims else f"missing: {missing_claims}",
    )
    check(
        "title boundary",
        "limited blood-replicated component" in title and "mechanism" not in title,
        manuscript.splitlines()[0].lstrip("# "),
    )

    forbidden = re.findall(
        r"\b(?:AUC|ROC|Mendelian|CellChat|NicheNet|TF-miRNA|nomogram|decision curve)\b",
        manuscript,
        flags=re.IGNORECASE,
    )
    check(
        "excluded submission modules",
        not forbidden,
        "no excluded-module narrative" if not forbidden else f"found={sorted(set(forbidden))}",
    )
    check(
        "nominal-only names absent",
        "KPNA2" not in manuscript and "PRKX" not in manuscript,
        "nominal-only genes not named in submission manuscript",
    )

    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
    check("DOCX embedded figures", len(media) == 12, f"embedded media={len(media)}")

    try:
        pages = pdf_pages(pdfinfo, pdf)
        check("PDF rendering", pages >= 25, f"pages={pages}")
    except Exception as exc:
        check("PDF rendering", False, str(exc))

    jwt_pattern = re.compile(rb"eyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}")
    secret_hits: list[str] = []
    for root in [project_root / "manuscript", out, project_root / "scripts", project_root / "R"]:
        for path in root.rglob("*"):
            if path.is_file() and path.stat().st_size <= 10_000_000:
                try:
                    if jwt_pattern.search(path.read_bytes()):
                        secret_hits.append(str(path.relative_to(project_root)))
                except OSError:
                    pass
    check("secret scan", not secret_hits, "no JWT-like token found" if not secret_hits else f"hits={secret_hits}")
    check(
        "one-command entry point",
        (project_root / "run_submission_v33.ps1").exists(),
        "run_submission_v33.ps1 present",
    )

    passed = all(item["passed"] for item in checks)
    return {"status": "PASS" if passed else "FAIL", "checks": checks}


def write_report(result: dict, report: Path, json_path: Path) -> None:
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# V3.3 submission audit",
        "",
        f"**Overall status: {result['status']}**",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    for item in result["checks"]:
        status = "PASS" if item["passed"] else "FAIL"
        detail = str(item["detail"]).replace("|", "\\|")
        lines.append(f"| {item['name']} | {status} | {detail} |")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--pdfinfo", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = audit(root, args.docx.resolve(), args.pdf.resolve(), args.pdfinfo.resolve())
    out = root / "results" / "submission_v33"
    write_report(result, out / "submission_audit_v33.md", out / "submission_audit_v33.json")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
