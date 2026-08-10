from __future__ import annotations

import argparse
import csv
import json
import re
import zipfile
from pathlib import Path


FIXED_CANDIDATES = {
    "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1",
    "CFI", "DDIT3", "DIRAS3", "EFEMP1", "HK2",
}

MAIN_FIGURES = [
    "Figure1_study_workflow",
    "Figure2_shared_transcriptomic_alterations",
    "Figure3_functional_characterization",
    "Figure4_transcriptional_divergence",
    "Figure5_candidate_gene_characterization",
    "Figure6_single_cell_localization",
    "Figure7_summary_model",
]

SUPPLEMENTARY_FIGURES = [
    "SupplementaryFigure1_core_sensitivity",
    "SupplementaryFigure2_candidate_set_sensitivity",
    "SupplementaryFigure3_exploratory_classification",
    "SupplementaryFigure4_complete_Hallmark_direction",
    "SupplementaryFigure5_all_single_cell_UMAPs",
    "SupplementaryFigure6_bulk_PCA_QC",
    "SupplementaryFigure7_candidate_stability",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def audit(project_root: Path, docx: Path, pdf: Path | None) -> dict:
    out = project_root / "results" / "submission_v32"
    figures = out / "figures"
    tables = out / "supplementary_tables"
    manuscript_path = project_root / "manuscript" / "OC_OA_manuscript_revision_v32.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    expected = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing = [f"{stem}.{ext}" for stem in expected for ext in ("png", "pdf") if not (figures / f"{stem}.{ext}").exists()]
    check("figure files", not missing, "28 PNG/PDF files present" if not missing else f"missing: {missing}")

    legends_text = (figures / "figure_legends.md").read_text(encoding="utf-8")
    legend_headings = re.findall(r"^## (?:Supplementary )?Figure \d+\.", legends_text, flags=re.MULTILINE)
    check("figure legends", len(legend_headings) == 14, f"{len(legend_headings)} unique legend headings")

    shared = read_csv(tables / "Table_S2_shared_differentially_expressed_genes.csv")
    concordant = sum((float(row["logFC_OA"]) * float(row["logFC_OC"])) > 0 for row in shared)
    quadrants = {
        "OA+/OC+": sum(float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) > 0 for r in shared),
        "OA-/OC-": sum(float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) < 0 for r in shared),
        "OA+/OC-": sum(float(r["logFC_OA"]) > 0 and float(r["logFC_OC"]) < 0 for r in shared),
        "OA-/OC+": sum(float(r["logFC_OA"]) < 0 and float(r["logFC_OC"]) > 0 for r in shared),
    }
    check("shared-gene counts", len(shared) == 286 and concordant == 146, f"n={len(shared)}; concordant={concordant}; discordant={len(shared)-concordant}")
    check("direction quadrants", quadrants == {"OA+/OC+": 112, "OA-/OC-": 34, "OA+/OC-": 86, "OA-/OC+": 54}, json.dumps(quadrants))

    candidate_rows = read_csv(tables / "Table_S16_candidate_prioritization_matrix.csv")
    candidates = {row["gene"] for row in candidate_rows}
    check("fixed candidate set", candidates == FIXED_CANDIDATES and len(candidate_rows) == 10, ", ".join(sorted(candidates)))

    hallmark = read_csv(tables / "Table_S18_Hallmark_pathway_direction_matrix.csv")
    joint = [r for r in hallmark if truth(r["both_significant"])]
    discordant_pathways = [r for r in joint if r["direction_class"].lower() == "discordant"]
    check("Hallmark direction", len(joint) == 10 and len(discordant_pathways) == 6, f"joint={len(joint)}; opposite={len(discordant_pathways)}")

    go = read_csv(tables / "Table_S11b_GO_shared_genes.csv")
    significant_go = [r for r in go if float(r["p.adjust"]) < 0.05]
    aging = [r for r in significant_go if "aging" in r["Description"].lower()]
    check("GO aging boundary", not aging, "no FDR-significant GO term containing 'aging'")

    kegg = read_csv(tables / "Table_S11c_KEGG_shared_genes.csv")
    significant_kegg = [r for r in kegg if float(r["p.adjust"]) < 0.05]
    check("KEGG result", len(significant_kegg) == 1 and significant_kegg[0]["Description"].lower() == "cell cycle", f"significant={len(significant_kegg)}")

    mapping = read_csv(tables / "Table_S25a_STRING_mapping_audit.csv")
    topology = read_csv(tables / "Table_S25c_STRING_node_topology.csv")
    physical = [r for r in topology if r["network_type"] == "high-confidence physical"]
    connected = [r for r in physical if float(r["degree"]) > 0]
    connected_candidates = sorted(r["gene"] for r in connected if truth(r["candidate"]))
    mapped = sum(truth(r.get("mapped", "")) for r in mapping)
    check("STRING mapping", mapped == 275 and len(mapping) == 286, f"mapped={mapped}/286")
    check("STRING candidate boundary", len(connected) == 46 and connected_candidates == ["JUNB"], f"connected products={len(connected)}; connected candidates={connected_candidates}")

    single_source = read_csv(figures / "source_data" / "Figure6_representative_UMAPs.csv")
    representative_ids = {r["dataset_id"] for r in single_source}
    labels = {r["label"] for r in single_source}
    check("representative single-cell atlases", representative_ids == {"GSE255460", "GSE154600"}, ", ".join(sorted(representative_ids)))
    check("exact OC labels", {"Fibroblast", "Ovarian.cancer.cell"}.issubset(labels), "exact source labels retained")

    scope = read_csv(out / "analysis" / "V32_scope_decisions.csv")
    scope_map = {r["module"]: r["v32_status"] for r in scope}
    expected_excluded = {"Mendelian randomization", "CellChat/NicheNet", "TF-miRNA", "TCGA/immune/HPA context"}
    scope_ok = all("excluded" in scope_map.get(module, "") for module in expected_excluded)
    check("strategic scope", scope_ok and scope_map.get("Exploratory classification") == "supplement only", "excluded modules documented; classification supplement only")

    main_names = " ".join(MAIN_FIGURES).lower()
    check("classification demotion", "classification" not in main_names and "roc" not in main_names, "no classification or ROC main-figure filename")

    title = manuscript.splitlines()[0]
    check("title positioning", "transcriptional divergence" in title.lower() and "conserved" not in title.lower(), title.lstrip("# "))
    required_phrases = [
        "286 shared genes", "146 (51.0%)", "140 (49.0%)", "six had opposite",
        "Exploratory classification analysis showed cohort-dependent discrimination performance",
        "Fibroblast was not relabeled", "Ovarian.cancer.cell was not relabeled",
    ]
    missing_phrases = [phrase for phrase in required_phrases if phrase not in manuscript]
    check("manuscript claims", not missing_phrases, "required bounded claims present" if not missing_phrases else f"missing: {missing_phrases}")

    named_excluded = re.findall(r"\b(?:CellChat|NicheNet|TF-miRNA|Mendelian randomization|\bMR\b)\b", manuscript, flags=re.IGNORECASE)
    check("excluded-module narrative", not named_excluded, "excluded modules absent from submitted manuscript" if not named_excluded else f"found: {sorted(set(named_excluded))}")

    with zipfile.ZipFile(docx) as archive:
        media = [name for name in archive.namelist() if name.startswith("word/media/")]
    check("DOCX embedded figures", len(media) == 14, f"embedded media={len(media)}")

    if pdf and pdf.exists():
        try:
            import fitz  # type: ignore
            pages = fitz.open(pdf).page_count
            check("PDF rendering", pages > 20, f"pages={pages}")
        except Exception as exc:
            check("PDF rendering", False, f"could not inspect PDF: {exc}")

    jwt_pattern = re.compile(rb"eyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}")
    secret_hits: list[str] = []
    scan_roots = [project_root / "manuscript", out, project_root / "scripts", project_root / "R"]
    for root in scan_roots:
        for path in root.rglob("*"):
            if path.is_file() and path.stat().st_size <= 10_000_000:
                try:
                    if jwt_pattern.search(path.read_bytes()):
                        secret_hits.append(str(path.relative_to(project_root)))
                except OSError:
                    pass
    check("secret scan", not secret_hits, "no JWT-like token found" if not secret_hits else f"hits: {secret_hits}")
    check("one-command entry point", (project_root / "run_submission_v32.ps1").exists(), "run_submission_v32.ps1 present")

    passed = all(item["passed"] for item in checks)
    return {"status": "PASS" if passed else "FAIL", "checks": checks}


def write_report(result: dict, report_path: Path, json_path: Path) -> None:
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# V3.2 submission audit", "", f"**Overall status: {result['status']}**", "", "| Check | Status | Detail |", "|---|---:|---|"]
    for item in result["checks"]:
        status = "PASS" if item["passed"] else "FAIL"
        detail = str(item["detail"]).replace("|", "\\|")
        lines.append(f"| {item['name']} | {status} | {detail} |")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result = audit(root, args.docx.resolve(), args.pdf.resolve() if args.pdf else None)
    out = root / "results" / "submission_v32"
    write_report(result, out / "submission_audit_v32.md", out / "submission_audit_v32.json")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
