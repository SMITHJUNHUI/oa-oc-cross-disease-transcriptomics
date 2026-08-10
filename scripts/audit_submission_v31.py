#!/usr/bin/env python3
"""Audit the V3.1 scientific claims, deliverables, dependencies, and credentials."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from zipfile import ZipFile

from pypdf import PdfReader


JWT = re.compile(rb"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}")
FIXED_PANEL = [
    "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1",
    "CFI", "DDIT3", "DIRAS3", "EFEMP1", "HK2",
]
RANKED_TOP15 = FIXED_PANEL + ["KIT", "MYZAP", "NOD2", "OGN", "RTN1"]
MAIN_FIGURES = [
    "Figure1_question_driven_framework.png",
    "Figure2_gene_and_direction_heterogeneity.png",
    "Figure3_pathway_direction.png",
    "Figure4_cellular_context.png",
    "Figure5_molecular_separability.png",
    "Figure6_integrated_context_model.png",
]
SUPPLEMENTARY_FIGURES = [
    "SupplementaryFigure1_sensitivity_details.png",
    "SupplementaryFigure2_negative_bidirectional_MR.png",
    "SupplementaryFigure3_single_cell_UMAPs.png",
    "SupplementaryFigure4_HPA_normal_tissue_context.png",
    "SupplementaryFigure5_TCGA_relative_context.png",
    "SupplementaryFigure6_cell_type_functional_annotation.png",
    "SupplementaryFigure7_complete_pathway_direction.png",
    "SupplementaryFigure8_candidate_evidence_stability.png",
    "SupplementaryFigure9_candidate_centered_Hallmark_context.png",
    "SupplementaryFigure10_design_and_reliability.png",
    "SupplementaryFigure11_upstream_regulatory_context.png",
    "SupplementaryFigure12_direction_aware_STRING_network.png",
    "SupplementaryFigure13_sample_consensus_CellChat.png",
    "SupplementaryFigure14_NicheNet_prior_overlay.png",
    "SupplementaryFigure15_panel_size_sensitivity.png",
    "SupplementaryFigure16_bulk_PCA_and_QC.png",
]
REQUIRED_TABLES = [
    "Table_S16_candidate_prioritization_matrix.csv",
    "Table_S25d_STRING_network_topology.csv",
    "Table_S26a_CellChat_sample_audit.csv",
    "Table_S26b_CellChat_sample_interactions.csv",
    "Table_S26c_CellChat_consensus_pathways.csv",
    "Table_S26d_shared_DEG_anchored_CellChat.csv",
    "Table_S27_NicheNet_prior_overlay.csv",
    "Table_S28_communication_feasibility_boundaries.csv",
    "Table_S29a_panel_size_composition.csv",
    "Table_S29b_panel_size_direction_sensitivity.csv",
    "Table_S29c_panel_size_Hallmark_sensitivity.csv",
    "Table_S29d_panel_size_pathway_profile_comparison.csv",
    "Table_S30a_extended_gene_cell_detection.csv",
    "Table_S30b_panel_size_localization_sensitivity.csv",
    "Table_S30c_panel_size_localization_profile_comparison.csv",
    "Table_S31a_bulk_unsupervised_PCA.csv",
    "Table_S31b_bulk_sample_correlation_QC.csv",
]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def close(value: str, expected: float, tolerance: float = 1e-9) -> bool:
    return math.isclose(float(value), expected, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    outputs = args.outputs.resolve()
    submission = project / "results" / "submission_v31"
    figures = submission / "figures"
    tables = submission / "supplementary_tables"
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(condition), "detail": detail})

    expected_figures = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing_figures = [name for name in expected_figures if not (figures / name).is_file()]
    check("22 publication figures", not missing_figures, f"missing={missing_figures}")
    missing_tables = [name for name in REQUIRED_TABLES if not (tables / name).is_file()]
    check("V3.1 supplementary tables", not missing_tables, f"missing={missing_tables}")

    legend_text = (figures / "figure_legends.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", legend_text, flags=re.MULTILINE)
    check(
        "figure legends unique",
        len(headings) == 22 and len(set(headings)) == 22,
        f"headings={len(headings)}, unique={len(set(headings))}",
    )

    candidate_rows = rows(tables / "Table_S16_candidate_prioritization_matrix.csv")
    candidate_order = [row["gene"] for row in sorted(candidate_rows, key=lambda row: int(row["prioritization_rank"]))]
    check("fixed ten-gene set and order", candidate_order == FIXED_PANEL, f"genes={candidate_order}")

    composition = rows(tables / "Table_S29a_panel_size_composition.csv")
    ranked = [row["gene"] for row in sorted(composition, key=lambda row: int(row["original_rank"]))]
    check("deterministic top-15 extension", ranked == RANKED_TOP15, f"genes={ranked}")

    direction = rows(tables / "Table_S29b_panel_size_direction_sensitivity.csv")
    discordant = {
        row["panel"]: (int(row["genes"]), float(row["proportion"]))
        for row in direction if row["direction_class"] == "discordant"
    }
    expected_direction = {
        "Top 5": (3, 0.6),
        "Top 10": (7, 0.7),
        "Top 15": (10, 2 / 3),
    }
    direction_ok = discordant.keys() == expected_direction.keys() and all(
        discordant[panel][0] == expected[0]
        and math.isclose(discordant[panel][1], expected[1], rel_tol=1e-12, abs_tol=1e-12)
        for panel, expected in expected_direction.items()
    )
    check("panel-size direction sensitivity", direction_ok, f"discordant={discordant}")

    hallmark = rows(tables / "Table_S29c_panel_size_Hallmark_sensitivity.csv")
    hallmark_fdr = [float(row["FDR"]) for row in hallmark]
    hallmark_ok = len(hallmark) == 150 and min(hallmark_fdr) >= 0.05
    check("small-panel Hallmark boundary", hallmark_ok, f"rows={len(hallmark)}, min_FDR={min(hallmark_fdr):.6g}")

    pathway_profile = rows(tables / "Table_S29d_panel_size_pathway_profile_comparison.csv")
    expected_pathway = {
        ("Top 5", "Top 10"): (0.628715584010892, 0.428571428571429),
        ("Top 5", "Top 15"): (0.540367424058794, 0.428571428571429),
        ("Top 10", "Top 15"): (0.883114821537604, 0.666666666666667),
    }
    observed_pathway = {
        (row["panel_left"], row["panel_right"]):
        (float(row["spearman_fold_enrichment"]), float(row["top10_jaccard"]))
        for row in pathway_profile
    }
    pathway_ok = observed_pathway.keys() == expected_pathway.keys() and all(
        math.isclose(observed_pathway[key][0], expected[0], rel_tol=1e-9, abs_tol=1e-9)
        and math.isclose(observed_pathway[key][1], expected[1], rel_tol=1e-9, abs_tol=1e-9)
        for key, expected in expected_pathway.items()
    )
    check("panel pathway-profile sensitivity", pathway_ok, f"comparisons={observed_pathway}")

    localization = rows(tables / "Table_S30c_panel_size_localization_profile_comparison.csv")
    localization_lookup = {
        (row["dataset_id"], row["panel_left"], row["panel_right"]): row
        for row in localization
    }
    localization_ok = (
        len(localization) == 6
        and localization_lookup[("GSE154600", "Top 5", "Top 10")]["top_cell_type_left"] == "Ovarian.cancer.cell"
        and localization_lookup[("GSE154600", "Top 5", "Top 10")]["top_cell_type_right"] == "Fibroblast"
        and localization_lookup[("GSE154600", "Top 10", "Top 15")]["top_cell_type_left"] == "Fibroblast"
        and localization_lookup[("GSE154600", "Top 10", "Top 15")]["top_cell_type_right"] == "Fibroblast"
        and all(
            row["top_cell_type_left"] == "HomC" and row["top_cell_type_right"] == "HomC"
            for row in localization if row["dataset_id"] == "GSE255460"
        )
    )
    check("exact-label localization sensitivity", localization_ok, f"rows={len(localization)}")

    pca = rows(tables / "Table_S31a_bulk_unsupervised_PCA.csv")
    correlations = rows(tables / "Table_S31b_bulk_sample_correlation_QC.csv")
    pca_counts = {dataset: sum(row["dataset"] == dataset for row in pca) for dataset in {row["dataset"] for row in pca}}
    correlation_counts = {
        dataset: sum(row["dataset"] == dataset for row in correlations)
        for dataset in {row["dataset"] for row in correlations}
    }
    pca_ok = pca_counts == {"OA discovery": 38, "OC discovery": 63} and correlation_counts == pca_counts
    check("bulk PCA and sample-correlation coverage", pca_ok, f"PCA={pca_counts}, correlation={correlation_counts}")

    topology = rows(tables / "Table_S25d_STRING_network_topology.csv")
    physical = [
        row for row in topology
        if row["network_type"] == "high-confidence physical" and row["subset"] == "all_mapped_shared_DEGs"
    ]
    topology_ok = (
        len(physical) == 1 and physical[0]["mapped_nodes"] == "275"
        and physical[0]["connected_nodes"] == "46" and physical[0]["edges"] == "62"
    )
    check("STRING primary topology", topology_ok, "expected 275 mapped, 46 connected, 62 edges")

    cellchat_audit = rows(tables / "Table_S26a_CellChat_sample_audit.csv")
    completed = [row for row in cellchat_audit if row["status"] == "completed"]
    interactions = rows(tables / "Table_S26b_CellChat_sample_interactions.csv")
    consensus = rows(tables / "Table_S26c_CellChat_consensus_pathways.csv")
    retained_consensus = [row for row in consensus if row.get("consensus_status", "").upper() == "TRUE"]
    anchored = rows(tables / "Table_S26d_shared_DEG_anchored_CellChat.csv")
    direct_candidates = [row for row in anchored if row.get("candidate_overlap", "").strip()]
    check("CellChat biological samples", len(completed) == 16, f"completed={len(completed)}")
    check("CellChat sample-level interactions", len(interactions) == 49647, f"rows={len(interactions)}")
    check("CellChat consensus records", len(retained_consensus) == 199, f"retained={len(retained_consensus)}")
    check("shared-DEG-anchored consensus", len(anchored) == 526, f"rows={len(anchored)}")
    check("no direct fixed-panel ligand/receptor", not direct_candidates, f"rows={len(direct_candidates)}")

    nichenet = rows(tables / "Table_S27_NicheNet_prior_overlay.csv")
    niche_targets = {row["target"] for row in nichenet}
    check(
        "NicheNet bounded prior overlay",
        len(nichenet) == 10240 and niche_targets == set(FIXED_PANEL),
        f"rows={len(nichenet)}, targets={sorted(niche_targets)}",
    )

    manuscript = project / "manuscript" / "OC_OA_manuscript_revision_v31.md"
    manuscript_text = manuscript.read_text(encoding="utf-8")
    required_phrases = [
        "146 (51.0%) were concordant and 140 (49.0%) discordant",
        "no panel-level Hallmark term survived FDR correction",
        "the highest label shifted from Ovarian.cancer.cell for top 5 to Fibroblast for top 10 and top 15",
        "none contained a fixed panel gene as a direct ligand or receptor",
        "ligand activity was not estimable under the available design",
        "not a universal diagnostic model",
        "Bidirectional MR found no evidence that genetic liability to OA affected OC risk, or that genetic liability to OC affected OA risk, under the selected datasets, instruments, and assumptions.",
    ]
    missing_phrases = [phrase for phrase in required_phrases if phrase not in manuscript_text]
    check("manuscript pivotal claims and boundaries", not missing_phrases, f"missing={missing_phrases}")
    forbidden_relabels = re.findall(r"\bCAF\b|tumou?r epithelial", manuscript_text, flags=re.IGNORECASE)
    check("no strengthened single-cell relabeling", not forbidden_relabels, f"hits={forbidden_relabels}")

    registry = submission / "claim_evidence_registry_v31.csv"
    checklist = submission / "reproducibility_checklist_v31.csv"
    check("claim-evidence registry", registry.is_file() and len(rows(registry)) >= 10, f"rows={len(rows(registry)) if registry.exists() else 0}")
    check("reproducibility checklist", checklist.is_file() and len(rows(checklist)) >= 5, f"rows={len(rows(checklist)) if checklist.exists() else 0}")

    cellchat_archive = project / "environment" / "vendor" / "CellChat_2.2.0.9001.tar.gz"
    nichenet_prior = project / "results" / "api_cache" / "nichenet_v2" / "ligand_target_matrix_nsga2r_final.rds"
    expected_cellchat = "DA6CEA9B0F8AD59A44F71733AD48148297EB10768C3A9F77918B458EC0DD8667"
    expected_nichenet = "699FCE17FF65E2511277359696306BB130AA833BD7F30DD9290AD9EFD9DC9C5D"
    check(
        "vendored CellChat integrity",
        cellchat_archive.is_file() and sha256(cellchat_archive) == expected_cellchat,
        expected_cellchat,
    )
    check(
        "NicheNet v2 prior integrity",
        nichenet_prior.is_file() and sha256(nichenet_prior) == expected_nichenet,
        expected_nichenet,
    )

    docx = outputs / "OC_OA_manuscript_revision_v31_with_figures.docx"
    pdf = outputs / "OC_OA_manuscript_revision_v31_with_figures.pdf"
    media_count = 0
    if docx.is_file():
        with ZipFile(docx) as archive:
            media_count = sum(name.startswith("word/media/") for name in archive.namelist())
    check("Word manuscript embeds 22 figures", media_count == 22, f"media={media_count}")
    page_count = len(PdfReader(pdf).pages) if pdf.is_file() else 0
    check("PDF export has 43 pages", page_count == 43, f"pages={page_count}")

    text_suffixes = {
        ".r", ".py", ".ps1", ".bat", ".md", ".txt", ".csv", ".tsv",
        ".json", ".yaml", ".yml", ".lock", ".log", ".env",
    }
    secret_hits: list[str] = []
    for path in project.rglob("*"):
        if (
            not path.is_file()
            or (path.suffix.lower() not in text_suffixes and path.name not in {".Renviron", ".Rprofile"})
            or path.stat().st_size > 12_000_000
        ):
            continue
        if JWT.search(path.read_bytes()):
            secret_hits.append(str(path.relative_to(project)))
    check("credential scan", not secret_hits, f"hits={secret_hits}")

    passed = all(item["passed"] for item in checks)
    payload = {"status": "PASS" if passed else "FAIL", "checks": checks}
    args.json.resolve().parent.mkdir(parents=True, exist_ok=True)
    args.json.resolve().write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# OC-OA V3.1 submission audit",
        "",
        f"**Overall status: {'PASS' if passed else 'FAIL'}**",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    for item in checks:
        detail = str(item["detail"]).replace("|", "\\|")
        lines.append(f"| {item['check']} | {'PASS' if item['passed'] else 'FAIL'} | {detail} |")
    lines.extend([
        "",
        "This audit verifies the V3.1 figure/table set, exact pivotal counts, fixed-panel identity, panel-size sensitivity, exact single-cell labels, bulk-QC coverage, dependency integrity, manuscript inference boundaries, embedded media, PDF pagination, and absence of JWT-like credentials.",
        "",
    ])
    args.report.resolve().write_text("\n".join(lines), encoding="utf-8")
    print(f"Audit: {payload['status']}")
    print(f"Checks: {len(checks)}")
    print(f"Report: {args.report.resolve()}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
