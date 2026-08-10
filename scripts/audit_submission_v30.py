#!/usr/bin/env python3
"""Audit the V3.0 manuscript, figures, tables, dependencies, and credentials."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from zipfile import ZipFile


JWT = re.compile(rb"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}")
CANDIDATES = {
    "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1",
    "CFI", "DDIT3", "DIRAS3", "EFEMP1", "HK2",
}
MAIN_FIGURES = [
    "Figure1_systems_framework.png",
    "Figure2_gene_direction.png",
    "Figure3_pathway_direction.png",
    "Figure4_cellular_context.png",
    "Figure5_molecular_separability.png",
    "Figure6_genetic_liability_boundary.png",
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
]
REQUIRED_TABLES = [
    "Table_S24a_dataset_context_CCSS.csv",
    "Table_S24b_disease_consensus_CCSS.csv",
    "Table_S24c_sample_aware_UCell_summary.csv",
    "Table_S25a_STRING_mapping_audit.csv",
    "Table_S25b_direction_aware_STRING_edges.csv",
    "Table_S25c_STRING_node_topology.csv",
    "Table_S25d_STRING_network_topology.csv",
    "Table_S25e_STRING_direction_label_permutation.csv",
    "Table_S26a_CellChat_sample_audit.csv",
    "Table_S26b_CellChat_sample_interactions.csv",
    "Table_S26c_CellChat_consensus_pathways.csv",
    "Table_S26d_shared_DEG_anchored_CellChat.csv",
    "Table_S27_NicheNet_prior_overlay.csv",
    "Table_S28_communication_feasibility_boundaries.csv",
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    outputs = args.outputs.resolve()
    submission = project / "results" / "submission_v30"
    figures = submission / "figures"
    tables = submission / "supplementary_tables"
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(condition), "detail": detail})

    expected_figures = MAIN_FIGURES + SUPPLEMENTARY_FIGURES
    missing_figures = [name for name in expected_figures if not (figures / name).is_file()]
    check("20 publication figures", not missing_figures, f"missing={missing_figures}")
    missing_tables = [name for name in REQUIRED_TABLES if not (tables / name).is_file()]
    check("V3.0 supplementary tables", not missing_tables, f"missing={missing_tables}")

    legend_text = (figures / "figure_legends.md").read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", legend_text, flags=re.MULTILINE)
    check("figure legends unique", len(headings) == 20 and len(set(headings)) == 20,
          f"headings={len(headings)}, unique={len(set(headings))}")

    topology = rows(tables / "Table_S25d_STRING_network_topology.csv")
    physical = [row for row in topology if row["network_type"] == "high-confidence physical"
                and row["subset"] == "all_mapped_shared_DEGs"]
    topology_ok = len(physical) == 1 and physical[0]["mapped_nodes"] == "275" \
        and physical[0]["connected_nodes"] == "46" and physical[0]["edges"] == "62"
    check("STRING primary topology", topology_ok, "expected 275 mapped, 46 connected, 62 edges")

    candidate_rows = rows(tables / "Table_S16_candidate_prioritization_matrix.csv")
    genes = {row["gene"] for row in candidate_rows}
    check("fixed ten-gene set", genes == CANDIDATES, f"genes={sorted(genes)}")

    audit = rows(tables / "Table_S26a_CellChat_sample_audit.csv")
    completed = [row for row in audit if row["status"] == "completed"]
    interactions = rows(tables / "Table_S26b_CellChat_sample_interactions.csv")
    consensus = rows(tables / "Table_S26c_CellChat_consensus_pathways.csv")
    consensus_retained = [
        row for row in consensus
        if row.get("consensus_status", "").strip().upper() == "TRUE"
    ]
    anchored = rows(tables / "Table_S26d_shared_DEG_anchored_CellChat.csv")
    check("CellChat biological samples", len(completed) == 16, f"completed={len(completed)}")
    check("CellChat sample-level interactions", len(interactions) == 49647,
          f"rows={len(interactions)}")
    check("CellChat consensus context-pathway records", len(consensus_retained) == 199,
          f"all_rows={len(consensus)}, consensus_rows={len(consensus_retained)}")
    check("shared-DEG-anchored consensus", len(anchored) == 526, f"rows={len(anchored)}")
    direct_candidates = [row for row in anchored if row.get("candidate_overlap", "").strip()]
    check("no direct fixed-candidate ligand/receptor", not direct_candidates,
          f"rows_with_candidate_overlap={len(direct_candidates)}")

    nichenet = rows(tables / "Table_S27_NicheNet_prior_overlay.csv")
    niche_targets = {row["target"] for row in nichenet}
    check("NicheNet bounded prior overlay", len(nichenet) == 10240 and niche_targets == CANDIDATES,
          f"rows={len(nichenet)}, targets={sorted(niche_targets)}")

    manuscript = project / "manuscript" / "OC_OA_manuscript_revision_v30.md"
    manuscript_text = manuscript.read_text(encoding="utf-8")
    required_phrases = [
        "146 (51.0%) were concordant and 140 (49.0%) discordant",
        "none contained a fixed candidate as a direct ligand or receptor",
        "ligand activity was not estimable under the available design",
        "not a universal diagnostic model",
    ]
    absent_phrases = [phrase for phrase in required_phrases if phrase not in manuscript_text]
    check("manuscript inference boundaries", not absent_phrases, f"missing={absent_phrases}")

    cellchat_archive = project / "environment" / "vendor" / "CellChat_2.2.0.9001.tar.gz"
    nichenet_prior = project / "results" / "api_cache" / "nichenet_v2" / "ligand_target_matrix_nsga2r_final.rds"
    expected_cellchat = "DA6CEA9B0F8AD59A44F71733AD48148297EB10768C3A9F77918B458EC0DD8667"
    expected_nichenet = "699FCE17FF65E2511277359696306BB130AA833BD7F30DD9290AD9EFD9DC9C5D"
    check("vendored CellChat integrity", cellchat_archive.is_file() and sha256(cellchat_archive) == expected_cellchat,
          expected_cellchat)
    check("NicheNet v2 prior integrity", nichenet_prior.is_file() and sha256(nichenet_prior) == expected_nichenet,
          expected_nichenet)

    docx = outputs / "OC_OA_manuscript_revision_v30_with_figures.docx"
    pdf = outputs / "OC_OA_manuscript_revision_v30_with_figures.pdf"
    media_count = 0
    if docx.is_file():
        with ZipFile(docx) as archive:
            media_count = sum(name.startswith("word/media/") for name in archive.namelist())
    check("Word manuscript embeds 20 figures", media_count == 20, f"media={media_count}")
    check("PDF export present", pdf.is_file() and pdf.stat().st_size > 1_000_000,
          f"bytes={pdf.stat().st_size if pdf.exists() else 0}")

    text_suffixes = {".r", ".py", ".ps1", ".bat", ".md", ".txt", ".csv", ".tsv", ".json", ".yaml", ".yml", ".lock", ".log", ".env"}
    secret_hits: list[str] = []
    for path in project.rglob("*"):
        if (not path.is_file() or
                (path.suffix.lower() not in text_suffixes and path.name not in {".Renviron", ".Rprofile"}) or
                path.stat().st_size > 12_000_000):
            continue
        if JWT.search(path.read_bytes()):
            secret_hits.append(str(path.relative_to(project)))
    check("credential scan", not secret_hits, f"hits={secret_hits}")

    passed = all(item["passed"] for item in checks)
    payload = {"status": "PASS" if passed else "FAIL", "checks": checks}
    args.json.resolve().parent.mkdir(parents=True, exist_ok=True)
    args.json.resolve().write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = ["# OC-OA V3.0 submission audit", "", f"**Overall status: {'PASS' if passed else 'FAIL'}**", "",
             "| Check | Status | Detail |", "|---|---|---|"]
    for item in checks:
        detail = str(item["detail"]).replace("|", "\\|")
        lines.append(f"| {item['check']} | {'PASS' if item['passed'] else 'FAIL'} | {detail} |")
    lines.extend(["", "The audit verifies file presence, exact pivotal record counts, fixed-gene identity, dependency integrity, embedded figures, inference-boundary wording, and absence of JWT-like credentials.", ""])
    args.report.resolve().write_text("\n".join(lines), encoding="utf-8")
    print(f"Audit: {payload['status']}")
    print(f"Checks: {len(checks)}")
    print(f"Report: {args.report.resolve()}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
