#!/usr/bin/env python3
"""Audit manuscript-facing claims, source data, figures, and disclosure boundaries."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, label: str, details: str, checks: list[dict[str, str]]) -> None:
    checks.append(
        {
            "check": label,
            "status": "PASS" if condition else "FAIL",
            "details": details,
        }
    )


def numeric_equal(value: float, expected: float, tolerance: float = 1e-9) -> bool:
    return math.isfinite(value) and abs(value - expected) <= tolerance


def stream_contains(path: Path, needle: bytes) -> bool:
    overlap = max(0, len(needle) - 1)
    previous = b""
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                return False
            block = previous + chunk
            if needle in block:
                return True
            previous = block[-overlap:] if overlap else b""


def zip_entry_contains(archive: ZipFile, name: str, needle: bytes) -> bool:
    overlap = max(0, len(needle) - 1)
    previous = b""
    with archive.open(name, "r") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                return False
            block = previous + chunk
            if needle in block:
                return True
            previous = block[-overlap:] if overlap else b""


def auditable_files(project: Path, outputs: Path) -> list[Path]:
    candidates: set[Path] = set()
    subtrees = [
        project / "R",
        project / "config",
        project / "docs",
        project / "scripts",
        project / "tests",
        project / "tools",
        project / "manuscript",
        project / "results" / "logs",
        project / "results" / "submission",
    ]
    for subtree in subtrees:
        if subtree.exists():
            candidates.update(path for path in subtree.rglob("*") if path.is_file())
    candidates.update(path for path in project.iterdir() if path.is_file())
    candidates.update(
        path
        for path in (
            outputs / "OC_OA_manuscript_revision_v2.docx",
            outputs / "OC_OA_manuscript_revision_v2.pdf",
            outputs / "OC_OA_manuscript_revision_v2.md",
        )
        if path.exists()
    )
    return sorted(candidates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--manuscript", required=True, type=Path)
    parser.add_argument("--output-md", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    outputs = project.parent.resolve()
    submission = project / "results" / "submission"
    sensitivity = submission / "sensitivity"
    figures = submission / "figures"
    supplementary = submission / "supplementary_tables"
    manuscript = args.manuscript.resolve()
    checks: list[dict[str, str]] = []

    # Primary discovery counts.
    deg = {row["disease"]: row for row in read_csv(project / "results" / "tables" / "DEG_selection_summary.csv")}
    require(int(deg["OA"]["significant_genes"]) == 2008, "OA primary DEG count", "2,008", checks)
    require(int(deg["OC"]["significant_genes"]) == 2310, "OC primary DEG count", "2,310", checks)
    require(int(deg["shared"]["significant_genes"]) == 286, "Primary shared DEG count", "286", checks)

    threshold = read_csv(sensitivity / "deg_threshold_sensitivity_summary.csv")
    primary_threshold = next(row for row in threshold if row["is_primary"].upper() == "TRUE")
    require(
        int(primary_threshold["directionally_concordant_count"]) == 146,
        "Directionally concordant shared genes",
        "146/286",
        checks,
    )
    require(
        int(primary_threshold["retained_hub_count"]) == 10,
        "Primary hub-gene retention",
        "10/10",
        checks,
    )

    # WGCNA stability.
    wgcna = {row["dataset_id"]: row for row in read_csv(sensitivity / "wgcna_module_trait_bootstrap.csv")}
    require(
        numeric_equal(float(wgcna["GSE114007"]["observed_correlation"]), -0.950676532217701, 1e-12),
        "OA WGCNA module–trait correlation",
        "r = -0.951",
        checks,
    )
    require(
        numeric_equal(float(wgcna["GSE18520"]["observed_correlation"]), -0.878533931987627, 1e-12),
        "OC WGCNA module–trait correlation",
        "r = -0.879",
        checks,
    )
    require(
        all(float(row["sign_stability"]) == 1.0 for row in wgcna.values()),
        "WGCNA bootstrap sign stability",
        "1.000 in both discovery datasets",
        checks,
    )

    # Strict nested ML audit boundary.
    ml = read_csv(sensitivity / "machine_learning_repeated_cv_summary.csv")
    require(
        len(ml) == 4 and all(float(row["auc_median"]) == 1.0 for row in ml),
        "Strict nested internal ML AUCs",
        "Four disease/model summaries; median AUC 1.000",
        checks,
    )
    require(
        all(row["candidate_space"] == "all measured genes" for row in ml)
        and all(
            row["feature_selection_scope"] == "outer training fold only"
            for row in ml
        ),
        "Strict nested feature-selection scope",
        "All measured genes; screening repeated inside each outer training fold",
        checks,
    )

    # Direction-fixed external composite.
    external = {
        row["dataset_id"]: row
        for row in read_csv(sensitivity / "external_validation_signed_composite_score.csv")
    }
    expected_auc = {
        "GSE117999": 0.520,
        "GSE82107": 0.628571428571429,
        "GSE54388": 1.000,
        "GSE12470": 0.97906976744186,
    }
    require(
        all(numeric_equal(float(external[key]["auc"]), expected, 1e-12) for key, expected in expected_auc.items()),
        "Direction-fixed external composite AUCs",
        "OA 0.520/0.629; OC 1.000/0.979",
        checks,
    )
    require(
        numeric_equal(
            float(external["GSE54388"]["permutation_empirical_p"]),
            1 / 1001,
            1e-12,
        )
        and numeric_equal(
            float(external["GSE12470"]["permutation_empirical_p"]),
            1 / 1001,
            1e-12,
        ),
        "OC external score permutation tests",
        "Both empirical P=0.001 with 1,000 permutations",
        checks,
    )
    require(
        float(external["GSE54388"]["leave_one_out_auc_minimum"]) == 1.0
        and float(external["GSE12470"]["leave_one_out_auc_minimum"]) >= 0.97,
        "OC leave-one-sample-out influence",
        "Minimum AUC 1.000 and 0.971",
        checks,
    )

    hpa = read_csv(sensitivity / "hpa_normal_tissue_context.csv")
    require(
        len(hpa) == 10
        and sum(
            row["ovary_listed_as_specific"].upper() == "TRUE" for row in hpa
        )
        == 1
        and all(
            row["cartilage_in_reference"].upper() == "FALSE" for row in hpa
        ),
        "HPA normal-tissue context boundary",
        "10 genes; only DIRAS3 lists ovary; cartilage absent",
        checks,
    )

    # TCGA.
    tcga = read_csv(sensitivity / "tcga_cox_model_sensitivity.csv")
    adjusted_risk = next(
        row
        for row in tcga
        if row["term"] == "risk_z" and row["model"] == "risk_adjusted_for_age_and_stage"
    )
    require(
        numeric_equal(float(adjusted_risk["hazard_ratio"]), 1.26160811049363, 1e-12),
        "TCGA age/stage-adjusted risk association",
        "HR 1.262 (95% CI 1.076–1.479)",
        checks,
    )
    tcga_optimism = read_csv(sensitivity / "tcga_optimism_bootstrap_summary.csv")
    corrected = float(tcga_optimism[0]["optimism_corrected_cindex"])
    require(
        numeric_equal(corrected, 0.586430939510612, 1e-12),
        "TCGA optimism-corrected discrimination",
        "C-index 0.586",
        checks,
    )

    # Single-cell gate totals.
    sc_rows = [
        row
        for row in read_csv(supplementary / "Table_S9_single_cell_QC_and_status.csv")
        if row["source_table"] == "qc_status"
    ]
    total_cells = sum(int(row["cells"]) for row in sc_rows)
    qc_pass = sum(int(row["qc_pass"]) for row in sc_rows)
    require(total_cells == 1187436, "Single-cell audited input total", "1,187,436 cells", checks)
    require(qc_pass == 1025361, "Single-cell QC-pass total", "1,025,361 cells", checks)
    require(
        len(sc_rows) == 5 and all(row["downstream_ready"].upper() == "TRUE" for row in sc_rows),
        "Single-cell downstream gate",
        "5/5 datasets downstream-ready",
        checks,
    )

    # MR is present and null in both directions.
    mr = read_csv(project / "results" / "mr" / "MR_combined_estimates.csv")
    ivw = [row for row in mr if row["method"] == "Inverse variance weighted"]
    require(
        len(ivw) == 2 and all(float(row["pval"]) > 0.05 for row in ivw),
        "Bidirectional IVW MR",
        "Both P values >0.05; reported as no detected causal evidence",
        checks,
    )

    # Manuscript-facing artifacts.
    main_figures = [f"Figure{i}" for i in range(1, 7)]
    supp_figures = [f"SupplementaryFigure{i}" for i in range(1, 5)]
    for prefix in main_figures + supp_figures:
        pdf_matches = list(figures.glob(f"{prefix}_*.pdf"))
        png_matches = list(figures.glob(f"{prefix}_*.png"))
        require(
            len(pdf_matches) == 1 and len(png_matches) == 1,
            f"{prefix} vector/raster pair",
            "one PDF and one 300-dpi PNG",
            checks,
        )
    source_count = len(list((figures / "source_data").glob("*.csv")))
    require(source_count == 32, "Exact figure source-data tables", "32 CSV files", checks)

    supplementary_csvs = list(supplementary.glob("Table_S*.csv"))
    require(
        len(supplementary_csvs) == 20,
        "Supplementary table files",
        "20 CSV files representing 17 indexed table entries",
        checks,
    )
    claim_registry = read_csv(submission / "claim_evidence_registry.csv")
    require(
        len(claim_registry) == 13
        and all(row["status"] == "verified against current outputs" for row in claim_registry),
        "Claim–evidence registry",
        "13/13 claims verified",
        checks,
    )

    package_zip = outputs / "OC_OA_submission_package_v2.zip"

    text = manuscript.read_text(encoding="utf-8")
    required_boundaries = [
        "feature screening was repeated",
        "not a clinically validated diagnostic panel",
        "not a clinically ready prognostic model",
        "not proof that a causal effect is absent",
        "does not establish a causal relationship",
        "OA and OC datasets were not integrated into a shared latent space",
        "cartilage is absent",
        "molecular separation rather than clinical diagnostic validation",
    ]
    for phrase in required_boundaries:
        require(
            phrase.casefold() in text.casefold(),
            f"Required interpretation boundary: {phrase}",
            "present in manuscript",
            checks,
        )

    forbidden_patterns = [
        r"\bOA causes OC\b",
        r"\bOC causes OA\b",
        r"\bMR proves\b",
        r"\buniversally validated signature\b",
        r"\breverse causality established\b",
        r"\bimmune mechanism proven\b",
        r"\bcell type causes disease\b",
        r"\bcausal modules?\b",
    ]
    forbidden_hits = [
        pattern
        for pattern in forbidden_patterns
        if re.search(pattern, text, flags=re.IGNORECASE)
    ]
    require(
        not forbidden_hits,
        "Affirmative causal/clinical overclaim scan",
        "no prohibited affirmative wording" if not forbidden_hits else "; ".join(forbidden_hits),
        checks,
    )

    # Credential audit. The user supplied a JWT in chat; it must never enter deliverables.
    jwt_prefix = b"eyJhbGciOi" + b"JSUzI1Ni"
    secret_hits = [
        str(path)
        for path in auditable_files(project, outputs)
        if stream_contains(path, jwt_prefix)
    ]
    if package_zip.exists():
        with ZipFile(package_zip, "r") as archive:
            secret_hits.extend(
                f"{package_zip}!{entry.filename}"
                for entry in archive.infolist()
                if not entry.is_dir() and zip_entry_contains(archive, entry.filename, jwt_prefix)
            )
    require(
        not secret_hits,
        "OpenGWAS JWT leakage scan",
        "0 credential-prefix matches" if not secret_hits else "; ".join(secret_hits),
        checks,
    )

    failed = [check for check in checks if check["status"] != "PASS"]
    report_lines = [
        "# Submission audit report",
        "",
        f"- Generated: {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %Z')}",
        f"- Overall status: **{'PASS' if not failed else 'FAIL'}**",
        f"- Checks passed: **{len(checks) - len(failed)}/{len(checks)}**",
        "- Scope: numerical claims, non-causal interpretation boundaries, figure/table completeness, and credential leakage.",
        "",
        "| Check | Status | Evidence |",
        "|---|---|---|",
    ]
    for check in checks:
        report_lines.append(
            f"| {check['check']} | {check['status']} | {check['details']} |"
        )
    report_lines.extend(
        [
            "",
            "## Human sign-off still required",
            "",
            "- Author names, affiliations, corresponding-author details, and CRediT roles.",
            "- Target journal and journal-specific formatting, word limits, and reporting checklist.",
            "- Funding, competing interests, and the submitting institution’s secondary-analysis ethics statement.",
            "- Final scientific and clinical review by all authors before submission.",
            "",
            "A passing computational audit verifies consistency with the current project outputs; it is not a substitute for author accountability or peer review.",
        ]
    )

    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    args.output_json.write_text(
        json.dumps(
            {
                "generated": datetime.now().astimezone().isoformat(),
                "overall_status": "PASS" if not failed else "FAIL",
                "passed": len(checks) - len(failed),
                "total": len(checks),
                "checks": checks,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"STATUS={'PASS' if not failed else 'FAIL'}")
    print(f"CHECKS={len(checks) - len(failed)}/{len(checks)}")
    print(f"REPORT={args.output_md.resolve()}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
