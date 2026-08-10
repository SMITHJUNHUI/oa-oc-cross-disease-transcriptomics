#!/usr/bin/env python3
"""Build and verify the curated V2.4 submission/reproducibility archive."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


JWT = re.compile(
    rb"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}"
)


def add_tree(files: dict[str, Path], root: Path, prefix: str) -> None:
    if not root.exists():
        return
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        files[f"{prefix}/{path.relative_to(root).as_posix()}"] = path


def add_file(files: dict[str, Path], path: Path, archive_name: str) -> None:
    if path.exists() and path.is_file():
        files[archive_name] = path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    outputs = args.outputs.resolve()
    submission = project / "results" / "submission_v24"
    files: dict[str, Path] = {}

    add_file(
        files,
        outputs / "OC_OA_manuscript_revision_v24_with_figures.docx",
        "manuscript/OC_OA_manuscript_revision_v24_with_figures.docx",
    )
    add_file(
        files,
        outputs / "OC_OA_manuscript_revision_v24_with_figures.pdf",
        "manuscript/OC_OA_manuscript_revision_v24_with_figures.pdf",
    )
    add_file(
        files,
        project / "manuscript" / "OC_OA_manuscript_revision_v24.docx",
        "manuscript/OC_OA_manuscript_revision_v24.docx",
    )
    add_file(
        files,
        project / "manuscript" / "OC_OA_manuscript_revision_v24.md",
        "manuscript/OC_OA_manuscript_revision_v24.md",
    )
    add_file(
        files,
        project / "manuscript" / "OC_OA_revision_response_matrix_v24.md",
        "manuscript/OC_OA_revision_response_matrix_v24.md",
    )
    add_file(
        files,
        outputs / "OC_OA_submission_audit_report_v24.md",
        "manuscript/OC_OA_submission_audit_report_v24.md",
    )

    add_tree(files, submission / "figures", "figures")
    add_tree(files, submission / "supplementary_tables", "supplementary_tables")
    add_tree(files, submission / "analysis", "analysis")
    add_tree(files, submission / "logs", "logs")
    for pattern in ("*.csv", "*.json", "*.md"):
        for path in sorted(submission.glob(pattern)):
            add_file(files, path, f"manifests/{path.name}")

    add_tree(files, project / "R", "rebuild/R")
    add_tree(files, project / "tests", "rebuild/tests")
    add_tree(files, project / "docs", "rebuild/docs")
    add_tree(files, project / "data" / "reference", "rebuild/data/reference")
    add_tree(
        files,
        project / "data" / "external" / "HPA",
        "rebuild/data/external/HPA",
    )
    add_file(
        files,
        project / "data" / "raw" / "regulatory" / "human_500diff.gmt",
        "rebuild/data/raw/regulatory/human_500diff.gmt",
    )

    for path in sorted((project / "config").glob("*.yml")):
        if path.name != "local.yml":
            add_file(files, path, f"rebuild/config/{path.name}")
    for path in sorted((project / "config").glob("*.yaml")):
        add_file(files, path, f"rebuild/config/{path.name}")
    for path in sorted((project / "scripts").glob("*")):
        if path.is_file():
            add_file(files, path, f"rebuild/scripts/{path.name}")
    for pattern in (
        "run_*.R",
        "run_*.ps1",
        "run_*.bat",
        "setup.R",
        "setup_project.ps1",
        "setup_project.bat",
        "README.md",
        "DESCRIPTION",
        "renv.lock",
        ".Rprofile",
        "VALIDATION_REPORT.md",
    ):
        for path in sorted(project.glob(pattern)):
            add_file(files, path, f"rebuild/{path.name}")

    required = {
        "manuscript/OC_OA_manuscript_revision_v24_with_figures.docx",
        "manuscript/OC_OA_manuscript_revision_v24_with_figures.pdf",
        "manuscript/OC_OA_revision_response_matrix_v24.md",
        "manuscript/OC_OA_submission_audit_report_v24.md",
        "figures/Figure1_study_design.pdf",
        "figures/Figure4_external_validation.pdf",
        "figures/SupplementaryFigure9_candidate_centered_Hallmark_context.pdf",
        "figures/SupplementaryFigure10_external_evaluation_context.pdf",
        "figures/SupplementaryFigure11_upstream_regulatory_context.pdf",
        "supplementary_tables/Table_S20_candidate_centered_Hallmark_context.csv",
        "supplementary_tables/Table_S21_cross_cohort_molecular_separability_context.csv",
        "supplementary_tables/Table_S22a_external_signed_score_effect_sizes.csv",
        "supplementary_tables/Table_S22b_cross_fitted_calibration_metrics.csv",
        "supplementary_tables/Table_S23a_KnockTF_candidate_regulatory_context.csv",
        "supplementary_tables/Table_S23b_miRTarBase_candidate_regulatory_context.csv",
        "manifests/submission_audit_v24.json",
        "rebuild/run_submission_v24.R",
        "rebuild/run_submission_v24.ps1",
        "rebuild/scripts/audit_submission_v24.py",
        "rebuild/scripts/build_manuscript_v24.py",
        "rebuild/scripts/build_submission_archive_v24.py",
        "rebuild/data/raw/regulatory/human_500diff.gmt",
    }
    missing = sorted(required.difference(files))
    if missing:
        raise SystemExit(f"Missing required V2.4 entries: {missing}")

    secret_hits = [
        archive_name
        for archive_name, path in files.items()
        if JWT.search(path.read_bytes())
    ]
    if secret_hits:
        raise SystemExit(f"Credential-like content found in: {secret_hits}")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(output, "w", ZIP_DEFLATED, compresslevel=9) as archive:
        for archive_name, path in sorted(files.items()):
            archive.write(path, archive_name)

    with ZipFile(output) as archive:
        names = {entry.filename for entry in archive.infolist() if not entry.is_dir()}
        if not required.issubset(names):
            raise SystemExit("V2.4 archive verification failed after writing.")

    checksum = hashlib.sha256(output.read_bytes()).hexdigest()
    checksum_path = output.with_suffix(".sha256.txt")
    checksum_path.write_text(f"{checksum}  {output.name}\n", encoding="ascii")
    print(f"Archive: {output}")
    print(f"Files: {len(files)}")
    print(f"Bytes: {output.stat().st_size}")
    print(f"SHA256: {checksum}")


if __name__ == "__main__":
    main()
