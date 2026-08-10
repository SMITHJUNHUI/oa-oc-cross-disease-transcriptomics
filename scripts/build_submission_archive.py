#!/usr/bin/env python3
"""Build a curated revision-v2 submission archive without raw matrices or caches."""

from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


JWT_PREFIX = b"eyJhbGciOi" + b"JSUzI1Ni"


def add_tree(files: dict[str, Path], root: Path, prefix: str) -> None:
    if not root.exists():
        return
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        files[f"{prefix}/{relative}"] = path


def add_file(
    files: dict[str, Path],
    path: Path,
    archive_name: str,
) -> None:
    if path.exists() and path.is_file():
        files[archive_name] = path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--docx", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    submission = project / "results" / "submission"
    files: dict[str, Path] = {}
    add_file(files, args.docx.resolve(), "manuscript/OC_OA_manuscript_revision_v2.docx")
    add_file(files, args.pdf.resolve(), "manuscript/OC_OA_manuscript_revision_v2.pdf")
    add_file(files, args.markdown.resolve(), "manuscript/OC_OA_manuscript_revision_v2.md")
    add_file(files, args.audit.resolve(), "manuscript/OC_OA_submission_audit_report_v2.md")

    add_tree(files, submission / "figures", "figures")
    add_tree(files, submission / "supplementary_tables", "supplementary_tables")
    add_tree(files, project / "tests", "rebuild/tests")
    add_tree(files, project / "docs", "rebuild/docs")
    add_tree(files, project / "data" / "external" / "HPA", "rebuild/data/external/HPA")

    for name in (
        "claim_evidence_registry.csv",
        "data_source_manifest.csv",
        "parameter_manifest.csv",
        "reproducibility_checklist.csv",
        "submission_audit.json",
    ):
        add_file(files, submission / name, f"manifests/{name}")
    for name in (
        "non_mr_sensitivity_report.md",
        "sensitivity_parameter_manifest.csv",
        "machine_learning_repeated_cv_summary.csv",
        "machine_learning_selection_frequency.csv",
        "machine_learning_screen_frequency.csv",
        "external_validation_signed_composite_score.csv",
        "external_validation_permutation_auc.csv",
        "external_validation_leave_one_out_auc.csv",
        "hpa_normal_tissue_context.csv",
        "wgcna_module_trait_bootstrap.csv",
        "tcga_optimism_bootstrap_summary.csv",
    ):
        add_file(
            files,
            submission / "sensitivity" / name,
            f"sensitivity/{name}",
        )

    add_tree(files, project / "R", "rebuild/R")
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
        "manuscript/OC_OA_manuscript_revision_v2.docx",
        "manuscript/OC_OA_manuscript_revision_v2.pdf",
        "figures/Figure1_study_design.pdf",
        "figures/SupplementaryFigure4_HPA_normal_tissue_context.pdf",
        "supplementary_tables/Table_S14_HPA_normal_tissue_context.csv",
        "manifests/claim_evidence_registry.csv",
        "rebuild/run_submission_package.ps1",
    }
    missing = sorted(required.difference(files))
    if missing:
        raise SystemExit(f"Missing required archive entries: {missing}")

    secret_hits = [
        archive_name
        for archive_name, path in files.items()
        if JWT_PREFIX in path.read_bytes()
    ]
    if secret_hits:
        raise SystemExit(f"Credential-like content found in: {secret_hits}")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(output, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for archive_name, path in sorted(files.items()):
            archive.write(path, archive_name)

    with ZipFile(output, "r") as archive:
        names = {entry.filename for entry in archive.infolist() if not entry.is_dir()}
        if not required.issubset(names):
            raise SystemExit("Archive verification failed after writing.")
    print(f"ARCHIVE={output}")
    print(f"FILES={len(files)}")
    print(f"BYTES={output.stat().st_size}")


if __name__ == "__main__":
    main()
