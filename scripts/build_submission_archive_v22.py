#!/usr/bin/env python3
"""Build the curated V2.2 submission package without raw matrices or caches."""

from __future__ import annotations

import argparse
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
    submission = project / "results" / "submission_v22"
    files: dict[str, Path] = {}

    deliverables = {
        "OC_OA_manuscript_revision_v22_with_figures.docx": "manuscript/OC_OA_manuscript_revision_v22_with_figures.docx",
        "OC_OA_manuscript_revision_v22_with_figures.pdf": "manuscript/OC_OA_manuscript_revision_v22_with_figures.pdf",
        "OC_OA_manuscript_revision_v22.docx": "manuscript/OC_OA_manuscript_revision_v22.docx",
        "OC_OA_manuscript_revision_v22.pdf": "manuscript/OC_OA_manuscript_revision_v22.pdf",
        "OC_OA_manuscript_revision_v22.md": "manuscript/OC_OA_manuscript_revision_v22.md",
        "OC_OA_revision_v22_response_matrix.md": "manuscript/OC_OA_revision_v22_response_matrix.md",
        "OC_OA_submission_audit_report_v22.md": "manuscript/OC_OA_submission_audit_report_v22.md",
    }
    for filename, archive_name in deliverables.items():
        add_file(files, outputs / filename, archive_name)

    add_tree(files, submission / "figures", "figures")
    add_tree(files, submission / "supplementary_tables", "supplementary_tables")
    add_tree(files, submission / "analysis", "analysis")
    add_tree(files, submission / "logs", "logs")
    for filename in (
        "claim_evidence_registry_v22.csv",
        "data_source_manifest.csv",
        "parameter_manifest.csv",
        "reproducibility_checklist_v22.csv",
        "submission_audit_v22.json",
    ):
        add_file(files, submission / filename, f"manifests/{filename}")

    add_tree(files, project / "R", "rebuild/R")
    add_tree(files, project / "tests", "rebuild/tests")
    add_tree(files, project / "docs", "rebuild/docs")
    add_tree(files, project / "data" / "reference", "rebuild/data/reference")
    add_tree(
        files,
        project / "data" / "external" / "HPA",
        "rebuild/data/external/HPA",
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
    add_file(
        files,
        project / "manuscript" / "OC_OA_manuscript_revision_v22.md",
        "rebuild/manuscript/OC_OA_manuscript_revision_v22.md",
    )

    required = {
        "manuscript/OC_OA_manuscript_revision_v22_with_figures.docx",
        "manuscript/OC_OA_manuscript_revision_v22_with_figures.pdf",
        "manuscript/OC_OA_revision_v22_response_matrix.md",
        "manuscript/OC_OA_submission_audit_report_v22.md",
        "figures/Figure2_bulk_discovery.pdf",
        "figures/Figure4_external_validation.pdf",
        "figures/Figure5_single_cell_localization.pdf",
        "figures/SupplementaryFigure5_TCGA_relative_context.pdf",
        "figures/SupplementaryFigure6_cell_type_functional_annotation.pdf",
        "supplementary_tables/Table_S12a_MR_estimates_and_provenance.csv",
        "supplementary_tables/Table_S16_candidate_prioritization_matrix.csv",
        "supplementary_tables/Table_S17_cell_type_marker_GO_annotation.csv",
        "manifests/claim_evidence_registry_v22.csv",
        "manifests/submission_audit_v22.json",
        "rebuild/run_submission_v22.ps1",
        "rebuild/scripts/audit_submission_v22.py",
    }
    missing = sorted(required.difference(files))
    if missing:
        raise SystemExit(f"Missing required V2.2 entries: {missing}")

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
            raise SystemExit("V2.2 archive verification failed after writing.")
    print(f"Archive: {output}")
    print(f"Files: {len(files)}")
    print(f"Bytes: {output.stat().st_size}")


if __name__ == "__main__":
    main()
