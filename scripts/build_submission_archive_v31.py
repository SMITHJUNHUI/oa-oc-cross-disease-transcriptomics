#!/usr/bin/env python3
"""Build and verify the curated offline-capable V3.1 reproducibility archive."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


JWT = re.compile(rb"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}")
TEXT_SUFFIXES = {
    ".r", ".py", ".ps1", ".bat", ".md", ".txt", ".csv", ".tsv",
    ".json", ".yaml", ".yml", ".lock", ".log", ".env",
}


def add_file(files: dict[str, Path], path: Path, name: str) -> None:
    if path.is_file():
        files[name] = path


def add_tree(files: dict[str, Path], root: Path, prefix: str) -> None:
    if not root.is_dir():
        return
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        files[f"{prefix}/{path.relative_to(root).as_posix()}"] = path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--outputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    outputs = args.outputs.resolve()
    submission = project / "results" / "submission_v31"
    files: dict[str, Path] = {}

    for filename in (
        "OC_OA_manuscript_revision_v31_with_figures.docx",
        "OC_OA_manuscript_revision_v31_with_figures.pdf",
        "OC_OA_submission_audit_v31.md",
    ):
        add_file(files, outputs / filename, f"manuscript/{filename}")
    for filename in (
        "OC_OA_manuscript_revision_v31.docx",
        "OC_OA_manuscript_revision_v31.md",
        "OC_OA_revision_response_matrix_v31.md",
    ):
        add_file(files, project / "manuscript" / filename, f"manuscript/{filename}")

    add_tree(files, submission / "figures", "results/submission_v31/figures")
    add_tree(files, submission / "supplementary_tables", "results/submission_v31/supplementary_tables")
    add_tree(files, submission / "analysis", "results/submission_v31/analysis")
    add_tree(files, submission / "logs", "results/submission_v31/logs")
    for path in sorted(submission.glob("*")):
        if path.is_file():
            add_file(files, path, f"results/submission_v31/{path.name}")

    add_tree(files, project / "R", "rebuild/R")
    add_tree(files, project / "scripts", "rebuild/scripts")
    add_tree(files, project / "config", "rebuild/config")
    add_tree(files, project / "docs", "rebuild/docs")
    add_tree(files, project / "tests", "rebuild/tests")
    add_tree(files, project / "environment" / "vendor", "rebuild/environment/vendor")
    add_file(
        files,
        project / "environment" / "V30_dependency_manifest.md",
        "rebuild/environment/V30_dependency_manifest.md",
    )
    add_tree(
        files,
        project / "results" / "api_cache" / "string_v12",
        "rebuild/results/api_cache/string_v12",
    )
    add_tree(
        files,
        project / "results" / "api_cache" / "nichenet_v2",
        "rebuild/results/api_cache/nichenet_v2",
    )

    for pattern in (
        "run_*.R", "run_*.ps1", "run_*.bat", "setup.R", "setup_project.ps1",
        "setup_project.bat", "README.md", "DESCRIPTION", "renv.lock", ".Rprofile",
        "VALIDATION_REPORT.md",
    ):
        for path in sorted(project.glob(pattern)):
            add_file(files, path, f"rebuild/{path.name}")

    required = {
        "manuscript/OC_OA_manuscript_revision_v31_with_figures.docx",
        "manuscript/OC_OA_manuscript_revision_v31_with_figures.pdf",
        "manuscript/OC_OA_submission_audit_v31.md",
        "manuscript/OC_OA_revision_response_matrix_v31.md",
        "results/submission_v31/figures/Figure1_question_driven_framework.png",
        "results/submission_v31/figures/Figure6_integrated_context_model.png",
        "results/submission_v31/figures/SupplementaryFigure13_sample_consensus_CellChat.png",
        "results/submission_v31/figures/SupplementaryFigure14_NicheNet_prior_overlay.png",
        "results/submission_v31/figures/SupplementaryFigure15_panel_size_sensitivity.png",
        "results/submission_v31/figures/SupplementaryFigure16_bulk_PCA_and_QC.png",
        "results/submission_v31/supplementary_tables/Table_S29a_panel_size_composition.csv",
        "results/submission_v31/supplementary_tables/Table_S31b_bulk_sample_correlation_QC.csv",
        "results/submission_v31/submission_audit_v31.json",
        "rebuild/R/21_reviewer_v31.R",
        "rebuild/run_submission_v31.R",
        "rebuild/run_submission_v31.ps1",
        "rebuild/scripts/audit_submission_v31.py",
        "rebuild/scripts/build_manuscript_v31.py",
        "rebuild/scripts/build_manuscript_with_figures_v31.py",
        "rebuild/scripts/build_submission_archive_v31.py",
        "rebuild/docs/V31_README.md",
        "rebuild/environment/vendor/CellChat_2.2.0.9001.tar.gz",
        "rebuild/results/api_cache/nichenet_v2/ligand_target_matrix_nsga2r_final.rds",
    }
    missing = sorted(required.difference(files))
    if missing:
        raise SystemExit(f"Missing required V3.1 entries: {missing}")

    secret_hits = []
    for name, path in files.items():
        if (
            (path.suffix.lower() in TEXT_SUFFIXES or path.name in {".Renviron", ".Rprofile"})
            and path.stat().st_size <= 12_000_000
            and JWT.search(path.read_bytes())
        ):
            secret_hits.append(name)
    if secret_hits:
        raise SystemExit(f"Credential-like content found in: {secret_hits}")

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(output, "w", ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        for name, path in sorted(files.items()):
            archive.write(path, name)

    with ZipFile(output) as archive:
        names = {entry.filename for entry in archive.infolist() if not entry.is_dir()}
        if not required.issubset(names):
            raise SystemExit("V3.1 archive verification failed after writing")

    checksum = sha256(output)
    checksum_path = output.with_suffix(".sha256.txt")
    checksum_path.write_text(f"{checksum}  {output.name}\n", encoding="ascii")
    print(f"Archive: {output}")
    print(f"Files: {len(files)}")
    print(f"Bytes: {output.stat().st_size}")
    print(f"SHA256: {checksum}")


if __name__ == "__main__":
    main()
