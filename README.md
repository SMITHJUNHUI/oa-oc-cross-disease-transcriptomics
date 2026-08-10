# OA-OC cross-disease transcriptomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21876011.svg)](https://doi.org/10.5281/zenodo.21876011)

This repository accompanies the manuscript *Shared molecular features between osteoarthritis and ovarian cancer revealed by multi-layer transcriptomic analyses*.

**Authors:** Junhui Shi, Mengxiang Liu, Repkat Inayatilla, Ke Li and Lei Chen. Junhui Shi and Mengxiang Liu contributed equally. Lei Chen is the corresponding author.

The archive contains the versioned R and Python code, configuration templates, dependency records, tests and execution instructions used for the reported analyses. Raw expression files are not redistributed because all study datasets are publicly available from the NCBI Gene Expression Omnibus (GEO); the required accessions are listed in `SOURCE_DATASETS.tsv` and in the manuscript.

## Contents

- `R/`: analysis functions for bulk transcriptomics, enrichment, external validation, single-cell quality control and localization, blood validation, plotting and reporting.
- `scripts/`: manuscript-focused build and audit utilities.
- `tools/`: streaming utilities used for large single-cell matrices.
- `config/`: portable configuration templates. `local.yml` is intentionally excluded because it contains machine-specific paths.
- `tests/`: unit and regression tests.
- `manifests/`: dependency and runtime records that do not contain access credentials.
- `renv.lock`: R package lockfile.
- `run_project.R` and `run_project.ps1`: main analysis entry points.
- `run_submission_v42.R` and `run_submission_v42.ps1`: manuscript-focused result and figure rebuild.

## Reproduction outline

1. Install R 4.5.x and Python 3.8 or later with NumPy.
2. Restore the R environment from `renv.lock` by running `Rscript setup.R`.
3. Download the public GEO inputs and place them under the relative paths defined in `config/config.yml`, or copy `config/config.example.yml` to `config/local.yml` and update only local paths.
4. Run a preflight check:

   `powershell -NoProfile -ExecutionPolicy Bypass -File .\run_project.ps1 -Mode preflight`

5. Run the configured analysis:

   `powershell -NoProfile -ExecutionPolicy Bypass -File .\run_project.ps1 -Mode full`

6. Rebuild the manuscript-focused outputs:

   `powershell -NoProfile -ExecutionPolicy Bypass -File .\run_submission_v42.ps1`

   The final R build stage regenerates all seven main figures, the five title-free panel composites for Supplementary Figures S1-S5, and the five-page `Additional_file_2_supplementary_figures.pdf`.

The fixed project seed is `20260726`. Large single-cell analyses can require substantial disk space and memory; cached expression matrices and raw public inputs are not included in this archive.

## Credential handling

No access token, password or private key is included. The optional OpenGWAS module reads `OPENGWAS_JWT` from the process environment and is disabled in the portable configuration. The submitted manuscript treats Mendelian-randomization results as non-core and does not require an OpenGWAS credential for the reported main analyses.

## Scope

The project contains versioned modules retained for provenance. Some optional modules produce exploratory outputs that are not reported in the submitted article. The article-specific figures, tables and claims should be interpreted against the final manuscript and its two data/figure supplements.

## Data availability

All expression datasets used by the workflow are public GEO records. Their accession numbers and analytical roles are listed in `SOURCE_DATASETS.tsv`. No participant-level private data or redistributed raw GEO files are included here.

## Citation and archival status

Use the metadata in `CITATION.cff` when citing this software. The versioned source repository is available at https://github.com/SMITHJUNHUI/oa-oc-cross-disease-transcriptomics. The immutable v1.0.0 release is archived at https://doi.org/10.5281/zenodo.21876012. The all-version concept DOI is https://doi.org/10.5281/zenodo.21876011.

## Licence

The analysis code is distributed under the MIT License. Reused public datasets remain subject to the terms of their source repositories.
