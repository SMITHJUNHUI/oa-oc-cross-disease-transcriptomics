#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/rebuild_bmc_figures.R <project_root>", call. = FALSE)
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
setwd(project_root)

source_files <- sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE))
standalone_scripts <- c("29_final_discovery_iteration.R", "30_final_figure_layout_polish.R")
source_files <- source_files[!basename(source_files) %in% standalone_scripts]
for (source_file in source_files) {
  sys.source(source_file, envir = .GlobalEnv)
}

suppressPackageStartupMessages(library(patchwork))

initialize_logging(file.path(project_root, "results", "submission_v34"))
run_blood_systemic_v33(project_root)
run_reviewer_v34(project_root)

initialize_logging(file.path(project_root, "results", "submission_v40"))
run_reviewer_v40(project_root)

initialize_logging(file.path(project_root, "results", "submission_v41"))
run_reviewer_v41(project_root)

initialize_logging(file.path(project_root, "results", "submission_v42"))
run_reviewer_v42(project_root)

status_29 <- system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(file.path(project_root, "R", "29_final_discovery_iteration.R"), project_root)
)
if (!identical(status_29, 0L)) stop("Final discovery-figure iteration failed.", call. = FALSE)

status_30 <- system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(file.path(project_root, "R", "30_final_figure_layout_polish.R"), project_root)
)
if (!identical(status_30, 0L)) stop("Final figure-layout polish failed.", call. = FALSE)

message("BMC figure package rebuilt at the journal's 170-mm double-column width.")
