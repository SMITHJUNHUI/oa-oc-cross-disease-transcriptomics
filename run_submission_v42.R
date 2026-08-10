#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_CTYPE", "English_United States.utf8"), silent = TRUE)
arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
project_root <- if (length(file_argument)) normalizePath(dirname(sub("^--file=", "", file_argument[[1L]])), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(project_root)
standalone_scripts <- c("29_final_discovery_iteration.R", "30_final_figure_layout_polish.R")
source_files <- sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE))
source_files <- source_files[!basename(source_files) %in% standalone_scripts]
for (file in source_files) sys.source(file, envir = .GlobalEnv)
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
