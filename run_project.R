#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  repos = c(
    MRCIEU = "https://mrcieu.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match) == 0L) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  normalizePath(
    dirname(sub("^--file=", "", match[[1L]])),
    winslash = "/",
    mustWork = TRUE
  )
}

project_root <- script_path()
setwd(project_root)

r_files <- sort(list.files(
  file.path(project_root, "R"),
  pattern = "\\.[Rr]$",
  full.names = TRUE
))

if (length(r_files) == 0L) {
  stop("No R modules were found in the R/ directory.", call. = FALSE)
}

for (file in r_files) {
  sys.source(file, envir = .GlobalEnv)
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))

if (identical(cli$mode, "tests")) {
  source(file.path(project_root, "tests", "run_tests.R"), local = .GlobalEnv)
} else {
  run_pipeline(
    project_root = project_root,
    config_path = cli$config,
    mode = cli$mode,
    force = cli$force,
    from = cli$from,
    to = cli$to
  )
}

