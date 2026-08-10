v40_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v40"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs")),
    reference_audit = ensure_dir(file.path(root, "reference_audit"))
  )
}

v40_copy_tree <- function(source, target) {
  if (!dir.exists(source)) stop("Missing V3.4 source directory: ", source, call. = FALSE)
  ensure_dir(target)
  directories <- list.dirs(source, recursive = TRUE, full.names = TRUE)
  for (directory in directories) {
    relative <- substring(directory, nchar(source) + 2L)
    if (nzchar(relative)) ensure_dir(file.path(target, relative))
  }
  files <- list.files(source, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  for (file in files) {
    relative <- substring(file, nchar(source) + 2L)
    destination <- file.path(target, relative)
    ensure_dir(dirname(destination))
    if (!file.copy(file, destination, overwrite = TRUE)) {
      stop("Could not copy V3.4 artifact into V4.0: ", file, call. = FALSE)
    }
  }
  invisible(target)
}

v40_figure4_builder <- function() {
  source_text <- paste(deparse(v34_build_figure4), collapse = "\n")
  source_text <- gsub(
    "Representative genes, not a predictive signature",
    "Candidate molecular features across evidence layers",
    source_text,
    fixed = TRUE
  )
  eval(parse(text = source_text), envir = .GlobalEnv)
}

v40_write_documentation <- function(project_root, paths) {
  legends_source <- file.path(project_root, "manuscript", "figure_legends_v40.md")
  if (!file.copy(legends_source, file.path(paths$figures, "figure_legends.md"), overwrite = TRUE)) {
    stop("Could not install V4.0 figure legends.", call. = FALSE)
  }

  terminology <- data.frame(
    concept = c(
      "cross-disease overlap", "functional recurrence", "five-gene summary",
      "blood result", "single-cell populations", "overall interpretation"
    ),
    canonical_term = c(
      "shared molecular features", "recurring biological themes", "representative molecular features",
      "candidate systemic molecular signal", "exact source-defined cell populations",
      "partial molecular convergence with context dependence"
    ),
    avoided_term = c(
      "shared disease mechanism", "activated mechanism", "diagnostic signature",
      "blood biomarker", "forced homologous cell labels", "causal relationship"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(terminology, file.path(paths$analysis, "V40_terminology_ledger.csv"))

  registry <- data.frame(
    claim = c(
      "286 shared tissue genes", "146 concordant and 140 discordant genes",
      "external tissue replication differs between OA and OC",
      "representative genes occupy distinct cellular contexts",
      "G0S2 is the sole independent dual-blood-FDR result"
    ),
    evidence = c(
      "Figure 2; Table S2", "Figure 2; Table S2", "Figure 4; Tables S4-S6",
      "Figure 5; Tables S8a-S8b", "Figure 6; Tables S9-S11"
    ),
    boundary = c(
      "overlap is threshold dependent", "direction is not a mechanistic assay",
      "OA replication is weaker and cohort dependent", "atlas labels are not equated across diseases",
      "one candidate signal requires prospective and protein-level confirmation"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v40.csv"))
}

run_reviewer_v40 <- function(project_root) {
  source_root <- file.path(project_root, "results", "submission_v34")
  paths <- v40_output_paths(project_root)
  v40_copy_tree(source_root, paths$root)

  external <- v34_external_validation(project_root, paths)
  genes <- c("G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3")
  single_cell <- v34_single_cell_candidates(project_root, paths, genes)
  candidates <- v34_candidate_evidence(paths, external, single_cell, genes)
  v40_figure4_builder()(paths, external, candidates)
  v40_write_documentation(project_root, paths)

  log_info("V4.0 figure package completed: data are unchanged from V3.4; Figure 4 wording and all legends were aligned with the converged manuscript terminology.")
  invisible(paths)
}
