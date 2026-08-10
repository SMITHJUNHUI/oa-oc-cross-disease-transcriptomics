sc_audit_dataset_id <- function(dataset, fallback) {
  id <- dataset$id %||% fallback
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    stop("Single-cell dataset id must be one non-empty string.", call. = FALSE)
  }
  id
}

sc_audit_output_dir <- function(output_dir, dataset_id) {
  ensure_dir(file.path(output_dir, "single_cell", clean_filename(dataset_id)))
}

sc_audit_capabilities <- function(dataset_id, capability, status, reason) {
  stopifnot(
    length(capability) == length(status),
    length(status) == length(reason)
  )
  data.frame(
    dataset_id = dataset_id,
    capability = capability,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

sc_audit_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt", encoding = "UTF-8")
  } else {
    file(path, open = "rt", encoding = "UTF-8")
  }
}

sc_count_text_rows <- function(path, chunk_size = 100000L) {
  if (!file.exists(path)) {
    return(list(
      total = NA_real_,
      nonempty = NA_real_,
      blank = NA_real_,
      error = paste0("File does not exist: ", path)
    ))
  }

  connection <- sc_audit_text_connection(path)
  on.exit(close(connection), add = TRUE)
  total <- 0
  nonempty <- 0

  repeat {
    lines <- readLines(
      connection,
      n = as.integer(chunk_size),
      warn = FALSE,
      skipNul = TRUE
    )
    if (length(lines) == 0L) break
    total <- total + length(lines)
    nonempty <- nonempty + sum(nzchar(trimws(lines)))
  }

  list(
    total = as.numeric(total),
    nonempty = as.numeric(nonempty),
    blank = as.numeric(total - nonempty),
    error = ""
  )
}

sc_read_matrix_market_header <- function(path, max_header_lines = 100L) {
  if (!file.exists(path)) {
    stop("Matrix file does not exist: ", path, call. = FALSE)
  }

  connection <- sc_audit_text_connection(path)
  on.exit(close(connection), add = TRUE)
  banner <- readLines(connection, n = 1L, warn = FALSE)
  if (length(banner) != 1L) {
    stop("Matrix Market file is empty: ", path, call. = FALSE)
  }

  tokens <- strsplit(trimws(banner), "[[:space:]]+")[[1L]]
  if (
    length(tokens) != 5L ||
      !identical(tokens[[1L]], "%%MatrixMarket") ||
      !identical(tolower(tokens[[2L]]), "matrix") ||
      !identical(tolower(tokens[[3L]]), "coordinate")
  ) {
    stop(
      "Expected a Matrix Market coordinate matrix banner; observed: ",
      banner,
      call. = FALSE
    )
  }

  dimension_line <- character()
  lines_read <- 1L
  for (index in seq_len(as.integer(max_header_lines))) {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    lines_read <- lines_read + 1L
    trimmed <- trimws(line)
    if (!nzchar(trimmed) || startsWith(trimmed, "%")) next
    dimension_line <- trimmed
    break
  }
  if (length(dimension_line) == 0L) {
    stop(
      "No Matrix Market dimension line was found in the first ",
      max_header_lines,
      " header lines.",
      call. = FALSE
    )
  }

  dimensions <- suppressWarnings(scan(
    text = dimension_line,
    what = double(),
    quiet = TRUE
  ))
  if (
    length(dimensions) != 3L ||
      any(!is.finite(dimensions)) ||
      any(dimensions < 0) ||
      any(dimensions != floor(dimensions))
  ) {
    stop(
      "Invalid Matrix Market dimensions: ",
      dimension_line,
      call. = FALSE
    )
  }

  list(
    valid = TRUE,
    banner = banner,
    object = tolower(tokens[[2L]]),
    representation = tolower(tokens[[3L]]),
    field = tolower(tokens[[4L]]),
    symmetry = tolower(tokens[[5L]]),
    n_rows = dimensions[[1L]],
    n_columns = dimensions[[2L]],
    n_nonzero = dimensions[[3L]],
    header_lines_read = lines_read,
    matrix_body_read = FALSE
  )
}

sc_format_audit_count <- function(value) {
  if (length(value) == 0L || is.na(value)) return("unavailable")
  format(value, scientific = FALSE, trim = TRUE, big.mark = ",")
}

sc_atomic_fwrite_gzip <- function(x, path) {
  require_namespace("data.table", "large single-cell QC audit tables")
  ensure_dir(dirname(path))
  temporary <- paste0(path, ".tmp-", Sys.getpid(), ".gz")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(x, temporary, na = "NA", compress = "gzip")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write compressed QC table: ", path, call. = FALSE)
  }
  invisible(path)
}

audit_gse169454_sc_adapter <- function(dataset, output_dir) {
  dataset_id <- sc_audit_dataset_id(dataset, "GSE169454")
  dataset_dir <- sc_audit_output_dir(output_dir, dataset_id)
  features_path <- dataset$features_path %||% ""
  matrix_path <- dataset$matrix_path %||% ""
  barcodes_path <- dataset$barcodes_path %||% ""

  feature_rows <- sc_count_text_rows(features_path)
  barcode_rows <- sc_count_text_rows(barcodes_path)
  matrix_header <- tryCatch(
    sc_read_matrix_market_header(matrix_path),
    error = function(error) {
      list(
        valid = FALSE,
        error = conditionMessage(error),
        n_rows = NA_real_,
        n_columns = NA_real_,
        n_nonzero = NA_real_,
        field = NA_character_,
        representation = NA_character_,
        symmetry = NA_character_,
        header_lines_read = NA_integer_,
        matrix_body_read = FALSE
      )
    }
  )

  feature_match <- isTRUE(matrix_header$valid) &&
    !is.na(feature_rows$nonempty) &&
    identical(
      as.numeric(feature_rows$nonempty),
      as.numeric(matrix_header$n_rows)
    )
  barcode_match <- isTRUE(matrix_header$valid) &&
    !is.na(barcode_rows$nonempty) &&
    identical(
      as.numeric(barcode_rows$nonempty),
      as.numeric(matrix_header$n_columns)
    )
  structure_consistent <- isTRUE(feature_match) && isTRUE(barcode_match)
  status <- if (structure_consistent) {
    "series_bundle_structure_valid_not_loaded"
  } else {
    "blocked_series_bundle"
  }

  file_size <- function(path) {
    if (!file.exists(path)) return(NA_real_)
    as.numeric(file.info(path)$size)
  }
  summary <- data.frame(
    dataset_id = dataset_id,
    status = status,
    features_path = features_path,
    features_size_bytes = file_size(features_path),
    features_rows = feature_rows$nonempty,
    features_blank_rows = feature_rows$blank,
    matrix_path = matrix_path,
    matrix_size_bytes = file_size(matrix_path),
    matrix_market_valid = isTRUE(matrix_header$valid),
    matrix_field = matrix_header$field %||% NA_character_,
    matrix_rows = matrix_header$n_rows,
    matrix_columns = matrix_header$n_columns,
    matrix_nonzero = matrix_header$n_nonzero,
    matrix_header_lines_read = matrix_header$header_lines_read,
    matrix_body_read = FALSE,
    barcodes_path = barcodes_path,
    barcodes_size_bytes = file_size(barcodes_path),
    barcodes_rows = barcode_rows$nonempty,
    barcodes_blank_rows = barcode_rows$blank,
    features_match_matrix_rows = feature_match,
    barcodes_match_matrix_columns = barcode_match,
    structure_consistent = structure_consistent,
    audited_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )

  structure_reason <- if (structure_consistent) {
    "Matrix dimensions match the feature and barcode row counts."
  } else {
    paste(
      "The Matrix Market dimensions do not match both companion files:",
      paste0(
        "matrix=",
        sc_format_audit_count(matrix_header$n_rows),
        "x",
        sc_format_audit_count(matrix_header$n_columns),
        ", features=",
        sc_format_audit_count(feature_rows$nonempty),
        ", barcodes=",
        sc_format_audit_count(barcode_rows$nonempty),
        "."
      )
    )
  }
  capability <- sc_audit_capabilities(
    dataset_id,
    capability = c(
      "matrix_market_header_audit",
      "feature_barcode_mapping",
      "raw_count_matrix_loading",
      "cell_qc",
      "scDblFinder",
      "ambient_rna",
      "expression_localization"
    ),
    status = c(
      if (isTRUE(matrix_header$valid)) "available" else "blocked",
      if (structure_consistent) "available" else "blocked",
      if (structure_consistent) "not_run_by_design" else "blocked",
      "blocked",
      "blocked",
      "blocked",
      "blocked"
    ),
    reason = c(
      if (isTRUE(matrix_header$valid)) {
        paste0(
          "Read ",
          matrix_header$header_lines_read,
          " header lines only; the matrix body was not read."
        )
      } else {
        matrix_header$error %||% "Matrix Market header could not be read."
      },
      structure_reason,
      if (structure_consistent) {
        "This audit adapter validates structure only and never loads the matrix body."
      } else {
        "Gene-by-cell counts cannot be mapped safely while dimensions disagree."
      },
      "Cell-level QC metrics cannot be recomputed from an inconsistent bundle.",
      "scDblFinder requires a valid feature-by-cell raw count matrix.",
      paste(
        "The series bundle does not provide a validated unfiltered-droplet and",
        "filtered-cell pair."
      ),
      "Hub-gene expression cannot be localized without valid gene-row mapping."
    )
  )

  summary_path <- file.path(dataset_dir, "structure_audit.csv")
  capability_path <- file.path(dataset_dir, "capabilities.csv")
  report_path <- file.path(dataset_dir, "blocker_report.md")
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  write_utf8(
    c(
      paste0("# ", dataset_id, " single-cell bundle audit"),
      "",
      paste0("- Status: `", status, "`."),
      paste0(
        "- Matrix header: ",
        sc_format_audit_count(matrix_header$n_rows),
        " rows x ",
        sc_format_audit_count(matrix_header$n_columns),
        " columns; ",
        sc_format_audit_count(matrix_header$n_nonzero),
        " non-zero entries."
      ),
      paste0(
        "- Companion files: ",
        sc_format_audit_count(feature_rows$nonempty),
        " feature rows and ",
        sc_format_audit_count(barcode_rows$nonempty),
        " barcode rows."
      ),
      "- The matrix body was deliberately not read.",
      "",
      "## Decision",
      "",
      structure_reason,
      paste(
        "Do not create a Seurat or SingleCellExperiment object from this series",
        "bundle while the structural blocker remains."
      ),
      paste(
        "Use the seven per-GSM raw/filtered 10x bundles, preserve integer counts,",
        "and attach sample/condition metadata before QC."
      )
    ),
    report_path
  )

  list(
    status = status,
    capability = capability,
    summary = summary,
    thresholds = data.frame(),
    cell_metrics = data.frame(),
    artifacts = list(
      summary = summary_path,
      capability = capability_path,
      report = report_path
    )
  )
}

sc_gse180661_blocked_result <- function(
    dataset_id,
    dataset_dir,
    cells_path,
    reason
) {
  summary <- data.frame(
    dataset_id = dataset_id,
    status = "metadata_only_blocked_matrix",
    cells_path = cells_path,
    n_cells = NA_real_,
    n_samples = NA_real_,
    n_patients = NA_real_,
    duplicated_cell_ids = NA_real_,
    cells_passing_qc_audit = NA_real_,
    cells_flagged_qc_audit = NA_real_,
    raw_count_matrix_available = FALSE,
    audited_at = as.character(Sys.time()),
    details = reason,
    stringsAsFactors = FALSE
  )
  capability <- sc_audit_capabilities(
    dataset_id,
    capability = c(
      "metadata_read",
      "existing_qc_audit",
      "raw_count_matrix",
      "scDblFinder",
      "ambient_rna",
      "expression_localization"
    ),
    status = rep("blocked", 6L),
    reason = c(
      reason,
      reason,
      "No validated local raw-count matrix is configured.",
      "scDblFinder requires the missing count matrix and sample partitions.",
      "Ambient-RNA estimation requires count data; empty droplets are also absent.",
      "Hub-gene expression cannot be localized from metadata alone."
    )
  )
  summary_path <- file.path(dataset_dir, "summary.csv")
  capability_path <- file.path(dataset_dir, "capabilities.csv")
  report_path <- file.path(dataset_dir, "blocker_report.md")
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  write_utf8(
    c(
      paste0("# ", dataset_id, " metadata audit"),
      "",
      "- Status: `metadata_only_blocked_matrix`.",
      paste0("- Blocker: ", reason)
    ),
    report_path
  )
  list(
    status = "metadata_only_blocked_matrix",
    capability = capability,
    summary = summary,
    thresholds = data.frame(),
    cell_metrics = data.frame(),
    artifacts = list(
      summary = summary_path,
      capability = capability_path,
      report = report_path
    )
  )
}

audit_gse180661_sc_adapter <- function(
    dataset,
    output_dir,
    lower_nmads = 3,
    upper_nmads = 5,
    mt_nmads = 3
) {
  require_namespace("data.table", "streamlined single-cell metadata audit")
  dataset_id <- sc_audit_dataset_id(dataset, "GSE180661")
  dataset_dir <- sc_audit_output_dir(output_dir, dataset_id)
  cells_path <- dataset$cells_path %||% dataset$metadata_path %||% ""
  if (!file.exists(cells_path)) {
    return(sc_gse180661_blocked_result(
      dataset_id,
      dataset_dir,
      cells_path,
      paste0("Cell metadata file does not exist: ", cells_path)
    ))
  }

  header <- names(data.table::fread(
    cells_path,
    nrows = 0L,
    showProgress = FALSE
  ))
  required <- c(
    "cell_id", "sample", "percent.mt", "nCount_RNA", "nFeature_RNA"
  )
  missing_required <- setdiff(required, header)
  if (length(missing_required) > 0L) {
    return(sc_gse180661_blocked_result(
      dataset_id,
      dataset_dir,
      cells_path,
      paste0(
        "Required metadata columns are missing: ",
        paste(missing_required, collapse = ", "),
        "."
      )
    ))
  }

  optional <- c("patient_id", "cell_type", "cell_type_super")
  selected <- c(required, intersect(optional, header))
  metrics <- data.table::fread(
    cells_path,
    select = selected,
    na.strings = c("", "NA"),
    showProgress = FALSE
  )
  metrics[, `:=`(
    batch = as.character(sample),
    nCount = suppressWarnings(as.numeric(nCount_RNA)),
    nFeature = suppressWarnings(as.numeric(nFeature_RNA)),
    percent_mt = suppressWarnings(as.numeric(percent.mt))
  )]
  metrics[, dataset_id := dataset_id]

  missing_core <- sum(
    is.na(metrics$cell_id) |
      is.na(metrics$batch) |
      is.na(metrics$nCount) |
      is.na(metrics$nFeature) |
      is.na(metrics$percent_mt)
  )
  if (missing_core > 0L) {
    return(sc_gse180661_blocked_result(
      dataset_id,
      dataset_dir,
      cells_path,
      paste0(
        missing_core,
        " metadata rows have missing or non-numeric core QC values."
      )
    ))
  }

  for (helper in c(
    "derive_sc_qc_thresholds",
    "apply_sc_qc_thresholds",
    "plot_sc_qc_metrics"
  )) {
    if (!exists(helper, mode = "function", inherits = TRUE)) {
      stop("Required single-cell QC helper is unavailable: ", helper, call. = FALSE)
    }
  }

  thresholds <- derive_sc_qc_thresholds(
    metrics,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt",
    lower_nmads = lower_nmads,
    upper_nmads = upper_nmads,
    mt_nmads = mt_nmads
  )
  audited <- apply_sc_qc_thresholds(
    metrics,
    thresholds,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt"
  )
  if (is.null(audited)) audited <- metrics
  if (nrow(audited) != nrow(metrics)) {
    stop(
      "QC audit must flag cells without dropping rows; expected ",
      nrow(metrics),
      " rows but received ",
      nrow(audited),
      ".",
      call. = FALSE
    )
  }
  if (!"dataset_id" %in% names(thresholds)) {
    thresholds$dataset_id <- dataset_id
  }

  plot_error <- ""
  plot_result <- tryCatch(
    plot_sc_qc_metrics(
      audited,
      thresholds = thresholds,
      output_dir = dataset_dir,
      dataset_id = dataset_id,
      batch_col = "batch",
      count_col = "nCount",
      feature_col = "nFeature",
      mt_col = "percent_mt"
    ),
    error = function(error) {
      plot_error <<- conditionMessage(error)
      log_warn(dataset_id, " QC audit plots failed: ", plot_error)
      NULL
    }
  )

  pass_candidates <- names(audited)[grepl(
    "(pass.*qc|qc.*pass)",
    names(audited),
    ignore.case = TRUE
  )]
  pass_column <- if (length(pass_candidates) > 0L) {
    pass_candidates[[1L]]
  } else {
    NA_character_
  }
  pass_values <- if (!is.na(pass_column)) {
    as.logical(audited[[pass_column]])
  } else {
    rep(NA, nrow(audited))
  }
  n_pass <- if (all(is.na(pass_values))) NA_real_ else sum(pass_values, na.rm = TRUE)
  n_flagged <- if (is.na(n_pass)) NA_real_ else nrow(audited) - n_pass
  n_patients <- if ("patient_id" %in% names(audited)) {
    data.table::uniqueN(audited$patient_id, na.rm = TRUE)
  } else {
    NA_integer_
  }

  matrix_path <- dataset$matrix_h5_path %||%
    dataset$h5_path %||%
    dataset$count_matrix_path %||%
    ""
  matrix_available <- nzchar(matrix_path) && file.exists(matrix_path)
  summary <- data.frame(
    dataset_id = dataset_id,
    status = "metadata_only_blocked_matrix",
    cells_path = cells_path,
    cells_size_bytes = as.numeric(file.info(cells_path)$size),
    n_cells = nrow(audited),
    n_samples = data.table::uniqueN(audited$batch),
    n_patients = n_patients,
    duplicated_cell_ids = sum(duplicated(audited$cell_id)),
    qc_pass_column = pass_column,
    cells_passing_qc_audit = n_pass,
    cells_flagged_qc_audit = n_flagged,
    raw_count_matrix_available = matrix_available,
    matrix_path = matrix_path,
    thresholds_are_audit_flags = TRUE,
    cells_removed = 0,
    audited_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  capability <- sc_audit_capabilities(
    dataset_id,
    capability = c(
      "metadata_read",
      "existing_qc_audit",
      "qc_threshold_plots",
      "raw_count_matrix",
      "scDblFinder",
      "ambient_rna",
      "expression_localization"
    ),
    status = c(
      "available",
      "available",
      if (nzchar(plot_error)) "failed" else "available",
      "blocked",
      "blocked",
      "blocked",
      "blocked"
    ),
    reason = c(
      paste0(nrow(audited), " metadata rows were read without dropping cells."),
      paste(
        "Per-sample MAD thresholds audit existing nCount_RNA, nFeature_RNA,",
        "and percent.mt values."
      ),
      if (nzchar(plot_error)) plot_error else "Threshold-justification plots were written.",
      if (matrix_available) {
        "A matrix path exists but this blocked adapter does not validate or load it."
      } else {
        "Only cells.tsv.gz is configured locally; the raw UMI H5 matrix is absent."
      },
      "scDblFinder is blocked until the raw UMI matrix is validated and joined by cell_id.",
      paste(
        "Ambient-RNA correction is blocked without count data and an",
        "unfiltered-droplet matrix."
      ),
      "Hub-gene expression cannot be localized from cell metadata alone."
    )
  )

  thresholds_path <- file.path(dataset_dir, "qc_thresholds.csv")
  metrics_path <- file.path(dataset_dir, "cell_qc_audit.csv.gz")
  summary_path <- file.path(dataset_dir, "summary.csv")
  capability_path <- file.path(dataset_dir, "capabilities.csv")
  report_path <- file.path(dataset_dir, "blocker_report.md")
  safe_write_csv(thresholds, thresholds_path)
  sc_atomic_fwrite_gzip(audited, metrics_path)
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  write_utf8(
    c(
      paste0("# ", dataset_id, " metadata-only QC audit"),
      "",
      "- Status: `metadata_only_blocked_matrix`.",
      paste0("- Cells audited: ", sc_format_audit_count(nrow(audited)), "."),
      paste0("- Sample partitions: ", data.table::uniqueN(audited$batch), "."),
      paste0("- Cells removed: 0; QC results are audit flags only."),
      "",
      "## Capabilities",
      "",
      paste(
        "Existing nCount_RNA, nFeature_RNA and percent.mt values were audited",
        "with data-driven MAD thresholds within each sample."
      ),
      paste(
        "scDblFinder, ambient-RNA correction and expression localization remain",
        "blocked because the raw UMI matrix is not validated by this adapter."
      )
    ),
    report_path
  )

  list(
    status = "metadata_only_blocked_matrix",
    capability = capability,
    summary = summary,
    thresholds = thresholds,
    cell_metrics = audited,
    artifacts = list(
      thresholds = thresholds_path,
      cell_metrics = metrics_path,
      summary = summary_path,
      capability = capability_path,
      report = report_path,
      plots = plot_result
    )
  )
}

audit_blocked_sc_adapter <- function(dataset, output_dir, ...) {
  dataset_id <- toupper(sc_audit_dataset_id(dataset, ""))
  switch(
    dataset_id,
    GSE169454 = audit_gse169454_sc_adapter(dataset, output_dir),
    GSE180661 = audit_gse180661_sc_adapter(dataset, output_dir, ...),
    stop(
      "No blocked single-cell adapter is implemented for dataset: ",
      dataset_id,
      call. = FALSE
    )
  )
}
