.oa_sc_require_function <- function(name) {
  if (!exists(name, mode = "function", inherits = TRUE)) {
    stop(
      "Single-cell adapter requires helper function: ", name, "().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.oa_sc_output_dir <- function(config, dataset_id) {
  ensure_dir(file.path(
    config$project$output_dir,
    "single_cell",
    clean_filename(dataset_id)
  ))
}

.oa_sc_capability_table <- function(dataset_id, capability, status, reason) {
  data.frame(
    dataset_id = dataset_id,
    capability = capability,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

.oa_sc_empty_thresholds <- function() {
  data.frame(
    batch = character(),
    n_cells = integer(),
    nCount_lower = numeric(),
    nCount_upper = numeric(),
    nFeature_lower = numeric(),
    nFeature_upper = numeric(),
    percent_mt_upper = numeric(),
    stringsAsFactors = FALSE
  )
}

.oa_sc_tar_members <- function(path) {
  if (!file.exists(path)) {
    stop("Archive does not exist: ", path, call. = FALSE)
  }
  members <- utils::untar(path, list = TRUE)
  members <- members[nzchar(members) & !grepl("/$", members)]
  if (length(members) == 0L) {
    stop("Archive contains no regular-file members: ", path, call. = FALSE)
  }
  members
}

.oa_sc_tar_executable <- function() {
  candidates <- c(Sys.which("tar"), Sys.which("bsdtar"))
  executable <- candidates[nzchar(candidates)][1L]
  if (is.na(executable) || !nzchar(executable)) {
    stop(
      "A tar/bsdtar executable is required for streaming archive members.",
      call. = FALSE
    )
  }
  executable
}

.oa_sc_tar_member_command <- function(archive, member) {
  executable <- .oa_sc_tar_executable()
  if (.Platform$OS.type == "windows") {
    # cmd.exe misparses a command whose first token is quoted. The system tar
    # path is converted to a short path, while archive/member remain arguments.
    executable <- utils::shortPathName(executable)
    return(paste(
      executable,
      "-xOf",
      shQuote(archive, type = "cmd"),
      shQuote(member, type = "cmd")
    ))
  }
  paste(shQuote(executable), "-xOf", shQuote(archive), shQuote(member))
}

.oa_sc_first_existing_column <- function(table, candidates) {
  match <- candidates[candidates %in% names(table)][1L]
  if (length(match) == 0L || is.na(match)) NULL else match
}

.oa_sc_read_gse255460_metadata <- function(path) {
  if (!file.exists(path)) {
    stop("GSE255460 metadata does not exist: ", path, call. = FALSE)
  }
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  metadata <- utils::read.csv(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required <- c("nCount_RNA", "nFeature_RNA", "ID")
  missing <- setdiff(required, names(metadata))
  if (length(missing) > 0L) {
    stop(
      "GSE255460 metadata is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  cell_column <- .oa_sc_first_existing_column(
    metadata,
    c("cell_id", "barcode", "X")
  )
  if (is.null(cell_column)) {
    stop(
      "GSE255460 metadata has no cell identifier column (cell_id/barcode/X).",
      call. = FALSE
    )
  }
  metadata$.cell_id <- as.character(metadata[[cell_column]])
  if (
    anyNA(metadata$.cell_id) ||
      any(!nzchar(metadata$.cell_id)) ||
      anyDuplicated(metadata$.cell_id)
  ) {
    stop("GSE255460 metadata cell identifiers are invalid or duplicated.", call. = FALSE)
  }
  metadata
}

.oa_sc_scan_gse255460_counts <- function(
    archive,
    member,
    metadata,
    spotcheck_cells = 5L,
    stream_chunk_lines = 25L,
    structure_check_every = 5000L
) {
  spotcheck_cells <- max(1L, min(as.integer(spotcheck_cells), nrow(metadata)))
  stream_chunk_lines <- max(1L, as.integer(stream_chunk_lines))
  structure_check_every <- max(1L, as.integer(structure_check_every))

  command <- .oa_sc_tar_member_command(archive, member)
  connection <- pipe(command, open = "rt", encoding = "UTF-8")
  on.exit(
    if (!is.null(connection)) try(close(connection), silent = TRUE),
    add = TRUE
  )

  header <- readLines(connection, n = 1L, warn = FALSE)
  if (length(header) != 1L || !nzchar(header)) {
    stop("GSE255460 count member has no header line.", call. = FALSE)
  }
  matrix_cells <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  expected_cells <- make.names(metadata$.cell_id, unique = FALSE)
  if (anyDuplicated(matrix_cells)) {
    stop("GSE255460 matrix header contains duplicated cell names.", call. = FALSE)
  }
  if (!identical(matrix_cells, expected_cells)) {
    stop(
      "GSE255460 matrix columns do not exactly match make.names(metadata cell IDs) ",
      "in the same order.",
      call. = FALSE
    )
  }

  selected <- seq_len(spotcheck_cells)
  observed_counts <- numeric(spotcheck_cells)
  observed_features <- integer(spotcheck_cells)
  integer_nonnegative <- TRUE
  genes <- character()
  row_index <- 0L
  structure_rows <- integer()
  structure_fields <- integer()
  last_line <- NULL

  repeat {
    lines <- readLines(connection, n = stream_chunk_lines, warn = FALSE)
    if (length(lines) == 0L) break
    for (line in lines) {
      row_index <- row_index + 1L
      last_line <- line
      tokens <- scan(
        text = line,
        what = character(),
        sep = "\t",
        quote = "",
        nmax = spotcheck_cells + 1L,
        quiet = TRUE
      )
      if (length(tokens) != spotcheck_cells + 1L) {
        stop(
          "GSE255460 count row ", row_index,
          " ended before the requested spot-check columns.",
          call. = FALSE
        )
      }
      gene <- tokens[[1L]]
      if (!nzchar(gene)) {
        stop("GSE255460 contains an empty gene name at row ", row_index, ".", call. = FALSE)
      }
      values_text <- tokens[-1L]
      values <- suppressWarnings(as.numeric(values_text))
      valid <- !is.na(values) &
        values >= 0 &
        abs(values - round(values)) <= sqrt(.Machine$double.eps)
      integer_nonnegative <- integer_nonnegative && all(valid)
      if (!all(valid)) {
        stop(
          "GSE255460 count spot-check found a negative, fractional, or nonnumeric ",
          "value at gene row ", row_index, ".",
          call. = FALSE
        )
      }
      observed_counts <- observed_counts + values
      observed_features <- observed_features + as.integer(values > 0)
      genes[[row_index]] <- gene

      if (row_index <= 3L || row_index %% structure_check_every == 0L) {
        structure_rows <- c(structure_rows, row_index)
        structure_fields <- c(
          structure_fields,
          length(strsplit(line, "\t", fixed = TRUE)[[1L]])
        )
      }
    }
  }
  close(connection)
  connection <- NULL

  if (
    row_index > 0L &&
      !row_index %in% structure_rows &&
      !is.null(last_line)
  ) {
    structure_rows <- c(structure_rows, row_index)
    structure_fields <- c(
      structure_fields,
      length(strsplit(last_line, "\t", fixed = TRUE)[[1L]])
    )
  }
  expected_fields <- length(matrix_cells) + 1L
  if (length(structure_fields) == 0L || any(structure_fields != expected_fields)) {
    stop(
      "GSE255460 sampled count rows do not have the expected ",
      expected_fields, " tab-delimited fields.",
      call. = FALSE
    )
  }
  if (anyDuplicated(genes)) {
    stop("GSE255460 matrix contains duplicated gene names.", call. = FALSE)
  }

  metadata_counts <- as.numeric(metadata$nCount_RNA[selected])
  metadata_features <- as.integer(metadata$nFeature_RNA[selected])
  spotcheck <- data.frame(
    cell_id = metadata$.cell_id[selected],
    matrix_cell_id = matrix_cells[selected],
    observed_nCount = observed_counts,
    metadata_nCount = metadata_counts,
    observed_nFeature = observed_features,
    metadata_nFeature = metadata_features,
    nCount_matches = observed_counts == metadata_counts,
    nFeature_matches = observed_features == metadata_features,
    stringsAsFactors = FALSE
  )
  if (!all(spotcheck$nCount_matches & spotcheck$nFeature_matches)) {
    stop(
      "GSE255460 raw-count spot check disagrees with metadata QC metrics.",
      call. = FALSE
    )
  }

  list(
    member = member,
    n_cells = length(matrix_cells),
    n_genes = row_index,
    unique_genes = length(unique(genes)),
    header_mapping = "identical(make.names(metadata_cell_id), matrix_header)",
    header_mapping_ok = TRUE,
    integer_nonnegative_spotcheck = integer_nonnegative,
    structure_rows_checked = structure_rows,
    structure_field_counts = structure_fields,
    spotcheck = spotcheck
  )
}

adapt_gse255460 <- function(
    dataset,
    config,
    spotcheck_cells = 5L,
    stream_chunk_lines = 25L,
    structure_check_every = 5000L,
    lower_nmads = 3,
    upper_nmads = 5,
    mt_nmads = 3
) {
  .oa_sc_require_function("derive_sc_qc_thresholds")
  .oa_sc_require_function("apply_sc_qc_thresholds")
  .oa_sc_require_function("plot_sc_qc_metrics")

  dataset_id <- dataset$id %||% "GSE255460"
  metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
  members <- .oa_sc_tar_members(dataset$counts_archive)
  count_members <- members[basename(members) == "sc_counts.txt"]
  if (length(count_members) != 1L || length(members) != 1L) {
    stop(
      "GSE255460 counts archive must contain exactly one file named sc_counts.txt.",
      call. = FALSE
    )
  }

  metrics <- data.frame(
    dataset_id = dataset_id,
    cell_id = metadata$.cell_id,
    matrix_cell_id = make.names(metadata$.cell_id, unique = FALSE),
    batch = as.character(metadata$ID),
    nCount = as.numeric(metadata$nCount_RNA),
    nFeature = as.numeric(metadata$nFeature_RNA),
    stringsAsFactors = FALSE
  )
  if (
    anyNA(metrics$batch) ||
      any(!nzchar(metrics$batch)) ||
      anyNA(metrics$nCount) ||
      anyNA(metrics$nFeature)
  ) {
    stop("GSE255460 metadata contains missing batch/count/feature metrics.", call. = FALSE)
  }

  mt_column <- .oa_sc_first_existing_column(
    metadata,
    c("percent.mt", "percent_mt", "pct_counts_mt")
  )
  if (is.null(mt_column)) {
    metrics$percent_mt <- NA_real_
    mt_metric_available <- FALSE
    mt_metric_source <- "absent_from_local_metadata"
  } else {
    metrics$percent_mt <- suppressWarnings(as.numeric(metadata[[mt_column]]))
    mt_metric_available <- any(is.finite(metrics$percent_mt))
    mt_metric_source <- mt_column
  }

  metadata_columns <- intersect(
    c(
      "sample", "trait", "group", "new_ID", "celltype",
      "seurat_clusters", "orig.ident"
    ),
    names(metadata)
  )
  for (column in metadata_columns) {
    metrics[[column]] <- metadata[[column]]
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
  metrics <- apply_sc_qc_thresholds(
    metrics,
    thresholds,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt"
  )

  audit <- .oa_sc_scan_gse255460_counts(
    archive = dataset$counts_archive,
    member = count_members[[1L]],
    metadata = metadata,
    spotcheck_cells = spotcheck_cells,
    stream_chunk_lines = stream_chunk_lines,
    structure_check_every = structure_check_every
  )

  output_dir <- .oa_sc_output_dir(config, dataset_id)
  summary <- data.frame(
    dataset_id = dataset_id,
    data_type = "raw_integer_counts",
    matrix_orientation = "genes_by_cells",
    archive_members = length(members),
    matrix_member = audit$member,
    n_genes = audit$n_genes,
    n_cells = audit$n_cells,
    n_batches = length(unique(metrics$batch)),
    mt_metric_available = mt_metric_available,
    mt_metric_source = mt_metric_source,
    header_mapping_ok = audit$header_mapping_ok,
    raw_count_spotcheck_cells = nrow(audit$spotcheck),
    raw_count_spotcheck_ok = all(
      audit$spotcheck$nCount_matches & audit$spotcheck$nFeature_matches
    ),
    structure_rows_checked = length(audit$structure_rows_checked),
    cells_passing_qc = sum(metrics$passes_QC %in% TRUE, na.rm = TRUE),
    cells_failing_qc = sum(metrics$passes_QC %in% FALSE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  capability <- rbind(
    .oa_sc_capability_table(
      dataset_id,
      "raw_counts",
      "available",
      "Integer count spot checks match metadata nCount_RNA/nFeature_RNA."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "mad_qc_by_capture",
      "available",
      paste0(
        "Count/feature MAD QC is run by ID; mitochondrial QC is ",
        if (mt_metric_available) "available." else "not applied because percent.mt is absent."
      )
    ),
    .oa_sc_capability_table(
      dataset_id,
      "scDblFinder",
      "available",
      "Run on raw counts with samples=ID; metadata sample must not be used as capture batch."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "ambient_rna",
      "limited",
      "Cell-only methods are possible by ID, but SoupX is blocked without unfiltered droplets."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "dense_in_memory_import",
      "blocked",
      "The 10.5 GB ultra-wide text member must be converted to a sparse/on-disk representation."
    )
  )

  summary_path <- file.path(output_dir, "summary.csv")
  capability_path <- file.path(output_dir, "capabilities.csv")
  threshold_path <- file.path(output_dir, "qc_thresholds.csv")
  metrics_path <- file.path(output_dir, "cell_metrics.csv")
  spotcheck_path <- file.path(output_dir, "raw_count_spotcheck.csv")
  structure_path <- file.path(output_dir, "archive_structure_checks.csv")
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  safe_write_csv(thresholds, threshold_path)
  safe_write_csv(metrics, metrics_path)
  safe_write_csv(audit$spotcheck, spotcheck_path)
  safe_write_csv(
    data.frame(
      gene_row = audit$structure_rows_checked,
      observed_fields = audit$structure_field_counts,
      expected_fields = audit$n_cells + 1L,
      stringsAsFactors = FALSE
    ),
    structure_path
  )
  plot_paths <- plot_sc_qc_metrics(
    metrics,
    thresholds = thresholds,
    output_dir = output_dir,
    dataset_id = dataset_id,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt"
  )

  list(
    status = "ok",
    capability = capability,
    summary = summary,
    thresholds = thresholds,
    cell_metrics = metrics,
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = threshold_path,
      cell_metrics = metrics_path,
      raw_count_spotcheck = spotcheck_path,
      archive_structure = structure_path,
      plots = plot_paths
    )
  )
}

.oa_sc_read_gse104782_member <- function(archive, member) {
  command <- .oa_sc_tar_member_command(archive, member)
  connection <- gzcon(pipe(command, open = "rb"), text = TRUE)
  on.exit(
    if (!is.null(connection)) try(close(connection), silent = TRUE),
    add = TRUE
  )
  table <- utils::read.delim(
    connection,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  close(connection)
  connection <- NULL

  if (ncol(table) < 2L || nrow(table) == 0L) {
    stop("GSE104782 member is empty or malformed: ", member, call. = FALSE)
  }
  gene_column <- names(table)[[1L]]
  genes <- trimws(as.character(table[[gene_column]]))
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("GSE104782 member has invalid/duplicated genes: ", member, call. = FALSE)
  }
  expression_frame <- table[-1L]
  numeric_columns <- vapply(expression_frame, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop(
      "GSE104782 member contains nonnumeric expression columns: ", member,
      call. = FALSE
    )
  }
  expression <- as.matrix(expression_frame)
  storage.mode(expression) <- "double"
  if (anyNA(expression) || any(!is.finite(expression)) || any(expression < 0)) {
    stop("GSE104782 member contains invalid TPM values: ", member, call. = FALSE)
  }

  cells <- colnames(expression)
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
    stop("GSE104782 member has invalid/duplicated cell names: ", member, call. = FALSE)
  }
  mitochondrial <- grepl("^MT-", genes, ignore.case = TRUE)
  total_tpm <- colSums(expression)
  mt_tpm <- if (any(mitochondrial)) {
    colSums(expression[mitochondrial, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(expression))
  }
  percent_mt <- ifelse(
    is.finite(total_tpm) & total_tpm > 0 & is.finite(mt_tpm),
    100 * mt_tpm / total_tpm,
    NA_real_
  )
  basename_member <- basename(member)
  gsm <- sub("_.*$", "", basename_member)
  plate_match <- regexec("_(P[0-9]+)_", basename_member, perl = TRUE)
  plate_parts <- regmatches(basename_member, plate_match)[[1L]]
  plate <- if (length(plate_parts) >= 2L) plate_parts[[2L]] else NA_character_
  donor <- sub("_.*$", "", cells)
  state_match <- regexec("^[^_]+_(S[0-9]+)\\.", cells, perl = TRUE)
  state_parts <- regmatches(cells, state_match)
  state_code <- vapply(
    state_parts,
    function(parts) if (length(parts) >= 2L) parts[[2L]] else NA_character_,
    character(1)
  )

  metrics <- data.frame(
    dataset_id = "GSE104782",
    cell_id = cells,
    batch = gsm,
    plate = plate,
    archive_member = member,
    donor = donor,
    state_code = state_code,
    total_tpm = total_tpm,
    nFeature = colSums(expression > 0),
    percent_mt = percent_mt,
    passes_QC = rep(NA, ncol(expression)),
    failure_reason = rep("not_evaluated_tpm_only", ncol(expression)),
    stringsAsFactors = FALSE
  )
  fractional <- any(
    abs(expression - round(expression)) > sqrt(.Machine$double.eps),
    na.rm = TRUE
  )
  member_summary <- data.frame(
    member = member,
    batch = gsm,
    plate = plate,
    n_cells = ncol(expression),
    n_genes = nrow(expression),
    fractional_values = fractional,
    minimum_total_tpm = min(total_tpm),
    median_total_tpm = stats::median(total_tpm),
    maximum_total_tpm = max(total_tpm),
    nonzero_fraction = sum(expression > 0) / length(expression),
    stringsAsFactors = FALSE
  )
  list(
    genes = genes,
    metrics = metrics,
    member_summary = member_summary
  )
}

adapt_gse104782 <- function(
    dataset,
    config,
    tpm_sum_target = 1e6,
    tpm_relative_tolerance = 1e-3
) {
  .oa_sc_require_function("plot_sc_qc_metrics")

  dataset_id <- dataset$id %||% "GSE104782"
  members <- .oa_sc_tar_members(dataset$archive_path)
  tpm_members <- members[grepl("_TPM[.]txt[.]gz$", members, ignore.case = TRUE)]
  if (length(tpm_members) == 0L || length(tpm_members) != length(members)) {
    stop(
      "GSE104782 archive must contain only *_TPM.txt.gz members.",
      call. = FALSE
    )
  }

  base_genes <- NULL
  metric_parts <- vector("list", length(tpm_members))
  member_parts <- vector("list", length(tpm_members))
  for (index in seq_along(tpm_members)) {
    result <- .oa_sc_read_gse104782_member(
      dataset$archive_path,
      tpm_members[[index]]
    )
    if (is.null(base_genes)) {
      base_genes <- result$genes
    } else if (!identical(base_genes, result$genes)) {
      stop(
        "GSE104782 TPM members do not share identical gene names/order: ",
        tpm_members[[index]],
        call. = FALSE
      )
    }
    result$metrics$dataset_id <- dataset_id
    metric_parts[[index]] <- result$metrics
    member_parts[[index]] <- result$member_summary
  }
  metrics <- do.call(rbind, metric_parts)
  member_summary <- do.call(rbind, member_parts)
  rownames(metrics) <- NULL
  rownames(member_summary) <- NULL
  if (anyDuplicated(metrics$cell_id)) {
    stop("GSE104782 has duplicated cell names across archive members.", call. = FALSE)
  }

  relative_error <- abs(metrics$total_tpm - tpm_sum_target) / tpm_sum_target
  fractional_confirmed <- all(member_summary$fractional_values)
  tpm_sum_confirmed <- all(
    is.finite(relative_error) & relative_error <= tpm_relative_tolerance
  )
  if (!fractional_confirmed || !tpm_sum_confirmed) {
    stop(
      "GSE104782 failed TPM validation: fractional values and column sums near ",
      format(tpm_sum_target, scientific = FALSE), " are required.",
      call. = FALSE
    )
  }

  thresholds <- .oa_sc_empty_thresholds()
  output_dir <- .oa_sc_output_dir(config, dataset_id)
  summary <- data.frame(
    dataset_id = dataset_id,
    data_type = "TPM_normalized_not_raw_counts",
    matrix_orientation = "genes_by_cells",
    archive_members = length(tpm_members),
    n_genes = length(base_genes),
    n_cells = nrow(metrics),
    n_batches = length(unique(metrics$batch)),
    n_plates = length(unique(metrics$plate)),
    n_donors = length(unique(metrics$donor)),
    n_state_codes = length(unique(metrics$state_code)),
    fractional_values_confirmed = fractional_confirmed,
    tpm_sum_target = tpm_sum_target,
    tpm_sum_relative_tolerance = tpm_relative_tolerance,
    minimum_total_tpm = min(metrics$total_tpm),
    median_total_tpm = stats::median(metrics$total_tpm),
    maximum_total_tpm = max(metrics$total_tpm),
    passes_qc_available = FALSE,
    stringsAsFactors = FALSE
  )
  capability <- rbind(
    .oa_sc_capability_table(
      dataset_id,
      "exploratory_TPM_expression",
      "available",
      "TPM matrices can be combined after gene-order validation."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "raw_counts",
      "blocked",
      "All local matrices contain fractional TPM values with column sums near one million."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "MAD_count_QC",
      "blocked",
      "Count-based thresholds are not valid for TPM; passes_QC is intentionally NA."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "scDblFinder",
      "blocked",
      "Raw counts are absent and each plate contains only 32-48 cells."
    ),
    .oa_sc_capability_table(
      dataset_id,
      "ambient_rna",
      "blocked",
      "Raw counts, unfiltered droplets/empty wells, and an ambient profile are absent."
    )
  )

  summary_path <- file.path(output_dir, "summary.csv")
  capability_path <- file.path(output_dir, "capabilities.csv")
  threshold_path <- file.path(output_dir, "qc_thresholds.csv")
  metrics_path <- file.path(output_dir, "cell_metrics.csv")
  member_path <- file.path(output_dir, "archive_member_summary.csv")
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  safe_write_csv(thresholds, threshold_path)
  safe_write_csv(metrics, metrics_path)
  safe_write_csv(member_summary, member_path)
  plot_paths <- plot_sc_qc_metrics(
    metrics,
    thresholds = NULL,
    output_dir = output_dir,
    dataset_id = dataset_id,
    batch_col = "batch",
    count_col = "total_tpm",
    feature_col = "nFeature",
    mt_col = "percent_mt"
  )

  list(
    status = "ok",
    capability = capability,
    summary = summary,
    thresholds = thresholds,
    cell_metrics = metrics,
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = threshold_path,
      cell_metrics = metrics_path,
      archive_members = member_path,
      plots = plot_paths
    )
  )
}

run_oa_sc_adapter <- function(dataset, config) {
  dataset_id <- dataset$id %||% ""
  switch(
    dataset_id,
    GSE255460 = adapt_gse255460(dataset, config),
    GSE104782 = adapt_gse104782(dataset, config),
    stop("No OA single-cell adapter is implemented for: ", dataset_id, call. = FALSE)
  )
}
