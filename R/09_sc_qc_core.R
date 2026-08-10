sc_metric_bounds <- function(
    values,
    lower_nmads = 3,
    upper_nmads = 5,
    percentage = FALSE
) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(list(
      median = NA_real_,
      mad = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    ))
  }

  centre <- stats::median(values)
  spread <- stats::mad(values, center = centre, constant = 1.4826)
  quantiles <- as.numeric(stats::quantile(
    values,
    probs = c(0.01, 0.99),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))
  if (!is.finite(spread) || spread <= sqrt(.Machine$double.eps)) {
    lower <- quantiles[[1L]]
    upper <- quantiles[[2L]]
  } else {
    lower <- centre - lower_nmads * spread
    upper <- centre + upper_nmads * spread
  }
  lower <- max(0, lower)
  if (percentage) {
    upper <- min(100, max(0, upper))
  }
  if (!is.finite(upper) || upper <= lower) {
    upper <- max(values)
  }

  list(
    median = centre,
    mad = spread,
    lower = lower,
    upper = upper
  )
}

# This definition intentionally replaces the slower pure-R implementation
# sourced from 09_sc_adapters_oa.R. It validates every matrix row while keeping
# memory bounded, then returns the same adapter contract.
.oa_sc_scan_gse255460_counts <- function(
    archive,
    member,
    metadata,
    spotcheck_cells = 5L,
    stream_chunk_lines = 25L,
    structure_check_every = 5000L
) {
  require_namespace("digest", "GSE255460 header validation")
  require_namespace("jsonlite", "GSE255460 streaming audit")
  spotcheck_cells <- max(1L, min(as.integer(spotcheck_cells), nrow(metadata)))
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  script <- file.path(project_root, "tools", "scan_gse255460.py")
  if (!file.exists(script)) {
    stop("GSE255460 scanner is missing: ", script, call. = FALSE)
  }
  python_candidates <- unique(c(
    Sys.getenv("OC_OA_PYTHON", ""),
    Sys.which("python3"),
    Sys.which("python")
  ))
  python_candidates <- python_candidates[
    nzchar(python_candidates) & file.exists(python_candidates)
  ]
  if (length(python_candidates) == 0L) {
    stop(
      "GSE255460 requires a Python executable for bounded-memory streaming. ",
      "Set OC_OA_PYTHON to Python 2.7+ or Python 3.",
      call. = FALSE
    )
  }
  python <- python_candidates[[1L]]
  output <- tempfile("gse255460-audit-", fileext = ".json")
  on.exit(unlink(output, force = TRUE), add = TRUE)
  status <- system2(
    python,
    args = c(
      shQuote(script),
      "--archive", shQuote(archive),
      "--member", shQuote(member),
      "--output", shQuote(output),
      "--spotcheck-cells", as.character(spotcheck_cells)
    ),
    stdout = "",
    stderr = ""
  )
  if (!identical(as.integer(status), 0L) || !file.exists(output)) {
    stop(
      "The GSE255460 bounded-memory scanner failed with exit status ",
      status, ".",
      call. = FALSE
    )
  }
  audit <- jsonlite::read_json(output, simplifyVector = TRUE)
  expected_cells <- make.names(metadata$.cell_id, unique = FALSE)
  expected_header <- paste(expected_cells, collapse = "\t")
  expected_md5 <- digest::digest(
    expected_header,
    algo = "md5",
    serialize = FALSE
  )
  if (
    !identical(as.numeric(audit$n_cells), as.numeric(length(expected_cells))) ||
      !identical(tolower(audit$header_md5), tolower(expected_md5))
  ) {
    stop(
      "GSE255460 matrix header does not exactly equal ",
      "make.names(metadata cell IDs) in order.",
      call. = FALSE
    )
  }
  if (as.numeric(audit$structure_field_mismatches) != 0) {
    stop(
      "GSE255460 has a malformed field count at gene row ",
      audit$first_structure_mismatch_row, ".",
      call. = FALSE
    )
  }

  selected <- seq_len(spotcheck_cells)
  observed_counts <- as.numeric(audit$observed_counts)
  observed_features <- as.integer(audit$observed_features)
  metadata_counts <- as.numeric(metadata$nCount_RNA[selected])
  metadata_features <- as.integer(metadata$nFeature_RNA[selected])
  spotcheck <- data.frame(
    cell_id = metadata$.cell_id[selected],
    matrix_cell_id = expected_cells[selected],
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
  n_genes <- as.integer(audit$n_genes)
  list(
    member = member,
    n_cells = as.integer(audit$n_cells),
    n_genes = n_genes,
    unique_genes = as.integer(audit$unique_genes),
    header_mapping = "identical(make.names(metadata_cell_id), matrix_header)",
    header_mapping_ok = TRUE,
    integer_nonnegative_spotcheck = isTRUE(
      audit$integer_nonnegative_spotcheck
    ),
    structure_rows_checked = seq_len(n_genes),
    structure_field_counts = rep(as.integer(audit$expected_fields), n_genes),
    spotcheck = spotcheck,
    scanner = basename(script),
    scanner_elapsed_seconds = as.numeric(audit$elapsed_seconds),
    full_structure_check = TRUE
  )
}

derive_sc_qc_thresholds <- function(
    metrics,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt",
    lower_nmads = 3,
    upper_nmads = 5,
    mt_nmads = 3
) {
  required <- c(batch_col, count_col, feature_col)
  missing <- setdiff(required, names(metrics))
  if (length(missing) > 0L) {
    stop(
      "QC metrics are missing required columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!mt_col %in% names(metrics)) {
    metrics[[mt_col]] <- NA_real_
  }

  batch <- as.character(metrics[[batch_col]])
  batch[is.na(batch) | !nzchar(batch)] <- "__missing__"
  batches <- unique(batch)
  rows <- lapply(batches, function(value) {
    keep <- batch == value
    count <- sc_metric_bounds(
      metrics[[count_col]][keep],
      lower_nmads = lower_nmads,
      upper_nmads = upper_nmads
    )
    feature <- sc_metric_bounds(
      metrics[[feature_col]][keep],
      lower_nmads = lower_nmads,
      upper_nmads = upper_nmads
    )
    mt <- sc_metric_bounds(
      metrics[[mt_col]][keep],
      lower_nmads = 0,
      upper_nmads = mt_nmads,
      percentage = TRUE
    )
    data.frame(
      batch = value,
      n_cells = sum(keep),
      nCount_median = count$median,
      nCount_mad = count$mad,
      nCount_lower = count$lower,
      nCount_upper = count$upper,
      nFeature_median = feature$median,
      nFeature_mad = feature$mad,
      nFeature_lower = feature$lower,
      nFeature_upper = feature$upper,
      percent_mt_median = mt$median,
      percent_mt_mad = mt$mad,
      percent_mt_upper = mt$upper,
      lower_nmads = lower_nmads,
      upper_nmads = upper_nmads,
      mt_nmads = mt_nmads,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

append_qc_reason <- function(existing, failed, reason) {
  existing <- as.character(existing)
  existing[is.na(existing)] <- ""
  indices <- which(failed %in% TRUE)
  if (length(indices) > 0L) {
    existing[indices] <- ifelse(
      nzchar(existing[indices]),
      paste(existing[indices], reason, sep = ";"),
      reason
    )
  }
  existing
}

apply_sc_qc_thresholds <- function(
    metrics,
    thresholds,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt"
) {
  if (!mt_col %in% names(metrics)) {
    metrics[[mt_col]] <- NA_real_
  }
  batch <- as.character(metrics[[batch_col]])
  batch[is.na(batch) | !nzchar(batch)] <- "__missing__"
  index <- match(batch, thresholds$batch)
  if (anyNA(index)) {
    stop("At least one QC batch has no threshold row.", call. = FALSE)
  }

  count <- as.numeric(metrics[[count_col]])
  feature <- as.numeric(metrics[[feature_col]])
  mt <- as.numeric(metrics[[mt_col]])
  count_low <- thresholds$nCount_lower[index]
  count_high <- thresholds$nCount_upper[index]
  feature_low <- thresholds$nFeature_lower[index]
  feature_high <- thresholds$nFeature_upper[index]
  mt_high <- thresholds$percent_mt_upper[index]

  fail_count_low <- !is.finite(count) | count < count_low
  fail_count_high <- is.finite(count_high) & count > count_high
  fail_feature_low <- !is.finite(feature) | feature < feature_low
  fail_feature_high <- is.finite(feature_high) & feature > feature_high
  fail_mt <- is.finite(mt_high) & is.finite(mt) & mt > mt_high

  reasons <- rep("", nrow(metrics))
  reasons <- append_qc_reason(reasons, fail_count_low, "low_counts")
  reasons <- append_qc_reason(reasons, fail_count_high, "high_counts")
  reasons <- append_qc_reason(reasons, fail_feature_low, "low_features")
  reasons <- append_qc_reason(reasons, fail_feature_high, "high_features")
  reasons <- append_qc_reason(reasons, fail_mt, "high_mitochondrial_fraction")
  metrics$passes_core_QC <- !(
    fail_count_low | fail_count_high | fail_feature_low |
      fail_feature_high | fail_mt
  )
  metrics$passes_QC <- metrics$passes_core_QC
  metrics$failure_reason <- reasons
  metrics
}

plot_sc_qc_metrics <- function(
    metrics,
    thresholds = NULL,
    output_dir,
    dataset_id,
    batch_col = "batch",
    count_col = "nCount",
    feature_col = "nFeature",
    mt_col = "percent_mt"
) {
  require_namespace("ggplot2", "single-cell QC figures")
  figure_dir <- ensure_dir(file.path(output_dir, "figures"))
  paths <- character()
  if (!mt_col %in% names(metrics)) {
    metrics[[mt_col]] <- NA_real_
  }

  max_points <- 100000L
  if (nrow(metrics) > max_points) {
    evenly_spaced <- unique(round(seq(1, nrow(metrics), length.out = max_points)))
    plot_data <- metrics[evenly_spaced, , drop = FALSE]
  } else {
    plot_data <- metrics
  }
  pass_label <- if ("passes_QC" %in% names(plot_data)) {
    ifelse(
      is.na(plot_data$passes_QC),
      "not_assessed",
      ifelse(plot_data$passes_QC, "pass", "fail")
    )
  } else {
    rep("not_assessed", nrow(plot_data))
  }
  scatter <- data.frame(
    nCount = as.numeric(plot_data[[count_col]]),
    nFeature = as.numeric(plot_data[[feature_col]]),
    pass = pass_label,
    stringsAsFactors = FALSE
  )
  scatter <- scatter[
    is.finite(scatter$nCount) & is.finite(scatter$nFeature),
    ,
    drop = FALSE
  ]
  if (nrow(scatter) > 0L) {
    scatter_path <- file.path(
      figure_dir,
      paste0(clean_filename(dataset_id), "_count_feature_scatter.png")
    )
    plot <- ggplot2::ggplot(
      scatter,
      ggplot2::aes(x = nCount, y = nFeature, colour = pass)
    ) +
      ggplot2::geom_point(alpha = 0.22, size = 0.35) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c(fail = "#D73027", pass = "#1A9850", not_assessed = "#4575B4")
      ) +
      ggplot2::labs(
        title = paste(dataset_id, "single-cell QC"),
        x = "Total counts (log10)",
        y = "Detected features (log10)",
        colour = "QC"
      ) +
      ggplot2::theme_bw(base_size = 10)
    ggplot2::ggsave(scatter_path, plot, width = 7, height = 5, dpi = 160)
    paths <- c(paths, scatter_path)
  }

  mt <- data.frame(
    percent_mt = as.numeric(metrics[[mt_col]]),
    stringsAsFactors = FALSE
  )
  mt <- mt[is.finite(mt$percent_mt), , drop = FALSE]
  if (nrow(mt) > 0L) {
    mt_path <- file.path(
      figure_dir,
      paste0(clean_filename(dataset_id), "_mitochondrial_distribution.png")
    )
    plot <- ggplot2::ggplot(mt, ggplot2::aes(x = percent_mt)) +
      ggplot2::geom_histogram(bins = 80, fill = "#4575B4", colour = "white") +
      ggplot2::labs(
        title = paste(dataset_id, "mitochondrial fraction"),
        x = "Mitochondrial counts (%)",
        y = "Cells"
      ) +
      ggplot2::theme_bw(base_size = 10)
    ggplot2::ggsave(mt_path, plot, width = 7, height = 4.5, dpi = 160)
    paths <- c(paths, mt_path)
  }

  if (
    batch_col %in% names(metrics) &&
      length(unique(metrics[[batch_col]])) > 1L
  ) {
    batch_values <- as.character(metrics[[batch_col]])
    split_indices <- split(seq_len(nrow(metrics)), batch_values)
    batch_summary <- do.call(rbind, lapply(names(split_indices), function(batch) {
      indices <- split_indices[[batch]]
      data.frame(
        batch = batch,
        cells = length(indices),
        median_counts = stats::median(
          as.numeric(metrics[[count_col]][indices]),
          na.rm = TRUE
        ),
        median_features = stats::median(
          as.numeric(metrics[[feature_col]][indices]),
          na.rm = TRUE
        ),
        median_percent_mt = stats::median(
          as.numeric(metrics[[mt_col]][indices]),
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }))
    batch_path <- file.path(
      figure_dir,
      paste0(clean_filename(dataset_id), "_batch_medians.png")
    )
    plot <- ggplot2::ggplot(
      batch_summary,
      ggplot2::aes(x = median_counts, y = median_features, size = cells)
    ) +
      ggplot2::geom_point(alpha = 0.7, colour = "#313695") +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = paste(dataset_id, "batch-level QC medians"),
        x = "Median counts (log10)",
        y = "Median features (log10)"
      ) +
      ggplot2::theme_bw(base_size = 10)
    ggplot2::ggsave(batch_path, plot, width = 7, height = 5, dpi = 160)
    paths <- c(paths, batch_path)
  }
  unname(paths)
}

write_sc_table <- function(x, path) {
  ensure_dir(dirname(path))
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(
      x,
      file = path,
      sep = "\t",
      quote = FALSE,
      na = "",
      compress = if (endsWith(tolower(path), ".gz")) "gzip" else "none"
    )
  } else {
    connection <- if (endsWith(tolower(path), ".gz")) {
      gzfile(path, open = "wt")
    } else {
      file(path, open = "wt")
    }
    on.exit(close(connection), add = TRUE)
    utils::write.table(
      x,
      file = connection,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      na = ""
    )
  }
  invisible(path)
}

safe_archive_members <- function(archive) {
  members <- utils::untar(archive, list = TRUE)
  unsafe <- grepl(
    "(^|[/\\\\])\\.\\.([/\\\\]|$)|^[/\\\\]|^[A-Za-z]:",
    members
  )
  if (any(unsafe)) {
    stop(
      "Archive contains unsafe member paths: ",
      paste(head(members[unsafe], 5L), collapse = ", "),
      call. = FALSE
    )
  }
  members
}

extract_sc_archive <- function(archive, destination, expected_members = NULL) {
  ensure_dir(destination)
  members <- safe_archive_members(archive)
  selected <- expected_members %||% members
  if (!all(selected %in% members)) {
    stop("Archive is missing one or more expected members.", call. = FALSE)
  }
  targets <- file.path(destination, selected)
  if (!all(file.exists(targets))) {
    utils::untar(archive, files = selected, exdir = destination)
  }
  if (!all(file.exists(targets))) {
    stop("Archive extraction did not produce every expected file.", call. = FALSE)
  }
  normalizePath(targets, winslash = "/", mustWork = TRUE)
}

read_gzip_lines <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  readLines(connection, warn = FALSE)
}

read_gse154600_sample <- function(paths, sample_id) {
  require_namespace("Matrix", "10x sparse matrix import")
  require_namespace("SingleCellExperiment", "single-cell object storage")
  require_namespace("S4Vectors", "single-cell metadata storage")

  gene_path <- paths[grepl("_genes\\.tsv\\.gz$", paths)]
  barcode_path <- paths[grepl("_barcodes\\.tsv\\.gz$", paths)]
  matrix_path <- paths[grepl("_matrix\\.mtx\\.gz$", paths)]
  if (length(gene_path) != 1L || length(barcode_path) != 1L ||
      length(matrix_path) != 1L) {
    stop("Incomplete 10x triplet for ", sample_id, ".", call. = FALSE)
  }

  genes <- data.table::fread(
    gene_path,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )
  barcodes <- data.table::fread(
    barcode_path,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )[[1L]]
  connection <- gzfile(matrix_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  counts <- methods::as(Matrix::readMM(connection), "CsparseMatrix")
  if (nrow(counts) != nrow(genes) || ncol(counts) != length(barcodes)) {
    stop(
      sample_id, " 10x dimensions do not match genes/barcodes.",
      call. = FALSE
    )
  }
  if (length(counts@x) > 0L &&
      (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8))) {
    stop(sample_id, " matrix is not non-negative integer raw counts.", call. = FALSE)
  }

  gene_id <- as.character(genes[[1L]])
  gene_symbol <- if (ncol(genes) >= 2L) {
    as.character(genes[[2L]])
  } else {
    gene_id
  }
  rownames(counts) <- make.unique(gene_id)
  cell_id <- paste(sample_id, as.character(barcodes), sep = "_")
  colnames(counts) <- cell_id
  row_data <- S4Vectors::DataFrame(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    row.names = rownames(counts)
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )
  list(
    sce = sce,
    sample_id = sample_id,
    gene_symbol = gene_symbol,
    cell_id = cell_id
  )
}

audit_gse154600 <- function(dataset, config) {
  require_namespace("Matrix", "single-cell sparse matrices")
  require_namespace("SingleCellExperiment", "single-cell object storage")
  require_namespace("S4Vectors", "single-cell metadata")
  require_namespace("scDblFinder", "mandatory doublet detection")
  require_namespace("BiocParallel", "serial doublet detection")
  require_namespace("data.table", "single-cell input tables")

  dataset_id <- dataset$id %||% "GSE154600"
  output_dir <- ensure_dir(file.path(
    config$project$output_dir,
    "single_cell",
    dataset_id
  ))
  cache_dir <- ensure_dir(file.path(
    config$project$cache_dir,
    "single_cell",
    dataset_id,
    "extracted"
  ))
  archive <- dataset$archive_path
  members <- safe_archive_members(archive)
  expected_suffixes <- c(
    "_barcodes.tsv.gz", "_genes.tsv.gz", "_matrix.mtx.gz"
  )
  sample_ids <- unique(sub(
    "(_barcodes\\.tsv\\.gz|_genes\\.tsv\\.gz|_matrix\\.mtx\\.gz)$",
    "",
    members
  ))
  triplet_counts <- vapply(
    sample_ids,
    function(sample_id) sum(startsWith(members, paste0(sample_id, "_"))),
    integer(1)
  )
  if (
    length(sample_ids) == 0L ||
      any(triplet_counts != length(expected_suffixes)) ||
      any(!vapply(
        sample_ids,
        function(sample_id) {
          all(paste0(sample_id, expected_suffixes) %in% members)
        },
        logical(1)
      ))
  ) {
    stop("GSE154600 archive is not a complete set of 10x triplets.", call. = FALSE)
  }
  extracted <- extract_sc_archive(archive, cache_dir, members)

  sample_summaries <- list()
  all_thresholds <- list()
  all_metric_paths <- character()
  object_paths <- character()
  for (index in seq_along(sample_ids)) {
    sample_id <- sample_ids[[index]]
    log_info(
      "Single-cell ", dataset_id, ": importing ", sample_id,
      " (", index, "/", length(sample_ids), ")."
    )
    paths <- extracted[startsWith(basename(extracted), paste0(sample_id, "_"))]
    imported <- read_gse154600_sample(paths, sample_id)
    sce <- imported$sce
    counts <- SummarizedExperiment::assay(sce, "counts")
    total_counts <- Matrix::colSums(counts)
    detected <- Matrix::colSums(counts > 0)
    mt <- grepl("^MT-", imported$gene_symbol, ignore.case = TRUE)
    mt_counts <- if (any(mt)) {
      Matrix::colSums(counts[mt, , drop = FALSE])
    } else {
      rep(NA_real_, ncol(counts))
    }
    percent_mt <- ifelse(
      is.finite(mt_counts) & total_counts > 0,
      100 * mt_counts / total_counts,
      NA_real_
    )
    metrics <- data.frame(
      cell_id = imported$cell_id,
      dataset_id = dataset_id,
      batch = sample_id,
      nCount = as.numeric(total_counts),
      nFeature = as.numeric(detected),
      percent_mt = as.numeric(percent_mt),
      doublet_score = NA_real_,
      doublet_class = NA_character_,
      stringsAsFactors = FALSE
    )
    thresholds <- derive_sc_qc_thresholds(
      metrics,
      lower_nmads = config$single_cell$lower_nmads %||% 3,
      upper_nmads = config$single_cell$upper_nmads %||% 5,
      mt_nmads = config$single_cell$mt_nmads %||% 3
    )
    metrics <- apply_sc_qc_thresholds(metrics, thresholds)
    core_pass <- which(metrics$passes_core_QC)
    if (length(core_pass) < 100L) {
      stop(
        sample_id,
        " has fewer than 100 cells after core QC; scDblFinder is not valid.",
        call. = FALSE
      )
    }

    set.seed(as.integer(config$project$seed) + index)
    doublet_sce <- scDblFinder::scDblFinder(
      sce[, core_pass, drop = FALSE],
      BPPARAM = BiocParallel::SerialParam(),
      verbose = TRUE
    )
    doublet_data <- as.data.frame(SummarizedExperiment::colData(doublet_sce))
    if (!all(c("scDblFinder.score", "scDblFinder.class") %in% names(doublet_data))) {
      stop(
        "scDblFinder did not return its mandatory score/class columns for ",
        sample_id, ".",
        call. = FALSE
      )
    }
    metrics$doublet_score[core_pass] <- doublet_data$scDblFinder.score
    metrics$doublet_class[core_pass] <- as.character(
      doublet_data$scDblFinder.class
    )
    predicted_doublet <- metrics$passes_core_QC &
      metrics$doublet_class != "singlet"
    predicted_doublet[is.na(predicted_doublet)] <- FALSE
    metrics$failure_reason <- append_qc_reason(
      metrics$failure_reason,
      predicted_doublet,
      "predicted_doublet"
    )
    metrics$passes_QC <- metrics$passes_core_QC &
      !is.na(metrics$doublet_class) &
      metrics$doublet_class == "singlet"

    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
      metrics,
      row.names = metrics$cell_id
    )
    object_path <- file.path(
      config$project$cache_dir,
      "single_cell",
      dataset_id,
      paste0(sample_id, "_qc_sce.rds")
    )
    atomic_save_rds(sce, object_path, compress = FALSE)
    metric_path <- file.path(
      output_dir,
      paste0(sample_id, "_cell_qc.tsv.gz")
    )
    write_sc_table(metrics, metric_path)
    plot_sc_qc_metrics(
      metrics,
      thresholds,
      output_dir = output_dir,
      dataset_id = paste(dataset_id, sample_id, sep = "_")
    )
    thresholds$dataset_id <- dataset_id
    thresholds$sample_id <- sample_id
    all_thresholds[[sample_id]] <- thresholds
    all_metric_paths <- c(all_metric_paths, metric_path)
    object_paths <- c(object_paths, object_path)
    sample_summaries[[sample_id]] <- data.frame(
      dataset_id = dataset_id,
      sample_id = sample_id,
      genes = nrow(sce),
      cells = ncol(sce),
      core_qc_pass = sum(metrics$passes_core_QC, na.rm = TRUE),
      predicted_doublets = sum(
        metrics$doublet_class == "doublet",
        na.rm = TRUE
      ),
      final_qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
      median_counts = stats::median(metrics$nCount),
      median_features = stats::median(metrics$nFeature),
      median_percent_mt = stats::median(metrics$percent_mt, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    rm(sce, doublet_sce, doublet_data, counts, metrics)
    invisible(gc())
  }

  summary <- do.call(rbind, sample_summaries)
  thresholds <- do.call(rbind, all_thresholds)
  safe_write_csv(summary, file.path(output_dir, "sample_qc_summary.csv"))
  safe_write_csv(thresholds, file.path(output_dir, "qc_thresholds.csv"))
  capability <- data.frame(
    dataset_id = dataset_id,
    capability = c(
      "adapter_validation", "raw_counts", "core_qc", "scDblFinder",
      "ambient_rna", "downstream_qc_object"
    ),
    status = c(
      "passed", "available", "completed", "completed",
      "blocked", "available"
    ),
    reason = c(
      "All five per-sample 10x triplets passed dimension and integer-count checks.",
      "Non-negative integer UMI counts are preserved in the counts assay.",
      "Per-sample median/MAD thresholds were applied without discarding source cells.",
      "scDblFinder completed independently for every technical sample.",
      "No unfiltered/empty droplets are present, so ambient-RNA correction is unsupported.",
      "Each cached SingleCellExperiment preserves counts and passes_QC metadata."
    ),
    stringsAsFactors = FALSE
  )
  list(
    id = dataset_id,
    status = "validated",
    reason = "",
    summary = summary,
    thresholds = thresholds,
    capability = capability,
    cell_metric_paths = all_metric_paths,
    object_paths = object_paths,
    cells = sum(summary$cells),
    qc_pass = sum(summary$final_qc_pass),
    downstream_ready = TRUE
  )
}

audit_gse255460 <- function(dataset, config) {
  adapt_gse255460_partitioned_csr(dataset, config)
}

audit_gse104782 <- function(dataset, config) {
  adapt_gse104782_umi_counts(dataset, config)
}

audit_gse169454 <- function(dataset, config) {
  adapt_gse169454_tenx_raw_tar(dataset, config)
}

audit_gse180661 <- function(dataset, config) {
  adapt_gse180661_h5_csr(dataset, config)
}

single_cell_error_result <- function(dataset_id, adapter, error) {
  list(
    id = dataset_id,
    status = "blocked_adapter_error",
    reason = compact_error(error),
    summary = data.frame(),
    thresholds = data.frame(),
    capability = data.frame(
      dataset_id = dataset_id,
      capability = "adapter_validation",
      status = "failed",
      reason = compact_error(error),
      stringsAsFactors = FALSE
    ),
    cells = NA_real_,
    qc_pass = NA_real_,
    downstream_ready = FALSE
  )
}

normalize_single_cell_result <- function(result, dataset, adapter) {
  result$id <- result$id %||% dataset$id
  result$status <- result$status %||% "blocked_unspecified"
  result$reason <- result$reason %||% ""
  result$summary <- result$summary %||% data.frame()
  result$thresholds <- result$thresholds %||% data.frame()
  result$cells <- as.numeric(result$cells %||% NA_real_)
  result$qc_pass <- as.numeric(result$qc_pass %||% NA_real_)
  result$downstream_ready <- isTRUE(result$downstream_ready)
  if (
    !is.null(result$cell_metrics) &&
      (is.data.frame(result$cell_metrics) ||
        data.table::is.data.table(result$cell_metrics))
  ) {
    result$cell_metrics_rows <- nrow(result$cell_metrics)
    result$cell_metrics <- NULL
  }
  if (is.null(result$capability) || nrow(result$capability) == 0L) {
    result$capability <- data.frame(
      dataset_id = result$id,
      capability = "adapter_validation",
      status = result$status,
      reason = result$reason,
      stringsAsFactors = FALSE
    )
  }
  required_capability_columns <- c(
    "dataset_id", "capability", "status", "reason"
  )
  missing_capability_columns <- setdiff(
    required_capability_columns,
    names(result$capability)
  )
  if (length(missing_capability_columns) > 0L) {
    stop(
      "Adapter ", result$id, " returned an invalid capability table: ",
      paste(missing_capability_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  result$capability <- result$capability[
    ,
    required_capability_columns,
    drop = FALSE
  ]
  result$capability$adapter <- adapter
  result$capability <- result$capability[
    ,
    c("dataset_id", "adapter", "capability", "status", "reason"),
    drop = FALSE
  ]
  result
}

write_single_cell_gate_report <- function(results, status_table, config) {
  report_path <- file.path(
    config$project$output_dir,
    "reports",
    "single_cell_qc_gate_report.md"
  )
  lines <- c(
    "# Single-cell adapter and QC gate report",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- Scope: import validation, cell-level QC, and doublet/ambient capability gating.",
    "- Cell annotation, integration, UMAP interpretation, and differential expression are outside this QC-only stage.",
    "",
    "## Dataset status",
    "",
    "| Dataset | Adapter status | Cells audited | Final QC pass | Downstream ready | Reason |",
    "|---|---|---:|---:|---|---|"
  )
  for (index in seq_len(nrow(status_table))) {
    row <- status_table[index, , drop = FALSE]
    reason <- gsub("\\|", "/", row$reason)
    lines <- c(lines, sprintf(
      "| %s | %s | %s | %s | %s | %s |",
      row$dataset_id,
      row$status,
      ifelse(is.na(row$cells), "NA", format(row$cells, big.mark = ",")),
      ifelse(is.na(row$qc_pass), "NA", format(row$qc_pass, big.mark = ",")),
      ifelse(row$downstream_ready, "yes", "no"),
      reason
    ))
  }
  lines <- c(
    lines,
    "",
    "## Gate interpretation",
    "",
    "- `validated` means the local source structure passed its adapter checks and mandatory analyses supported by that source completed.",
    "- A dataset with normalized-only values, a malformed bundle, or a missing count matrix remains blocked even when descriptive metadata QC is available.",
    "- Ambient-RNA correction is never claimed when unfiltered droplets or empty droplets are absent.",
    "- The module is enabled as a QC gate. Downstream single-cell inference remains closed unless every selected dataset required for that inference is validated.",
    ""
  )
  write_utf8(lines, report_path)
  report_path
}

run_single_cell_stage <- function(config) {
  if (!isTRUE(config$modules$single_cell)) {
    return(list(
      status = "disabled",
      reason = "The single-cell QC gate is disabled by configuration."
    ))
  }

  python_runtime <- find_single_cell_python(config)
  if (!identical(python_runtime$status[[1L]], "available")) {
    stop(
      "The single-cell QC gate requires Python >= 3.8 for GSE255460.",
      call. = FALSE
    )
  }
  Sys.setenv(OC_OA_PYTHON = python_runtime$executable[[1L]])
  datasets <- config$single_cell$datasets %||% list()
  adapter_functions <- c(
    gse255460_wide_counts = "audit_gse255460",
    gse104782_umi_counts = "audit_gse104782",
    gse169454_tenx_raw_tar = "audit_gse169454",
    gse180661_h5_csr = "audit_gse180661",
    tenx_tar_gse154600 = "audit_gse154600"
  )
  results <- list()
  for (dataset in datasets) {
    dataset_id <- dataset$id
    adapter <- dataset$adapter
    function_name <- unname(adapter_functions[[adapter]])
    log_info(
      "Single-cell QC gate: ", dataset_id,
      " using adapter ", adapter, "."
    )
    if (is.null(function_name) || !exists(function_name, mode = "function")) {
      result <- single_cell_error_result(
        dataset_id,
        adapter,
        simpleError(paste0("No registered adapter function for ", adapter, "."))
      )
    } else {
      result <- tryCatch(
        get(function_name, mode = "function")(dataset, config),
        error = function(error) {
          log_error(
            "Single-cell adapter ", dataset_id, " blocked: ",
            conditionMessage(error)
          )
          single_cell_error_result(dataset_id, adapter, error)
        }
      )
    }
    results[[dataset_id]] <- normalize_single_cell_result(
      result,
      dataset,
      adapter
    )
  }

  status_table <- do.call(rbind, lapply(results, function(result) {
    data.frame(
      dataset_id = result$id,
      status = result$status,
      cells = result$cells,
      qc_pass = result$qc_pass,
      downstream_ready = result$downstream_ready,
      reason = result$reason,
      stringsAsFactors = FALSE
    )
  }))
  capability <- do.call(rbind, lapply(results, `[[`, "capability"))
  safe_write_csv(
    status_table,
    file.path(
      config$project$output_dir,
      "tables",
      "single_cell_dataset_status.csv"
    )
  )
  safe_write_csv(
    capability,
    file.path(
      config$project$output_dir,
      "tables",
      "single_cell_capability_matrix.csv"
    )
  )
  report_path <- write_single_cell_gate_report(
    results,
    status_table,
    config
  )
  yaml::write_yaml(
    list(
      seed = config$project$seed,
      scope = "qc_gate_only",
      threshold_method = "per-batch median/MAD with quantile fallback",
      lower_nmads = config$single_cell$lower_nmads %||% 3,
      upper_nmads = config$single_cell$upper_nmads %||% 5,
      mitochondrial_nmads = config$single_cell$mt_nmads %||% 3,
      doublet_method = "scDblFinder per technical batch when raw counts exist",
      ambient_requirement = "unfiltered droplets or empty-droplet profile",
      gse255460_scanner = "tools/scan_gse255460.py",
      gse255460_sparse_converter = "tools/convert_gse255460_sparse.py",
      gse255460_sparse_format = "partitioned little-endian int32 CSR",
      python_executable = python_runtime$executable[[1L]],
      python_version = python_runtime$version[[1L]],
      numpy_version = python_runtime$numpy_version[[1L]] %||% NA_character_,
      allow_partial = isTRUE(config$single_cell$allow_partial)
    ),
    file.path(
      config$project$output_dir,
      "manifests",
      "single_cell_qc_parameters.yml"
    )
  )

  ready <- status_table$downstream_ready
  minimum_ready <- as.integer(config$single_cell$minimum_validated_datasets %||% 1L)
  stage_status <- if (all(ready)) {
    "validated"
  } else if (sum(ready) >= minimum_ready) {
    "partial"
  } else {
    "blocked"
  }
  if (!all(ready) && !isTRUE(config$single_cell$allow_partial)) {
    stop(
      "Single-cell QC gate did not validate every configured dataset. See ",
      report_path, ".",
      call. = FALSE
    )
  }
  list(
    status = stage_status,
    scope = "qc_gate_only",
    downstream_enabled = all(ready),
    validated_datasets = status_table$dataset_id[ready],
    blocked_datasets = status_table$dataset_id[!ready],
    dataset_status = status_table,
    capability = capability,
    report = report_path,
    datasets = results
  )
}
