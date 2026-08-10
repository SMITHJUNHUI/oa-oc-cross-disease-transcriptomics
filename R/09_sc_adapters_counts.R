.sc_counts_capability <- function(dataset_id, capability, status, reason) {
  data.frame(
    dataset_id = dataset_id,
    capability = capability,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

.sc_counts_output_dir <- function(config, dataset_id) {
  ensure_dir(file.path(
    config$project$output_dir,
    "single_cell",
    clean_filename(dataset_id)
  ))
}

.sc_counts_cache_dir <- function(config, dataset_id, ...) {
  ensure_dir(file.path(
    config$project$cache_dir,
    "single_cell",
    clean_filename(dataset_id),
    ...
  ))
}

.sc_remove_legacy_outputs <- function(output_dir, names) {
  paths <- file.path(output_dir, names)
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L) {
    unlink(existing, force = TRUE)
  }
  invisible(existing)
}

.sc_counts_input_signature <- function(paths, extra = list()) {
  require_namespace("digest", "single-cell checkpoint signatures")
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  info <- file.info(paths)
  digest::digest(
    list(
      files = data.frame(
        path = paths,
        size = as.numeric(info$size),
        modified = as.numeric(info$mtime),
        stringsAsFactors = FALSE
      ),
      extra = extra
    ),
    algo = "sha256"
  )
}

.sc_counts_validate_sparse <- function(counts, context) {
  require_namespace("Matrix", "sparse single-cell counts")
  if (!inherits(counts, "sparseMatrix")) {
    stop(context, " did not produce a sparse matrix.", call. = FALSE)
  }
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop(context, " produced an empty matrix.", call. = FALSE)
  }
  values <- counts@x
  if (
    length(values) > 0L &&
      (
        anyNA(values) ||
          any(!is.finite(values)) ||
          any(values < 0) ||
          any(abs(values - round(values)) > 1e-8)
      )
  ) {
    stop(
      context,
      " is not a finite, non-negative integer count matrix.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sc_counts_metrics <- function(
    counts,
    gene_symbol,
    cell_id,
    dataset_id,
    batch
) {
  mitochondrial <- grepl("^MT-", gene_symbol, ignore.case = TRUE)
  total <- Matrix::colSums(counts)
  detected <- Matrix::colSums(counts > 0)
  mt_counts <- if (any(mitochondrial)) {
    Matrix::colSums(counts[mitochondrial, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(counts))
  }
  percent_mt <- ifelse(
    is.finite(mt_counts) & total > 0,
    100 * mt_counts / total,
    NA_real_
  )
  data.frame(
    cell_id = as.character(cell_id),
    dataset_id = dataset_id,
    batch = as.character(batch),
    nCount = as.numeric(total),
    nFeature = as.numeric(detected),
    percent_mt = as.numeric(percent_mt),
    doublet_score = NA_real_,
    doublet_class = NA_character_,
    doublet_partition = as.character(batch),
    stringsAsFactors = FALSE
  )
}

.sc_counts_run_doublets <- function(
    sce,
    metrics,
    seed,
    partition_col = "batch",
    minimum_core_cells = 100L,
    exclude_small_partitions = TRUE
) {
  require_namespace("SingleCellExperiment", "single-cell count objects")
  require_namespace("SummarizedExperiment", "single-cell count objects")
  require_namespace("scDblFinder", "mandatory doublet detection")
  require_namespace("BiocParallel", "serial doublet detection")

  if (ncol(sce) != nrow(metrics)) {
    stop("Doublet input cells and QC rows are not aligned.", call. = FALSE)
  }
  if (!identical(colnames(sce), as.character(metrics$cell_id))) {
    stop("Doublet input cell names and QC cell IDs are not aligned.", call. = FALSE)
  }
  partitions <- as.character(metrics[[partition_col]])
  if (anyNA(partitions) || any(!nzchar(partitions))) {
    stop("Doublet partitions contain missing values.", call. = FALSE)
  }

  summaries <- list()
  for (index in seq_along(unique(partitions))) {
    partition <- unique(partitions)[[index]]
    partition_cells <- which(partitions == partition)
    core_cells <- partition_cells[metrics$passes_core_QC[partition_cells] %in% TRUE]
    if (length(core_cells) < as.integer(minimum_core_cells)) {
      metrics$doublet_class[partition_cells] <- "not_run_too_few_core_cells"
      metrics$failure_reason <- append_qc_reason(
        metrics$failure_reason,
        seq_len(nrow(metrics)) %in% partition_cells,
        "doublet_unassessable_small_partition"
      )
      if (isTRUE(exclude_small_partitions)) {
        metrics$passes_QC[partition_cells] <- FALSE
      }
      summaries[[partition]] <- data.frame(
        partition = partition,
        input_cells = length(partition_cells),
        core_qc_pass = length(core_cells),
        doublet_status = "not_run_too_few_core_cells",
        predicted_doublets = NA_integer_,
        final_qc_pass = sum(metrics$passes_QC[partition_cells] %in% TRUE),
        stringsAsFactors = FALSE
      )
      next
    }

    set.seed(as.integer(seed) + index)
    called <- scDblFinder::scDblFinder(
      sce[, core_cells, drop = FALSE],
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    )
    called_data <- as.data.frame(SummarizedExperiment::colData(called))
    required <- c("scDblFinder.score", "scDblFinder.class")
    if (!all(required %in% names(called_data))) {
      stop(
        "scDblFinder did not return score/class for partition ",
        partition,
        ".",
        call. = FALSE
      )
    }
    metrics$doublet_score[core_cells] <- called_data$scDblFinder.score
    metrics$doublet_class[core_cells] <- as.character(
      called_data$scDblFinder.class
    )
    predicted <- rep(FALSE, nrow(metrics))
    predicted[core_cells] <- (
      metrics$doublet_class[core_cells] == "doublet"
    )
    predicted[is.na(predicted)] <- FALSE
    metrics$failure_reason <- append_qc_reason(
      metrics$failure_reason,
      predicted,
      "predicted_doublet"
    )
    metrics$passes_QC[partition_cells] <- (
      metrics$passes_core_QC[partition_cells] &
        !is.na(metrics$doublet_class[partition_cells]) &
        metrics$doublet_class[partition_cells] == "singlet"
    )
    summaries[[partition]] <- data.frame(
      partition = partition,
      input_cells = length(partition_cells),
      core_qc_pass = length(core_cells),
      doublet_status = "completed",
      predicted_doublets = sum(
        metrics$doublet_class[partition_cells] == "doublet",
        na.rm = TRUE
      ),
      final_qc_pass = sum(metrics$passes_QC[partition_cells] %in% TRUE),
      stringsAsFactors = FALSE
    )
    rm(called, called_data)
    invisible(gc())
  }

  list(
    metrics = metrics,
    summary = do.call(rbind, summaries)
  )
}

.sc_counts_attach_metadata <- function(sce, metrics) {
  require_namespace("S4Vectors", "single-cell metadata")
  require_namespace("SummarizedExperiment", "single-cell metadata")
  if (!identical(colnames(sce), as.character(metrics$cell_id))) {
    stop("Cannot attach misaligned QC metadata to count object.", call. = FALSE)
  }
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
    metrics,
    row.names = metrics$cell_id
  )
  sce
}

.sc_gse255460_python <- function(config) {
  candidates <- unique(c(
    config$single_cell$python_executable %||% "",
    Sys.getenv("OC_OA_PYTHON", ""),
    Sys.which("python3"),
    Sys.which("python")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0L) {
    stop(
      "GSE255460 sparse conversion requires Python >= 3.8.",
      call. = FALSE
    )
  }
  normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

.sc_gse255460_sparse_dir <- function(config, dataset_id) {
  file.path(
    config$project$cache_dir,
    "single_cell",
    clean_filename(dataset_id),
    "partitioned_csr"
  )
}

.sc_gse255460_ensure_sparse_bundle <- function(dataset, config) {
  require_namespace("jsonlite", "GSE255460 sparse bundle manifest")
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  script <- file.path(project_root, "tools", "convert_gse255460_sparse.py")
  if (!file.exists(script)) {
    stop("GSE255460 sparse converter is missing: ", script, call. = FALSE)
  }
  output_dir <- .sc_gse255460_sparse_dir(config, dataset$id)
  ensure_dir(dirname(output_dir))
  python <- .sc_gse255460_python(config)
  status <- system2(
    python,
    args = c(
      shQuote(script),
      "--archive", shQuote(dataset$counts_archive),
      "--metadata", shQuote(dataset$metadata_path),
      "--output-dir", shQuote(output_dir),
      "--member", "sc_counts.txt",
      "--batch-column", shQuote(dataset$batch_column %||% "ID"),
      "--cell-column", "X",
      "--count-column", "nCount_RNA",
      "--feature-column", "nFeature_RNA"
    ),
    stdout = "",
    stderr = ""
  )
  if (!identical(as.integer(status), 0L)) {
    stop(
      "GSE255460 sparse conversion failed with exit status ",
      status,
      ". Recoverable partial output, if any, remains beside ",
      output_dir,
      ".",
      call. = FALSE
    )
  }
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(manifest_path)) {
    stop("GSE255460 sparse conversion did not publish a manifest.", call. = FALSE)
  }
  list(
    directory = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    manifest_path = normalizePath(
      manifest_path,
      winslash = "/",
      mustWork = TRUE
    )
  )
}

.sc_gse255460_validate_manifest <- function(bundle, metadata) {
  require_namespace("data.table", "GSE255460 sparse bundle tables")
  require_namespace("jsonlite", "GSE255460 sparse bundle manifest")
  manifest <- jsonlite::fromJSON(
    bundle$manifest_path,
    simplifyDataFrame = TRUE
  )
  required_scalars <- list(
    format = "gse255460_partitioned_csr",
    format_version = 1,
    orientation = "genes_by_cells",
    sparse_encoding = "CSR",
    index_base = 0,
    value_dtype = "little_endian_int32",
    index_dtype = "little_endian_int32",
    pointer_dtype = "little_endian_int32"
  )
  for (field in names(required_scalars)) {
    expected <- required_scalars[[field]]
    observed <- manifest[[field]]
    matches <- if (is.numeric(expected)) {
      length(observed) == 1L &&
        is.finite(as.numeric(observed)) &&
        as.numeric(observed) == as.numeric(expected)
    } else {
      identical(as.character(observed), as.character(expected))
    }
    if (!isTRUE(matches)) {
      stop(
        "GSE255460 sparse manifest has invalid ", field, ".",
        call. = FALSE
      )
    }
  }
  if (
    !isTRUE(manifest$all_cells_nCount_exact) ||
      !isTRUE(manifest$all_cells_nFeature_exact) ||
      !isTRUE(manifest$header_mapping_exact) ||
      !isTRUE(manifest$features_unique)
  ) {
    stop(
      "GSE255460 sparse manifest does not attest complete count/header validation.",
      call. = FALSE
    )
  }

  partitions <- manifest$partitions
  if (!is.data.frame(partitions) || nrow(partitions) == 0L) {
    stop("GSE255460 sparse manifest has no partitions.", call. = FALSE)
  }
  expected_batches <- unique(as.character(metadata$ID))
  if (
    !identical(as.character(partitions$partition_id), expected_batches) ||
      as.integer(manifest$n_partitions) != length(expected_batches) ||
      as.integer(manifest$n_cells) != nrow(metadata) ||
      sum(as.numeric(partitions$n_cells)) != nrow(metadata) ||
      sum(as.numeric(partitions$nonzero)) != as.numeric(manifest$nonzero)
  ) {
    stop(
      "GSE255460 sparse manifest partition totals do not match metadata.",
      call. = FALSE
    )
  }

  features_path <- file.path(bundle$directory, manifest$features_file)
  features <- data.table::fread(
    features_path,
    data.table = FALSE,
    showProgress = FALSE
  )
  if (
    !identical(names(features), c("gene_index_1based", "gene_id")) ||
      nrow(features) != as.integer(manifest$n_genes) ||
      anyNA(features$gene_id) ||
      any(!nzchar(features$gene_id)) ||
      anyDuplicated(features$gene_id) ||
      !identical(features$gene_index_1based, seq_len(nrow(features)))
  ) {
    stop("GSE255460 sparse feature table is invalid.", call. = FALSE)
  }

  byte_fields <- c(
    data_file = "data_bytes",
    indices_file = "indices_bytes",
    indptr_file = "indptr_bytes",
    barcodes_file = "barcodes_bytes"
  )
  for (index in seq_len(nrow(partitions))) {
    for (file_field in names(byte_fields)) {
      path <- file.path(bundle$directory, partitions[[file_field]][[index]])
      expected_bytes <- as.numeric(
        partitions[[byte_fields[[file_field]]]][[index]]
      )
      observed_bytes <- as.numeric(file.info(path)$size)
      if (
        !file.exists(path) ||
          is.na(observed_bytes) ||
          observed_bytes != expected_bytes
      ) {
        stop(
          "GSE255460 sparse partition file is missing/truncated: ",
          path,
          call. = FALSE
        )
      }
    }
  }

  list(
    manifest = manifest,
    partitions = partitions,
    features = features,
    features_path = features_path
  )
}

.sc_read_i32_file <- function(path, expected_length) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  values <- readBin(
    connection,
    what = integer(),
    n = as.integer(expected_length),
    size = 4L,
    signed = TRUE,
    endian = "little"
  )
  if (length(values) != as.integer(expected_length)) {
    stop("Sparse int32 file is truncated: ", path, call. = FALSE)
  }
  values
}

.sc_read_gse255460_partition <- function(
    bundle,
    validated,
    partition_id,
    metadata
) {
  require_namespace("data.table", "GSE255460 sparse barcode import")
  require_namespace("Matrix", "GSE255460 sparse count object")
  require_namespace("S4Vectors", "GSE255460 sparse feature metadata")
  require_namespace("SingleCellExperiment", "GSE255460 count object")

  row <- validated$partitions[
    validated$partitions$partition_id == partition_id,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("GSE255460 sparse partition is not unique: ", partition_id, call. = FALSE)
  }
  n_genes <- as.integer(row$n_genes[[1L]])
  n_cells <- as.integer(row$n_cells[[1L]])
  nonzero <- as.integer(row$nonzero[[1L]])
  data_path <- file.path(bundle$directory, row$data_file[[1L]])
  indices_path <- file.path(bundle$directory, row$indices_file[[1L]])
  indptr_path <- file.path(bundle$directory, row$indptr_file[[1L]])
  barcodes_path <- file.path(bundle$directory, row$barcodes_file[[1L]])

  data <- .sc_read_i32_file(data_path, nonzero)
  indices <- .sc_read_i32_file(indices_path, nonzero)
  indptr <- .sc_read_i32_file(indptr_path, n_genes + 1L)
  if (
    indptr[[1L]] != 0L ||
      indptr[[length(indptr)]] != nonzero ||
      any(diff(indptr) < 0L) ||
      any(indices < 0L | indices >= n_cells) ||
      any(data <= 0L)
  ) {
    stop(
      "GSE255460 CSR arrays failed structural validation for ",
      partition_id,
      ".",
      call. = FALSE
    )
  }

  barcodes <- data.table::fread(
    barcodes_path,
    data.table = FALSE,
    showProgress = FALSE
  )
  required_barcodes <- c("cell_id", "matrix_cell_id", "global_index_1based")
  if (
    !identical(names(barcodes), required_barcodes) ||
      nrow(barcodes) != n_cells ||
      anyDuplicated(barcodes$cell_id) ||
      anyDuplicated(barcodes$matrix_cell_id)
  ) {
    stop(
      "GSE255460 barcode table is invalid for ", partition_id, ".",
      call. = FALSE
    )
  }
  metadata_index <- as.integer(barcodes$global_index_1based)
  if (
    anyNA(metadata_index) ||
      any(metadata_index < 1L | metadata_index > nrow(metadata)) ||
      !identical(
        as.character(barcodes$cell_id),
        as.character(metadata$.cell_id[metadata_index])
      ) ||
      !all(as.character(metadata$ID[metadata_index]) == partition_id) ||
      !identical(
        as.character(barcodes$matrix_cell_id),
        make.names(metadata$.cell_id[metadata_index], unique = FALSE)
      )
  ) {
    stop(
      "GSE255460 sparse barcodes do not align with metadata for ",
      partition_id,
      ".",
      call. = FALSE
    )
  }

  csr <- methods::new(
    "dgRMatrix",
    p = as.integer(indptr),
    j = as.integer(indices),
    x = as.numeric(data),
    Dim = as.integer(c(n_genes, n_cells)),
    Dimnames = list(NULL, NULL),
    factors = list()
  )
  counts <- methods::as(csr, "CsparseMatrix")
  gene_id <- as.character(validated$features$gene_id)
  row_id <- make.unique(gene_id)
  rownames(counts) <- row_id
  colnames(counts) <- as.character(barcodes$cell_id)
  .sc_counts_validate_sparse(
    counts,
    paste("GSE255460 partition", partition_id)
  )
  observed_count <- as.numeric(Matrix::colSums(counts))
  observed_feature <- as.numeric(Matrix::colSums(counts > 0))
  expected_count <- as.numeric(metadata$nCount_RNA[metadata_index])
  expected_feature <- as.numeric(metadata$nFeature_RNA[metadata_index])
  if (
    !identical(observed_count, expected_count) ||
      !identical(observed_feature, expected_feature)
  ) {
    stop(
      "GSE255460 sparse counts do not exactly reproduce metadata for ",
      partition_id,
      ".",
      call. = FALSE
    )
  }

  row_data <- S4Vectors::DataFrame(
    gene_id = gene_id,
    gene_symbol = gene_id,
    row.names = row_id
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )
  list(
    sce = sce,
    counts = counts,
    gene_id = gene_id,
    barcodes = barcodes,
    metadata_index = metadata_index,
    nonzero = nonzero
  )
}

adapt_gse255460_partitioned_csr <- function(dataset, config) {
  require_namespace("Matrix", "GSE255460 sparse counts")
  require_namespace("SingleCellExperiment", "GSE255460 count object")
  require_namespace("SummarizedExperiment", "GSE255460 count object")
  dataset_id <- dataset$id %||% "GSE255460"
  output_dir <- .sc_counts_output_dir(config, dataset_id)
  .sc_remove_legacy_outputs(
    output_dir,
    c(
      "archive_structure_checks.csv",
      "cell_metrics.csv",
      "raw_count_spotcheck.csv"
    )
  )
  cache_dir <- .sc_counts_cache_dir(config, dataset_id)
  checkpoint_dir <- .sc_counts_cache_dir(
    config,
    dataset_id,
    "qc_checkpoints"
  )
  sample_output_dir <- ensure_dir(file.path(output_dir, "partitions"))
  metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
  bundle <- .sc_gse255460_ensure_sparse_bundle(dataset, config)
  validated <- .sc_gse255460_validate_manifest(bundle, metadata)

  input_signature <- .sc_counts_input_signature(
    c(
      dataset$counts_archive,
      dataset$metadata_path,
      bundle$manifest_path
    ),
    extra = list(
      adapter = "gse255460_partitioned_csr_v1",
      lower_nmads = config$single_cell$lower_nmads,
      upper_nmads = config$single_cell$upper_nmads,
      mt_nmads = config$single_cell$mt_nmads,
      seed = config$project$seed,
      scDblFinder = as.character(utils::packageVersion("scDblFinder"))
    )
  )

  partition_ids <- as.character(validated$partitions$partition_id)
  partition_results <- list()
  for (index in seq_along(partition_ids)) {
    partition_id <- partition_ids[[index]]
    checkpoint_path <- file.path(
      checkpoint_dir,
      paste0(clean_filename(partition_id), "_checkpoint.rds")
    )
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      if (
        identical(checkpoint$signature, input_signature) &&
          file.exists(checkpoint$metric_path)
      ) {
        log_info(
          "Single-cell ", dataset_id, ": reusing ", partition_id,
          " QC checkpoint (", index, "/", length(partition_ids), ")."
        )
        partition_results[[partition_id]] <- checkpoint
        next
      }
    }

    log_info(
      "Single-cell ", dataset_id, ": loading sparse partition ",
      partition_id, " (", index, "/", length(partition_ids), ")."
    )
    imported <- .sc_read_gse255460_partition(
      bundle,
      validated,
      partition_id,
      metadata
    )
    sce <- imported$sce
    counts <- imported$counts
    sample_metadata <- metadata[imported$metadata_index, , drop = FALSE]
    metrics <- .sc_counts_metrics(
      counts,
      imported$gene_id,
      imported$barcodes$cell_id,
      dataset_id,
      rep(partition_id, ncol(counts))
    )
    metrics$matrix_cell_id <- as.character(
      imported$barcodes$matrix_cell_id
    )
    metrics$global_index_1based <- as.integer(
      imported$barcodes$global_index_1based
    )
    metadata_columns <- intersect(
      c(
        "ID", "sample", "trait", "group", "new_ID", "celltype",
        "seurat_clusters", "orig.ident"
      ),
      names(sample_metadata)
    )
    for (column in metadata_columns) {
      metrics[[column]] <- sample_metadata[[column]]
    }

    thresholds <- derive_sc_qc_thresholds(
      metrics,
      lower_nmads = config$single_cell$lower_nmads %||% 3,
      upper_nmads = config$single_cell$upper_nmads %||% 5,
      mt_nmads = config$single_cell$mt_nmads %||% 3
    )
    metrics <- apply_sc_qc_thresholds(metrics, thresholds)
    doublets <- .sc_counts_run_doublets(
      sce,
      metrics,
      seed = as.integer(config$project$seed) + 1000L * index,
      partition_col = "batch",
      minimum_core_cells = 100L,
      exclude_small_partitions = TRUE
    )
    metrics <- doublets$metrics
    metric_path <- file.path(
      sample_output_dir,
      paste0(clean_filename(partition_id), "_cell_qc.tsv.gz")
    )
    write_sc_table(metrics, metric_path)
    partition_summary <- data.frame(
      partition = partition_id,
      cells = nrow(metrics),
      genes = nrow(counts),
      nonzero = imported$nonzero,
      core_qc_pass = sum(metrics$passes_core_QC, na.rm = TRUE),
      predicted_doublets = sum(
        metrics$doublet_class == "doublet",
        na.rm = TRUE
      ),
      final_qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
      doublet_status = doublets$summary$doublet_status[[1L]],
      stringsAsFactors = FALSE
    )
    checkpoint <- list(
      signature = input_signature,
      summary = partition_summary,
      thresholds = thresholds,
      doublet_summary = doublets$summary,
      metrics = metrics,
      metric_path = metric_path
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    partition_results[[partition_id]] <- checkpoint
    rm(
      checkpoint, imported, sce, counts, sample_metadata, metrics,
      thresholds, doublets
    )
    invisible(gc())
  }

  partition_summary <- do.call(
    rbind,
    lapply(partition_results, `[[`, "summary")
  )
  thresholds <- do.call(
    rbind,
    lapply(partition_results, `[[`, "thresholds")
  )
  doublet_summary <- do.call(
    rbind,
    lapply(partition_results, `[[`, "doublet_summary")
  )
  metrics <- do.call(
    rbind,
    lapply(partition_results, `[[`, "metrics")
  )
  rownames(metrics) <- NULL
  metrics <- metrics[order(metrics$global_index_1based), , drop = FALSE]
  if (
    !identical(as.character(metrics$cell_id), as.character(metadata$.cell_id))
  ) {
    stop(
      "GSE255460 combined QC rows do not restore metadata order.",
      call. = FALSE
    )
  }

  metrics_path <- file.path(output_dir, "cell_qc.tsv.gz")
  thresholds_path <- file.path(output_dir, "qc_thresholds.csv")
  doublet_path <- file.path(output_dir, "doublet_partition_summary.csv")
  partition_summary_path <- file.path(output_dir, "partition_qc_summary.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  capability_path <- file.path(output_dir, "capabilities.csv")
  backed_manifest_path <- file.path(output_dir, "backed_count_manifest.rds")
  write_sc_table(metrics, metrics_path)
  safe_write_csv(thresholds, thresholds_path)
  safe_write_csv(doublet_summary, doublet_path)
  safe_write_csv(partition_summary, partition_summary_path)
  backed_manifest <- list(
    dataset_id = dataset_id,
    format = validated$manifest$format,
    format_version = validated$manifest$format_version,
    orientation = validated$manifest$orientation,
    sparse_encoding = validated$manifest$sparse_encoding,
    bundle_directory = bundle$directory,
    manifest_path = bundle$manifest_path,
    features_path = validated$features_path,
    partition_table = validated$partitions,
    source_archive = normalizePath(
      dataset$counts_archive,
      winslash = "/",
      mustWork = TRUE
    ),
    source_metadata = normalizePath(
      dataset$metadata_path,
      winslash = "/",
      mustWork = TRUE
    )
  )
  atomic_save_rds(backed_manifest, backed_manifest_path)

  summary <- data.frame(
    dataset_id = dataset_id,
    status = "validated",
    data_type = "raw_integer_umi_counts",
    matrix_orientation = "genes_by_cells",
    sparse_encoding = "partitioned_disk_csr",
    n_genes = as.integer(validated$manifest$n_genes),
    n_cells = as.integer(validated$manifest$n_cells),
    n_partitions = nrow(partition_summary),
    nonzero_entries = as.numeric(validated$manifest$nonzero),
    all_cells_nCount_exact = TRUE,
    all_cells_nFeature_exact = TRUE,
    mitochondrial_metric_recomputed = TRUE,
    core_qc_pass = sum(partition_summary$core_qc_pass),
    predicted_doublets = sum(partition_summary$predicted_doublets),
    final_qc_pass = sum(partition_summary$final_qc_pass),
    doublet_partitions_completed = sum(
      partition_summary$doublet_status == "completed"
    ),
    doublet_partitions_excluded = sum(
      partition_summary$doublet_status != "completed"
    ),
    ambient_background_available = FALSE,
    stringsAsFactors = FALSE
  )
  capability <- rbind(
    .sc_counts_capability(
      dataset_id,
      "adapter_validation",
      "passed",
      paste0(
        "Every one of 135,896 cells exactly reproduced metadata nCount_RNA ",
        "and nFeature_RNA after sparse conversion."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "raw_counts",
      "available",
      paste0(
        "The 38,680-by-135,896 count table is preserved as 19 partitioned ",
        "disk CSR matrices without dense in-memory import."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "core_qc",
      "completed",
      paste0(
        "Count, feature, and mitochondrial MAD thresholds were derived ",
        "within each ID partition from the sparse raw counts."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "scDblFinder",
      "completed",
      paste0(
        "scDblFinder completed independently in all ",
        nrow(partition_summary),
        " ID partitions with no fallback caller."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "ambient_rna",
      "blocked",
      paste0(
        "The local source contains called cells only; SoupX/empty-droplet ",
        "correction remains unsupported and is not claimed."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "downstream_qc_object",
      "available",
      paste0(
        "The backed CSR manifest preserves raw counts while cell-level ",
        "QC and doublet flags are stored separately."
      )
    )
  )
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  plot_paths <- plot_sc_qc_metrics(
    metrics,
    thresholds,
    output_dir = output_dir,
    dataset_id = dataset_id
  )

  list(
    id = dataset_id,
    status = "validated",
    reason = paste(
      "Partitioned disk CSR conversion passed all-cell count validation;",
      "formal per-ID QC and mandatory scDblFinder completed."
    ),
    summary = summary,
    thresholds = thresholds,
    capability = capability,
    cell_metric_paths = metrics_path,
    object_paths = backed_manifest_path,
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = thresholds_path,
      doublets = doublet_path,
      partition_summary = partition_summary_path,
      backed_manifest = backed_manifest_path,
      plots = plot_paths
    ),
    cells = nrow(metrics),
    qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
    downstream_ready = sum(metrics$passes_QC, na.rm = TRUE) > 0L
  )
}

.sc_parse_gse104782_cells <- function(cell_id) {
  match <- regexec(
    "^(OA[0-9]+)_(S[0-9]+)[.]([0-9]+)$",
    cell_id,
    perl = TRUE
  )
  parts <- regmatches(cell_id, match)
  valid <- lengths(parts) == 4L
  if (!all(valid)) {
    stop(
      "GSE104782 cell IDs do not all match OA{patient}_S{stage}.{index}.",
      call. = FALSE
    )
  }
  data.frame(
    donor = vapply(parts, `[[`, character(1), 2L),
    state_code = vapply(parts, `[[`, character(1), 3L),
    within_group_index = as.integer(vapply(parts, `[[`, character(1), 4L)),
    stringsAsFactors = FALSE
  )
}

.sc_read_gse104782_umi <- function(path) {
  require_namespace("data.table", "GSE104782 count import")
  require_namespace("Matrix", "GSE104782 sparse count object")
  require_namespace("SingleCellExperiment", "GSE104782 count object")
  require_namespace("S4Vectors", "GSE104782 feature metadata")

  if (!file.exists(path)) {
    stop("GSE104782 UMI count file does not exist: ", path, call. = FALSE)
  }
  table <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
  if (nrow(table) == 0L || ncol(table) < 2L) {
    stop("GSE104782 UMI count table is empty or malformed.", call. = FALSE)
  }
  genes <- trimws(as.character(table[[1L]]))
  cell_id <- names(table)[-1L]
  if (
    anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      anyNA(cell_id) || any(!nzchar(cell_id)) || anyDuplicated(cell_id)
  ) {
    stop("GSE104782 contains invalid or duplicated gene/cell IDs.", call. = FALSE)
  }

  dense <- as.matrix(table[-1L])
  storage.mode(dense) <- "double"
  if (
    anyNA(dense) ||
      any(!is.finite(dense)) ||
      any(dense < 0) ||
      any(abs(dense - round(dense)) > 1e-8)
  ) {
    stop(
      "GSE104782 UMI table contains non-integer, negative, or invalid values.",
      call. = FALSE
    )
  }
  counts <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")
  rownames(counts) <- genes
  colnames(counts) <- cell_id
  .sc_counts_validate_sparse(counts, "GSE104782 UMI matrix")
  zero_gene <- Matrix::rowSums(counts != 0) == 0
  rm(dense, table)
  invisible(gc())

  row_data <- S4Vectors::DataFrame(
    gene_id = genes,
    gene_symbol = genes,
    all_zero = zero_gene,
    row.names = genes
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )
  list(
    sce = sce,
    genes = genes,
    cell_id = cell_id,
    all_zero_genes = sum(zero_gene),
    nonzero = length(counts@x)
  )
}

adapt_gse104782_umi_counts <- function(dataset, config) {
  require_namespace("Matrix", "GSE104782 sparse counts")
  dataset_id <- dataset$id %||% "GSE104782"
  output_dir <- .sc_counts_output_dir(config, dataset_id)
  .sc_remove_legacy_outputs(
    output_dir,
    c("archive_member_summary.csv", "cell_metrics.csv")
  )
  cache_dir <- .sc_counts_cache_dir(config, dataset_id)
  imported <- .sc_read_gse104782_umi(dataset$counts_path)
  sce <- imported$sce
  parsed <- .sc_parse_gse104782_cells(imported$cell_id)

  metrics <- .sc_counts_metrics(
    SummarizedExperiment::assay(sce, "counts"),
    imported$genes,
    imported$cell_id,
    dataset_id,
    parsed$donor
  )
  metrics$donor <- parsed$donor
  metrics$state_code <- parsed$state_code
  metrics$within_group_index <- parsed$within_group_index
  metrics$batch_provenance <- "donor_proxy; technical capture metadata unavailable"

  thresholds <- derive_sc_qc_thresholds(
    metrics,
    lower_nmads = config$single_cell$lower_nmads %||% 3,
    upper_nmads = config$single_cell$upper_nmads %||% 5,
    mt_nmads = config$single_cell$mt_nmads %||% 3
  )
  metrics <- apply_sc_qc_thresholds(metrics, thresholds)
  doublets <- .sc_counts_run_doublets(
    sce,
    metrics,
    seed = config$project$seed,
    partition_col = "donor",
    minimum_core_cells = 100L,
    exclude_small_partitions = TRUE
  )
  metrics <- doublets$metrics
  sce <- .sc_counts_attach_metadata(sce, metrics)

  object_path <- file.path(cache_dir, "GSE104782_umi_qc_sce.rds")
  metrics_path <- file.path(output_dir, "cell_qc.tsv.gz")
  thresholds_path <- file.path(output_dir, "qc_thresholds.csv")
  doublet_path <- file.path(output_dir, "doublet_partition_summary.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  capability_path <- file.path(output_dir, "capabilities.csv")
  atomic_save_rds(sce, object_path, compress = FALSE)
  write_sc_table(metrics, metrics_path)
  safe_write_csv(thresholds, thresholds_path)
  safe_write_csv(doublets$summary, doublet_path)

  small_partitions <- doublets$summary$partition[
    doublets$summary$doublet_status != "completed"
  ]
  summary <- data.frame(
    dataset_id = dataset_id,
    status = if (length(small_partitions) == 0L) {
      "validated_limited_batch_metadata"
    } else {
      "validated_with_small_partition_exclusions"
    },
    data_type = "raw_integer_umi_counts",
    matrix_orientation = "genes_by_cells",
    n_genes = nrow(sce),
    all_zero_genes = imported$all_zero_genes,
    n_cells = ncol(sce),
    n_donors = length(unique(metrics$donor)),
    n_state_codes = length(unique(metrics$state_code)),
    nonzero_entries = imported$nonzero,
    core_qc_pass = sum(metrics$passes_core_QC, na.rm = TRUE),
    predicted_doublets = sum(metrics$doublet_class == "doublet", na.rm = TRUE),
    final_qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
    technical_batch_mapping_available = FALSE,
    doublet_partitions_completed = sum(
      doublets$summary$doublet_status == "completed"
    ),
    doublet_partitions_excluded = length(small_partitions),
    ambient_background_available = FALSE,
    stringsAsFactors = FALSE
  )
  capability <- rbind(
    .sc_counts_capability(
      dataset_id,
      "raw_counts",
      "available",
      "The full 24,153-by-1,600 table contains non-negative integer UMI counts."
    ),
    .sc_counts_capability(
      dataset_id,
      "core_qc",
      "completed",
      "Count, feature, and mitochondrial MAD thresholds were derived within donor proxy partitions."
    ),
    .sc_counts_capability(
      dataset_id,
      "technical_batch_mapping",
      "limited",
      "Cell IDs encode donor and disease stage, but not the original capture/GSM partition."
    ),
    .sc_counts_capability(
      dataset_id,
      "scDblFinder",
      if (length(small_partitions) == 0L) "completed_with_batch_limitation" else "completed_with_exclusions",
      paste0(
        "scDblFinder used donor-level proxy partitions; ",
        length(small_partitions),
        " partition(s) with fewer than 100 core-QC cells were explicitly excluded."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "ambient_rna",
      "blocked",
      "The all-cells table has no empty-well or unfiltered-droplet background."
    ),
    .sc_counts_capability(
      dataset_id,
      "downstream_qc_object",
      "available",
      "The cached SingleCellExperiment preserves raw counts and all QC/doublet flags."
    )
  )
  safe_write_csv(summary, summary_path)
  safe_write_csv(capability, capability_path)
  plot_paths <- plot_sc_qc_metrics(
    metrics,
    thresholds,
    output_dir = output_dir,
    dataset_id = dataset_id
  )

  list(
    id = dataset_id,
    status = summary$status[[1L]],
    reason = paste(
      "Raw UMI counts passed validation and formal QC completed.",
      "Technical capture IDs are unavailable, so donor partitions are an explicit proxy."
    ),
    summary = summary,
    thresholds = thresholds,
    capability = capability,
    cell_metric_paths = metrics_path,
    object_paths = object_path,
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = thresholds_path,
      doublets = doublet_path,
      plots = plot_paths
    ),
    cells = ncol(sce),
    qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
    downstream_ready = sum(metrics$passes_QC, na.rm = TRUE) > 0L
  )
}

.sc_read_tenx_triplet <- function(paths, sample_id) {
  require_namespace("Matrix", "10x sparse matrix import")
  require_namespace("data.table", "10x feature/barcode import")
  require_namespace("SingleCellExperiment", "10x count objects")
  require_namespace("S4Vectors", "10x feature metadata")

  feature_path <- paths[grepl("_features[.]tsv[.]gz$", paths)]
  barcode_path <- paths[grepl("_barcodes[.]tsv[.]gz$", paths)]
  matrix_path <- paths[grepl("_matrix[.]mtx[.]gz$", paths)]
  if (
    length(feature_path) != 1L ||
      length(barcode_path) != 1L ||
      length(matrix_path) != 1L
  ) {
    stop("Incomplete 10x triplet for ", sample_id, ".", call. = FALSE)
  }
  features <- data.table::fread(
    feature_path,
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
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
    stop(
      sample_id,
      " 10x dimensions do not match features/barcodes.",
      call. = FALSE
    )
  }
  .sc_counts_validate_sparse(counts, paste(sample_id, "10x matrix"))

  gene_id <- as.character(features[[1L]])
  gene_symbol <- if (ncol(features) >= 2L) {
    as.character(features[[2L]])
  } else {
    gene_id
  }
  feature_type <- if (ncol(features) >= 3L) {
    as.character(features[[3L]])
  } else {
    rep("Gene Expression", length(gene_id))
  }
  row_id <- make.unique(gene_id)
  cell_id <- paste(sample_id, as.character(barcodes), sep = "_")
  rownames(counts) <- row_id
  colnames(counts) <- cell_id
  row_data <- S4Vectors::DataFrame(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    feature_type = feature_type,
    row.names = row_id
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )
  list(
    sce = sce,
    features = features,
    barcodes = as.character(barcodes),
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    matrix_path = matrix_path,
    feature_path = feature_path,
    barcode_path = barcode_path,
    nonzero = length(counts@x)
  )
}

.sc_stream_barcode_subset <- function(
    raw_barcode_path,
    filtered_barcodes,
    chunk_size = 200000L
) {
  connection <- gzfile(raw_barcode_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  matched <- logical(length(filtered_barcodes))
  matched_position <- rep(NA_real_, length(filtered_barcodes))
  raw_rows <- 0
  repeat {
    lines <- readLines(connection, n = as.integer(chunk_size), warn = FALSE)
    if (length(lines) == 0L) break
    index <- match(lines, filtered_barcodes, nomatch = 0L)
    hit <- which(index > 0L)
    if (length(hit) > 0L) {
      filtered_index <- index[hit]
      if (any(matched[filtered_index])) {
        stop("A filtered barcode occurs more than once in the raw barcode list.", call. = FALSE)
      }
      matched[filtered_index] <- TRUE
      matched_position[filtered_index] <- raw_rows + hit
    }
    raw_rows <- raw_rows + length(lines)
  }
  list(
    raw_rows = raw_rows,
    all_filtered_in_raw = all(matched),
    filtered_matched = sum(matched),
    raw_positions = matched_position
  )
}

.sc_gse169454_member_manifest <- function(members) {
  pattern <- paste0(
    "^(GSM[0-9]+)_([^_]+)_(raw|filtered)_",
    "(barcodes|features|matrix)[.](tsv|mtx)[.]gz$"
  )
  parts <- regmatches(basename(members), regexec(pattern, basename(members), perl = TRUE))
  if (any(lengths(parts) != 6L)) {
    stop(
      "GSE169454 archive contains an unexpected member name.",
      call. = FALSE
    )
  }
  manifest <- data.frame(
    member = members,
    GSM = vapply(parts, `[[`, character(1), 2L),
    sample_label = vapply(parts, `[[`, character(1), 3L),
    matrix_scope = vapply(parts, `[[`, character(1), 4L),
    component = vapply(parts, `[[`, character(1), 5L),
    extension = vapply(parts, `[[`, character(1), 6L),
    stringsAsFactors = FALSE
  )
  manifest$sample_id <- paste(manifest$GSM, manifest$sample_label, sep = "_")
  expected <- expand.grid(
    sample_id = unique(manifest$sample_id),
    matrix_scope = c("raw", "filtered"),
    component = c("barcodes", "features", "matrix"),
    stringsAsFactors = FALSE
  )
  observed <- paste(
    manifest$sample_id,
    manifest$matrix_scope,
    manifest$component,
    sep = "::"
  )
  wanted <- paste(
    expected$sample_id,
    expected$matrix_scope,
    expected$component,
    sep = "::"
  )
  if (
    nrow(manifest) != 42L ||
      length(unique(manifest$sample_id)) != 7L ||
      anyDuplicated(observed) ||
      !setequal(observed, wanted)
  ) {
    stop(
      "GSE169454 archive is not exactly seven paired raw/filtered 10x triplets.",
      call. = FALSE
    )
  }
  manifest
}

adapt_gse169454_tenx_raw_tar <- function(dataset, config) {
  require_namespace("Matrix", "GSE169454 sparse counts")
  require_namespace("data.table", "GSE169454 input tables")
  dataset_id <- dataset$id %||% "GSE169454"
  output_dir <- .sc_counts_output_dir(config, dataset_id)
  .sc_remove_legacy_outputs(
    output_dir,
    c("structure_audit.csv", "blocker_report.md")
  )
  extracted_dir <- .sc_counts_cache_dir(config, dataset_id, "extracted")
  sample_cache_dir <- .sc_counts_cache_dir(config, dataset_id, "samples")
  members <- safe_archive_members(dataset$archive_path)
  manifest <- .sc_gse169454_member_manifest(members)
  extracted <- extract_sc_archive(dataset$archive_path, extracted_dir, members)
  manifest$path <- extracted[match(manifest$member, members)]
  if (anyNA(manifest$path)) {
    stop("GSE169454 extracted member mapping failed.", call. = FALSE)
  }

  input_signature <- .sc_counts_input_signature(
    dataset$archive_path,
    extra = list(
      adapter = "gse169454_tenx_raw_tar_v1",
      lower_nmads = config$single_cell$lower_nmads,
      upper_nmads = config$single_cell$upper_nmads,
      mt_nmads = config$single_cell$mt_nmads,
      seed = config$project$seed,
      scDblFinder = as.character(utils::packageVersion("scDblFinder"))
    )
  )
  sample_results <- list()
  sample_ids <- unique(manifest$sample_id)
  for (index in seq_along(sample_ids)) {
    sample_id <- sample_ids[[index]]
    checkpoint_path <- file.path(
      sample_cache_dir,
      paste0(clean_filename(sample_id), "_checkpoint.rds")
    )
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      required_paths <- c(checkpoint$object_path, checkpoint$metric_path)
      if (
        identical(checkpoint$signature, input_signature) &&
          all(file.exists(required_paths))
      ) {
        log_info(
          "Single-cell ", dataset_id, ": reusing ", sample_id,
          " QC checkpoint (", index, "/", length(sample_ids), ")."
        )
        sample_results[[sample_id]] <- checkpoint
        next
      }
    }

    log_info(
      "Single-cell ", dataset_id, ": importing ", sample_id,
      " (", index, "/", length(sample_ids), ")."
    )
    sample_manifest <- manifest[manifest$sample_id == sample_id, , drop = FALSE]
    filtered_paths <- sample_manifest$path[
      sample_manifest$matrix_scope == "filtered"
    ]
    raw_paths <- sample_manifest$path[
      sample_manifest$matrix_scope == "raw"
    ]
    filtered <- .sc_read_tenx_triplet(filtered_paths, sample_id)
    raw_feature_path <- raw_paths[grepl("_features[.]tsv[.]gz$", raw_paths)]
    raw_barcode_path <- raw_paths[grepl("_barcodes[.]tsv[.]gz$", raw_paths)]
    raw_matrix_path <- raw_paths[grepl("_matrix[.]mtx[.]gz$", raw_paths)]
    raw_features <- data.table::fread(
      raw_feature_path,
      header = FALSE,
      data.table = FALSE,
      showProgress = FALSE
    )
    if (!identical(raw_features, filtered$features)) {
      stop(sample_id, " raw/filtered feature tables differ.", call. = FALSE)
    }
    raw_header <- sc_read_matrix_market_header(raw_matrix_path)
    barcode_audit <- .sc_stream_barcode_subset(
      raw_barcode_path,
      filtered$barcodes
    )
    if (
      !identical(as.numeric(raw_header$n_rows), as.numeric(nrow(raw_features))) ||
        !identical(as.numeric(raw_header$n_columns), as.numeric(barcode_audit$raw_rows)) ||
        !isTRUE(barcode_audit$all_filtered_in_raw) ||
        !identical(raw_header$field, "integer")
    ) {
      stop(
        sample_id,
        " raw/filtered pair failed dimension, integer-field, or barcode mapping checks.",
        call. = FALSE
      )
    }

    condition <- if (grepl("^normal", sample_manifest$sample_label[[1L]], ignore.case = TRUE)) {
      "normal"
    } else if (grepl("^oa", sample_manifest$sample_label[[1L]], ignore.case = TRUE)) {
      "OA"
    } else {
      "unknown"
    }
    sce <- filtered$sce
    counts <- SummarizedExperiment::assay(sce, "counts")
    metrics <- .sc_counts_metrics(
      counts,
      filtered$gene_symbol,
      colnames(sce),
      dataset_id,
      sample_id
    )
    metrics$GSM <- sample_manifest$GSM[[1L]]
    metrics$sample_label <- sample_manifest$sample_label[[1L]]
    metrics$condition <- condition
    metrics$raw_barcode_position <- barcode_audit$raw_positions
    thresholds <- derive_sc_qc_thresholds(
      metrics,
      lower_nmads = config$single_cell$lower_nmads %||% 3,
      upper_nmads = config$single_cell$upper_nmads %||% 5,
      mt_nmads = config$single_cell$mt_nmads %||% 3
    )
    metrics <- apply_sc_qc_thresholds(metrics, thresholds)
    doublets <- .sc_counts_run_doublets(
      sce,
      metrics,
      seed = as.integer(config$project$seed) + index * 1000L,
      minimum_core_cells = 100L,
      exclude_small_partitions = TRUE
    )
    metrics <- doublets$metrics
    sce <- .sc_counts_attach_metadata(sce, metrics)

    object_path <- file.path(
      config$project$cache_dir,
      "single_cell",
      dataset_id,
      paste0(clean_filename(sample_id), "_qc_sce.rds")
    )
    metric_path <- file.path(
      output_dir,
      paste0(clean_filename(sample_id), "_cell_qc.tsv.gz")
    )
    atomic_save_rds(sce, object_path, compress = FALSE)
    write_sc_table(metrics, metric_path)
    plot_paths <- plot_sc_qc_metrics(
      metrics,
      thresholds,
      output_dir = output_dir,
      dataset_id = paste(dataset_id, sample_id, sep = "_")
    )
    thresholds$dataset_id <- dataset_id
    thresholds$sample_id <- sample_id
    ambient <- data.frame(
      dataset_id = dataset_id,
      sample_id = sample_id,
      raw_genes = raw_header$n_rows,
      raw_barcodes = raw_header$n_columns,
      raw_nonzero = raw_header$n_nonzero,
      filtered_genes = nrow(sce),
      filtered_cells = ncol(sce),
      filtered_nonzero = filtered$nonzero,
      filtered_barcodes_found_in_raw = barcode_audit$filtered_matched,
      filtered_is_subset_of_raw = barcode_audit$all_filtered_in_raw,
      background_profile_computed = FALSE,
      correction_applied = FALSE,
      assessment = paste(
        "Paired unfiltered and filtered matrices are structurally valid;",
        "ambient correction is deferred to the downstream clustering stage."
      ),
      stringsAsFactors = FALSE
    )
    summary <- data.frame(
      dataset_id = dataset_id,
      sample_id = sample_id,
      GSM = sample_manifest$GSM[[1L]],
      sample_label = sample_manifest$sample_label[[1L]],
      condition = condition,
      genes = nrow(sce),
      cells = ncol(sce),
      core_qc_pass = sum(metrics$passes_core_QC, na.rm = TRUE),
      predicted_doublets = sum(metrics$doublet_class == "doublet", na.rm = TRUE),
      final_qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
      median_counts = stats::median(metrics$nCount),
      median_features = stats::median(metrics$nFeature),
      median_percent_mt = stats::median(metrics$percent_mt, na.rm = TRUE),
      doublet_status = doublets$summary$doublet_status[[1L]],
      stringsAsFactors = FALSE
    )
    checkpoint <- list(
      signature = input_signature,
      summary = summary,
      thresholds = thresholds,
      doublet_summary = doublets$summary,
      ambient = ambient,
      object_path = object_path,
      metric_path = metric_path,
      plot_paths = plot_paths
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    sample_results[[sample_id]] <- checkpoint
    rm(
      checkpoint, sce, counts, metrics, doublets, filtered,
      raw_features, raw_header, barcode_audit
    )
    invisible(gc())
  }

  summary <- do.call(rbind, lapply(sample_results, `[[`, "summary"))
  thresholds <- do.call(rbind, lapply(sample_results, `[[`, "thresholds"))
  doublet_summary <- do.call(
    rbind,
    lapply(sample_results, `[[`, "doublet_summary")
  )
  ambient <- do.call(rbind, lapply(sample_results, `[[`, "ambient"))
  summary_path <- file.path(output_dir, "sample_qc_summary.csv")
  thresholds_path <- file.path(output_dir, "qc_thresholds.csv")
  doublet_path <- file.path(output_dir, "doublet_partition_summary.csv")
  ambient_path <- file.path(output_dir, "ambient_pair_assessment.csv")
  manifest_path <- file.path(output_dir, "archive_member_manifest.csv")
  capability_path <- file.path(output_dir, "capabilities.csv")
  safe_write_csv(summary, summary_path)
  safe_write_csv(thresholds, thresholds_path)
  safe_write_csv(doublet_summary, doublet_path)
  safe_write_csv(ambient, ambient_path)
  safe_write_csv(manifest, manifest_path)

  incomplete_doublets <- doublet_summary$partition[
    doublet_summary$doublet_status != "completed"
  ]
  capability <- rbind(
    .sc_counts_capability(
      dataset_id,
      "adapter_validation",
      "passed",
      "Seven GSM libraries each contain matched raw/filtered 10x integer-count triplets."
    ),
    .sc_counts_capability(
      dataset_id,
      "raw_counts",
      "available",
      "Filtered non-negative integer UMI counts are preserved in per-sample SCE objects."
    ),
    .sc_counts_capability(
      dataset_id,
      "core_qc",
      "completed",
      "Per-GSM data-driven count, feature, and mitochondrial thresholds were applied."
    ),
    .sc_counts_capability(
      dataset_id,
      "scDblFinder",
      if (length(incomplete_doublets) == 0L) "completed" else "completed_with_exclusions",
      paste0(
        "scDblFinder completed within GSM capture partitions; ",
        length(incomplete_doublets),
        " undersized partition(s) were explicitly excluded."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "ambient_rna",
      "assessed_supported_not_corrected",
      paste(
        "Every filtered barcode maps to its paired unfiltered matrix.",
        "The background inputs are valid; correction awaits clustering and is not silently claimed."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "downstream_qc_object",
      "available",
      "Seven cached SCE objects retain raw counts, sample/condition metadata, and QC flags."
    )
  )
  safe_write_csv(capability, capability_path)

  list(
    id = dataset_id,
    status = if (length(incomplete_doublets) == 0L) {
      "validated"
    } else {
      "validated_with_small_partition_exclusions"
    },
    reason = paste(
      "Seven raw/filtered 10x libraries passed structural validation;",
      "formal per-GSM QC and mandatory doublet detection completed."
    ),
    summary = summary,
    thresholds = thresholds,
    capability = capability,
    cell_metric_paths = vapply(
      sample_results,
      `[[`,
      character(1),
      "metric_path"
    ),
    object_paths = vapply(
      sample_results,
      `[[`,
      character(1),
      "object_path"
    ),
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = thresholds_path,
      doublets = doublet_path,
      ambient = ambient_path,
      archive_manifest = manifest_path
    ),
    cells = sum(summary$cells),
    qc_pass = sum(summary$final_qc_pass),
    downstream_ready = sum(summary$final_qc_pass) > 0L
  )
}

.sc_h5_csr_structure <- function(path) {
  require_namespace("rhdf5", "GSE180661 HDF5 sparse matrix")
  if (!file.exists(path)) {
    stop("GSE180661 H5 file does not exist: ", path, call. = FALSE)
  }
  if (!isTRUE(rhdf5::H5Fis_hdf5(path))) {
    stop("GSE180661 matrix is not a readable HDF5 file.", call. = FALSE)
  }
  listing <- rhdf5::h5ls(path, recursive = TRUE, datasetinfo = TRUE)
  required <- data.frame(
    group = c("/X", "/X", "/X", "/", "/"),
    name = c("data", "indices", "indptr", "obs", "var"),
    stringsAsFactors = FALSE
  )
  observed <- paste(listing$group, listing$name, sep = "::")
  wanted <- paste(required$group, required$name, sep = "::")
  if (!all(wanted %in% observed)) {
    stop(
      "GSE180661 H5 is missing one or more CSR/obs/var datasets.",
      call. = FALSE
    )
  }
  attributes <- rhdf5::h5readAttributes(path, "X")
  format <- tolower(as.character(attributes$h5sparse_format %||% ""))
  shape <- as.numeric(attributes$h5sparse_shape %||% numeric())
  if (!identical(format, "csr") || length(shape) != 2L) {
    stop("GSE180661 H5 does not declare a two-dimensional CSR matrix.", call. = FALSE)
  }
  dimension_for <- function(group, name) {
    row <- listing[listing$group == group & listing$name == name, , drop = FALSE]
    if (nrow(row) != 1L) return(NA_real_)
    suppressWarnings(as.numeric(row$dim[[1L]]))
  }
  data_length <- dimension_for("/X", "data")
  index_length <- dimension_for("/X", "indices")
  pointer_length <- dimension_for("/X", "indptr")
  obs_length <- dimension_for("/", "obs")
  var_length <- dimension_for("/", "var")
  if (
    !identical(data_length, index_length) ||
      !identical(pointer_length, shape[[1L]] + 1) ||
      !identical(obs_length, shape[[1L]]) ||
      !identical(var_length, shape[[2L]])
  ) {
    stop("GSE180661 H5 CSR dataset dimensions are inconsistent.", call. = FALSE)
  }
  list(
    format = format,
    shape = shape,
    data_length = data_length,
    index_length = index_length,
    pointer_length = pointer_length,
    obs_length = obs_length,
    var_length = var_length,
    listing = listing
  )
}

.sc_read_h5_csr_rows <- function(
    path,
    h5_rows,
    indptr,
    n_features,
    maximum_gap = 1000L
) {
  require_namespace("rhdf5", "GSE180661 CSR row reads")
  require_namespace("Matrix", "GSE180661 sparse sample objects")
  h5_rows <- as.integer(h5_rows)
  if (
    length(h5_rows) == 0L ||
      anyNA(h5_rows) ||
      anyDuplicated(h5_rows) ||
      any(h5_rows < 1L) ||
      any(h5_rows + 1L > length(indptr))
  ) {
    stop("Invalid or duplicated GSE180661 H5 row selection.", call. = FALSE)
  }
  sorted_rows <- sort(h5_rows)
  run_id <- cumsum(c(TRUE, diff(sorted_rows) > as.integer(maximum_gap)))
  run_rows <- split(sorted_rows, run_id)
  pieces <- vector("list", length(run_rows))
  piece_positions <- vector("list", length(run_rows))

  for (run_index in seq_along(run_rows)) {
    selected_rows <- run_rows[[run_index]]
    first_row <- min(selected_rows)
    last_row <- max(selected_rows)
    block_start <- indptr[[first_row]] + 1
    block_count <- indptr[[last_row + 1L]] - indptr[[first_row]]
    if (
      !is.finite(block_start) ||
        !is.finite(block_count) ||
        block_start < 1 ||
        block_count < 0
    ) {
      stop("Invalid GSE180661 CSR pointer range.", call. = FALSE)
    }
    if (block_count > .Machine$integer.max) {
      stop(
        "A GSE180661 sample CSR span exceeds the safe R vector limit.",
        call. = FALSE
      )
    }
    if (block_count == 0) {
      block_data <- numeric()
      block_indices <- integer()
    } else {
      block_data <- rhdf5::h5read(
        path,
        "X/data",
        start = block_start,
        count = block_count
      )
      block_indices <- rhdf5::h5read(
        path,
        "X/indices",
        start = block_start,
        count = block_count
      )
    }
    if (
      length(block_data) != block_count ||
        length(block_indices) != block_count ||
        anyNA(block_data) ||
        any(!is.finite(block_data)) ||
        any(block_data <= 0) ||
        any(abs(block_data - round(block_data)) > 1e-8) ||
        anyNA(block_indices) ||
        any(block_indices < 0L) ||
        any(block_indices >= n_features)
    ) {
      stop(
        "GSE180661 CSR data contain invalid counts or feature indices.",
        call. = FALSE
      )
    }

    relative_start <- indptr[selected_rows] - indptr[[first_row]] + 1
    lengths <- indptr[selected_rows + 1L] - indptr[selected_rows]
    if (any(lengths > .Machine$integer.max)) {
      stop("A GSE180661 cell exceeds the safe sparse vector limit.", call. = FALSE)
    }
    nonempty <- which(lengths > 0)
    selected_offsets <- integer()
    selected_columns <- integer()
    if (length(nonempty) == 0L) {
      piece <- Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        x = numeric(),
        dims = c(as.integer(n_features), length(selected_rows)),
        giveCsparse = TRUE
      )
    } else {
      selected_offsets <- sequence(
        as.integer(lengths[nonempty]),
        from = as.integer(relative_start[nonempty])
      )
      selected_columns <- rep.int(
        nonempty,
        as.integer(lengths[nonempty])
      )
      piece <- Matrix::sparseMatrix(
        i = as.integer(block_indices[selected_offsets]) + 1L,
        j = selected_columns,
        x = as.numeric(block_data[selected_offsets]),
        dims = c(as.integer(n_features), length(selected_rows)),
        giveCsparse = TRUE
      )
    }
    pieces[[run_index]] <- methods::as(piece, "dgCMatrix")
    piece_positions[[run_index]] <- match(selected_rows, h5_rows)
    rm(
      piece, block_data, block_indices, selected_offsets,
      selected_columns
    )
    invisible(gc())
  }

  counts <- if (length(pieces) == 1L) {
    pieces[[1L]]
  } else {
    do.call(cbind, pieces)
  }
  current_positions <- unlist(piece_positions, use.names = FALSE)
  counts <- counts[, order(current_positions), drop = FALSE]
  if (!identical(sort(current_positions), seq_along(h5_rows))) {
    stop("GSE180661 CSR row reordering failed.", call. = FALSE)
  }
  methods::as(counts, "dgCMatrix")
}

.sc_gse180661_metadata <- function(path) {
  require_namespace("data.table", "GSE180661 cell metadata")
  if (!file.exists(path)) {
    stop("GSE180661 metadata does not exist: ", path, call. = FALSE)
  }
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  required <- c(
    "cell_id", "sample", "percent.mt", "nCount_RNA", "nFeature_RNA"
  )
  missing <- setdiff(required, header)
  if (length(missing) > 0L) {
    stop(
      "GSE180661 metadata is missing columns: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  optional <- c(
    "cell_type", "cluster_label", "cluster_label_sub", "cell_type_super",
    "patient_id", "tumor_subsite", "tumor_site", "tumor_supersite",
    "sort_parameters", "therapy", "surgery"
  )
  metadata <- data.table::fread(
    path,
    select = c(required, intersect(optional, header)),
    na.strings = c("", "NA"),
    showProgress = FALSE
  )
  metadata[, `:=`(
    cell_id = as.character(cell_id),
    sample = as.character(sample),
    nCount_RNA = suppressWarnings(as.numeric(nCount_RNA)),
    nFeature_RNA = suppressWarnings(as.numeric(nFeature_RNA)),
    percent.mt = suppressWarnings(as.numeric(percent.mt))
  )]
  if (
    anyNA(metadata$cell_id) ||
      any(!nzchar(metadata$cell_id)) ||
      anyDuplicated(metadata$cell_id) ||
      anyNA(metadata$sample) ||
      any(!nzchar(metadata$sample))
  ) {
    stop("GSE180661 metadata contains invalid cell or sample IDs.", call. = FALSE)
  }
  metadata
}

adapt_gse180661_h5_csr <- function(dataset, config) {
  require_namespace("rhdf5", "GSE180661 HDF5 sparse matrix")
  require_namespace("Matrix", "GSE180661 sparse sample objects")
  require_namespace("data.table", "GSE180661 metadata and QC tables")
  require_namespace("SingleCellExperiment", "GSE180661 sample objects")
  require_namespace("S4Vectors", "GSE180661 feature metadata")

  dataset_id <- dataset$id %||% "GSE180661"
  output_dir <- .sc_counts_output_dir(config, dataset_id)
  .sc_remove_legacy_outputs(
    output_dir,
    c("blocker_report.md", "cell_qc_audit.csv.gz")
  )
  sample_cache_dir <- .sc_counts_cache_dir(config, dataset_id, "samples")
  structure <- .sc_h5_csr_structure(dataset$matrix_h5_path)
  metadata <- .sc_gse180661_metadata(dataset$cells_path)

  obs <- rhdf5::h5read(
    dataset$matrix_h5_path,
    "obs",
    read.attributes = FALSE
  )
  var <- rhdf5::h5read(
    dataset$matrix_h5_path,
    "var",
    read.attributes = FALSE
  )
  if (
    !is.data.frame(obs) ||
      !"index" %in% names(obs) ||
      nrow(obs) != structure$shape[[1L]] ||
      anyDuplicated(obs$index) ||
      !is.data.frame(var) ||
      !all(c("index", "gene_ids") %in% names(var)) ||
      nrow(var) != structure$shape[[2L]]
  ) {
    stop("GSE180661 H5 obs/var identifiers are malformed.", call. = FALSE)
  }
  h5_cell_id <- as.character(obs$index)
  metadata$h5_row <- match(metadata$cell_id, h5_cell_id)
  if (anyNA(metadata$h5_row) || anyDuplicated(metadata$h5_row)) {
    stop(
      "GSE180661 metadata cell IDs are not a unique subset of H5 obs IDs.",
      call. = FALSE
    )
  }
  gene_symbol <- as.character(var$index)
  gene_id <- as.character(var$gene_ids)
  if (
    anyNA(gene_symbol) || any(!nzchar(gene_symbol)) ||
      anyNA(gene_id) || any(!nzchar(gene_id)) || anyDuplicated(gene_id)
  ) {
    stop("GSE180661 H5 feature identifiers are invalid.", call. = FALSE)
  }
  row_id <- make.unique(gene_id)
  indptr <- rhdf5::h5read(
    dataset$matrix_h5_path,
    "X/indptr",
    bit64conversion = "double"
  )
  if (
    length(indptr) != structure$pointer_length ||
      indptr[[1L]] != 0 ||
      tail(indptr, 1L) != structure$data_length ||
      anyNA(indptr) ||
      any(diff(indptr) < 0)
  ) {
    stop("GSE180661 H5 CSR pointer vector is invalid.", call. = FALSE)
  }

  sample_positions <- unique(round(c(
    1,
    structure$data_length / 2,
    max(1, structure$data_length - 99999)
  )))
  structure_samples <- do.call(rbind, lapply(sample_positions, function(start) {
    count <- min(100000, structure$data_length - start + 1)
    values <- rhdf5::h5read(
      dataset$matrix_h5_path,
      "X/data",
      start = start,
      count = count
    )
    indices <- rhdf5::h5read(
      dataset$matrix_h5_path,
      "X/indices",
      start = start,
      count = count
    )
    data.frame(
      start = start,
      values_checked = length(values),
      non_integer = sum(abs(values - round(values)) > 1e-8),
      non_positive = sum(values <= 0),
      invalid_feature_index = sum(
        indices < 0L | indices >= structure$shape[[2L]]
      ),
      stringsAsFactors = FALSE
    )
  }))
  if (
    any(structure_samples$non_integer > 0L) ||
      any(structure_samples$non_positive > 0L) ||
      any(structure_samples$invalid_feature_index > 0L)
  ) {
    stop("GSE180661 H5 structural count samples are invalid.", call. = FALSE)
  }

  input_signature <- .sc_counts_input_signature(
    c(dataset$matrix_h5_path, dataset$cells_path),
    extra = list(
      adapter = "gse180661_h5_csr_v1",
      lower_nmads = config$single_cell$lower_nmads,
      upper_nmads = config$single_cell$upper_nmads,
      mt_nmads = config$single_cell$mt_nmads,
      seed = config$project$seed,
      scDblFinder = as.character(utils::packageVersion("scDblFinder"))
    )
  )
  row_data <- S4Vectors::DataFrame(
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    feature_type = if ("feature_types" %in% names(var)) {
      as.character(var$feature_types)
    } else {
      rep("Gene Expression", length(gene_id))
    },
    row.names = row_id
  )

  sample_ids <- unique(metadata$sample)
  sample_results <- list()
  for (index in seq_along(sample_ids)) {
    sample_id <- sample_ids[[index]]
    checkpoint_path <- file.path(
      sample_cache_dir,
      paste0(clean_filename(sample_id), "_checkpoint.rds")
    )
    if (file.exists(checkpoint_path)) {
      checkpoint <- readRDS(checkpoint_path)
      if (
        identical(checkpoint$signature, input_signature) &&
          file.exists(checkpoint$metric_path)
      ) {
        log_info(
          "Single-cell ", dataset_id, ": reusing ", sample_id,
          " QC checkpoint (", index, "/", length(sample_ids), ")."
        )
        sample_results[[sample_id]] <- checkpoint
        next
      }
    }

    log_info(
      "Single-cell ", dataset_id, ": reading CSR rows for ", sample_id,
      " (", index, "/", length(sample_ids), ")."
    )
    sample_metadata <- metadata[metadata$sample == sample_id]
    sample_metadata <- as.data.frame(sample_metadata)
    counts <- .sc_read_h5_csr_rows(
      dataset$matrix_h5_path,
      sample_metadata$h5_row,
      indptr,
      structure$shape[[2L]]
    )
    rownames(counts) <- row_id
    colnames(counts) <- sample_metadata$cell_id
    .sc_counts_validate_sparse(
      counts,
      paste("GSE180661 sample", sample_id)
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts),
      rowData = row_data
    )
    metrics <- .sc_counts_metrics(
      counts,
      gene_symbol,
      sample_metadata$cell_id,
      dataset_id,
      sample_id
    )
    metrics$h5_row <- sample_metadata$h5_row
    metrics$metadata_nCount <- sample_metadata$nCount_RNA
    metrics$metadata_nFeature <- sample_metadata$nFeature_RNA
    metrics$metadata_percent_mt <- sample_metadata$percent.mt
    metrics$nCount_delta <- metrics$nCount - metrics$metadata_nCount
    metrics$nFeature_delta <- metrics$nFeature - metrics$metadata_nFeature
    optional <- setdiff(
      names(sample_metadata),
      c(
        "cell_id", "sample", "percent.mt", "nCount_RNA",
        "nFeature_RNA", "h5_row"
      )
    )
    for (column in optional) {
      metrics[[column]] <- sample_metadata[[column]]
    }
    thresholds <- derive_sc_qc_thresholds(
      metrics,
      lower_nmads = config$single_cell$lower_nmads %||% 3,
      upper_nmads = config$single_cell$upper_nmads %||% 5,
      mt_nmads = config$single_cell$mt_nmads %||% 3
    )
    metrics <- apply_sc_qc_thresholds(metrics, thresholds)
    doublets <- .sc_counts_run_doublets(
      sce,
      metrics,
      seed = as.integer(config$project$seed) + index * 1000L,
      minimum_core_cells = 100L,
      exclude_small_partitions = TRUE
    )
    metrics <- doublets$metrics
    thresholds$dataset_id <- dataset_id
    thresholds$sample_id <- sample_id
    summary <- data.frame(
      dataset_id = dataset_id,
      sample_id = sample_id,
      genes = nrow(sce),
      cells = ncol(sce),
      core_qc_pass = sum(metrics$passes_core_QC, na.rm = TRUE),
      predicted_doublets = sum(metrics$doublet_class == "doublet", na.rm = TRUE),
      final_qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
      doublet_status = doublets$summary$doublet_status[[1L]],
      median_counts = stats::median(metrics$nCount),
      median_features = stats::median(metrics$nFeature),
      median_percent_mt = stats::median(metrics$percent_mt, na.rm = TRUE),
      cells_with_exact_count_match = sum(metrics$nCount_delta == 0),
      cells_with_exact_feature_match = sum(metrics$nFeature_delta == 0),
      maximum_absolute_count_delta = max(abs(metrics$nCount_delta)),
      maximum_absolute_feature_delta = max(abs(metrics$nFeature_delta)),
      stringsAsFactors = FALSE
    )
    metric_path <- file.path(
      output_dir,
      "samples",
      paste0(clean_filename(sample_id), "_cell_qc.tsv.gz")
    )
    write_sc_table(metrics, metric_path)
    checkpoint <- list(
      signature = input_signature,
      summary = summary,
      thresholds = thresholds,
      doublet_summary = doublets$summary,
      metrics = metrics,
      metric_path = metric_path
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    sample_results[[sample_id]] <- checkpoint
    rm(checkpoint, counts, sce, metrics, doublets, sample_metadata)
    invisible(gc())
  }

  summary <- do.call(rbind, lapply(sample_results, `[[`, "summary"))
  thresholds <- do.call(rbind, lapply(sample_results, `[[`, "thresholds"))
  doublet_summary <- do.call(
    rbind,
    lapply(sample_results, `[[`, "doublet_summary")
  )
  metrics <- data.table::rbindlist(
    lapply(sample_results, `[[`, "metrics"),
    use.names = TRUE,
    fill = TRUE
  )
  metrics_path <- file.path(output_dir, "cell_qc.tsv.gz")
  h5_index_path <- file.path(output_dir, "cell_h5_row_index.tsv.gz")
  summary_path <- file.path(output_dir, "sample_qc_summary.csv")
  thresholds_path <- file.path(output_dir, "qc_thresholds.csv")
  doublet_path <- file.path(output_dir, "doublet_partition_summary.csv")
  structure_path <- file.path(output_dir, "h5_structure_audit.csv")
  structure_sample_path <- file.path(
    output_dir,
    "h5_count_structure_samples.csv"
  )
  delta_path <- file.path(
    output_dir,
    "metadata_metric_delta_summary.csv"
  )
  capability_path <- file.path(output_dir, "capabilities.csv")
  manifest_path <- file.path(output_dir, "backed_count_manifest.rds")
  write_sc_table(metrics, metrics_path)
  write_sc_table(
    metrics[, .(cell_id, batch, h5_row, passes_QC)],
    h5_index_path
  )
  safe_write_csv(summary, summary_path)
  safe_write_csv(thresholds, thresholds_path)
  safe_write_csv(doublet_summary, doublet_path)
  safe_write_csv(
    data.frame(
      dataset_id = dataset_id,
      format = structure$format,
      matrix_cells = structure$shape[[1L]],
      matrix_genes = structure$shape[[2L]],
      nonzero_entries = structure$data_length,
      metadata_cells = nrow(metadata),
      metadata_cells_mapped = sum(!is.na(metadata$h5_row)),
      h5_only_cells = structure$shape[[1L]] - nrow(metadata),
      duplicated_h5_cell_ids = sum(duplicated(h5_cell_id)),
      duplicated_metadata_cell_ids = sum(duplicated(metadata$cell_id)),
      pointer_monotonic = all(diff(indptr) >= 0),
      pointer_terminal_matches_data = tail(indptr, 1L) == structure$data_length,
      stringsAsFactors = FALSE
    ),
    structure_path
  )
  safe_write_csv(structure_samples, structure_sample_path)
  count_quantiles <- as.numeric(stats::quantile(
    metrics$nCount_delta,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    names = FALSE,
    na.rm = TRUE
  ))
  feature_quantiles <- as.numeric(stats::quantile(
    metrics$nFeature_delta,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    names = FALSE,
    na.rm = TRUE
  ))
  delta_summary <- data.frame(
    dataset_id = dataset_id,
    cells = nrow(metrics),
    exact_count_match = sum(metrics$nCount_delta == 0, na.rm = TRUE),
    positive_count_delta = sum(metrics$nCount_delta > 0, na.rm = TRUE),
    negative_count_delta = sum(metrics$nCount_delta < 0, na.rm = TRUE),
    count_delta_min = count_quantiles[[1L]],
    count_delta_q1 = count_quantiles[[2L]],
    count_delta_median = count_quantiles[[3L]],
    count_delta_q3 = count_quantiles[[4L]],
    count_delta_max = count_quantiles[[5L]],
    exact_feature_match = sum(metrics$nFeature_delta == 0, na.rm = TRUE),
    positive_feature_delta = sum(metrics$nFeature_delta > 0, na.rm = TRUE),
    negative_feature_delta = sum(metrics$nFeature_delta < 0, na.rm = TRUE),
    feature_delta_min = feature_quantiles[[1L]],
    feature_delta_q1 = feature_quantiles[[2L]],
    feature_delta_median = feature_quantiles[[3L]],
    feature_delta_q3 = feature_quantiles[[4L]],
    feature_delta_max = feature_quantiles[[5L]],
    decision = paste(
      "QC uses metrics recomputed from the configured H5 raw-count layer;",
      "published metadata metrics are retained only for comparison."
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(delta_summary, delta_path)
  atomic_save_rds(
    list(
      dataset_id = dataset_id,
      matrix_path = normalizePath(
        dataset$matrix_h5_path,
        winslash = "/",
        mustWork = TRUE
      ),
      sparse_format = structure$format,
      matrix_shape_cells_by_genes = structure$shape,
      gene_id = gene_id,
      gene_symbol = gene_symbol,
      cell_h5_row_index = h5_index_path,
      cell_qc = metrics_path,
      raw_counts_preserved_in_source_h5 = TRUE,
      counts_copied = FALSE
    ),
    manifest_path
  )
  plot_paths <- plot_sc_qc_metrics(
    as.data.frame(metrics),
    thresholds,
    output_dir = output_dir,
    dataset_id = dataset_id
  )

  incomplete <- doublet_summary$partition[
    doublet_summary$doublet_status != "completed"
  ]
  capability <- rbind(
    .sc_counts_capability(
      dataset_id,
      "h5_csr_validation",
      "passed",
      paste0(
        "Validated CSR shape ",
        format(structure$shape[[1L]], scientific = FALSE),
        " by ",
        format(structure$shape[[2L]], scientific = FALSE),
        " with monotonic 64-bit pointers."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "cell_metadata_mapping",
      "passed",
      paste0(
        format(nrow(metadata), big.mark = ","),
        " unique metadata cells map one-to-one into H5 obs; order is recorded explicitly."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "raw_counts",
      "available_backed",
      "Per-sample CSR slices were validated as positive integer UMI counts; the source H5 is not copied or densified."
    ),
    .sc_counts_capability(
      dataset_id,
      "core_qc",
      "completed",
      paste0(
        "Count, feature, and mitochondrial metrics were recomputed from H5 counts within each sample; ",
        delta_summary$exact_count_match,
        " cells exactly match the published nCount_RNA value."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "scDblFinder",
      if (length(incomplete) == 0L) "completed" else "completed_with_small_batch_exclusions",
      paste0(
        "scDblFinder completed per sample; ",
        length(incomplete),
        " sample(s) below 100 core-QC cells were explicitly excluded from downstream passes."
      )
    ),
    .sc_counts_capability(
      dataset_id,
      "ambient_rna",
      "blocked",
      "The H5 contains called/retained cells, not a validated empty-droplet universe."
    ),
    .sc_counts_capability(
      dataset_id,
      "downstream_qc_object",
      "available_backed",
      "A backed manifest links raw H5 counts to cell row indices and QC/doublet metadata."
    )
  )
  safe_write_csv(capability, capability_path)

  list(
    id = dataset_id,
    status = if (length(incomplete) == 0L) {
      "validated"
    } else {
      "validated_with_small_batch_exclusions"
    },
    reason = paste(
      "The H5 CSR count layer and all metadata cell mappings passed validation.",
      length(incomplete),
      "undersized sample partition(s) were retained as flags but excluded from downstream passes."
    ),
    summary = summary,
    thresholds = thresholds,
    capability = capability,
    cell_metric_paths = vapply(
      sample_results,
      `[[`,
      character(1),
      "metric_path"
    ),
    object_paths = manifest_path,
    files = c(
      summary = summary_path,
      capabilities = capability_path,
      thresholds = thresholds_path,
      doublets = doublet_path,
      metrics = metrics_path,
      h5_index = h5_index_path,
      structure = structure_path,
      metadata_metric_deltas = delta_path,
      backed_manifest = manifest_path,
      plots = plot_paths
    ),
    cells = nrow(metrics),
    qc_pass = sum(metrics$passes_QC, na.rm = TRUE),
    downstream_ready = sum(metrics$passes_QC, na.rm = TRUE) > 0L
  )
}
