.scd_require_packages <- function() {
  packages <- c(
    "BiocSingular", "bluster", "data.table", "edgeR", "ggplot2",
    "Matrix", "readxl", "rhdf5", "scater", "scran", "scuttle",
    "SingleCellExperiment", "SingleR"
  )
  for (package in packages) {
    require_namespace(package, "single-cell downstream analysis")
  }
  invisible(packages)
}

.scd_output_dir <- function(config, dataset_id = NULL) {
  root <- ensure_dir(file.path(config$project$output_dir, "single_cell_downstream"))
  if (is.null(dataset_id)) root else ensure_dir(file.path(root, dataset_id))
}

.scd_cache_dir <- function(config, dataset_id = NULL) {
  root <- ensure_dir(file.path(
    config$project$cache_dir,
    "single_cell_downstream"
  ))
  if (is.null(dataset_id)) root else ensure_dir(file.path(root, dataset_id))
}

.scd_dataset_config <- function(config, dataset_id) {
  matches <- Filter(
    function(dataset) identical(as.character(dataset$id), dataset_id),
    config$single_cell$datasets
  )
  if (length(matches) != 1L) {
    stop(
      "Single-cell downstream dataset is not uniquely configured: ",
      dataset_id,
      call. = FALSE
    )
  }
  matches[[1L]]
}

.scd_hub_genes <- function(config) {
  path <- file.path(
    config$project$output_dir,
    "tables",
    "final_hub_genes.txt"
  )
  if (!file.exists(path)) {
    stop("Final hub-gene file is missing: ", path, call. = FALSE)
  }
  genes <- unique(toupper(trimws(readLines(path, warn = FALSE))))
  genes[nzchar(genes)]
}

.scd_even_indices <- function(n, maximum) {
  n <- as.integer(n)
  maximum <- as.integer(maximum)
  if (n <= maximum) {
    return(seq_len(n))
  }
  unique(as.integer(round(seq(1, n, length.out = maximum))))
}

.scd_collapse_gene_symbols <- function(counts, symbols, fallback = NULL) {
  require_namespace("Matrix", "sparse gene-symbol collapsing")
  symbols <- as.character(symbols)
  if (is.null(fallback)) {
    fallback <- rownames(counts)
  }
  fallback <- as.character(fallback)
  invalid <- is.na(symbols) | !nzchar(trimws(symbols))
  symbols[invalid] <- fallback[invalid]
  symbols <- toupper(trimws(symbols))
  keep <- !is.na(symbols) & nzchar(symbols)
  counts <- counts[keep, , drop = FALSE]
  symbols <- symbols[keep]
  unique_symbols <- unique(symbols)

  if (!anyDuplicated(symbols)) {
    rownames(counts) <- symbols
    return(methods::as(counts, "dgCMatrix"))
  }

  mapping <- Matrix::sparseMatrix(
    i = seq_along(symbols),
    j = match(symbols, unique_symbols),
    x = 1,
    dims = c(length(symbols), length(unique_symbols))
  )
  collapsed <- Matrix::t(mapping) %*% counts
  rownames(collapsed) <- unique_symbols
  colnames(collapsed) <- colnames(counts)
  methods::as(collapsed, "dgCMatrix")
}

.scd_sce_from_counts <- function(counts, cell_data = NULL) {
  if (is.null(cell_data)) {
    cell_data <- S4Vectors::DataFrame(
      row.names = colnames(counts)
    )
  } else {
    cell_data <- S4Vectors::DataFrame(
      cell_data,
      row.names = colnames(counts)
    )
  }
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = methods::as(counts, "dgCMatrix")),
    colData = cell_data
  )
}

.scd_merge_qc_sce <- function(paths) {
  if (length(paths) == 0L || any(!file.exists(paths))) {
    stop("One or more single-cell QC checkpoints are missing.", call. = FALSE)
  }
  count_pieces <- vector("list", length(paths))
  metadata_pieces <- vector("list", length(paths))
  reference_symbols <- NULL
  reference_ids <- NULL

  for (index in seq_along(paths)) {
    log_info(
      "Single-cell downstream: loading ", basename(paths[[index]]),
      " (", index, "/", length(paths), ")."
    )
    sce <- readRDS(paths[[index]])
    if (!"counts" %in% SummarizedExperiment::assayNames(sce)) {
      stop("QC checkpoint has no counts assay: ", paths[[index]], call. = FALSE)
    }
    pass <- as.logical(SummarizedExperiment::colData(sce)$passes_QC)
    pass[is.na(pass)] <- FALSE
    sce <- sce[, pass, drop = FALSE]
    symbols <- as.character(
      SummarizedExperiment::rowData(sce)$gene_symbol %||% rownames(sce)
    )
    ids <- as.character(
      SummarizedExperiment::rowData(sce)$gene_id %||% rownames(sce)
    )
    if (is.null(reference_symbols)) {
      reference_symbols <- symbols
      reference_ids <- ids
    } else if (
      !identical(symbols, reference_symbols) ||
        !identical(ids, reference_ids)
    ) {
      stop(
        "QC checkpoints do not use an identical feature order.",
        call. = FALSE
      )
    }
    counts <- SummarizedExperiment::assay(sce, "counts")
    count_pieces[[index]] <- methods::as(counts, "dgCMatrix")
    metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    rownames(metadata) <- colnames(counts)
    metadata_pieces[[index]] <- metadata
    rm(sce, counts)
    invisible(gc())
  }

  counts <- do.call(cbind, count_pieces)
  metadata <- data.table::rbindlist(
    metadata_pieces,
    use.names = TRUE,
    fill = TRUE
  )
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- colnames(counts)
  counts <- .scd_collapse_gene_symbols(
    counts,
    reference_symbols,
    reference_ids
  )
  list(
    sce = .scd_sce_from_counts(counts, metadata),
    source_paths = normalizePath(
      paths,
      winslash = "/",
      mustWork = TRUE
    )
  )
}

.scd_select_hvgs <- function(sce, block, n) {
  block <- as.factor(block)
  if (length(levels(block)) >= 2L && all(table(block) >= 2L)) {
    variance_fit <- scran::modelGeneVar(
      sce,
      block = block,
      BPPARAM = BiocParallel::SerialParam()
    )
  } else {
    variance_fit <- scran::modelGeneVar(
      sce,
      BPPARAM = BiocParallel::SerialParam()
    )
  }
  selected <- scran::getTopHVGs(
    variance_fit,
    n = min(as.integer(n), nrow(sce))
  )
  selected[
    !grepl(
      "^(MT-|RPS|RPL|HBA[12]$|HBB$|MALAT1$)",
      selected,
      ignore.case = TRUE
    )
  ]
}

.scd_embed_and_cluster <- function(
    sce,
    block,
    config,
    correct_batch = FALSE
) {
  seed <- as.integer(config$project$seed)
  parameters <- config$single_cell_downstream
  set.seed(seed)
  sce <- scuttle::logNormCounts(sce)
  hvg <- .scd_select_hvgs(sce, block, parameters$hvg_n)
  if (length(hvg) < 50L) {
    stop("Fewer than 50 informative HVGs were selected.", call. = FALSE)
  }
  components <- min(
    as.integer(parameters$pca_n),
    length(hvg) - 1L,
    ncol(sce) - 1L
  )
  sce <- scater::runPCA(
    sce,
    subset_row = hvg,
    ncomponents = components,
    BSPARAM = BiocSingular::IrlbaParam(),
    name = "PCA"
  )
  coordinates <- SingleCellExperiment::reducedDim(sce, "PCA")
  latent_method <- "PCA on batch-aware HVGs; no batch regression"
  if (isTRUE(correct_batch) && length(unique(block)) > 1L) {
    coordinates <- t(limma::removeBatchEffect(
      t(coordinates),
      batch = as.factor(block)
    ))
    latent_method <- paste0(
      "PCA on batch-aware HVGs; sample effect removed in PC space"
    )
  }
  SingleCellExperiment::reducedDim(sce, "PCA_analysis") <- coordinates
  clusters <- bluster::clusterRows(
    coordinates,
    BLUSPARAM = bluster::NNGraphParam(
      k = min(as.integer(parameters$neighbors_k), nrow(coordinates) - 1L),
      cluster.fun = "louvain"
    )
  )
  set.seed(seed)
  umap <- uwot::umap(
    coordinates,
    n_neighbors = min(30L, nrow(coordinates) - 1L),
    min_dist = 0.3,
    n_components = 2L,
    metric = "cosine",
    n_threads = min(4L, parallel::detectCores(logical = TRUE)),
    verbose = FALSE,
    ret_model = FALSE
  )
  colnames(umap) <- c("UMAP1", "UMAP2")
  rownames(umap) <- colnames(sce)
  SingleCellExperiment::reducedDim(sce, "UMAP") <- umap
  SummarizedExperiment::colData(sce)$analysis_cluster <- as.character(clusters)
  list(
    sce = sce,
    hvg = hvg,
    latent_method = latent_method,
    cluster_method = paste0(
      "Louvain SNN, k=", as.integer(parameters$neighbors_k)
    ),
    umap_method = "uwot cosine UMAP; n_neighbors=30; min_dist=0.3"
  )
}

.scd_cluster_markers <- function(sce, clusters, top_n = 25L) {
  clusters <- as.character(clusters)
  levels <- sort(unique(clusters))
  counts <- SummarizedExperiment::assay(sce, "counts")
  logcounts <- SummarizedExperiment::assay(sce, "logcounts")
  means <- vapply(
    levels,
    function(level) {
      Matrix::rowMeans(logcounts[, clusters == level, drop = FALSE])
    },
    numeric(nrow(sce))
  )
  detected <- vapply(
    levels,
    function(level) {
      Matrix::rowMeans(counts[, clusters == level, drop = FALSE] > 0)
    },
    numeric(nrow(sce))
  )
  rownames(means) <- rownames(sce)
  rownames(detected) <- rownames(sce)
  results <- vector("list", length(levels))
  for (index in seq_along(levels)) {
    other <- setdiff(seq_along(levels), index)
    contrast <- if (length(other) == 0L) {
      means[, index]
    } else {
      means[, index] - apply(means[, other, drop = FALSE], 1L, max)
    }
    order_index <- order(
      contrast,
      detected[, index],
      decreasing = TRUE,
      na.last = NA
    )
    order_index <- head(order_index[detected[order_index, index] >= 0.05], top_n)
    results[[index]] <- data.frame(
      cluster = levels[[index]],
      rank = seq_along(order_index),
      gene = rownames(sce)[order_index],
      mean_logcounts = means[order_index, index],
      delta_vs_best_other_cluster = contrast[order_index],
      fraction_detected = detected[order_index, index],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}

.scd_singleR_clusters <- function(test, reference, labels, clusters) {
  cluster_levels <- sort(unique(as.character(clusters)))
  cluster_ids <- factor(as.character(clusters), levels = cluster_levels)
  aggregated_test <- scuttle::aggregateAcrossCells(
    test,
    ids = cluster_ids,
    use.assay.type = "logcounts",
    statistics = "mean",
    use.dimred = FALSE
  )
  colnames(aggregated_test) <- cluster_levels
  prediction <- SingleR::SingleR(
    test = aggregated_test,
    ref = reference,
    labels = as.character(labels),
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts",
    de.method = "classic",
    BPPARAM = BiocParallel::SerialParam()
  )
  prediction <- as.data.frame(prediction)
  prediction$cluster <- rownames(prediction)
  prediction$cell_type <- ifelse(
    is.na(prediction$pruned.labels) | !nzchar(prediction$pruned.labels),
    "Ambiguous",
    prediction$pruned.labels
  )
  prediction$annotation_confidence <- ifelse(
    prediction$cell_type == "Ambiguous",
    "ambiguous",
    "reference_supported"
  )
  prediction
}

.scd_candidate_expression <- function(sce, groups, genes) {
  groups <- as.character(groups)
  present <- intersect(toupper(genes), rownames(sce))
  if (length(present) == 0L) {
    return(data.frame())
  }
  counts <- SummarizedExperiment::assay(sce, "counts")[present, , drop = FALSE]
  logcounts <- SummarizedExperiment::assay(
    sce,
    "logcounts"
  )[present, , drop = FALSE]
  levels <- sort(unique(groups))
  rows <- lapply(levels, function(level) {
    selected <- groups == level
    data.frame(
      group = level,
      gene = present,
      cells = sum(selected),
      mean_logcounts = Matrix::rowMeans(
        logcounts[, selected, drop = FALSE]
      ),
      fraction_detected = Matrix::rowMeans(
        counts[, selected, drop = FALSE] > 0
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.scd_composition <- function(metadata, sample_column, label_column, extras = NULL) {
  columns <- unique(c(sample_column, extras, label_column))
  table_data <- data.table::as.data.table(metadata)[
    ,
    .N,
    by = columns
  ]
  denominator_columns <- unique(c(sample_column, extras))
  table_data[, proportion := N / sum(N), by = denominator_columns]
  data.frame(table_data)
}

.scd_aggregate_counts <- function(counts, groups) {
  groups <- as.character(groups)
  group_levels <- unique(groups)
  mapping <- Matrix::sparseMatrix(
    i = seq_along(groups),
    j = match(groups, group_levels),
    x = 1,
    dims = c(length(groups), length(group_levels))
  )
  aggregated <- counts %*% mapping
  colnames(aggregated) <- group_levels
  rownames(aggregated) <- rownames(counts)
  methods::as(aggregated, "dgCMatrix")
}

.scd_edger_by_cell_type <- function(
    pseudobulk,
    metadata,
    condition_column,
    reference_level,
    target_level,
    contrast_name,
    minimum_cells,
    minimum_replicates,
    pair_column = NULL
) {
  metadata <- as.data.frame(metadata)
  if (!identical(colnames(pseudobulk), rownames(metadata))) {
    stop("Pseudobulk columns and metadata rows are not aligned.", call. = FALSE)
  }
  cell_types <- sort(unique(as.character(metadata$cell_type)))
  result_rows <- list()
  audit_rows <- list()

  for (cell_type in cell_types) {
    selected <- metadata$cell_type == cell_type &
      metadata$n_cells >= as.integer(minimum_cells) &
      metadata[[condition_column]] %in% c(reference_level, target_level)
    cell_metadata <- metadata[selected, , drop = FALSE]
    cell_counts <- pseudobulk[, selected, drop = FALSE]
    if (!is.null(pair_column)) {
      pair_counts <- table(
        cell_metadata[[pair_column]],
        cell_metadata[[condition_column]]
      )
      complete_pairs <- rownames(pair_counts)[
        pair_counts[, reference_level] > 0 &
          pair_counts[, target_level] > 0
      ]
      retained <- cell_metadata[[pair_column]] %in% complete_pairs
      cell_metadata <- cell_metadata[retained, , drop = FALSE]
      cell_counts <- cell_counts[, retained, drop = FALSE]
    }

    condition <- factor(
      cell_metadata[[condition_column]],
      levels = c(reference_level, target_level)
    )
    replicate_counts <- table(condition)
    sufficient <- length(replicate_counts) == 2L &&
      all(replicate_counts >= as.integer(minimum_replicates))
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      contrast = contrast_name,
      cell_type = cell_type,
      reference = reference_level,
      target = target_level,
      reference_replicates = unname(replicate_counts[reference_level] %||% 0L),
      target_replicates = unname(replicate_counts[target_level] %||% 0L),
      paired = !is.null(pair_column),
      status = if (sufficient) "tested" else "insufficient_replicates",
      stringsAsFactors = FALSE
    )
    if (!sufficient) next

    design <- if (is.null(pair_column)) {
      stats::model.matrix(~condition)
    } else {
      pair <- factor(cell_metadata[[pair_column]])
      stats::model.matrix(~pair + condition)
    }
    dge <- edgeR::DGEList(counts = cell_counts)
    keep <- edgeR::filterByExpr(dge, design = design)
    if (sum(keep) < 10L) {
      audit_rows[[length(audit_rows)]]$status <- "insufficient_expressed_genes"
      next
    }
    dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    coefficient <- grep(
      paste0("^condition", make.names(target_level), "$"),
      colnames(design)
    )
    if (length(coefficient) != 1L) {
      stop(
        "Could not identify edgeR coefficient for ", contrast_name,
        " in ", cell_type, ".",
        call. = FALSE
      )
    }
    test <- edgeR::glmQLFTest(fit, coef = coefficient)
    result <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table
    result$gene <- rownames(result)
    result$cell_type <- cell_type
    result$contrast <- contrast_name
    result$reference <- reference_level
    result$target <- target_level
    result_rows[[length(result_rows) + 1L]] <- result
  }

  list(
    results = if (length(result_rows) == 0L) {
      data.frame()
    } else {
      data.table::rbindlist(result_rows, use.names = TRUE, fill = TRUE)
    },
    audit = data.table::rbindlist(
      audit_rows,
      use.names = TRUE,
      fill = TRUE
    )
  )
}

.scd_plot_umap <- function(metadata, label_column, path, title, maximum_cells) {
  plot_data <- data.table::as.data.table(metadata)
  labels <- as.character(plot_data[[label_column]])
  labels[is.na(labels) | !nzchar(labels)] <- "Unassigned"
  plot_data[, .plot_label := labels]
  if (nrow(plot_data) > maximum_cells) {
    selected <- plot_data[
      ,
      .I[.scd_even_indices(.N, max(
        1L,
        floor(maximum_cells / data.table::uniqueN(.plot_label))
      ))],
      by = .plot_label
    ]$V1
    plot_data <- plot_data[selected]
  }
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = UMAP1, y = UMAP2, colour = .plot_label)
  ) +
    ggplot2::geom_point(size = 0.12, alpha = 0.45) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      colour = "Cell type",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.32, "cm")
    )
  ggplot2::ggsave(path, plot, width = 9, height = 6.5, units = "in")
  invisible(path)
}

.scd_plot_composition <- function(
    composition,
    sample_column,
    label_column,
    path,
    title
) {
  plot_data <- as.data.frame(composition)
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[sample_column]],
      y = proportion,
      fill = .data[[label_column]]
    )
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Cell proportion",
      fill = "Cell type"
    ) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 70,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major.x = ggplot2::element_blank()
    )
  ggplot2::ggsave(path, plot, width = 11, height = 6.5, units = "in")
  invisible(path)
}

.scd_analyse_gse104782 <- function(config, hub_genes) {
  dataset_id <- "GSE104782"
  output_dir <- .scd_output_dir(config, dataset_id)
  cache_path <- file.path(
    config$project$cache_dir,
    "single_cell",
    dataset_id,
    "GSE104782_umi_qc_sce.rds"
  )
  if (!file.exists(cache_path)) {
    stop("GSE104782 QC checkpoint is missing.", call. = FALSE)
  }
  published <- as.data.frame(readxl::read_excel(
    config$single_cell_downstream$gse104782_cluster_path
  ))
  if (!all(c("Cell", "Cluster") %in% names(published))) {
    stop("GSE104782 published cluster table is malformed.", call. = FALSE)
  }
  if (anyDuplicated(published$Cell)) {
    stop("GSE104782 published cell IDs are duplicated.", call. = FALSE)
  }

  source <- readRDS(cache_path)
  pass <- as.logical(SummarizedExperiment::colData(source)$passes_QC)
  pass[is.na(pass)] <- FALSE
  source <- source[, pass, drop = FALSE]
  source_metadata <- as.data.frame(SummarizedExperiment::colData(source))
  source_metadata$published_cell_type <- as.character(
    published$Cluster[match(source_metadata$cell_id, published$Cell)]
  )
  source_metadata$published_cell_type[
    is.na(source_metadata$published_cell_type) |
      !nzchar(source_metadata$published_cell_type)
  ] <- "Unassigned"
  counts <- .scd_collapse_gene_symbols(
    SummarizedExperiment::assay(source, "counts"),
    SummarizedExperiment::rowData(source)$gene_symbol,
    SummarizedExperiment::rowData(source)$gene_id
  )
  sce <- .scd_sce_from_counts(counts, source_metadata)
  embedded <- .scd_embed_and_cluster(
    sce,
    block = source_metadata$donor,
    config = config,
    correct_batch = FALSE
  )
  sce <- embedded$sce
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  umap <- SingleCellExperiment::reducedDim(sce, "UMAP")
  metadata$UMAP1 <- umap[, 1L]
  metadata$UMAP2 <- umap[, 2L]
  rownames(metadata) <- colnames(sce)

  composition <- .scd_composition(
    metadata,
    sample_column = "donor",
    label_column = "published_cell_type"
  )
  markers <- .scd_cluster_markers(
    sce,
    metadata$analysis_cluster
  )
  candidate_expression <- .scd_candidate_expression(
    sce,
    metadata$published_cell_type,
    hub_genes
  )
  write_sc_table(metadata, file.path(output_dir, "cell_annotations.tsv.gz"))
  safe_write_csv(composition, file.path(output_dir, "cell_composition.csv"))
  safe_write_csv(markers, file.path(output_dir, "cluster_markers.csv"))
  safe_write_csv(
    candidate_expression,
    file.path(output_dir, "hub_gene_expression_by_cell_type.csv")
  )
  .scd_plot_umap(
    metadata,
    "published_cell_type",
    file.path(output_dir, "UMAP_published_cell_types.pdf"),
    "GSE104782 OA cartilage: published cell types",
    config$single_cell_downstream$visualization_max_cells
  )
  .scd_plot_composition(
    composition,
    "donor",
    "published_cell_type",
    file.path(output_dir, "cell_composition_by_donor.pdf"),
    "GSE104782 cell-type composition by donor"
  )

  reference_keep <- metadata$published_cell_type != "Unassigned"
  reference <- sce[, reference_keep, drop = FALSE]
  reference$reference_label <- metadata$published_cell_type[reference_keep]
  checkpoint <- list(
    dataset_id = dataset_id,
    raw_count_checkpoint = normalizePath(
      cache_path,
      winslash = "/",
      mustWork = TRUE
    ),
    raw_counts_preserved = TRUE,
    published_annotation_path = normalizePath(
      config$single_cell_downstream$gse104782_cluster_path,
      winslash = "/",
      mustWork = TRUE
    ),
    hvg = embedded$hvg,
    latent_method = embedded$latent_method,
    cluster_method = embedded$cluster_method,
    umap_method = embedded$umap_method,
    cell_metadata = metadata[, c(
      "cell_id", "donor", "state_code", "published_cell_type",
      "analysis_cluster", "UMAP1", "UMAP2"
    )]
  )
  checkpoint_path <- file.path(
    .scd_cache_dir(config, dataset_id),
    "downstream_reduced_checkpoint.rds"
  )
  atomic_save_rds(checkpoint, checkpoint_path)
  summary <- data.frame(
    dataset_id = dataset_id,
    disease = "OA",
    analysis_status = "completed",
    cells = ncol(sce),
    annotated_cells = sum(metadata$published_cell_type != "Unassigned"),
    cell_types = length(unique(
      metadata$published_cell_type[
        metadata$published_cell_type != "Unassigned"
      ]
    )),
    annotation_source = "author-published GSE104782 cluster table",
    umap_cells = ncol(sce),
    pseudobulk_contrasts = 0L,
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(output_dir, "summary.csv"))
  list(
    summary = summary,
    reference = reference,
    output_dir = output_dir,
    checkpoint = checkpoint_path
  )
}

.scd_prepare_gse180661_reference <- function(config, hub_genes) {
  dataset_id <- "GSE180661"
  output_dir <- .scd_output_dir(config, dataset_id)
  reference_cache_path <- file.path(
    config$project$cache_dir,
    "single_cell_downstream",
    dataset_id,
    "balanced_reference_sce.rds"
  )
  summary_path <- file.path(output_dir, "summary.csv")
  qc_path <- file.path(
    config$project$output_dir,
    "single_cell",
    dataset_id,
    "cell_qc.tsv.gz"
  )
  manifest_path <- file.path(
    config$project$output_dir,
    "single_cell",
    dataset_id,
    "backed_count_manifest.rds"
  )
  dataset <- .scd_dataset_config(config, dataset_id)
  if (
    file.exists(reference_cache_path) &&
      file.exists(summary_path) &&
      file.info(reference_cache_path)$mtime >=
        max(file.info(c(qc_path, manifest_path))$mtime)
  ) {
    return(list(
      summary = utils::read.csv(summary_path, stringsAsFactors = FALSE),
      reference = readRDS(reference_cache_path),
      output_dir = output_dir,
      source_count_manifest = manifest_path,
      reference_checkpoint = reference_cache_path
    ))
  }
  qc <- data.table::fread(qc_path, showProgress = FALSE)
  qc <- qc[
    passes_QC == TRUE &
      !is.na(cell_type) &
      nzchar(cell_type) &
      !is.na(h5_row)
  ]
  if (nrow(qc) == 0L) {
    stop("GSE180661 has no annotated QC-passing cells.", call. = FALSE)
  }

  source_metadata <- data.table::fread(
    dataset$cells_path,
    select = c(
      "cell_id", "umap50_1", "umap50_2"
    ),
    showProgress = FALSE
  )
  source_index <- match(qc$cell_id, source_metadata$cell_id)
  if (anyNA(source_index)) {
    stop("GSE180661 UMAP metadata do not cover all QC cells.", call. = FALSE)
  }
  annotation_metadata <- data.frame(
    cell_id = qc$cell_id,
    batch = qc$batch,
    cell_type = qc$cell_type,
    cluster_label = qc$cluster_label,
    cluster_label_sub = qc$cluster_label_sub,
    cell_type_super = qc$cell_type_super,
    patient_id = qc$patient_id,
    tumor_subsite = qc$tumor_subsite,
    tumor_site = qc$tumor_site,
    tumor_supersite = qc$tumor_supersite,
    sort_parameters = qc$sort_parameters,
    therapy = qc$therapy,
    surgery = qc$surgery,
    UMAP1 = source_metadata$umap50_1[source_index],
    UMAP2 = source_metadata$umap50_2[source_index],
    stringsAsFactors = FALSE
  )
  write_sc_table(
    annotation_metadata,
    file.path(output_dir, "cell_annotations.tsv.gz")
  )
  composition <- .scd_composition(
    annotation_metadata,
    sample_column = "patient_id",
    label_column = "cell_type",
    extras = "tumor_supersite"
  )
  safe_write_csv(composition, file.path(output_dir, "cell_composition.csv"))
  .scd_plot_umap(
    annotation_metadata,
    "cell_type",
    file.path(output_dir, "UMAP_published_cell_types.pdf"),
    "GSE180661 HGSOC: published Harmony UMAP and cell types",
    config$single_cell_downstream$visualization_max_cells
  )

  reference_cells <- qc[
    order(cell_type, patient_id, h5_row),
    .SD[.scd_even_indices(.N, min(.N, 5L))],
    by = .(cell_type, patient_id)
  ]
  reference_cells <- reference_cells[
    ,
    .SD[.scd_even_indices(.N, min(.N, 200L))],
    by = cell_type
  ]
  if (anyDuplicated(reference_cells$h5_row)) {
    stop("GSE180661 reference sampling selected duplicated H5 rows.", call. = FALSE)
  }
  manifest <- readRDS(manifest_path)
  indptr <- rhdf5::h5read(
    manifest$matrix_path,
    "X/indptr",
    bit64conversion = "double"
  )
  if (
    anyNA(indptr) ||
      any(!is.finite(indptr)) ||
      any(indptr < 0) ||
      any(diff(indptr) < 0) ||
      any(abs(indptr - round(indptr)) > 0)
  ) {
    stop("GSE180661 CSR indptr is invalid after 64-bit-safe read.", call. = FALSE)
  }
  reference_counts <- .sc_read_h5_csr_rows(
    manifest$matrix_path,
    reference_cells$h5_row,
    indptr,
    n_features = length(manifest$gene_symbol),
    maximum_gap = 1000L
  )
  colnames(reference_counts) <- reference_cells$cell_id
  reference_counts <- .scd_collapse_gene_symbols(
    reference_counts,
    manifest$gene_symbol,
    manifest$gene_id
  )
  reference <- .scd_sce_from_counts(
    reference_counts,
    data.frame(
      cell_id = reference_cells$cell_id,
      reference_label = reference_cells$cell_type,
      patient_id = reference_cells$patient_id,
      batch = reference_cells$batch,
      stringsAsFactors = FALSE
    )
  )
  reference <- scuttle::logNormCounts(reference)
  atomic_save_rds(reference, reference_cache_path)
  reference_candidate_expression <- .scd_candidate_expression(
    reference,
    reference$reference_label,
    hub_genes
  )
  safe_write_csv(
    reference_candidate_expression,
    file.path(
      output_dir,
      "hub_gene_expression_in_balanced_reference_sample.csv"
    )
  )
  sampling_manifest <- data.frame(
    cell_id = reference_cells$cell_id,
    h5_row = reference_cells$h5_row,
    cell_type = reference_cells$cell_type,
    patient_id = reference_cells$patient_id,
    batch = reference_cells$batch,
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    sampling_manifest,
    file.path(output_dir, "annotation_reference_sampling.csv")
  )
  summary <- data.frame(
    dataset_id = dataset_id,
    disease = "OC",
    analysis_status = "completed_published_embedding",
    cells = nrow(annotation_metadata),
    annotated_cells = nrow(annotation_metadata),
    cell_types = length(unique(annotation_metadata$cell_type)),
    annotation_source = "author-published GSE180661 metadata",
    umap_cells = nrow(annotation_metadata),
    pseudobulk_contrasts = 0L,
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(output_dir, "summary.csv"))
  list(
    summary = summary,
    reference = reference,
    output_dir = output_dir,
    source_count_manifest = manifest_path,
    reference_checkpoint = reference_cache_path
  )
}

.scd_analyse_gse169454 <- function(config, hub_genes, reference) {
  dataset_id <- "GSE169454"
  output_dir <- .scd_output_dir(config, dataset_id)
  checkpoint_dir <- file.path(
    config$project$cache_dir,
    "single_cell",
    dataset_id
  )
  paths <- sort(list.files(
    checkpoint_dir,
    pattern = "_qc_sce\\.rds$",
    full.names = TRUE
  ))
  merged <- .scd_merge_qc_sce(paths)
  sce <- merged$sce
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  embedded <- .scd_embed_and_cluster(
    sce,
    block = metadata$batch,
    config = config,
    correct_batch = FALSE
  )
  sce <- embedded$sce
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  prediction <- .scd_singleR_clusters(
    sce,
    reference,
    labels = reference$reference_label,
    clusters = metadata$analysis_cluster
  )
  metadata$cell_type <- prediction$cell_type[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  metadata$annotation_delta <- prediction$delta.next[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  metadata$annotation_confidence <- prediction$annotation_confidence[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  umap <- SingleCellExperiment::reducedDim(sce, "UMAP")
  metadata$UMAP1 <- umap[, 1L]
  metadata$UMAP2 <- umap[, 2L]
  rownames(metadata) <- colnames(sce)

  composition <- .scd_composition(
    metadata,
    sample_column = "batch",
    label_column = "cell_type",
    extras = "condition"
  )
  markers <- .scd_cluster_markers(sce, metadata$analysis_cluster)
  candidate_expression <- .scd_candidate_expression(
    sce,
    metadata$cell_type,
    hub_genes
  )
  pseudobulk_key <- paste(metadata$batch, metadata$cell_type, sep = "||")
  pseudobulk <- .scd_aggregate_counts(
    SummarizedExperiment::assay(sce, "counts"),
    pseudobulk_key
  )
  pseudobulk_metadata <- unique(data.frame(
    key = pseudobulk_key,
    sample = metadata$batch,
    condition = metadata$condition,
    cell_type = metadata$cell_type,
    stringsAsFactors = FALSE
  ))
  pseudobulk_metadata$n_cells <- as.integer(table(pseudobulk_key)[
    pseudobulk_metadata$key
  ])
  rownames(pseudobulk_metadata) <- pseudobulk_metadata$key
  pseudobulk_metadata <- pseudobulk_metadata[colnames(pseudobulk), , drop = FALSE]
  differential <- .scd_edger_by_cell_type(
    pseudobulk,
    pseudobulk_metadata,
    condition_column = "condition",
    reference_level = "normal",
    target_level = "OA",
    contrast_name = "OA_vs_normal",
    minimum_cells = config$single_cell_downstream$minimum_pseudobulk_cells,
    minimum_replicates = config$single_cell_downstream$minimum_group_replicates
  )
  significant <- if (nrow(differential$results) == 0L) {
    data.frame()
  } else {
    differential$results[
      differential$results$FDR < 0.05,
      ,
      drop = FALSE
    ]
  }

  write_sc_table(metadata, file.path(output_dir, "cell_annotations.tsv.gz"))
  safe_write_csv(composition, file.path(output_dir, "cell_composition.csv"))
  safe_write_csv(prediction, file.path(output_dir, "cluster_annotation.csv"))
  safe_write_csv(markers, file.path(output_dir, "cluster_markers.csv"))
  safe_write_csv(
    candidate_expression,
    file.path(output_dir, "hub_gene_expression_by_cell_type.csv")
  )
  write_sc_table(
    differential$results,
    file.path(output_dir, "pseudobulk_OA_vs_normal_all.tsv.gz")
  )
  safe_write_csv(
    significant,
    file.path(output_dir, "pseudobulk_OA_vs_normal_FDR05.csv")
  )
  safe_write_csv(
    differential$audit,
    file.path(output_dir, "pseudobulk_contrast_audit.csv")
  )
  atomic_save_rds(
    list(counts = pseudobulk, metadata = pseudobulk_metadata),
    file.path(
      .scd_cache_dir(config, dataset_id),
      "pseudobulk_counts.rds"
    ),
    compress = FALSE
  )
  .scd_plot_umap(
    metadata,
    "cell_type",
    file.path(output_dir, "UMAP_reference_cell_types.pdf"),
    "GSE169454 OA and normal cartilage: reference-supported cell types",
    config$single_cell_downstream$visualization_max_cells
  )
  .scd_plot_composition(
    composition,
    "batch",
    "cell_type",
    file.path(output_dir, "cell_composition_by_sample.pdf"),
    "GSE169454 cell-type composition by sample"
  )
  checkpoint <- list(
    dataset_id = dataset_id,
    source_raw_count_checkpoints = merged$source_paths,
    raw_counts_preserved = TRUE,
    reference_dataset = "GSE104782",
    hvg = embedded$hvg,
    latent_method = embedded$latent_method,
    cluster_method = embedded$cluster_method,
    umap_method = embedded$umap_method,
    cell_metadata = metadata[, c(
      "cell_id", "batch", "condition", "analysis_cluster", "cell_type",
      "annotation_delta", "annotation_confidence", "UMAP1", "UMAP2"
    )]
  )
  checkpoint_path <- file.path(
    .scd_cache_dir(config, dataset_id),
    "downstream_reduced_checkpoint.rds"
  )
  atomic_save_rds(checkpoint, checkpoint_path)
  summary <- data.frame(
    dataset_id = dataset_id,
    disease = "OA",
    analysis_status = "completed",
    cells = ncol(sce),
    annotated_cells = sum(metadata$cell_type != "Ambiguous"),
    cell_types = length(unique(
      metadata$cell_type[metadata$cell_type != "Ambiguous"]
    )),
    annotation_source = "SingleR cluster transfer from GSE104782",
    umap_cells = ncol(sce),
    pseudobulk_contrasts = sum(differential$audit$status == "tested"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(output_dir, "summary.csv"))
  rm(sce, pseudobulk)
  invisible(gc())
  list(
    summary = summary,
    output_dir = output_dir,
    checkpoint = checkpoint_path
  )
}

.scd_analyse_gse154600 <- function(config, hub_genes, reference) {
  dataset_id <- "GSE154600"
  output_dir <- .scd_output_dir(config, dataset_id)
  checkpoint_dir <- file.path(
    config$project$cache_dir,
    "single_cell",
    dataset_id
  )
  paths <- sort(list.files(
    checkpoint_dir,
    pattern = "_qc_sce\\.rds$",
    full.names = TRUE
  ))
  merged <- .scd_merge_qc_sce(paths)
  sce <- merged$sce
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  embedded <- .scd_embed_and_cluster(
    sce,
    block = metadata$batch,
    config = config,
    correct_batch = FALSE
  )
  sce <- embedded$sce
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  prediction <- .scd_singleR_clusters(
    sce,
    reference,
    labels = reference$reference_label,
    clusters = metadata$analysis_cluster
  )
  metadata$cell_type <- prediction$cell_type[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  metadata$annotation_delta <- prediction$delta.next[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  metadata$annotation_confidence <- prediction$annotation_confidence[
    match(metadata$analysis_cluster, prediction$cluster)
  ]
  umap <- SingleCellExperiment::reducedDim(sce, "UMAP")
  metadata$UMAP1 <- umap[, 1L]
  metadata$UMAP2 <- umap[, 2L]
  rownames(metadata) <- colnames(sce)

  composition <- .scd_composition(
    metadata,
    sample_column = "batch",
    label_column = "cell_type"
  )
  markers <- .scd_cluster_markers(sce, metadata$analysis_cluster)
  candidate_expression <- .scd_candidate_expression(
    sce,
    metadata$cell_type,
    hub_genes
  )
  pseudobulk_key <- paste(metadata$batch, metadata$cell_type, sep = "||")
  pseudobulk <- .scd_aggregate_counts(
    SummarizedExperiment::assay(sce, "counts"),
    pseudobulk_key
  )
  pseudobulk_metadata <- unique(data.frame(
    key = pseudobulk_key,
    sample = metadata$batch,
    cell_type = metadata$cell_type,
    stringsAsFactors = FALSE
  ))
  pseudobulk_metadata$n_cells <- as.integer(table(pseudobulk_key)[
    pseudobulk_metadata$key
  ])
  rownames(pseudobulk_metadata) <- pseudobulk_metadata$key
  pseudobulk_metadata <- pseudobulk_metadata[colnames(pseudobulk), , drop = FALSE]

  write_sc_table(metadata, file.path(output_dir, "cell_annotations.tsv.gz"))
  safe_write_csv(composition, file.path(output_dir, "cell_composition.csv"))
  safe_write_csv(prediction, file.path(output_dir, "cluster_annotation.csv"))
  safe_write_csv(markers, file.path(output_dir, "cluster_markers.csv"))
  safe_write_csv(
    candidate_expression,
    file.path(output_dir, "hub_gene_expression_by_cell_type.csv")
  )
  atomic_save_rds(
    list(counts = pseudobulk, metadata = pseudobulk_metadata),
    file.path(
      .scd_cache_dir(config, dataset_id),
      "pseudobulk_counts.rds"
    ),
    compress = FALSE
  )
  .scd_plot_umap(
    metadata,
    "cell_type",
    file.path(output_dir, "UMAP_reference_cell_types.pdf"),
    "GSE154600 HGSOC: GSE180661 reference-supported cell types",
    config$single_cell_downstream$visualization_max_cells
  )
  .scd_plot_composition(
    composition,
    "batch",
    "cell_type",
    file.path(output_dir, "cell_composition_by_tumor.pdf"),
    "GSE154600 cell-type composition by tumor"
  )
  checkpoint <- list(
    dataset_id = dataset_id,
    source_raw_count_checkpoints = merged$source_paths,
    raw_counts_preserved = TRUE,
    reference_dataset = "GSE180661",
    hvg = embedded$hvg,
    latent_method = embedded$latent_method,
    cluster_method = embedded$cluster_method,
    umap_method = embedded$umap_method,
    cell_metadata = metadata[, c(
      "cell_id", "batch", "analysis_cluster", "cell_type",
      "annotation_delta", "annotation_confidence", "UMAP1", "UMAP2"
    )]
  )
  checkpoint_path <- file.path(
    .scd_cache_dir(config, dataset_id),
    "downstream_reduced_checkpoint.rds"
  )
  atomic_save_rds(checkpoint, checkpoint_path)
  summary <- data.frame(
    dataset_id = dataset_id,
    disease = "OC",
    analysis_status = "completed",
    cells = ncol(sce),
    annotated_cells = sum(metadata$cell_type != "Ambiguous"),
    cell_types = length(unique(
      metadata$cell_type[metadata$cell_type != "Ambiguous"]
    )),
    annotation_source = "SingleR cluster transfer from GSE180661",
    umap_cells = ncol(sce),
    pseudobulk_contrasts = 0L,
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(output_dir, "summary.csv"))
  rm(sce, pseudobulk)
  invisible(gc())
  list(
    summary = summary,
    output_dir = output_dir,
    checkpoint = checkpoint_path
  )
}

.scd_merge_pseudobulk_pieces <- function(pieces) {
  if (length(pieces) == 0L) {
    stop("No pseudobulk pieces were produced.", call. = FALSE)
  }
  combined <- do.call(cbind, pieces)
  keys <- colnames(combined)
  unique_keys <- unique(keys)
  if (!anyDuplicated(keys)) {
    return(methods::as(combined, "dgCMatrix"))
  }
  mapping <- Matrix::sparseMatrix(
    i = seq_along(keys),
    j = match(keys, unique_keys),
    x = 1,
    dims = c(length(keys), length(unique_keys))
  )
  aggregated <- combined %*% mapping
  rownames(aggregated) <- rownames(combined)
  colnames(aggregated) <- unique_keys
  methods::as(aggregated, "dgCMatrix")
}

.scd_analyse_gse255460 <- function(config, hub_genes) {
  dataset_id <- "GSE255460"
  output_dir <- .scd_output_dir(config, dataset_id)
  dataset <- .scd_dataset_config(config, dataset_id)
  qc_path <- file.path(
    config$project$output_dir,
    "single_cell",
    dataset_id,
    "cell_qc.tsv.gz"
  )
  qc <- data.table::fread(qc_path, showProgress = FALSE)
  qc <- qc[passes_QC == TRUE]
  qc[, donor := ifelse(trait == "Control", ID, sample)]
  if (
    nrow(qc) == 0L ||
      anyNA(qc$cell_id) ||
      anyNA(qc$celltype) ||
      anyDuplicated(qc$cell_id)
  ) {
    stop("GSE255460 QC-passing annotation table is invalid.", call. = FALSE)
  }

  composition <- .scd_composition(
    qc,
    sample_column = "ID",
    label_column = "celltype",
    extras = c("donor", "trait", "group")
  )
  safe_write_csv(composition, file.path(output_dir, "cell_composition.csv"))

  maximum_cells <- as.integer(
    config$single_cell_downstream$visualization_max_cells
  )
  strata <- qc[, .N, by = .(ID, celltype)]
  per_stratum <- max(25L, floor(maximum_cells / nrow(strata)))
  visualization_cells <- qc[
    order(ID, celltype, cell_id),
    .SD[.scd_even_indices(.N, min(.N, per_stratum))],
    by = .(ID, celltype)
  ]$cell_id
  if (length(visualization_cells) > maximum_cells) {
    visualization_cells <- visualization_cells[
      .scd_even_indices(length(visualization_cells), maximum_cells)
    ]
  }

  metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
  bundle <- .sc_gse255460_ensure_sparse_bundle(dataset, config)
  validated <- .sc_gse255460_validate_manifest(bundle, metadata)
  partition_ids <- as.character(validated$partitions$partition_id)
  visualization_pieces <- list()
  visualization_metadata <- list()
  region_pseudobulk_pieces <- list()
  region_metadata_rows <- list()
  candidate_rows <- list()

  for (index in seq_along(partition_ids)) {
    partition_id <- partition_ids[[index]]
    log_info(
      "Single-cell downstream GSE255460: reading ", partition_id,
      " (", index, "/", length(partition_ids), ")."
    )
    imported <- .sc_read_gse255460_partition(
      bundle,
      validated,
      partition_id,
      metadata
    )
    counts <- imported$counts
    partition_qc <- qc[ID == partition_id]
    cell_index <- match(partition_qc$cell_id, colnames(counts))
    if (anyNA(cell_index)) {
      stop(
        "GSE255460 QC cells do not map to partition ", partition_id,
        ".",
        call. = FALSE
      )
    }
    selected_counts <- counts[, cell_index, drop = FALSE]
    candidate_counts <- NULL
    group_key <- paste(
      partition_qc$donor,
      partition_qc$group,
      partition_qc$celltype,
      sep = "||"
    )
    region_piece <- .scd_aggregate_counts(selected_counts, group_key)
    region_pseudobulk_pieces[[length(region_pseudobulk_pieces) + 1L]] <-
      region_piece
    group_table <- data.table::data.table(
      group_id = group_key,
      donor = partition_qc$donor,
      condition = partition_qc$trait,
      region = partition_qc$group,
      cell_type = partition_qc$celltype
    )[
      ,
      .(n_cells = .N),
      by = .(group_id, donor, condition, region, cell_type)
    ]
    region_metadata_rows[[length(region_metadata_rows) + 1L]] <- group_table

    candidate_present <- intersect(toupper(hub_genes), rownames(selected_counts))
    if (length(candidate_present) > 0L) {
      candidate_counts <- selected_counts[
        candidate_present,
        ,
        drop = FALSE
      ]
      candidate_group <- paste(
        partition_qc$trait,
        partition_qc$group,
        partition_qc$celltype,
        sep = "||"
      )
      for (candidate_level in unique(candidate_group)) {
        selected <- candidate_group == candidate_level
        parts <- strsplit(candidate_level, "\\|\\|", fixed = FALSE)[[1L]]
        candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
          trait = parts[[1L]],
          region = parts[[2L]],
          cell_type = parts[[3L]],
          gene = candidate_present,
          cells = sum(selected),
          total_counts = as.numeric(Matrix::rowSums(
            candidate_counts[, selected, drop = FALSE]
          )),
          detected_cells = as.numeric(Matrix::rowSums(
            candidate_counts[, selected, drop = FALSE] > 0
          )),
          stringsAsFactors = FALSE
        )
      }
    }

    visualization_index <- match(
      intersect(partition_qc$cell_id, visualization_cells),
      colnames(counts)
    )
    visualization_index <- visualization_index[!is.na(visualization_index)]
    if (length(visualization_index) > 0L) {
      visualization_pieces[[length(visualization_pieces) + 1L]] <-
        counts[, visualization_index, drop = FALSE]
      visualization_metadata[[length(visualization_metadata) + 1L]] <-
        partition_qc[
          match(
            colnames(counts)[visualization_index],
            partition_qc$cell_id
          )
        ]
    }
    rm(
      imported, counts, selected_counts, region_piece,
      candidate_counts
    )
    invisible(gc())
  }

  region_pseudobulk <- .scd_merge_pseudobulk_pieces(
    region_pseudobulk_pieces
  )
  region_metadata <- data.table::rbindlist(
    region_metadata_rows,
    use.names = TRUE,
    fill = TRUE
  )[
    ,
    .(n_cells = sum(n_cells)),
    by = .(group_id, donor, condition, region, cell_type)
  ]
  region_metadata <- as.data.frame(region_metadata)
  region_metadata <- region_metadata[
    match(colnames(region_pseudobulk), region_metadata$group_id),
    ,
    drop = FALSE
  ]
  rownames(region_metadata) <- region_metadata$group_id

  donor_key <- paste(
    region_metadata$donor,
    region_metadata$cell_type,
    sep = "||"
  )
  donor_pseudobulk <- .scd_aggregate_counts(
    region_pseudobulk,
    donor_key
  )
  donor_metadata <- data.table::as.data.table(region_metadata)[
    ,
    .(
      condition = unique(condition),
      n_cells = sum(n_cells)
    ),
    by = .(donor, cell_type)
  ]
  donor_metadata[, key := paste(donor, cell_type, sep = "||")]
  donor_metadata <- as.data.frame(donor_metadata)
  donor_metadata <- donor_metadata[
    match(colnames(donor_pseudobulk), donor_metadata$key),
    ,
    drop = FALSE
  ]
  rownames(donor_metadata) <- donor_metadata$key

  oa_vs_control <- .scd_edger_by_cell_type(
    donor_pseudobulk,
    donor_metadata,
    condition_column = "condition",
    reference_level = "Control",
    target_level = "OA",
    contrast_name = "OA_vs_Control",
    minimum_cells = config$single_cell_downstream$minimum_pseudobulk_cells,
    minimum_replicates = config$single_cell_downstream$minimum_group_replicates
  )
  wb_vs_nwb <- .scd_edger_by_cell_type(
    region_pseudobulk,
    region_metadata,
    condition_column = "region",
    reference_level = "NWB",
    target_level = "WB",
    contrast_name = "WB_vs_NWB_paired",
    minimum_cells = config$single_cell_downstream$minimum_pseudobulk_cells,
    minimum_replicates = config$single_cell_downstream$minimum_group_replicates,
    pair_column = "donor"
  )
  differential <- data.table::rbindlist(
    list(oa_vs_control$results, wb_vs_nwb$results),
    use.names = TRUE,
    fill = TRUE
  )
  differential_audit <- data.table::rbindlist(
    list(oa_vs_control$audit, wb_vs_nwb$audit),
    use.names = TRUE,
    fill = TRUE
  )
  significant <- if (nrow(differential) == 0L) {
    data.frame()
  } else {
    differential[differential$FDR < 0.05, , drop = FALSE]
  }

  candidate_expression <- data.table::rbindlist(
    candidate_rows,
    use.names = TRUE,
    fill = TRUE
  )[
    ,
    .(
      cells = sum(cells),
      total_counts = sum(total_counts),
      detected_cells = sum(detected_cells)
    ),
    by = .(trait, region, cell_type, gene)
  ]
  candidate_expression[, mean_umi_per_cell := total_counts / cells]
  candidate_expression[, fraction_detected := detected_cells / cells]

  visualization_counts <- do.call(cbind, visualization_pieces)
  visualization_metadata <- data.table::rbindlist(
    visualization_metadata,
    use.names = TRUE,
    fill = TRUE
  )
  visualization_metadata <- as.data.frame(visualization_metadata)
  rownames(visualization_metadata) <- visualization_metadata$cell_id
  visualization_metadata <- visualization_metadata[
    colnames(visualization_counts),
    ,
    drop = FALSE
  ]
  visualization_sce <- .scd_sce_from_counts(
    visualization_counts,
    visualization_metadata
  )
  embedded <- .scd_embed_and_cluster(
    visualization_sce,
    block = visualization_metadata$ID,
    config = config,
    correct_batch = FALSE
  )
  visualization_sce <- embedded$sce
  visualization_metadata <- as.data.frame(
    SummarizedExperiment::colData(visualization_sce)
  )
  visualization_umap <- SingleCellExperiment::reducedDim(
    visualization_sce,
    "UMAP"
  )
  visualization_metadata$UMAP1 <- visualization_umap[, 1L]
  visualization_metadata$UMAP2 <- visualization_umap[, 2L]

  write_sc_table(
    qc,
    file.path(output_dir, "cell_annotations_all_QC_pass.tsv.gz")
  )
  write_sc_table(
    visualization_metadata,
    file.path(output_dir, "visualization_subsample.tsv.gz")
  )
  safe_write_csv(
    candidate_expression,
    file.path(output_dir, "hub_gene_expression_by_cell_type.csv")
  )
  write_sc_table(
    differential,
    file.path(output_dir, "pseudobulk_differential_all.tsv.gz")
  )
  safe_write_csv(
    significant,
    file.path(output_dir, "pseudobulk_differential_FDR05.csv")
  )
  safe_write_csv(
    differential_audit,
    file.path(output_dir, "pseudobulk_contrast_audit.csv")
  )
  atomic_save_rds(
    list(counts = donor_pseudobulk, metadata = donor_metadata),
    file.path(
      .scd_cache_dir(config, dataset_id),
      "pseudobulk_donor_counts.rds"
    ),
    compress = FALSE
  )
  atomic_save_rds(
    list(counts = region_pseudobulk, metadata = region_metadata),
    file.path(
      .scd_cache_dir(config, dataset_id),
      "pseudobulk_region_counts.rds"
    ),
    compress = FALSE
  )
  .scd_plot_umap(
    visualization_metadata,
    "celltype",
    file.path(output_dir, "UMAP_published_cell_types_subsample.pdf"),
    paste0(
      "GSE255460 OA cartilage: published cell types (n=",
      nrow(visualization_metadata), " stratified cells)"
    ),
    nrow(visualization_metadata)
  )
  .scd_plot_composition(
    composition,
    "ID",
    "celltype",
    file.path(output_dir, "cell_composition_by_sample.pdf"),
    "GSE255460 cell-type composition by cartilage sample"
  )
  checkpoint <- list(
    dataset_id = dataset_id,
    raw_count_manifest = file.path(
      config$project$output_dir,
      "single_cell",
      dataset_id,
      "backed_count_manifest.rds"
    ),
    raw_counts_preserved = TRUE,
    annotation_source = "author-published GSE255460 metadata",
    visualization_sampling = paste0(
      "deterministic stratification by ID and celltype; maximum ",
      maximum_cells, " cells"
    ),
    hvg = embedded$hvg,
    latent_method = embedded$latent_method,
    cluster_method = embedded$cluster_method,
    umap_method = embedded$umap_method,
    visualization_metadata = visualization_metadata[, c(
      "cell_id", "ID", "donor", "trait", "group", "celltype",
      "analysis_cluster", "UMAP1", "UMAP2"
    )]
  )
  checkpoint_path <- file.path(
    .scd_cache_dir(config, dataset_id),
    "downstream_reduced_checkpoint.rds"
  )
  atomic_save_rds(checkpoint, checkpoint_path)
  summary <- data.frame(
    dataset_id = dataset_id,
    disease = "OA",
    analysis_status = "completed_published_annotation_subsampled_umap",
    cells = nrow(qc),
    annotated_cells = nrow(qc),
    cell_types = length(unique(qc$celltype)),
    annotation_source = "author-published GSE255460 metadata",
    umap_cells = nrow(visualization_metadata),
    pseudobulk_contrasts = sum(differential_audit$status == "tested"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(output_dir, "summary.csv"))
  rm(
    visualization_sce, visualization_counts, donor_pseudobulk,
    region_pseudobulk
  )
  invisible(gc())
  list(
    summary = summary,
    output_dir = output_dir,
    checkpoint = checkpoint_path
  )
}

.scd_collect_hub_evidence <- function(config) {
  root <- .scd_output_dir(config)
  specifications <- list(
    list(
      id = "GSE104782", disease = "OA",
      file = "hub_gene_expression_by_cell_type.csv",
      metric = "mean_logcounts", context = "published_cell_type"
    ),
    list(
      id = "GSE169454", disease = "OA",
      file = "hub_gene_expression_by_cell_type.csv",
      metric = "mean_logcounts", context = "transferred_cell_type"
    ),
    list(
      id = "GSE255460", disease = "OA",
      file = "hub_gene_expression_by_cell_type.csv",
      metric = "mean_umi_per_cell", context = "trait_region_cell_type"
    ),
    list(
      id = "GSE154600", disease = "OC",
      file = "hub_gene_expression_by_cell_type.csv",
      metric = "mean_logcounts", context = "transferred_cell_type"
    ),
    list(
      id = "GSE180661", disease = "OC",
      file = "hub_gene_expression_in_balanced_reference_sample.csv",
      metric = "mean_logcounts", context = "balanced_reference_cell_type"
    )
  )
  rows <- lapply(specifications, function(specification) {
    path <- file.path(root, specification$id, specification$file)
    if (!file.exists(path)) return(NULL)
    table_data <- data.table::fread(path, showProgress = FALSE)
    if (specification$id == "GSE255460") {
      table_data[, context_value := paste(trait, region, sep = ":")]
      table_data[, cell_type_value := cell_type]
    } else {
      table_data[, context_value := specification$context]
      table_data[, cell_type_value := group]
    }
    data.frame(
      dataset_id = specification$id,
      disease = specification$disease,
      context = table_data$context_value,
      cell_type = table_data$cell_type_value,
      gene = table_data$gene,
      cells = table_data$cells,
      expression_metric = specification$metric,
      expression_value = table_data[[specification$metric]],
      fraction_detected = table_data$fraction_detected,
      stringsAsFactors = FALSE
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.scd_write_report <- function(summaries, manifest, config) {
  report_path <- file.path(
    config$project$output_dir,
    "reports",
    "single_cell_downstream_report.md"
  )
  lines <- c(
    "# Single-cell downstream analysis report",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0(
      "- Scope: clustering, conservative annotation, UMAP, composition, ",
      "hub-gene localization, and replicate-aware pseudobulk statistics."
    ),
    paste0(
      "- OA and OC were analysed separately; no cross-disease latent ",
      "integration was performed."
    ),
    "",
    "## Dataset completion",
    "",
    paste0(
      "| Dataset | Disease | Status | QC-pass cells | Annotated | ",
      "Cell types | UMAP cells | Pseudobulk tests |"
    ),
    "|---|---|---|---:|---:|---:|---:|---:|"
  )
  for (index in seq_len(nrow(summaries))) {
    row <- summaries[index, , drop = FALSE]
    lines <- c(lines, paste0(
      "| ", row$dataset_id,
      " | ", row$disease,
      " | ", row$analysis_status,
      " | ", format(row$cells, big.mark = ",", scientific = FALSE),
      " | ", format(row$annotated_cells, big.mark = ",", scientific = FALSE),
      " | ", row$cell_types,
      " | ", format(row$umap_cells, big.mark = ",", scientific = FALSE),
      " | ", row$pseudobulk_contrasts,
      " |"
    ))
  }
  lines <- c(
    lines,
    "",
    "## Analysis decisions",
    "",
    "- Raw counts remain in the existing per-sample QC checkpoints or backed CSR/H5 manifests; downstream checkpoints store reduced coordinates and provenance instead of copying all count matrices.",
    "- GSE104782 and GSE255460 use author-published cartilage labels.",
    "- GSE169454 clusters are annotated against the same-tissue GSE104782 reference with SingleR; pruned calls remain `Ambiguous`.",
    "- GSE180661 uses author-published HGSOC labels and Harmony UMAP coordinates.",
    "- GSE154600 clusters are annotated against a patient-balanced GSE180661 HGSOC reference with SingleR.",
    "- Standard sparse log-normalization, batch-aware HVGs, randomized SVD PCA, Louvain SNN clustering, and cosine UMAP were selected. scVI was skipped because cross-disease integration is not the biological target and the 32 GB local system should avoid duplicating million-cell objects.",
    "- Pseudobulk differential expression uses biological sample/donor replicates with edgeR quasi-likelihood models. GSE255460 WB versus NWB is paired by OA donor.",
    "- FDR values are Benjamini-Hochberg adjusted within each cell-type contrast. Descriptive cluster-marker tables are not inferential tests.",
    "",
    "## Provenance",
    "",
    paste0(
      "- Parameter manifest: `results/single_cell_downstream/",
      "analysis_manifest.csv`"
    ),
    paste0(
      "- Combined hub-gene evidence: `results/tables/",
      "single_cell_hub_gene_evidence.csv`"
    ),
    paste0(
      "- Dataset outputs: `results/single_cell_downstream/<GSE_ID>/`"
    )
  )
  write_utf8(lines, report_path)
  report_path
}

run_single_cell_downstream_stage <- function(
    single_cell_gate,
    machine_learning,
    config
) {
  .scd_require_packages()
  gate_status_path <- file.path(
    config$project$output_dir,
    "tables",
    "single_cell_dataset_status.csv"
  )
  gate_status <- data.table::fread(gate_status_path, showProgress = FALSE)
  if (
    nrow(gate_status) != 5L ||
      any(as.character(gate_status$downstream_ready) != "TRUE")
  ) {
    stop(
      "Single-cell downstream analysis requires all five QC gates.",
      call. = FALSE
    )
  }
  hub_genes <- .scd_hub_genes(config)
  log_info(
    "Single-cell downstream: starting five-dataset analysis with ",
    length(hub_genes), " pre-specified hub genes."
  )

  gse104782 <- .scd_analyse_gse104782(config, hub_genes)
  gse169454 <- .scd_analyse_gse169454(
    config,
    hub_genes,
    gse104782$reference
  )
  gse104782$reference <- NULL
  invisible(gc())

  gse180661 <- .scd_prepare_gse180661_reference(config, hub_genes)
  gse154600 <- .scd_analyse_gse154600(
    config,
    hub_genes,
    gse180661$reference
  )
  gse180661$reference <- NULL
  invisible(gc())

  gse255460 <- .scd_analyse_gse255460(config, hub_genes)
  summaries <- data.table::rbindlist(
    list(
      gse104782$summary,
      gse169454$summary,
      gse255460$summary,
      gse154600$summary,
      gse180661$summary
    ),
    use.names = TRUE,
    fill = TRUE
  )
  summary_path <- file.path(
    .scd_output_dir(config),
    "single_cell_downstream_summary.csv"
  )
  safe_write_csv(summaries, summary_path)

  manifest <- data.frame(
    parameter = c(
      "analysis_scope",
      "disease_integration",
      "raw_counts",
      "normalization",
      "hvg_selection",
      "pca",
      "clustering",
      "umap",
      "annotation_oa_reference",
      "annotation_oc_reference",
      "annotation_policy",
      "pseudobulk_method",
      "pseudobulk_minimum_cells",
      "pseudobulk_minimum_group_replicates",
      "multiple_testing",
      "scvi",
      "random_seed"
    ),
    value = c(
      "post-QC clustering, annotation, UMAP, composition, hub genes, pseudobulk",
      "OA and OC analysed separately",
      "preserved in QC checkpoints or backed CSR/H5 manifests",
      "scuttle::logNormCounts",
      paste0(
        "scran::modelGeneVar by technical partition; top ",
        config$single_cell_downstream$hvg_n
      ),
      paste0(
        "BiocSingular randomized SVD; ",
        config$single_cell_downstream$pca_n, " components"
      ),
      paste0(
        "bluster Louvain SNN; k=",
        config$single_cell_downstream$neighbors_k
      ),
      "uwot cosine; n_neighbors=30; min_dist=0.3",
      "GSE104782 author-published cartilage cell types",
      "GSE180661 author-published HGSOC cell types",
      "SingleR cluster-level transfer; pruned calls are Ambiguous",
      "edgeR quasi-likelihood on donor/sample pseudobulk counts",
      config$single_cell_downstream$minimum_pseudobulk_cells,
      config$single_cell_downstream$minimum_group_replicates,
      "Benjamini-Hochberg FDR within each cell-type contrast",
      "skipped: no cross-disease integration target; local memory constraint",
      config$project$seed
    ),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(
    .scd_output_dir(config),
    "analysis_manifest.csv"
  )
  safe_write_csv(manifest, manifest_path)
  hub_evidence <- .scd_collect_hub_evidence(config)
  hub_evidence_path <- file.path(
    config$project$output_dir,
    "tables",
    "single_cell_hub_gene_evidence.csv"
  )
  safe_write_csv(hub_evidence, hub_evidence_path)
  report_path <- .scd_write_report(summaries, manifest, config)

  list(
    status = "completed",
    datasets_completed = as.character(summaries$dataset_id),
    cells_analysed = sum(summaries$cells),
    summary = summary_path,
    manifest = manifest_path,
    hub_evidence = hub_evidence_path,
    report = report_path
  )
}
