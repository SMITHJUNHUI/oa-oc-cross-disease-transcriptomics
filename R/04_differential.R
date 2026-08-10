run_limma_differential <- function(expression, group) {
  require_namespace("limma", "differential expression")
  group <- droplevels(factor(group, levels = c("Normal", "Disease")))
  if (!identical(colnames(expression), names(group))) {
    group <- group[colnames(expression)]
  }
  if (anyNA(group) || length(levels(group)) != 2L) {
    stop("Differential analysis requires aligned Normal and Disease groups.", call. = FALSE)
  }

  finite_rows <- rowSums(is.finite(expression)) == ncol(expression)
  variable_rows <- apply(expression, 1L, stats::var, na.rm = TRUE) > 0
  expression <- expression[finite_rows & variable_rows, , drop = FALSE]
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrast <- limma::makeContrasts(Disease - Normal, levels = design)
  fit <- limma::lmFit(expression, design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, robust = TRUE)
  result <- limma::topTable(
    fit,
    coef = 1L,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  result$gene <- rownames(result)
  result <- result[!is.na(result$adj.P.Val), , drop = FALSE]
  result[, c("gene", setdiff(names(result), "gene")), drop = FALSE]
}

plot_volcano <- function(result, fdr, log2_fc, title, path) {
  require_namespace("ggplot2", "volcano plots")
  plot_data <- result
  plot_data$status <- "Not significant"
  significant <- plot_data$adj.P.Val < fdr & abs(plot_data$logFC) >= log2_fc
  plot_data$status[significant & plot_data$logFC > 0] <- "Up"
  plot_data$status[significant & plot_data$logFC < 0] <- "Down"
  plot_data$minus_log10_fdr <- -log10(
    pmax(plot_data$adj.P.Val, .Machine$double.xmin)
  )

  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = logFC, y = minus_log10_fdr, colour = status)
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.0) +
    ggplot2::geom_vline(
      xintercept = c(-log2_fc, log2_fc),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::scale_colour_manual(values = c(
      Down = "#2166AC",
      `Not significant` = "#BDBDBD",
      Up = "#B2182B"
    )) +
    ggplot2::labs(
      title = title,
      x = "log2 fold change",
      y = "-log10 adjusted P value",
      colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
  ggplot2::ggsave(path, figure, width = 7.5, height = 6)
  invisible(path)
}

run_differential_stage <- function(bulk_datasets, config) {
  training <- bulk_datasets[
    vapply(bulk_datasets, function(x) identical(x$role, "train"), logical(1))
  ]
  if (length(training) < 2L) {
    stop("At least two training datasets are required.", call. = FALSE)
  }

  results <- lapply(
    training,
    function(dataset) {
      log_info("Running limma for ", dataset$id, ".")
      result <- run_limma_differential(dataset$expression, dataset$group)
      prefix <- paste(dataset$disease, dataset$id, sep = "_")
      safe_write_csv(
        result,
        file.path(
          config$project$output_dir,
          "tables",
          paste0("DEG_", prefix, "_all.csv")
        )
      )
      plot_volcano(
        result,
        fdr = config$differential$fdr,
        log2_fc = config$differential$log2_fc,
        title = paste(dataset$disease, dataset$id, "Disease vs Normal"),
        path = file.path(
          config$project$output_dir,
          "figures",
          paste0("volcano_", prefix, ".pdf")
        )
      )
      list(
        id = dataset$id,
        disease = dataset$disease,
        table = result
      )
    }
  )
  names(results) <- names(training)
  results
}

select_significant_genes <- function(result, fdr, log2_fc) {
  result$gene[
    result$adj.P.Val < fdr &
      abs(result$logFC) >= log2_fc
  ]
}

derive_shared_candidates <- function(differential_results, config) {
  by_disease <- split(
    differential_results,
    vapply(differential_results, `[[`, character(1), "disease")
  )
  if (!all(c("OA", "OC") %in% names(by_disease))) {
    stop("Both OA and OC differential results are required.", call. = FALSE)
  }
  oa <- by_disease$OA[[1L]]$table
  oc <- by_disease$OC[[1L]]$table

  fdr <- config$differential$fdr
  threshold <- config$differential$log2_fc
  oa_genes <- select_significant_genes(oa, fdr, threshold)
  oc_genes <- select_significant_genes(oc, fdr, threshold)
  shared <- intersect(oa_genes, oc_genes)
  threshold_used <- threshold
  selection_rule <- "primary"

  if (length(shared) == 0L) {
    threshold_used <- config$differential$fallback_log2_fc
    selection_rule <- "predeclared_fallback"
    log_warn(
      "No shared DEGs at |log2FC| >= ", threshold,
      "; applying configured fallback |log2FC| >= ", threshold_used, "."
    )
    oa_genes <- select_significant_genes(oa, fdr, threshold_used)
    oc_genes <- select_significant_genes(oc, fdr, threshold_used)
    shared <- intersect(oa_genes, oc_genes)
  }
  if (length(shared) == 0L) {
    stop("No shared OA-OC DEGs were found under either threshold.", call. = FALSE)
  }

  oa_columns <- oa[match(shared, oa$gene), c(
    "gene", "logFC", "AveExpr", "P.Value", "adj.P.Val"
  )]
  oc_columns <- oc[match(shared, oc$gene), c(
    "gene", "logFC", "AveExpr", "P.Value", "adj.P.Val"
  )]
  names(oa_columns)[-1L] <- paste0(names(oa_columns)[-1L], "_OA")
  names(oc_columns)[-1L] <- paste0(names(oc_columns)[-1L], "_OC")
  table <- merge(oa_columns, oc_columns, by = "gene", sort = FALSE)
  table$directionally_concordant <- sign(table$logFC_OA) == sign(table$logFC_OC)
  table$selection_rule <- selection_rule
  table$log2_fc_threshold <- threshold_used

  if (isTRUE(config$differential$require_directional_concordance)) {
    table <- table[table$directionally_concordant, , drop = FALSE]
    if (nrow(table) == 0L) {
      stop(
        "Shared DEGs exist, but none are directionally concordant.",
        call. = FALSE
      )
    }
  }

  safe_write_csv(
    table,
    file.path(config$project$output_dir, "tables", "shared_OA_OC_DEGs.csv")
  )
  write_utf8(
    table$gene,
    file.path(config$project$output_dir, "tables", "shared_OA_OC_genes.txt")
  )

  summary <- data.frame(
    disease = c("OA", "OC", "shared"),
    significant_genes = c(length(oa_genes), length(oc_genes), nrow(table)),
    fdr = fdr,
    log2_fc_threshold = threshold_used,
    selection_rule = selection_rule,
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    summary,
    file.path(config$project$output_dir, "tables", "DEG_selection_summary.csv")
  )

  log_info(
    "Shared OA-OC candidates: ", nrow(table),
    " (", sum(table$directionally_concordant), " directionally concordant)."
  )
  list(
    genes = table$gene,
    table = table,
    threshold = threshold_used,
    rule = selection_rule
  )
}

