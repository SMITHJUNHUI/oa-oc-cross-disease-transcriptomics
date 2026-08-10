write_enrichment_result <- function(result, path) {
  table <- as.data.frame(result)
  safe_write_csv(table, path)
  table
}

run_gsea_for_dataset <- function(
    differential,
    gene_sets,
    settings,
    output_dir
) {
  ranked <- differential$table$logFC
  names(ranked) <- differential$table$gene
  ranked <- sort(ranked[!duplicated(names(ranked))], decreasing = TRUE)
  outputs <- list()

  for (name in names(gene_sets)) {
    path <- gene_sets[[name]]
    term_to_gene <- clusterProfiler::read.gmt(path)
    result <- clusterProfiler::GSEA(
      geneList = ranked,
      TERM2GENE = term_to_gene,
      pvalueCutoff = as.numeric(settings$gsea_pvalue_cutoff %||% 0.25),
      minGSSize = as.integer(settings$gsea_min_size %||% 10L),
      maxGSSize = as.integer(settings$gsea_max_size %||% 500L),
      pAdjustMethod = settings$p_adjust_method %||% "BH",
      verbose = FALSE,
      seed = TRUE
    )
    prefix <- paste(differential$disease, differential$id, name, sep = "_")
    table <- write_enrichment_result(
      result,
      file.path(output_dir, "tables", paste0("GSEA_", prefix, ".csv"))
    )
    if (nrow(table) > 0L) {
      figure <- enrichplot::dotplot(result, showCategory = min(15L, nrow(table))) +
        ggplot2::labs(title = paste("GSEA", differential$disease, name)) +
        ggplot2::theme_minimal(base_size = 10)
      ggplot2::ggsave(
        file.path(output_dir, "figures", paste0("GSEA_", prefix, "_dotplot.pdf")),
        figure,
        width = 9,
        height = 7
      )
    }
    outputs[[name]] <- result
  }
  outputs
}

run_enrichment_stage <- function(shared, differential_results, config) {
  if (!isTRUE(config$modules$enrichment)) {
    return(list(status = "disabled"))
  }
  require_namespace("clusterProfiler", "GO, KEGG and GSEA enrichment")
  require_namespace("org.Hs.eg.db", "human gene identifiers")
  require_namespace("enrichplot", "enrichment plots")

  mapping <- clusterProfiler::bitr(
    shared$genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  mapping <- mapping[!duplicated(mapping$SYMBOL), , drop = FALSE]
  safe_write_csv(
    mapping,
    file.path(config$project$output_dir, "tables", "shared_gene_entrez_mapping.csv")
  )
  if (nrow(mapping) == 0L) {
    stop("None of the shared genes mapped to Entrez IDs.", call. = FALSE)
  }

  settings <- config$enrichment
  go <- clusterProfiler::enrichGO(
    gene = mapping$ENTREZID,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = settings$p_adjust_method %||% "BH",
    pvalueCutoff = as.numeric(settings$pvalue_cutoff %||% 0.05),
    qvalueCutoff = as.numeric(settings$qvalue_cutoff %||% 0.20),
    readable = TRUE
  )
  kegg <- clusterProfiler::enrichKEGG(
    gene = mapping$ENTREZID,
    organism = "hsa",
    pAdjustMethod = settings$p_adjust_method %||% "BH",
    pvalueCutoff = as.numeric(settings$pvalue_cutoff %||% 0.05),
    qvalueCutoff = as.numeric(settings$qvalue_cutoff %||% 0.20)
  )
  go_table <- write_enrichment_result(
    go,
    file.path(config$project$output_dir, "tables", "GO_shared_genes.csv")
  )
  kegg_table <- write_enrichment_result(
    kegg,
    file.path(config$project$output_dir, "tables", "KEGG_shared_genes.csv")
  )

  if (nrow(go_table) > 0L) {
    figure <- enrichplot::dotplot(go, showCategory = min(20L, nrow(go_table))) +
      ggplot2::theme_minimal(base_size = 10)
    ggplot2::ggsave(
      file.path(config$project$output_dir, "figures", "GO_shared_genes_dotplot.pdf"),
      figure,
      width = 10,
      height = 8
    )
  }
  if (nrow(kegg_table) > 0L) {
    figure <- enrichplot::dotplot(kegg, showCategory = min(20L, nrow(kegg_table))) +
      ggplot2::theme_minimal(base_size = 10)
    ggplot2::ggsave(
      file.path(config$project$output_dir, "figures", "KEGG_shared_genes_dotplot.pdf"),
      figure,
      width = 10,
      height = 8
    )
  }

  gsea <- lapply(
    differential_results,
    run_gsea_for_dataset,
    gene_sets = config$gene_sets,
    settings = settings,
    output_dir = config$project$output_dir
  )
  names(gsea) <- names(differential_results)

  list(go = go, kegg = kegg, gsea = gsea, mapping = mapping)
}

evaluate_gene_roc <- function(expression, group, gene, training_logfc) {
  if (
    length(training_logfc) != 1L ||
      !is.finite(training_logfc) ||
      training_logfc == 0
  ) {
    stop(
      "A finite, non-zero training logFC is required for direction-fixed ROC.",
      call. = FALSE
    )
  }
  expected_sign <- sign(training_logfc)
  predictor <- expected_sign * as.numeric(expression[gene, ])
  roc <- pROC::roc(
    response = group,
    predictor = predictor,
    levels = c("Normal", "Disease"),
    direction = "<",
    quiet = TRUE
  )
  interval <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  list(
    roc = roc,
    row = data.frame(
      gene = gene,
      auc = as.numeric(pROC::auc(roc)),
      ci_lower = interval[[1L]],
      ci_median = interval[[2L]],
      ci_upper = interval[[3L]],
      direction = roc$direction,
      training_logfc = training_logfc,
      expected_expression = ifelse(
        training_logfc > 0,
        "higher_in_disease",
        "lower_in_disease"
      ),
      stringsAsFactors = FALSE
    )
  )
}

validate_dataset_roc <- function(dataset, genes, training_logfc, output_dir) {
  require_namespace("pROC", "external diagnostic validation")
  present <- intersect(genes, rownames(dataset$expression))
  present <- intersect(present, names(training_logfc))
  if (length(present) == 0L) {
    log_warn(
      dataset$id,
      " contains no selected hub genes with a finite training direction."
    )
    return(list(
      dataset_id = dataset$id,
      disease = dataset$disease,
      table = data.frame(),
      rocs = list()
    ))
  }

  evaluations <- lapply(
    present,
    function(gene) {
      evaluate_gene_roc(
        dataset$expression,
        dataset$group,
        gene,
        training_logfc[[gene]]
      )
    }
  )
  result_table <- do.call(rbind, lapply(evaluations, `[[`, "row"))
  result_table$dataset_id <- dataset$id
  result_table$disease <- dataset$disease
  result_table$n_normal <- unname(base::table(dataset$group)["Normal"])
  result_table$n_disease <- unname(base::table(dataset$group)["Disease"])
  result_table <- result_table[, c(
    "dataset_id", "disease", "gene", "auc", "ci_lower", "ci_median",
    "ci_upper", "direction", "training_logfc", "expected_expression",
    "n_normal", "n_disease"
  )]
  rocs <- lapply(evaluations, `[[`, "roc")
  names(rocs) <- present

  safe_write_csv(
    result_table,
    file.path(
      output_dir,
      "tables",
      paste0("ROC_", dataset$disease, "_", dataset$id, ".csv")
    )
  )
  figure <- pROC::ggroc(rocs, legacy.axes = TRUE, linewidth = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = paste(
        "Direction-fixed external validation:",
        dataset$disease,
        dataset$id
      ),
      x = "1 - specificity",
      y = "sensitivity",
      colour = "Gene"
    ) +
    ggplot2::theme_minimal(base_size = 10)
  ggplot2::ggsave(
    file.path(
      output_dir,
      "figures",
      paste0("ROC_", dataset$disease, "_", dataset$id, ".pdf")
    ),
    figure,
    width = 7,
    height = 6
  )

  list(
    dataset_id = dataset$id,
    disease = dataset$disease,
    table = result_table,
    rocs = rocs
  )
}

run_validation_stage <- function(bulk_datasets, ml_results, shared, config) {
  if (!isTRUE(config$modules$validation)) {
    return(list(status = "disabled"))
  }
  genes <- ml_results$final_genes %||% shared$genes
  validation <- bulk_datasets[
    vapply(
      bulk_datasets,
      function(x) identical(x$role, "validation"),
      logical(1)
    )
  ]
  results <- lapply(
    validation,
    function(dataset) {
      logfc_column <- paste0("logFC_", dataset$disease)
      if (!logfc_column %in% names(shared$table)) {
        stop(
          "Shared DEG table is missing ", logfc_column,
          " for direction-fixed validation.",
          call. = FALSE
        )
      }
      training_logfc <- shared$table[[logfc_column]]
      names(training_logfc) <- shared$table$gene
      training_logfc <- training_logfc[
        is.finite(training_logfc) & training_logfc != 0
      ]
      validate_dataset_roc(
        dataset,
        genes = genes,
        training_logfc = training_logfc,
        output_dir = config$project$output_dir
      )
    }
  )
  names(results) <- names(validation)
  combined <- do.call(
    rbind,
    lapply(results, function(x) x$table)
  )
  if (is.null(combined) || nrow(combined) == 0L) {
    stop("No external validation AUC values could be calculated.", call. = FALSE)
  }
  safe_write_csv(
    combined,
    file.path(config$project$output_dir, "tables", "external_validation_AUC.csv")
  )

  auc_matrix <- stats::xtabs(auc ~ gene + dataset_id, data = combined)
  grDevices::pdf(
    file.path(config$project$output_dir, "figures", "external_validation_AUC_heatmap.pdf"),
    width = 8,
    height = max(5, 0.35 * nrow(auc_matrix) + 2)
  )
  pheatmap::pheatmap(
    as.matrix(auc_matrix),
    cluster_rows = nrow(auc_matrix) > 1L,
    cluster_cols = ncol(auc_matrix) > 1L,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = "External validation AUC",
    color = grDevices::colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(50)
  )
  grDevices::dev.off()

  list(genes = genes, datasets = results, table = combined)
}
