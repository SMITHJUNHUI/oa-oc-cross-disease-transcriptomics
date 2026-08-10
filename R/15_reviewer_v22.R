v22_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v22"))
  figures <- ensure_dir(file.path(root, "figures"))
  source <- ensure_dir(file.path(figures, "source_data"))
  tables <- ensure_dir(file.path(root, "supplementary_tables"))
  analysis <- ensure_dir(file.path(root, "analysis"))
  list(
    root = root,
    figures = figures,
    source = source,
    tables = tables,
    analysis = analysis
  )
}

v22_copy_tree <- function(source, target, skip = character()) {
  source <- normalizePath(source, winslash = "/", mustWork = TRUE)
  target <- ensure_dir(target)
  files <- list.files(
    source,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  files <- files[file.info(files)$isdir %in% FALSE]
  if (length(skip) > 0L) {
    files <- files[!basename(files) %in% skip]
  }
  relative <- substring(files, nchar(source) + 2L)
  destinations <- file.path(target, relative)
  invisible(lapply(unique(dirname(destinations)), ensure_dir))
  copied <- file.copy(files, destinations, overwrite = TRUE, copy.date = TRUE)
  if (!all(copied)) {
    stop("V2.2 could not copy every V2.1 baseline output.", call. = FALSE)
  }
  invisible(destinations)
}

v22_peak_frequency <- function(frequency, disease, genes) {
  subset <- frequency[
    frequency$disease == disease & frequency$gene %in% genes,
    ,
    drop = FALSE
  ]
  do.call(rbind, lapply(genes, function(gene) {
    table <- subset[subset$gene == gene, , drop = FALSE]
    if (nrow(table) == 0L) {
      return(data.frame(
        gene = gene,
        frequency = 0,
        peak_model = "not selected",
        stringsAsFactors = FALSE
      ))
    }
    maximum <- max(table$selection_frequency, na.rm = TRUE)
    models <- paste(
      sort(unique(table$model[table$selection_frequency == maximum])),
      collapse = "; "
    )
    data.frame(
      gene = gene,
      frequency = maximum,
      peak_model = models,
      stringsAsFactors = FALSE
    )
  }))
}

v22_prepare_candidate_matrix <- function(project_root, paths) {
  baseline <- file.path(project_root, "results", "submission_v21")
  shared <- utils::read.csv(file.path(
    baseline,
    "supplementary_tables",
    "Table_S2_shared_differentially_expressed_genes.csv"
  ))
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  genes <- ml$final_genes
  shared <- shared[match(genes, shared$gene), , drop = FALSE]
  if (anyNA(shared$gene)) {
    stop("Candidate genes do not align with the shared DEG table.", call. = FALSE)
  }

  oa_modules <- utils::read.csv(file.path(
    project_root,
    "results",
    "tables",
    "WGCNA_OA_GSE114007_gene_modules.csv"
  ))
  oc_modules <- utils::read.csv(file.path(
    project_root,
    "results",
    "tables",
    "WGCNA_OC_GSE18520_gene_modules.csv"
  ))
  original_evidence <- ml$evidence[match(genes, ml$evidence$gene), , drop = FALSE]

  disease_consensus <- split(
    ml$datasets,
    vapply(ml$datasets, `[[`, character(1), "disease")
  )
  disease_consensus <- lapply(
    disease_consensus,
    function(results) unique(unlist(lapply(results, `[[`, "consensus")))
  )
  cross_disease_consensus <- Reduce(intersect, disease_consensus)
  vote_count <- vapply(genes, function(gene) {
    sum(vapply(ml$datasets, function(result) {
      as.integer(gene %in% result$lasso_genes) +
        as.integer(gene %in% result$rf_genes)
    }, integer(1)))
  }, integer(1))

  frequency <- utils::read.csv(file.path(
    baseline,
    "figures",
    "source_data",
    "Figure3_ML_selection_frequency.csv"
  ))
  oa_frequency <- v22_peak_frequency(frequency, "OA", genes)
  oc_frequency <- v22_peak_frequency(frequency, "OC", genes)

  contexts <- utils::read.csv(file.path(
    baseline,
    "analysis",
    "single_cell_gene_disease_context_matrix.csv"
  ))
  context_row <- function(gene, disease) {
    row <- contexts[
      contexts$gene == gene & contexts$disease == disease,
      ,
      drop = FALSE
    ]
    if (nrow(row) != 1L) {
      stop("Missing unique single-cell context for ", disease, " ", gene, ".")
    }
    row
  }

  pseudobulk <- utils::read.csv(file.path(
    baseline,
    "figures",
    "source_data",
    "Figure5_hub_pseudobulk_evidence.csv"
  ))
  pseudobulk_summary <- lapply(genes, function(gene) {
    table <- pseudobulk[pseudobulk$gene == gene, , drop = FALSE]
    data.frame(
      gene = gene,
      significant_oa_pseudobulk_contrasts = nrow(table),
      minimum_oa_pseudobulk_fdr = if (nrow(table) > 0L) {
        min(table$FDR, na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  pseudobulk_summary <- do.call(rbind, pseudobulk_summary)

  matrix <- data.frame(
    prioritization_rank = seq_along(genes),
    gene = genes,
    shared_primary_DEG_in_OA = TRUE,
    shared_primary_DEG_in_OC = TRUE,
    log2FC_OA = shared$logFC_OA,
    log2FC_OC = shared$logFC_OC,
    direction_class = ifelse(
      shared$directionally_concordant,
      "concordant",
      "discordant"
    ),
    combined_DEG_evidence_score =
      -log10(pmax(shared$adj.P.Val_OA, .Machine$double.xmin)) +
      -log10(pmax(shared$adj.P.Val_OC, .Machine$double.xmin)) +
      abs(shared$logFC_OA) +
      abs(shared$logFC_OC),
    OA_primary_WGCNA_module = genes %in%
      oa_modules$gene[oa_modules$module == "green"],
    OC_primary_WGCNA_module = genes %in%
      oc_modules$gene[oc_modules$module == "brown"],
    original_OA_model_consensus = original_evidence$selected_by_OA,
    original_OC_model_consensus = original_evidence$selected_by_OC,
    original_model_vote_count = vote_count,
    selection_stage = ifelse(
      genes %in% cross_disease_consensus,
      "cross-disease model consensus",
      "ranked multi-model vote completion"
    ),
    strict_nested_OA_max_selection_frequency =
      oa_frequency$frequency[match(genes, oa_frequency$gene)],
    strict_nested_OA_peak_model =
      oa_frequency$peak_model[match(genes, oa_frequency$gene)],
    strict_nested_OC_max_selection_frequency =
      oc_frequency$frequency[match(genes, oc_frequency$gene)],
    strict_nested_OC_peak_model =
      oc_frequency$peak_model[match(genes, oc_frequency$gene)],
    stringsAsFactors = FALSE
  )

  matrix$OA_top_cell_context <- vapply(
    genes,
    function(gene) context_row(gene, "OA")$top_cell_context,
    character(1)
  )
  matrix$OA_top_context_detection_fraction <- vapply(
    genes,
    function(gene) context_row(gene, "OA")$fraction_detected,
    numeric(1)
  )
  matrix$OC_top_cell_context <- vapply(
    genes,
    function(gene) context_row(gene, "OC")$top_cell_context,
    character(1)
  )
  matrix$OC_top_context_detection_fraction <- vapply(
    genes,
    function(gene) context_row(gene, "OC")$fraction_detected,
    numeric(1)
  )
  matrix <- merge(
    matrix,
    pseudobulk_summary,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  matrix <- matrix[match(genes, matrix$gene), , drop = FALSE]
  matrix$fixed_score_member <- TRUE
  matrix$interpretation <- paste0(
    "Prioritized molecular candidate; not a validated diagnostic, ",
    "prognostic, therapeutic, or mechanistic marker"
  )

  safe_write_csv(
    matrix,
    file.path(
      paths$tables,
      "Table_S16_candidate_prioritization_matrix.csv"
    )
  )
  safe_write_csv(
    matrix,
    file.path(
      paths$source,
      "SupplementaryFigure6_candidate_prioritization_matrix.csv"
    )
  )
  matrix
}

v22_read_annotations <- function(project_root, dataset_id) {
  path <- file.path(
    project_root,
    "results",
    "single_cell_downstream",
    dataset_id,
    "cell_annotations.tsv.gz"
  )
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  utils::read.delim(
    connection,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

v22_prepare_marker_sets <- function(
    project_root,
    dataset_id,
    disease,
    label_column
) {
  directory <- file.path(
    project_root,
    "results",
    "single_cell_downstream",
    dataset_id
  )
  markers <- utils::read.csv(
    file.path(directory, "cluster_markers.csv"),
    stringsAsFactors = FALSE
  )
  annotations <- v22_read_annotations(project_root, dataset_id)
  mapping <- aggregate(
    rep(1L, nrow(annotations)),
    by = list(
      cluster = as.character(annotations$analysis_cluster),
      cell_type = annotations[[label_column]]
    ),
    FUN = sum
  )
  names(mapping)[[3L]] <- "cells"
  mapping <- mapping[
    order(mapping$cluster, -mapping$cells, mapping$cell_type),
    ,
    drop = FALSE
  ]
  mapping <- mapping[!duplicated(mapping$cluster), , drop = FALSE]
  markers$cluster <- as.character(markers$cluster)
  markers <- merge(markers, mapping, by = "cluster", all.x = TRUE)
  markers$dataset_id <- dataset_id
  markers$disease <- disease
  markers
}

v22_cell_type_go <- function(project_root, paths) {
  markers <- rbind(
    v22_prepare_marker_sets(
      project_root,
      "GSE104782",
      "OA",
      "published_cell_type"
    ),
    v22_prepare_marker_sets(
      project_root,
      "GSE154600",
      "OC",
      "cell_type"
    )
  )
  marker_sets <- split(
    markers$gene,
    paste(markers$disease, markers$cell_type, sep = " | ")
  )
  universe <- AnnotationDbi::keys(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "SYMBOL"
  )
  results <- lapply(names(marker_sets), function(label) {
    genes <- unique(marker_sets[[label]])
    result <- suppressMessages(clusterProfiler::enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    ))
    table <- as.data.frame(result)
    if (nrow(table) == 0L) {
      return(NULL)
    }
    parts <- strsplit(label, " | ", fixed = TRUE)[[1L]]
    table$dataset_id <- if (parts[[1L]] == "OA") {
      "GSE104782"
    } else {
      "GSE154600"
    }
    table$disease <- parts[[1L]]
    table$cell_type <- parts[[2L]]
    table$marker_genes <- length(genes)
    table$marker_rule <- paste0(
      "Top 25 mean-expression contrast genes per analysis cluster; ",
      "clusters mapped to majority dataset-specific label"
    )
    table$analysis_scope <- "exploratory descriptive functional annotation"
    table$inference_boundary <- paste0(
      "Cluster-marker ranks are descriptive; enrichment is not evidence ",
      "of a conserved mechanism"
    )
    table
  })
  results <- do.call(rbind, results)
  results <- results[
    order(results$disease, results$cell_type, results$p.adjust),
    ,
    drop = FALSE
  ]
  results$fdr_significant <- results$p.adjust < 0.05
  safe_write_csv(
    results,
    file.path(
      paths$tables,
      "Table_S17_cell_type_marker_GO_annotation.csv"
    )
  )

  plotted <- do.call(rbind, lapply(
    split(results, paste(results$disease, results$cell_type)),
    function(table) head(table, 2L)
  ))
  plotted$minus_log10_fdr <- -log10(pmax(
    plotted$p.adjust,
    .Machine$double.xmin
  ))
  plotted$significance <- ifelse(
    plotted$fdr_significant,
    "FDR < 0.05",
    "FDR >= 0.05"
  )
  short_term <- function(value, width = 58L) {
    vapply(value, function(term) {
      if (nchar(term) <= width) {
        term
      } else {
        paste0(substr(term, 1L, width - 3L), "...")
      }
    }, character(1))
  }
  plotted$display_label <- paste0(
    gsub("\\.", " ", plotted$cell_type),
    " | ",
    short_term(plotted$Description)
  )
  safe_write_csv(
    plotted,
    file.path(paths$source, "SupplementaryFigure6_cell_type_GO.csv")
  )

  panel <- function(table, title) {
    table <- table[order(table$cell_type, table$p.adjust), , drop = FALSE]
    table$display_label <- factor(
      table$display_label,
      levels = rev(unique(table$display_label))
    )
    ggplot2::ggplot(
      table,
      ggplot2::aes(
        x = minus_log10_fdr,
        y = display_label,
        fill = significance,
        size = Count
      )
    ) +
      ggplot2::geom_vline(
        xintercept = -log10(0.05),
        linetype = "dashed",
        colour = "#9CA3AF"
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = minus_log10_fdr, yend = display_label),
        colour = "#D1D5DB",
        linewidth = 0.45
      ) +
      ggplot2::geom_point(
        shape = 21,
        colour = "#374151",
        stroke = 0.25
      ) +
      ggplot2::scale_fill_manual(values = c(
        "FDR < 0.05" = submission_palette[["shared"]],
        "FDR >= 0.05" = "#D1D5DB"
      )) +
      ggplot2::scale_size(range = c(2.0, 5.5)) +
      ggplot2::labs(
        title = title,
        x = "-log10(BH-adjusted P)",
        y = NULL,
        fill = NULL,
        size = "Genes"
      ) +
      submission_theme(7.3) +
      ggplot2::theme(
        legend.position = "top",
        axis.text.y = ggplot2::element_text(size = 6.6)
      )
  }
  oa <- panel(
    plotted[plotted$disease == "OA", , drop = FALSE],
    "OA cartilage: top cluster-marker processes"
  )
  oc <- panel(
    plotted[plotted$disease == "OC", , drop = FALSE],
    "Ovarian tumor atlas: top cluster-marker processes"
  )
  note <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = -0.97,
      y = 0,
      hjust = 0,
      size = 2.25,
      lineheight = 1.08,
      colour = "#374151",
      label = paste0(
        "Interpretation boundary: top 25 genes per cluster were mapped to\n",
        "majority dataset-specific labels and annotated by GO BP.\n",
        "Cluster-marker ranks are descriptive and localize functional themes.\n",
        "They do not demonstrate conserved function or mechanism."
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-0.5, 0.5)) +
    ggplot2::theme_void(base_family = "Arial")
  oa <- oa + ggplot2::labs(tag = "A")
  oc <- oc + ggplot2::labs(tag = "B")
  figure <- oa / oc / note +
    patchwork::plot_layout(heights = c(1.0, 1.35, 0.30)) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 14)
    )
  submission_save_plot(
    figure,
    "SupplementaryFigure6_cell_type_functional_annotation",
    paths$figures,
    height_mm = 190
  )
  results
}

v22_update_figure_materials <- function(paths) {
  legend_path <- file.path(paths$figures, "figure_legends.md")
  legends <- readLines(legend_path, warn = FALSE, encoding = "UTF-8")
  legends[legends ==
    "## Figure 6. Functional and exploratory prognostic context"] <-
    "## Figure 6. Functional context and exploratory survival association"
  legends <- c(
    legends,
    "",
    "## Supplementary Figure 6. Exploratory cell-type functional annotation",
    "",
    paste0(
      "**A,** Gene Ontology Biological Process over-representation among ",
      "top cluster markers from the GSE104782 OA cartilage atlas. ",
      "**B,** Corresponding annotation in GSE154600 ovarian tumors. ",
      "Clusters were mapped to their majority dataset-specific cell label; ",
      "the top two terms per label are shown. The dashed line marks BH FDR ",
      "0.05. Cluster-marker ranks are descriptive, so these panels localize ",
      "functional themes but do not demonstrate conserved function or mechanism."
    )
  )
  writeLines(enc2utf8(legends), legend_path, useBytes = TRUE)

  style_path <- file.path(paths$figures, "figure_style_manifest.csv")
  style <- utils::read.csv(style_path, stringsAsFactors = FALSE)
  style$value[style$setting == "revision"] <- "V2.2"
  safe_write_csv(style, style_path)
}

v22_update_table_index <- function(paths) {
  index_path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  additions <- data.frame(
    table_id = c("Table S16", "Table S17"),
    filename = c(
      "Table_S16_candidate_prioritization_matrix.csv",
      "Table_S17_cell_type_marker_GO_annotation.csv"
    ),
    title = c(
      "Transparent candidate prioritization matrix",
      "Exploratory cell-type marker functional annotation"
    ),
    contents = c(
      paste0(
        "Gene-level shared-DEG effects, WGCNA membership, original model ",
        "votes, strict nested selection frequency, and single-cell evidence."
      ),
      paste0(
        "GO Biological Process over-representation among descriptive top ",
        "cluster markers mapped to dataset-specific cell labels."
      )
    ),
    source = c(
      "results/cache/07_machine_learning.rds and linked source tables",
      "GSE104782/GSE154600 cluster_markers.csv and cell annotations"
    ),
    stringsAsFactors = FALSE
  )
  index <- index[!index$table_id %in% additions$table_id, , drop = FALSE]
  index <- rbind(index, additions)
  safe_write_csv(index, index_path)

  readme <- c(
    "# Supplementary table index",
    "",
    paste0(
      "All tables are UTF-8 CSV files. Interaction tables remain ",
      "hypothesis-generating and are not treatment recommendations."
    ),
    "",
    "| Table | File | Title | Contents |",
    "|---|---|---|---|"
  )
  readme <- c(readme, vapply(seq_len(nrow(index)), function(row) {
    paste0(
      "| ", index$table_id[[row]],
      " | `", index$filename[[row]],
      "` | ", index$title[[row]],
      " | ", index$contents[[row]], " |"
    )
  }, character(1)))
  writeLines(
    enc2utf8(readme),
    file.path(paths$tables, "README.md"),
    useBytes = TRUE
  )
}

v22_update_registries <- function(project_root, paths) {
  baseline <- file.path(project_root, "results", "submission_v21")
  claims <- utils::read.csv(file.path(
    baseline,
    "claim_evidence_registry_v21.csv"
  ))
  additions <- data.frame(
    claim_id = c("C18", "C19"),
    manuscript_claim = c(
      paste0(
        "The ten-gene set follows a transparent hierarchy spanning shared ",
        "DEGs, WGCNA membership, model votes, strict nested feature ",
        "frequency, and single-cell localization."
      ),
      paste0(
        "Exploratory top-cluster-marker annotation places OA and OC ",
        "candidates in distinct functional cell contexts."
      )
    ),
    primary_data = c(
      "results/submission_v22/supplementary_tables/Table_S16_candidate_prioritization_matrix.csv",
      "results/single_cell_downstream/GSE104782 and GSE154600 cluster_markers.csv"
    ),
    figure_or_table = c("Table S16", "Figure S6; Table S17"),
    allowed_wording = c(
      "transparent prioritization; heterogeneous feature stability",
      "exploratory descriptive functional annotation"
    ),
    prohibited_wording = c(
      "uniformly stable diagnostic or mechanistic panel",
      "functional proof or conserved cell mechanism"
    ),
    status = c("verified", "verified with explicit inference boundary"),
    stringsAsFactors = FALSE
  )
  claims <- claims[!claims$claim_id %in% additions$claim_id, , drop = FALSE]
  claims <- rbind(claims, additions)
  claims <- lapply(claims, function(column) {
    if (is.character(column)) {
      gsub("submission_v21", "submission_v22", column, fixed = TRUE)
    } else {
      column
    }
  })
  claims <- as.data.frame(claims, stringsAsFactors = FALSE)
  safe_write_csv(
    claims,
    file.path(paths$root, "claim_evidence_registry_v22.csv")
  )

  checklist <- utils::read.csv(file.path(
    baseline,
    "reproducibility_checklist_v21.csv"
  ))
  checklist$evidence <- gsub(
    "submission_v21",
    "submission_v22",
    checklist$evidence,
    fixed = TRUE
  )
  checklist$item[checklist$item_id == "R19"] <-
    "V2.2 reviewer-strengthening analyses are one-command reproducible."
  checklist$evidence[checklist$item_id == "R19"] <- "run_submission_v22.ps1"
  checklist$item[checklist$item_id == "R21"] <-
    "Six main and six supplementary figures have paired PDF/PNG outputs."
  checklist$evidence[checklist$item_id == "R21"] <-
    "results/submission_v22/figures/"
  additions <- data.frame(
    item_id = c("R23", "R24"),
    domain = c("candidate transparency", "single-cell interpretation"),
    item = c(
      "Candidate prioritization hierarchy is enumerated gene by gene.",
      "Cell-type functional annotation is reproducible and explicitly descriptive."
    ),
    status = c("complete", "complete"),
    evidence = c("Table S16", "Figure S6; Table S17"),
    stringsAsFactors = FALSE
  )
  checklist <- checklist[
    !checklist$item_id %in% additions$item_id,
    ,
    drop = FALSE
  ]
  checklist <- rbind(checklist, additions)
  safe_write_csv(
    checklist,
    file.path(paths$root, "reproducibility_checklist_v22.csv")
  )
}

run_reviewer_v22 <- function(project_root) {
  for (package in c(
    "ggplot2",
    "patchwork",
    "ragg",
    "clusterProfiler",
    "org.Hs.eg.db",
    "AnnotationDbi"
  )) {
    require_namespace(package, "V2.2 reviewer revision")
  }
  baseline <- file.path(project_root, "results", "submission_v21")
  if (!dir.exists(baseline)) {
    stop("V2.1 baseline outputs are required before V2.2.", call. = FALSE)
  }
  paths <- v22_output_paths(project_root)
  v22_copy_tree(
    baseline,
    paths$root,
    skip = c(
      "claim_evidence_registry_v21.csv",
      "reproducibility_checklist_v21.csv",
      "submission_audit_v21.json"
    )
  )
  log_info("Building V2.2 candidate-prioritization matrix.")
  candidate_matrix <- v22_prepare_candidate_matrix(project_root, paths)
  log_info("Building exploratory cell-type functional annotation.")
  cell_go <- v22_cell_type_go(project_root, paths)
  v22_update_figure_materials(paths)
  v22_update_table_index(paths)
  v22_update_registries(project_root, paths)
  log_info(
    "V2.2 reviewer revision completed: ",
    nrow(candidate_matrix),
    " candidates and ",
    nrow(cell_go),
    " cell-type GO rows."
  )
  invisible(paths)
}
