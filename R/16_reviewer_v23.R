v23_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v23"))
  figures <- ensure_dir(file.path(root, "figures"))
  source <- ensure_dir(file.path(figures, "source_data"))
  tables <- ensure_dir(file.path(root, "supplementary_tables"))
  analysis <- ensure_dir(file.path(root, "analysis"))
  logs <- ensure_dir(file.path(root, "logs"))
  list(
    root = root,
    figures = figures,
    source = source,
    tables = tables,
    analysis = analysis,
    logs = logs
  )
}

v23_support_label <- function(oa, oc) {
  ifelse(
    oa & oc,
    "OA and OC",
    ifelse(oa, "OA only", ifelse(oc, "OC only", "not selected"))
  )
}

v23_upgrade_candidate_matrix <- function(project_root, paths) {
  input <- utils::read.csv(file.path(
    paths$tables,
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  genes <- input$gene
  oa <- ml$datasets$oa_train
  oc <- ml$datasets$oc_train

  lasso_oa <- genes %in% oa$lasso_genes
  lasso_oc <- genes %in% oc$lasso_genes
  rf_oa <- genes %in% oa$rf_genes
  rf_oc <- genes %in% oc$rf_genes

  input$shared_DEG <- "yes, primary threshold in both diseases"
  input$direction <- input$direction_class
  input$WGCNA_support <- ifelse(
    input$OA_primary_WGCNA_module & input$OC_primary_WGCNA_module,
    "OA and OC primary modules",
    ifelse(
      input$OA_primary_WGCNA_module,
      "OA primary module",
      ifelse(input$OC_primary_WGCNA_module, "OC primary module", "none")
    )
  )
  input$original_LASSO_support_OA <- lasso_oa
  input$original_LASSO_support_OC <- lasso_oc
  input$LASSO_support <- v23_support_label(lasso_oa, lasso_oc)
  input$original_random_forest_support_OA <- rf_oa
  input$original_random_forest_support_OC <- rf_oc
  input$random_forest_support <- v23_support_label(rf_oa, rf_oc)
  input$strict_nested_frequency <- sprintf(
    "OA %.3f; OC %.3f",
    input$strict_nested_OA_max_selection_frequency,
    input$strict_nested_OC_max_selection_frequency
  )
  input$single_cell_context <- sprintf(
    "OA: %s (detect %.3f); OC: %s (detect %.3f)",
    input$OA_top_cell_context,
    input$OA_top_context_detection_fraction,
    input$OC_top_cell_context,
    input$OC_top_context_detection_fraction
  )
  input$ten_gene_set_role <- paste0(
    "Interpretable evidence summary; not an optimized predictive signature"
  )

  leading <- c(
    "gene",
    "prioritization_rank",
    "shared_DEG",
    "direction",
    "WGCNA_support",
    "LASSO_support",
    "random_forest_support",
    "strict_nested_frequency",
    "single_cell_context",
    "ten_gene_set_role"
  )
  input <- input[, c(leading, setdiff(names(input), leading)), drop = FALSE]
  safe_write_csv(
    input,
    file.path(
      paths$tables,
      "Table_S16_candidate_prioritization_matrix.csv"
    )
  )
  safe_write_csv(
    input,
    file.path(
      paths$source,
      "SupplementaryFigure6_candidate_prioritization_matrix.csv"
    )
  )
  input
}

v23_complete_hallmark_gsea <- function(object, seed) {
  gene_sets <- object@geneSets
  term_to_gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unname(unlist(gene_sets, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  result <- suppressMessages(clusterProfiler::GSEA(
    geneList = object@geneList,
    TERM2GENE = term_to_gene,
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE
  ))
  as.data.frame(result)
}

v23_title_case_pathway <- function(ids) {
  vapply(ids, function(id) {
    value <- gsub("_", " ", sub("^HALLMARK_", "", id))
    tools::toTitleCase(tolower(value))
  }, character(1))
}

v23_pathway_direction_analysis <- function(project_root, paths) {
  enrichment <- submission_load_cache(project_root, "05_enrichment.rds")
  oa <- v23_complete_hallmark_gsea(
    enrichment$gsea$oa_train$hallmark,
    20260726L
  )
  oc <- v23_complete_hallmark_gsea(
    enrichment$gsea$oc_train$hallmark,
    20260727L
  )
  oa <- oa[, c(
    "ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment"
  )]
  oc <- oc[, c(
    "ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment"
  )]
  names(oa) <- c(
    "pathway_id",
    "OA_set_size",
    "OA_NES",
    "OA_P",
    "OA_FDR",
    "OA_leading_edge"
  )
  names(oc) <- c(
    "pathway_id",
    "OC_set_size",
    "OC_NES",
    "OC_P",
    "OC_FDR",
    "OC_leading_edge"
  )
  paired <- merge(oa, oc, by = "pathway_id", all = TRUE, sort = FALSE)
  paired$pathway <- v23_title_case_pathway(paired$pathway_id)
  paired$OA_significant <- paired$OA_FDR < 0.05
  paired$OC_significant <- paired$OC_FDR < 0.05
  paired$both_significant <- paired$OA_significant & paired$OC_significant
  paired$direction_class <- ifelse(
    sign(paired$OA_NES) == sign(paired$OC_NES),
    "concordant",
    "discordant"
  )
  paired$paired_direction_index <- paired$OA_NES * paired$OC_NES
  paired$evidence_class <- ifelse(
    paired$both_significant,
    paste("both FDR < 0.05", paired$direction_class),
    ifelse(
      paired$OA_significant | paired$OC_significant,
      "one disease FDR < 0.05",
      "neither disease FDR < 0.05"
    )
  )
  paired$analysis_scope <- paste0(
    "Secondary descriptive comparison of independently estimated Hallmark NES"
  )
  paired$inference_boundary <- paste0(
    "NES agreement or opposition does not establish a shared mechanism"
  )
  paired <- paired[
    order(
      !paired$both_significant,
      -abs(paired$paired_direction_index),
      paired$pathway
    ),
    ,
    drop = FALSE
  ]
  rownames(paired) <- NULL
  safe_write_csv(
    paired,
    file.path(
      paths$tables,
      "Table_S18_Hallmark_pathway_direction_matrix.csv"
    )
  )
  safe_write_csv(
    paired,
    file.path(paths$source, "SupplementaryFigure7_pathway_direction.csv")
  )

  paired$display_class <- factor(
    ifelse(
      paired$both_significant,
      paste0("Both sig: ", paired$direction_class),
      ifelse(
        paired$OA_significant | paired$OC_significant,
        "One sig",
        "Neither sig"
      )
    ),
    levels = c(
      "Both sig: concordant",
      "Both sig: discordant",
      "One sig",
      "Neither sig"
    )
  )
  label_rows <- paired[paired$both_significant, , drop = FALSE]
  label_rows <- head(
    label_rows[order(-abs(label_rows$paired_direction_index)), , drop = FALSE],
    8L
  )
  colors <- c(
    "Both sig: concordant" = "#009E73",
    "Both sig: discordant" = "#D55E00",
    "One sig" = "#6B7280",
    "Neither sig" = "#D1D5DB"
  )
  shapes <- c(
    "Both sig: concordant" = 16,
    "Both sig: discordant" = 17,
    "One sig" = 1,
    "Neither sig" = 3
  )
  p1 <- ggplot2::ggplot(
    paired,
    ggplot2::aes(
      x = OA_NES,
      y = OC_NES,
      colour = display_class,
      shape = display_class
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.35) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "#D1D5DB",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_rows,
      ggplot2::aes(label = pathway),
      size = 2.05,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(
      title = "Hallmark direction across diseases",
      subtitle = "Same sign = concordant; opposite signs = discordant",
      x = "OA normalized enrichment score",
      y = "OC normalized enrichment score",
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.3) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 5.8)
    )

  selected <- paired[paired$both_significant, , drop = FALSE]
  selected <- head(
    selected[order(-abs(selected$paired_direction_index)), , drop = FALSE],
    12L
  )
  selected$pathway_label <- vapply(selected$pathway, function(value) {
    if (nchar(value) <= 39L) value else paste0(substr(value, 1L, 36L), "...")
  }, character(1))
  selected$pathway_label <- factor(
    selected$pathway_label,
    levels = rev(selected$pathway_label)
  )
  long <- rbind(
    data.frame(
      pathway_label = selected$pathway_label,
      disease = "OA",
      NES = selected$OA_NES,
      stringsAsFactors = FALSE
    ),
    data.frame(
      pathway_label = selected$pathway_label,
      disease = "OC",
      NES = selected$OC_NES,
      stringsAsFactors = FALSE
    )
  )
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = 0,
      colour = "#9CA3AF",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    ggplot2::geom_segment(
      data = selected,
      ggplot2::aes(
        x = OA_NES,
        xend = OC_NES,
        y = pathway_label,
        yend = pathway_label
      ),
      colour = "#CBD5E1",
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = NES, y = pathway_label, colour = disease, shape = disease),
      size = 2.5
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17)) +
    ggplot2::labs(
      title = "Strongest jointly significant pathways",
      subtitle = "Independent OA and OC NES; FDR < 0.05 in both",
      x = "Normalized enrichment score",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.3) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 6.5)
    )
  figure <- submission_panel_tag(p1 | p2)
  submission_save_plot(
    figure,
    "SupplementaryFigure7_pathway_direction",
    paths$figures,
    height_mm = 125
  )
  paired
}

v23_gene_cell_function_matrix <- function(paths, candidates) {
  go <- utils::read.csv(file.path(
    paths$tables,
    "Table_S17_cell_type_marker_GO_annotation.csv"
  ))
  pick_term <- function(disease, cell_type) {
    rows <- go[
      go$disease == disease & go$cell_type == cell_type,
      ,
      drop = FALSE
    ]
    if (nrow(rows) == 0L) {
      return(list(
        term = "Not available for this dataset-specific cell label",
        fdr = NA_real_,
        status = "unmatched dataset-specific label"
      ))
    }
    rows <- rows[order(rows$p.adjust, -rows$Count), , drop = FALSE]
    list(
      term = rows$Description[[1L]],
      fdr = rows$p.adjust[[1L]],
      status = "matched descriptive cell-type enrichment"
    )
  }
  rows <- lapply(seq_len(nrow(candidates)), function(index) {
    row <- candidates[index, , drop = FALSE]
    oa <- pick_term("OA", row$OA_top_cell_context[[1L]])
    oc <- pick_term("OC", row$OC_top_cell_context[[1L]])
    data.frame(
      gene = row$gene,
      direction = row$direction,
      OA_cell_context = row$OA_top_cell_context,
      OA_detection_fraction = row$OA_top_context_detection_fraction,
      OA_functional_theme = oa$term,
      OA_functional_theme_FDR = oa$fdr,
      OA_functional_theme_status = oa$status,
      OC_cell_context = row$OC_top_cell_context,
      OC_detection_fraction = row$OC_top_context_detection_fraction,
      OC_functional_theme = oc$term,
      OC_functional_theme_FDR = oc$fdr,
      OC_functional_theme_status = oc$status,
      interpretation = paste0(
        "Gene-cell-function mapping is descriptive; same gene does not imply ",
        "the same biological meaning across diseases"
      ),
      stringsAsFactors = FALSE
    )
  })
  matrix <- do.call(rbind, rows)
  safe_write_csv(
    matrix,
    file.path(
      paths$tables,
      "Table_S19_gene_cell_function_context_matrix.csv"
    )
  )
  matrix
}

v23_build_figure2 <- function(project_root, paths) {
  common <- utils::read.csv(file.path(
    paths$source,
    "Figure2_common_gene_effects_quadrants.csv"
  ))
  counts <- utils::read.csv(file.path(
    paths$source,
    "Figure2_quadrant_counts.csv"
  ))
  threshold <- utils::read.csv(file.path(
    paths$source,
    "Figure2_DEG_threshold_sensitivity.csv"
  ))
  oa <- utils::read.csv(file.path(paths$source, "Figure2_OA_DEG.csv"))
  oc <- utils::read.csv(file.path(paths$source, "Figure2_OC_DEG.csv"))
  hub_genes <- submission_load_cache(project_root, "07_machine_learning.rds")$final_genes

  p1 <- submission_volcano_plot(oa, "OA", hub_genes)
  p2 <- submission_volcano_plot(oc, "OC", hub_genes)
  quadrant_levels <- c(
    "OA higher / OC higher",
    "OA lower / OC lower",
    "OA higher / OC lower",
    "OA lower / OC higher"
  )
  common$direction_quadrant <- factor(
    common$direction_quadrant,
    levels = c("Not primary shared", quadrant_levels)
  )
  shared_only <- common[common$primary_shared, , drop = FALSE]
  candidates <- common[common$hub, , drop = FALSE]
  quadrant_colors <- c(
    "OA higher / OC higher" = "#009E73",
    "OA lower / OC lower" = "#6A3D9A",
    "OA higher / OC lower" = "#0072B2",
    "OA lower / OC higher" = "#D55E00"
  )
  p3_base <- ggplot2::ggplot(
    common,
    ggplot2::aes(x = logFC_OA, y = logFC_OC)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_point(colour = "#D4D7DB", size = 0.31, alpha = 0.30) +
    ggplot2::geom_point(
      data = shared_only,
      ggplot2::aes(colour = direction_quadrant),
      size = 0.82,
      alpha = 0.72
    ) +
    ggplot2::scale_colour_manual(values = quadrant_colors) +
    ggplot2::labs(
      title = "Shared membership does not imply shared direction",
      subtitle = "All measured genes shown; 146/286 concordant and 140/286 discordant",
      x = "OA log2 fold change",
      y = "OC log2 fold change",
      colour = NULL
    ) +
    submission_theme(7.5) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(2.8, "mm"),
      legend.text = ggplot2::element_text(size = 6.0)
    ) +
    ggplot2::labs(tag = "C")
  inset <- ggplot2::ggplot(
    candidates,
    ggplot2::aes(
      x = logFC_OA,
      y = logFC_OC,
      colour = direction_quadrant
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.25) +
    ggplot2::geom_point(size = 1.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = gene),
      size = 2.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.18,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = quadrant_colors) +
    ggplot2::labs(title = "Ten-gene evidence summary", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 6.5, base_family = "Arial") +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(
        fill = "white",
        colour = "#6B7280",
        linewidth = 0.35
      ),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 7.2),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
  p3 <- p3_base + patchwork::inset_element(
    inset,
    left = 0.50,
    bottom = 0.04,
    right = 0.99,
    top = 0.55
  )

  threshold$fdr_label <- paste0("FDR < ", threshold$fdr_threshold)
  p4 <- ggplot2::ggplot(
    threshold,
    ggplot2::aes(
      x = factor(absolute_log2fc_threshold),
      y = shared_count,
      colour = fdr_label,
      group = fdr_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(
      ggplot2::aes(label = shared_count),
      vjust = -0.7,
      size = 2.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(
      "FDR < 0.01" = submission_palette[["OA"]],
      "FDR < 0.05" = submission_palette[["OC"]]
    )) +
    ggplot2::expand_limits(y = max(threshold$shared_count) * 1.13) +
    ggplot2::labs(
      title = "Prespecified threshold grid",
      x = "Absolute log2 fold-change threshold",
      y = "Shared genes",
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  figure <- (
    (p1 + ggplot2::labs(tag = "A")) |
      (p2 + ggplot2::labs(tag = "B"))
  ) / (
    p3 |
      (p4 + ggplot2::labs(tag = "D"))
  ) & ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 14)
  )
  submission_save_plot(
    figure,
    "Figure2_bulk_discovery",
    paths$figures,
    height_mm = 155
  )
}

v23_feature_summary <- function(grid) {
  groups <- split(grid, paste(grid$disease, grid$model, sep = " | "))
  do.call(rbind, lapply(groups, function(table) {
    data.frame(
      disease = table$disease[[1L]],
      model = table$model[[1L]],
      candidates = nrow(table),
      median_frequency = stats::median(table$selection_frequency),
      maximum_frequency = max(table$selection_frequency),
      candidates_frequency_ge_0_25 = sum(table$selection_frequency >= 0.25),
      candidates_frequency_ge_0_50 = sum(table$selection_frequency >= 0.50),
      stringsAsFactors = FALSE
    )
  }))
}

v23_build_figure3 <- function(project_root, paths) {
  powers <- utils::read.csv(file.path(paths$source, "Figure3_WGCNA_powers.csv"))
  grid <- utils::read.csv(file.path(
    paths$source,
    "Figure3_ML_selection_frequency.csv"
  ))
  hub_effects <- utils::read.csv(file.path(
    paths$source,
    "Figure3_hub_gene_effects.csv"
  ))
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  summary <- v23_feature_summary(grid)
  safe_write_csv(
    summary,
    file.path(paths$source, "Figure3_feature_stability_summary.csv")
  )

  p1 <- ggplot2::ggplot(
    powers,
    ggplot2::aes(
      x = soft_power,
      y = absolute_module_trait_correlation,
      colour = disease,
      group = dataset_id
    )
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(ggplot2::aes(shape = is_primary), size = 2.2) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16)) +
    ggplot2::coord_cartesian(ylim = c(0.8, 1)) +
    ggplot2::labs(
      title = "Module-trait stability across powers",
      x = "Soft power",
      y = expression("|r|"),
      colour = NULL,
      shape = "Primary"
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p2 <- ggplot2::ggplot(
    powers,
    ggplot2::aes(
      x = soft_power,
      y = primary_gene_retention,
      colour = disease,
      group = dataset_id
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0.8,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(ggplot2::aes(shape = is_primary), size = 2.2) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), labels = scales::percent) +
    ggplot2::labs(
      title = "Primary-module gene retention",
      x = "Soft power",
      y = "Retention",
      colour = NULL,
      shape = "Primary"
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "none")

  grid$model <- factor(grid$model, levels = c("LASSO", "RandomForest"))
  p3 <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = model,
      y = selection_frequency,
      colour = disease,
      shape = model
    )
  ) +
    ggplot2::geom_jitter(
      width = 0.10,
      height = 0,
      size = 1.55,
      alpha = 0.72
    ) +
    ggplot2::stat_summary(
      fun = max,
      geom = "point",
      shape = 18,
      size = 3.5,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~disease, nrow = 1) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(LASSO = 16, RandomForest = 17)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Feature stability summary",
      subtitle = "Ten candidates; points are genes and diamonds mark maxima",
      x = NULL,
      y = "Strict nested selection frequency",
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.6) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )

  hub_effects$gene <- factor(hub_effects$gene, levels = rev(ml$final_genes))
  p4 <- ggplot2::ggplot(
    hub_effects,
    ggplot2::aes(x = disease, y = gene, fill = logFC)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", logFC)),
      size = 2.1
    ) +
    ggplot2::scale_fill_gradient2(
      low = submission_palette[["negative"]],
      mid = "white",
      high = submission_palette[["positive"]],
      midpoint = 0
    ) +
    ggplot2::labs(
      title = "Candidate discovery effects",
      x = NULL,
      y = NULL,
      fill = "log2 FC"
    ) +
    submission_theme()
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure3_network_and_ml_stability",
    paths$figures
  )

  grid$gene <- factor(grid$gene, levels = rev(ml$final_genes))
  heatmap_panel <- function(table, disease) {
    table <- table[table$disease == disease, , drop = FALSE]
    ggplot2::ggplot(
      table,
      ggplot2::aes(x = model, y = gene, fill = selection_frequency)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.0f%%", 100 * selection_frequency)),
        size = 2.35
      ) +
      ggplot2::scale_fill_gradientn(
        colours = c("#F3F4F6", "#9ECAE1", submission_palette[[disease]]),
        limits = c(0, 1),
        labels = scales::percent
      ) +
      ggplot2::labs(
        title = paste(disease, "models"),
        x = NULL,
        y = NULL,
        fill = "Frequency"
      ) +
      submission_theme(8.0) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
        legend.position = "right"
      )
  }
  detailed <- submission_panel_tag(
    heatmap_panel(grid, "OA") | heatmap_panel(grid, "OC")
  )
  submission_save_plot(
    detailed,
    "SupplementaryFigure8_detailed_feature_stability",
    paths$figures,
    height_mm = 120
  )
}

v23_build_figure4 <- function(project_root, paths) {
  pca <- v21_prepare_gse54388_pca(project_root, paths)
  curves <- v21_prepare_roc_curves(project_root, paths)
  sensitivity_dir <- file.path(
    project_root,
    "results",
    "submission",
    "sensitivity"
  )
  composite <- utils::read.csv(file.path(
    sensitivity_dir,
    "external_validation_signed_composite_score.csv"
  ))
  permutations <- utils::read.csv(file.path(
    sensitivity_dir,
    "external_validation_permutation_auc.csv"
  ))
  leave_one_out <- utils::read.csv(file.path(
    sensitivity_dir,
    "external_validation_leave_one_out_auc.csv"
  ))
  levels_dataset <- c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
  for (object_name in c("composite", "permutations", "leave_one_out")) {
    object <- get(object_name)
    object$dataset_id <- factor(object$dataset_id, levels = levels_dataset)
    assign(object_name, object)
  }
  curves$dataset_id <- factor(curves$dataset_id, levels = levels_dataset)
  curve_labels <- unique(curves[, c("dataset_id", "disease", "auc")])
  curve_labels$legend_label <- sprintf(
    "%s (AUC %.3f)",
    curve_labels$dataset_id,
    curve_labels$auc
  )
  curves <- merge(
    curves,
    curve_labels[, c("dataset_id", "legend_label")],
    by = "dataset_id",
    all.x = TRUE,
    sort = FALSE
  )
  curves$legend_label <- factor(
    curves$legend_label,
    levels = sprintf(
      "%s (AUC %.3f)",
      levels_dataset,
      c(0.520, 0.629, 1.000, 0.979)
    )
  )
  safe_write_csv(
    pca$scores,
    file.path(paths$source, "Figure4_GSE54388_unsupervised_PCA.csv")
  )
  safe_write_csv(
    curves,
    file.path(paths$source, "Figure4_direction_fixed_ROC_curves.csv")
  )
  safe_write_csv(
    permutations,
    file.path(paths$source, "Figure4_permutation_AUC.csv")
  )
  safe_write_csv(
    leave_one_out,
    file.path(paths$source, "Figure4_leave_one_out_AUC.csv")
  )

  pc1_variance <- 100 * pca$variance$variance_fraction[[1L]]
  pc2_variance <- 100 * pca$variance$variance_fraction[[2L]]
  p1 <- ggplot2::ggplot(
    pca$scores,
    ggplot2::aes(x = PC1, y = PC2, colour = group, shape = group)
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = group),
      type = "norm",
      linewidth = 0.55,
      linetype = 2,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 2.4, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = c(
      Normal = "#7A8793",
      Disease = submission_palette[["OC"]]
    )) +
    ggplot2::scale_shape_manual(values = c(Normal = 17, Disease = 16)) +
    ggplot2::labs(
      title = "GSE54388 unsupervised PCA",
      subtitle = "Top 2,000 variable genes; labels used only for display",
      x = sprintf("PC1 (%.1f%%)", pc1_variance),
      y = sprintf("PC2 (%.1f%%)", pc2_variance),
      colour = NULL,
      shape = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")

  p2 <- ggplot2::ggplot(
    permutations,
    ggplot2::aes(x = auc, fill = disease)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      alpha = 0.75,
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::geom_vline(
      data = composite,
      ggplot2::aes(xintercept = auc, colour = disease),
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~dataset_id, ncol = 2) +
    ggplot2::scale_fill_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Label-permutation null distributions",
      subtitle = "Observed fixed-score AUC shown by vertical line",
      x = "Permuted-label AUC",
      y = "Frequency"
    ) +
    submission_theme(7.4) +
    ggplot2::theme(legend.position = "none")

  p3 <- ggplot2::ggplot(
    leave_one_out,
    ggplot2::aes(
      x = leave_one_out_auc,
      y = dataset_id,
      colour = disease
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0.5,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_jitter(height = 0.12, width = 0, alpha = 0.65, size = 1.4) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 18,
      size = 3.1
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Leave-one-sample-out influence",
      subtitle = "Diamonds mark cohort medians",
      x = "AUC after omitting one sample",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")

  p4 <- ggplot2::ggplot(
    curves,
    ggplot2::aes(
      x = false_positive_rate,
      y = sensitivity,
      colour = legend_label,
      linetype = legend_label
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "#B8BDC4",
      linetype = "dashed"
    ) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c(
      "GSE117999 (AUC 0.520)" = "#0072B2",
      "GSE82107 (AUC 0.629)" = "#56B4E9",
      "GSE54388 (AUC 1.000)" = "#D55E00",
      "GSE12470 (AUC 0.979)" = "#E69F00"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "GSE117999 (AUC 0.520)" = "solid",
      "GSE82107 (AUC 0.629)" = "dashed",
      "GSE54388 (AUC 1.000)" = "solid",
      "GSE12470 (AUC 0.979)" = "dashed"
    )) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Secondary direction-fixed ROC display",
      subtitle = "Retrospective molecular separation",
      x = "1 - specificity",
      y = "Sensitivity",
      colour = NULL,
      linetype = NULL
    ) +
    submission_theme(7.2) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = "none"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.2)
    )

  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure4_external_validation",
    paths$figures,
    height_mm = 160
  )
}

v23_replace_legend <- function(lines, heading, body) {
  target <- paste0("## ", heading)
  start <- which(lines == target)
  if (length(start) != 1L) {
    stop("Expected one figure legend heading: ", heading, call. = FALSE)
  }
  next_headings <- which(grepl("^## ", lines) & seq_along(lines) > start)
  stop_at <- if (length(next_headings) > 0L) next_headings[[1L]] - 1L else length(lines)
  before <- if (start > 1L) lines[seq_len(start - 1L)] else character()
  after <- if (stop_at < length(lines)) lines[(stop_at + 1L):length(lines)] else character()
  c(before, target, "", body, "", after)
}

v23_update_figure_materials <- function(paths) {
  legend_path <- file.path(paths$figures, "figure_legends.md")
  legends <- readLines(legend_path, warn = FALSE, encoding = "UTF-8")
  legends <- v23_replace_legend(
    legends,
    "Figure 2. Direction-aware cross-disease transcriptomic discovery",
    paste0(
      "**A-B,** OA and OC differential-expression volcano plots. ",
      "**C,** All commonly measured genes with the 286 primary shared DEGs ",
      "colored by direction; the inset isolates the ten-gene interpretable ",
      "evidence summary. **D,** Shared-gene counts under six prespecified ",
      "thresholds."
    )
  )
  legends <- v23_replace_legend(
    legends,
    "Figure 3. Robust identification of molecular candidates",
    paste0(
      "**A,** WGCNA module-trait association under soft-power perturbation. ",
      "**B,** Primary-module gene retention. **C,** Summary of strict nested ",
      "feature-selection frequency across the ten candidates; points are genes ",
      "and diamonds are within-model maxima. **D,** Disease-specific discovery ",
      "effects. Detailed OA and OC model heatmaps are in Figure S8."
    )
  )
  legends <- v23_replace_legend(
    legends,
    "Figure 4. Direction-fixed cross-cohort molecular reproducibility",
    paste0(
      "**A,** Unsupervised PCA of GSE54388. **B,** Null AUC distributions ",
      "from 1,000 label permutations. **C,** Leave-one-sample-out AUCs. ",
      "**D,** ROC curves for the fixed signed ten-gene score, deliberately ",
      "shown last as a secondary display. These panels show retrospective ",
      "molecular separation, not clinical diagnostic performance."
    )
  )
  additions <- c(
    "## Supplementary Figure 7. Hallmark pathway direction across diseases",
    "",
    paste0(
      "**A,** OA and OC normalized enrichment scores (NES) for the complete ",
      "Hallmark collection. Same-sign estimates indicate directional ",
      "concordance and opposite signs indicate discordance. **B,** Paired NES ",
      "for the strongest pathways significant in both diseases. The paired ",
      "direction index is descriptive and does not establish a shared mechanism."
    ),
    "",
    "## Supplementary Figure 8. Detailed strict nested feature stability",
    "",
    paste0(
      "**A,** OA LASSO and random-forest selection frequencies for the ten ",
      "candidates. **B,** Corresponding OC frequencies. These detailed ",
      "heatmaps are separated from the main feature-stability summary."
    )
  )
  for (heading in c(
    "## Supplementary Figure 7. Hallmark pathway direction across diseases",
    "## Supplementary Figure 8. Detailed strict nested feature stability"
  )) {
    if (heading %in% legends) {
      start <- which(legends == heading)
      next_headings <- which(grepl("^## ", legends) & seq_along(legends) > start)
      stop_at <- if (length(next_headings) > 0L) {
        next_headings[[1L]] - 1L
      } else {
        length(legends)
      }
      legends <- legends[-(start:stop_at)]
    }
  }
  legends <- c(legends, "", additions, "")
  writeLines(enc2utf8(legends), legend_path, useBytes = TRUE)

  style_path <- file.path(paths$figures, "figure_style_manifest.csv")
  style <- utils::read.csv(style_path, stringsAsFactors = FALSE)
  style$value[style$setting == "revision"] <- "V2.3"
  safe_write_csv(style, style_path)
}

v23_update_table_index <- function(paths) {
  index_path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  additions <- data.frame(
    table_id = c("Table S18", "Table S19"),
    filename = c(
      "Table_S18_Hallmark_pathway_direction_matrix.csv",
      "Table_S19_gene_cell_function_context_matrix.csv"
    ),
    title = c(
      "Hallmark pathway direction matrix",
      "Candidate gene-cell-function context matrix"
    ),
    contents = c(
      paste0(
        "Complete paired OA/OC Hallmark NES, FDR, direction class, and ",
        "descriptive paired direction index."
      ),
      paste0(
        "Candidate-specific OA and OC cell contexts linked to descriptive ",
        "cell-type functional themes."
      )
    ),
    source = c(
      "Complete Hallmark GSEA using the original full ranked gene lists",
      "Table S16 candidate contexts joined to Table S17 cell-type annotation"
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

v23_update_registries <- function(project_root, paths) {
  baseline <- file.path(project_root, "results", "submission_v22")
  claims <- utils::read.csv(file.path(
    baseline,
    "claim_evidence_registry_v22.csv"
  ))
  additions <- data.frame(
    claim_id = c("C20", "C21", "C22"),
    manuscript_claim = c(
      paste0(
        "The ten-gene set is an interpretable evidence summary rather than ",
        "an optimized predictive signature."
      ),
      paste0(
        "Complete Hallmark NES comparison identifies both concordant and ",
        "discordant pathway directions across OA and OC."
      ),
      paste0(
        "Candidate gene-cell-function mapping shows disease-specific cellular ",
        "and functional contexts."
      )
    ),
    primary_data = c(
      "results/submission_v23/supplementary_tables/Table_S16_candidate_prioritization_matrix.csv",
      "results/submission_v23/supplementary_tables/Table_S18_Hallmark_pathway_direction_matrix.csv",
      "results/submission_v23/supplementary_tables/Table_S19_gene_cell_function_context_matrix.csv"
    ),
    figure_or_table = c(
      "Table S16",
      "Figure S7; Table S18",
      "Table S19"
    ),
    allowed_wording = c(
      "interpretable evidence summary; not optimized signature",
      "directional concordance or discordance of independently estimated NES",
      "descriptive disease-specific gene-cell-function context"
    ),
    prohibited_wording = c(
      "optimized diagnostic signature",
      "shared pathway mechanism",
      "same biological function across diseases"
    ),
    status = c(
      "verified",
      "verified with descriptive metric",
      "verified with explicit inference boundary"
    ),
    stringsAsFactors = FALSE
  )
  claims <- claims[!claims$claim_id %in% additions$claim_id, , drop = FALSE]
  claims <- rbind(claims, additions)
  claims <- lapply(claims, function(column) {
    if (is.character(column)) {
      gsub("submission_v22", "submission_v23", column, fixed = TRUE)
    } else {
      column
    }
  })
  claims <- as.data.frame(claims, stringsAsFactors = FALSE)
  safe_write_csv(
    claims,
    file.path(paths$root, "claim_evidence_registry_v23.csv")
  )

  checklist <- utils::read.csv(file.path(
    baseline,
    "reproducibility_checklist_v22.csv"
  ))
  checklist$evidence <- gsub(
    "submission_v22",
    "submission_v23",
    checklist$evidence,
    fixed = TRUE
  )
  checklist$item[checklist$item_id == "R19"] <-
    "V2.3 pre-submission refinement is one-command reproducible."
  checklist$evidence[checklist$item_id == "R19"] <- "run_submission_v23.ps1"
  checklist$item[checklist$item_id == "R21"] <-
    "Six main and eight supplementary figures have paired PDF/PNG outputs."
  checklist$evidence[checklist$item_id == "R21"] <-
    "results/submission_v23/figures/"
  additions <- data.frame(
    item_id = c("R25", "R26", "R27"),
    domain = c(
      "candidate transparency",
      "pathway direction",
      "figure hierarchy"
    ),
    item = c(
      "The ten-gene role is explicitly an evidence summary, not an optimized signature.",
      "Complete Hallmark direction analysis is reproducible from full ranked lists.",
      "ROC is visually secondary and detailed feature heatmaps are supplementary."
    ),
    status = c("complete", "complete", "complete"),
    evidence = c(
      "Table S16",
      "Figure S7; Table S18",
      "Figures 3-4; Figure S8"
    ),
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
    file.path(paths$root, "reproducibility_checklist_v23.csv")
  )
}

run_reviewer_v23 <- function(project_root) {
  for (package in c(
    "ggplot2",
    "patchwork",
    "ragg",
    "ggrepel",
    "clusterProfiler",
    "DOSE"
  )) {
    require_namespace(package, "V2.3 pre-submission refinement")
  }
  baseline <- file.path(project_root, "results", "submission_v22")
  if (!dir.exists(baseline)) {
    stop("V2.2 baseline outputs are required before V2.3.", call. = FALSE)
  }
  paths <- v23_output_paths(project_root)
  v22_copy_tree(
    baseline,
    paths$root,
    skip = c(
      "claim_evidence_registry_v22.csv",
      "reproducibility_checklist_v22.csv",
      "submission_audit_v22.json"
    )
  )
  log_info("Upgrading V2.3 candidate evidence matrix.")
  candidates <- v23_upgrade_candidate_matrix(project_root, paths)
  log_info("Running complete Hallmark pathway-direction analysis.")
  pathways <- v23_pathway_direction_analysis(project_root, paths)
  log_info("Building candidate gene-cell-function context matrix.")
  gene_context <- v23_gene_cell_function_matrix(paths, candidates)
  log_info("Rebuilding Figures 2, 3, and 4 with revised visual hierarchy.")
  v23_build_figure2(project_root, paths)
  v23_build_figure3(project_root, paths)
  v23_build_figure4(project_root, paths)
  v23_update_figure_materials(paths)
  v23_update_table_index(paths)
  v23_update_registries(project_root, paths)
  log_info(
    "V2.3 refinement completed: ",
    nrow(candidates),
    " candidates, ",
    nrow(pathways),
    " pathway pairs, and ",
    nrow(gene_context),
    " gene-cell-function rows."
  )
  invisible(paths)
}
