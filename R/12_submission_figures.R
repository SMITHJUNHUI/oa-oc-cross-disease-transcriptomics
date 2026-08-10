submission_palette <- c(
  OA = "#0072B2",
  OC = "#D55E00",
  shared = "#6A3D9A",
  positive = "#009E73",
  negative = "#CC79A7",
  neutral = "#6B7280",
  light = "#E5E7EB"
)

submission_theme <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 0.5,
        hjust = 0
      ),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(colour = "#202124"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(3.5, "mm"),
      panel.grid.major.y = ggplot2::element_line(
        colour = "#ECEFF1",
        linewidth = 0.25
      ),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

submission_save_plot <- function(
    plot,
    filename,
    output_dir,
    width_mm = 180,
    height_mm = 145
) {
  ensure_dir(output_dir)
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  ggplot2::ggsave(
    pdf_path,
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  ggplot2::ggsave(
    png_path,
    plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

submission_panel_tag <- function(plot) {
  plot +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 10,
        family = "Arial"
      ),
      plot.tag.position = c(0, 1)
    )
}

submission_workflow_plot <- function() {
  nodes <- data.frame(
    x = 1:6,
    label = c(
      "Bulk discovery\nOA + OC",
      "Shared DEGs\n& pathways",
      "WGCNA +\nML selection",
      "External\nvalidation",
      "TCGA +\nsingle-cell",
      "Evidence-\nbounded synthesis"
    ),
    fill = c(
      submission_palette[["OA"]],
      submission_palette[["shared"]],
      submission_palette[["shared"]],
      submission_palette[["positive"]],
      submission_palette[["OC"]],
      "#374151"
    )
  )
  arrows <- data.frame(x = 1:5, xend = 2:6)
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x + 0.42, xend = xend - 0.42, y = 1, yend = 1),
      colour = "#9CA3AF",
      linewidth = 0.7,
      arrow = grid::arrow(
        length = grid::unit(2.2, "mm"),
        type = "closed"
      )
    ) +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = 1, label = label, fill = fill),
      colour = "white",
      size = 2.7,
      fontface = "bold",
      linewidth = 0,
      label.padding = grid::unit(2.4, "mm"),
      lineheight = 0.95
    ) +
    ggplot2::annotate(
      "label",
      x = 3.5,
      y = 0.47,
      label = "Bidirectional MR: null supplementary evidence",
      fill = "white",
      colour = submission_palette[["neutral"]],
      size = 2.4,
      linewidth = 0.25,
      label.r = grid::unit(1.5, "mm")
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.45, 6.55), ylim = c(0.25, 1.35)) +
    ggplot2::theme_void(base_family = "Arial")
}

submission_build_figure1 <- function(project_root, figure_dir, source_dir) {
  train <- utils::read.csv(file.path(
    project_root, "results", "tables", "dataset_qc_summary_train.csv"
  ))
  validation <- utils::read.csv(file.path(
    project_root, "results", "tables", "dataset_qc_summary_validation.csv"
  ))
  bulk <- rbind(train, validation)
  bulk_long <- rbind(
    transform(
      bulk[, c("dataset_id", "disease", "role", "normal")],
      group = "Normal",
      samples = normal
    )[, c("dataset_id", "disease", "role", "group", "samples")],
    transform(
      bulk[, c(
        "dataset_id", "disease", "role", "disease_samples"
      )],
      group = "Disease",
      samples = disease_samples
    )[, c("dataset_id", "disease", "role", "group", "samples")]
  )
  sc <- utils::read.csv(file.path(
    project_root, "results", "tables", "single_cell_dataset_status.csv"
  ))
  sc$disease <- ifelse(
    sc$dataset_id %in% c("GSE104782", "GSE169454", "GSE255460"),
    "OA",
    "OC"
  )
  sc$qc_pass_fraction <- sc$qc_pass / sc$cells
  safe_write_csv(bulk_long, file.path(source_dir, "Figure1_bulk_cohorts.csv"))
  safe_write_csv(sc, file.path(source_dir, "Figure1_single_cell_qc.csv"))

  p1 <- submission_workflow_plot()
  p2 <- ggplot2::ggplot(
    bulk_long,
    ggplot2::aes(
      x = dataset_id,
      y = samples,
      fill = group,
      alpha = role
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(~disease, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c(
      Normal = "#B9C0C8",
      Disease = "#374151"
    )) +
    ggplot2::scale_alpha_manual(values = c(train = 1, validation = 0.62)) +
    ggplot2::labs(
      title = "Bulk cohorts",
      x = NULL,
      y = "Samples",
      fill = NULL,
      alpha = "Role"
    ) +
    submission_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
  p3 <- ggplot2::ggplot(
    sc,
    ggplot2::aes(
      y = stats::reorder(dataset_id, cells),
      colour = disease
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1, xend = cells, yend = dataset_id),
      colour = "#D1D5DB",
      linewidth = 1.4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = cells),
      shape = 21,
      fill = "white",
      size = 2.4,
      stroke = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = qc_pass),
      size = 2.6
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::labs(
      title = "Single-cell count-level audit",
      x = "Cells (log scale)",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")

  figure <- p1 / (p2 | p3) +
    patchwork::plot_layout(heights = c(0.75, 1.5))
  figure <- submission_panel_tag(figure)
  submission_save_plot(figure, "Figure1_study_design", figure_dir)
}

submission_volcano_plot <- function(table, disease, hub_genes) {
  table$category <- "Other"
  table$category[
    table$adj.P.Val < 0.05 & table$logFC >= 1
  ] <- "Higher in disease"
  table$category[
    table$adj.P.Val < 0.05 & table$logFC <= -1
  ] <- "Lower in disease"
  table$minus_log10_fdr <- pmin(
    -log10(pmax(table$adj.P.Val, .Machine$double.xmin)),
    50
  )
  labels <- table[table$gene %in% hub_genes, , drop = FALSE]
  ggplot2::ggplot(
    table,
    ggplot2::aes(x = logFC, y = minus_log10_fdr, colour = category)
  ) +
    ggplot2::geom_point(alpha = 0.42, size = 0.45) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      colour = "#9CA3AF",
      linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      colour = "#9CA3AF",
      linewidth = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = gene),
      size = 2.2,
      colour = "#111827",
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.2
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Higher in disease" = submission_palette[["positive"]],
      "Lower in disease" = submission_palette[["negative"]],
      Other = "#C9CDD2"
    )) +
    ggplot2::labs(
      title = paste0(disease, " differential expression"),
      x = "log2 fold change",
      y = expression(-log[10]("FDR") ~ "(capped at 50)"),
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "none")
}

submission_build_figure2 <- function(project_root, figure_dir, source_dir) {
  differential <- submission_load_cache(project_root, "03_differential.rds")
  shared <- submission_load_cache(project_root, "04_shared.rds")
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  sensitivity_dir <- file.path(
    project_root, "results", "submission", "sensitivity"
  )
  threshold <- utils::read.csv(file.path(
    sensitivity_dir, "deg_threshold_sensitivity_summary.csv"
  ))
  oa <- differential$oa_train$table
  oc <- differential$oc_train$table
  common <- merge(
    oa[, c("gene", "logFC", "adj.P.Val")],
    oc[, c("gene", "logFC", "adj.P.Val")],
    by = "gene",
    suffixes = c("_OA", "_OC")
  )
  common$primary_shared <- common$gene %in% shared$genes
  common$hub <- common$gene %in% ml$final_genes
  safe_write_csv(oa, file.path(source_dir, "Figure2_OA_DEG.csv"))
  safe_write_csv(oc, file.path(source_dir, "Figure2_OC_DEG.csv"))
  safe_write_csv(common, file.path(source_dir, "Figure2_common_gene_effects.csv"))
  safe_write_csv(
    threshold,
    file.path(source_dir, "Figure2_DEG_threshold_sensitivity.csv")
  )

  p1 <- submission_volcano_plot(oa, "OA", ml$final_genes)
  p2 <- submission_volcano_plot(oc, "OC", ml$final_genes)
  p3 <- ggplot2::ggplot(
    common,
    ggplot2::aes(x = logFC_OA, y = logFC_OC)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B8BDC4", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8BDC4", linewidth = 0.3) +
    ggplot2::geom_point(colour = "#D4D7DB", size = 0.35, alpha = 0.38) +
    ggplot2::geom_point(
      data = common[common$primary_shared, ],
      colour = submission_palette[["shared"]],
      size = 0.75,
      alpha = 0.62
    ) +
    ggrepel::geom_text_repel(
      data = common[common$hub, ],
      ggplot2::aes(label = gene),
      size = 2.2,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.2
    ) +
    ggplot2::labs(
      title = "Cross-disease effect directions",
      x = "OA log2 fold change",
      y = "OC log2 fold change"
    ) +
    submission_theme()
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
      title = "Predeclared DEG threshold grid",
      x = "Absolute log2 fold-change threshold",
      y = "Shared genes",
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(figure, "Figure2_bulk_discovery", figure_dir)
}

submission_complete_selection_grid <- function(selection, hub_genes) {
  labels <- unique(selection[, c("dataset_id", "disease", "model")])
  grid <- merge(
    labels,
    data.frame(gene = hub_genes, stringsAsFactors = FALSE),
    all = TRUE
  )
  grid <- merge(
    grid,
    selection[, c(
      "dataset_id", "disease", "model", "gene", "selection_frequency"
    )],
    by = c("dataset_id", "disease", "model", "gene"),
    all.x = TRUE
  )
  grid$selection_frequency[is.na(grid$selection_frequency)] <- 0
  grid$model_label <- paste(grid$disease, grid$model, sep = " · ")
  grid
}

submission_build_figure3 <- function(project_root, figure_dir, source_dir) {
  sensitivity_dir <- file.path(
    project_root, "results", "submission", "sensitivity"
  )
  powers <- utils::read.csv(file.path(
    sensitivity_dir, "wgcna_soft_power_perturbation.csv"
  ))
  bootstrap <- utils::read.csv(file.path(
    sensitivity_dir, "wgcna_module_trait_bootstrap.csv"
  ))
  selection <- utils::read.csv(file.path(
    sensitivity_dir, "machine_learning_selection_frequency.csv"
  ))
  shared <- submission_load_cache(project_root, "04_shared.rds")
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  grid <- submission_complete_selection_grid(selection, ml$final_genes)
  hub_effects <- shared$table[
    match(ml$final_genes, shared$table$gene),
    c("gene", "logFC_OA", "logFC_OC")
  ]
  hub_effects <- reshape(
    hub_effects,
    varying = c("logFC_OA", "logFC_OC"),
    v.names = "logFC",
    timevar = "disease",
    times = c("OA", "OC"),
    direction = "long"
  )
  rownames(hub_effects) <- NULL
  safe_write_csv(powers, file.path(source_dir, "Figure3_WGCNA_powers.csv"))
  safe_write_csv(
    bootstrap,
    file.path(source_dir, "Figure3_WGCNA_bootstrap.csv")
  )
  safe_write_csv(
    grid,
    file.path(source_dir, "Figure3_ML_selection_frequency.csv")
  )
  safe_write_csv(
    hub_effects,
    file.path(source_dir, "Figure3_hub_gene_effects.csv")
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
    ggplot2::geom_point(
      ggplot2::aes(shape = is_primary),
      size = 2.2
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16)) +
    ggplot2::coord_cartesian(ylim = c(0.8, 1)) +
    ggplot2::labs(
      title = "Module–trait stability across powers",
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
    ggplot2::geom_point(
      ggplot2::aes(shape = is_primary),
      size = 2.2
    ) +
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
  grid$gene <- factor(grid$gene, levels = rev(ml$final_genes))
  p3 <- ggplot2::ggplot(
    grid,
    ggplot2::aes(x = model_label, y = gene, fill = selection_frequency)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%", 100 * selection_frequency)),
      size = 2
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#F3F4F6", "#9ECAE1", submission_palette[["OA"]]),
      limits = c(0, 1),
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Strict nested-CV feature-selection frequency",
      x = NULL,
      y = NULL,
      fill = "Frequency"
    ) +
    submission_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
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
      title = "Hub-gene discovery effects",
      x = NULL,
      y = NULL,
      fill = "log2 FC"
    ) +
    submission_theme()
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(figure, "Figure3_network_and_ml_stability", figure_dir)
}

submission_build_figure4 <- function(project_root, figure_dir, source_dir) {
  sensitivity_dir <- file.path(
    project_root, "results", "submission", "sensitivity"
  )
  fixed <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_direction_fixed_auc.csv"
  ))
  composite <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_signed_composite_score.csv"
  ))
  bias <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_direction_bias_audit.csv"
  ))
  consistency <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_cross_cohort_consistency.csv"
  ))
  hub_genes <- submission_load_cache(
    project_root, "07_machine_learning.rds"
  )$final_genes
  fixed$gene <- factor(fixed$gene, levels = rev(hub_genes))
  fixed$dataset_id <- factor(
    fixed$dataset_id,
    levels = c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
  )
  fixed$ci_support <- fixed$ci_lower > 0.5
  consistency$gene <- factor(consistency$gene, levels = rev(hub_genes))
  safe_write_csv(
    fixed,
    file.path(source_dir, "Figure4_direction_fixed_external_AUC.csv")
  )
  safe_write_csv(
    composite,
    file.path(source_dir, "Figure4_signed_composite_AUC.csv")
  )
  safe_write_csv(
    bias,
    file.path(source_dir, "Figure4_direction_bias_audit.csv")
  )
  safe_write_csv(
    consistency,
    file.path(source_dir, "Figure4_cross_cohort_consistency.csv")
  )

  p1 <- ggplot2::ggplot(
    fixed,
    ggplot2::aes(x = dataset_id, y = gene, fill = auc)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(sprintf("%.2f", auc), ifelse(ci_support, "*", ""))
      ),
      size = 2
    ) +
    ggplot2::scale_fill_gradient2(
      low = submission_palette[["negative"]],
      mid = "white",
      high = submission_palette[["positive"]],
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = "Direction-fixed gene-level validation",
      subtitle = "* 95% CI excludes AUC = 0.5",
      x = NULL,
      y = NULL,
      fill = "AUC"
    ) +
    submission_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  p2 <- ggplot2::ggplot(
    composite,
    ggplot2::aes(
      x = auc,
      y = stats::reorder(dataset_id, auc),
      colour = disease
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0.5,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      width = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      title = "Signed 10-gene composite score",
      x = "AUC (95% CI)",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p3 <- ggplot2::ggplot(
    bias,
    ggplot2::aes(
      x = fixed_direction_auc,
      y = legacy_auto_direction_auc,
      colour = disease
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.8) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Audit of automatic ROC direction",
      x = "Direction-fixed AUC",
      y = "Legacy auto-oriented AUC",
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p4 <- ggplot2::ggplot(
    consistency,
    ggplot2::aes(
      x = minimum_fixed_auc,
      y = gene,
      colour = disease,
      shape = all_at_or_above_robust_threshold
    )
  ) +
    ggplot2::geom_vline(
      xintercept = c(0.5, 0.6),
      linetype = c("solid", "dashed"),
      colour = "#B8BDC4",
      linewidth = 0.35
    ) +
    ggplot2::geom_point(size = 2.3, position = ggplot2::position_dodge(0.35)) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "Yes", `FALSE` = "No")
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Minimum AUC across two cohorts",
      x = "Minimum direction-fixed AUC",
      y = NULL,
      colour = NULL,
      shape = "Both ≥ 0.60"
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(figure, "Figure4_external_validation", figure_dir)
}

submission_read_cell_hub_evidence <- function(project_root) {
  utils::read.csv(file.path(
    project_root, "results", "tables", "single_cell_hub_gene_evidence.csv"
  ))
}

submission_build_figure5 <- function(project_root, figure_dir, source_dir) {
  evidence <- submission_read_cell_hub_evidence(project_root)
  status <- utils::read.csv(file.path(
    project_root, "results", "tables", "single_cell_dataset_status.csv"
  ))
  summary <- utils::read.csv(file.path(
    project_root,
    "results",
    "single_cell_downstream",
    "single_cell_downstream_summary.csv"
  ))
  hub_genes <- submission_load_cache(
    project_root, "07_machine_learning.rds"
  )$final_genes
  evidence$weighted_detected <- evidence$fraction_detected * evidence$cells
  weighted <- aggregate(
    cbind(weighted_detected, cells) ~ dataset_id + disease + gene,
    data = evidence,
    FUN = sum
  )
  weighted$fraction_detected <- weighted$weighted_detected / weighted$cells
  coverage <- aggregate(
    fraction_detected ~ dataset_id + disease + gene,
    data = transform(
      evidence,
      fraction_detected = as.integer(fraction_detected >= 0.25)
    ),
    FUN = sum
  )
  names(coverage)[[4L]] <- "cell_types_fraction_ge_0_25"

  pseudobulk_files <- c(
    GSE169454 = file.path(
      project_root,
      "results",
      "single_cell_downstream",
      "GSE169454",
      "pseudobulk_OA_vs_normal_FDR05.csv"
    ),
    GSE255460 = file.path(
      project_root,
      "results",
      "single_cell_downstream",
      "GSE255460",
      "pseudobulk_differential_FDR05.csv"
    )
  )
  pseudobulk <- do.call(rbind, lapply(names(pseudobulk_files), function(id) {
    table <- utils::read.csv(pseudobulk_files[[id]])
    table$dataset_id <- id
    table[table$gene %in% hub_genes, , drop = FALSE]
  }))
  status <- merge(
    status,
    summary[, c("dataset_id", "disease", "cell_types")],
    by = "dataset_id",
    all.x = TRUE
  )
  status$qc_pass_fraction <- status$qc_pass / status$cells
  weighted$gene <- factor(weighted$gene, levels = hub_genes)
  coverage$gene <- factor(coverage$gene, levels = hub_genes)
  safe_write_csv(
    weighted,
    file.path(source_dir, "Figure5_hub_detection_by_dataset.csv")
  )
  safe_write_csv(
    coverage,
    file.path(source_dir, "Figure5_hub_cell_type_coverage.csv")
  )
  safe_write_csv(
    pseudobulk,
    file.path(source_dir, "Figure5_hub_pseudobulk_evidence.csv")
  )
  safe_write_csv(status, file.path(source_dir, "Figure5_single_cell_status.csv"))

  p1 <- ggplot2::ggplot(
    status,
    ggplot2::aes(
      x = qc_pass_fraction,
      y = stats::reorder(dataset_id, qc_pass_fraction),
      fill = disease
    )
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          scales::comma(qc_pass),
          " / ",
          scales::comma(cells)
        )
      ),
      hjust = 1.05,
      colour = "white",
      size = 2.1
    ) +
    ggplot2::scale_fill_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Count-level QC retention",
      x = "QC-pass fraction",
      y = NULL,
      fill = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p2 <- ggplot2::ggplot(
    weighted,
    ggplot2::aes(
      x = gene,
      y = dataset_id,
      size = fraction_detected,
      fill = fraction_detected
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "#374151", stroke = 0.25) +
    ggplot2::scale_size(range = c(0.5, 5), labels = scales::percent) +
    ggplot2::scale_fill_gradient(
      low = "#F3F4F6",
      high = submission_palette[["shared"]],
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Hub-gene detection across cell atlases",
      subtitle = "Cell-count-weighted fraction detected",
      x = NULL,
      y = NULL,
      size = "Detected",
      fill = "Detected"
    ) +
    submission_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  p3 <- ggplot2::ggplot(
    coverage,
    ggplot2::aes(
      x = gene,
      y = dataset_id,
      fill = cell_types_fraction_ge_0_25
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = cell_types_fraction_ge_0_25),
      size = 2.1
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F3F4F6",
      high = submission_palette[["OA"]]
    ) +
    ggplot2::labs(
      title = "Cell types with ≥25% detection",
      x = NULL,
      y = NULL,
      fill = "Cell types"
    ) +
    submission_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )
  if (nrow(pseudobulk) > 0L) {
    pseudobulk$label <- paste(
      pseudobulk$dataset_id,
      pseudobulk$cell_type,
      pseudobulk$contrast,
      sep = " · "
    )
    p4 <- ggplot2::ggplot(
      pseudobulk,
      ggplot2::aes(
        x = logFC,
        y = stats::reorder(label, logFC),
        colour = gene
      )
    ) +
      ggplot2::geom_vline(xintercept = 0, colour = "#B8BDC4") +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
    title = "Significant pseudobulk hub effects",
        x = "log2 fold change",
        y = NULL,
        colour = "Gene"
      ) +
      submission_theme() +
      ggplot2::theme(legend.position = "right")
  } else {
    p4 <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No final hub gene met FDR < 0.05\nin eligible pseudobulk contrasts",
        size = 3
      ) +
      ggplot2::labs(title = "Pseudobulk hub-gene effects") +
      ggplot2::theme_void(base_family = "Arial")
  }
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure5_single_cell_localization",
    figure_dir,
    height_mm = 155
  )
}

submission_hallmark_label <- function(value) {
  value <- sub("^HALLMARK_", "", value)
  tools::toTitleCase(tolower(gsub("_", " ", value)))
}

submission_build_figure6 <- function(project_root, figure_dir, source_dir) {
  table_dir <- file.path(project_root, "results", "tables")
  oa_gsea <- utils::read.csv(file.path(
    table_dir, "GSEA_OA_GSE114007_hallmark.csv"
  ))
  oc_gsea <- utils::read.csv(file.path(
    table_dir, "GSEA_OC_GSE18520_hallmark.csv"
  ))
  oa_gsea$disease <- "OA"
  oc_gsea$disease <- "OC"
  gsea <- rbind(oa_gsea, oc_gsea)
  top_ids <- names(head(sort(
    tapply(abs(gsea$NES), gsea$ID, max, na.rm = TRUE),
    decreasing = TRUE
  ), 12))
  gsea <- gsea[gsea$ID %in% top_ids, ]
  gsea$pathway <- submission_hallmark_label(gsea$ID)
  order_pathways <- names(sort(
    tapply(abs(gsea$NES), gsea$pathway, max),
    decreasing = FALSE
  ))
  gsea$pathway <- factor(gsea$pathway, levels = order_pathways)

  immune_oa <- utils::read.csv(file.path(
    table_dir, "immune_OA_GSE114007_statistics.csv"
  ))
  immune_oc <- utils::read.csv(file.path(
    table_dir, "immune_OC_GSE18520_statistics.csv"
  ))
  immune <- rbind(immune_oa, immune_oc)
  top_signatures <- names(head(sort(
    tapply(immune$fdr, immune$signature, min, na.rm = TRUE)
  ), 10))
  immune <- immune[immune$signature %in% top_signatures, ]
  immune$signature <- gsub("_", " ", immune$signature)
  immune$significant <- immune$fdr < 0.05

  sensitivity_dir <- file.path(
    project_root, "results", "submission", "sensitivity"
  )
  cox <- utils::read.csv(file.path(
    sensitivity_dir, "tcga_cox_model_sensitivity.csv"
  ))
  cox <- cox[
    cox$model == "selected_gene_multivariable" |
      (cox$model == "risk_adjusted_for_age_and_stage" &
        cox$term == "risk_z"),
    ,
    drop = FALSE
  ]
  cox$label <- c(
    gsub("^z_", "", cox$term[cox$model == "selected_gene_multivariable"]),
    rep(
      "Risk score (age/stage adjusted)",
      sum(cox$model == "risk_adjusted_for_age_and_stage")
    )
  )
  bootstrap <- utils::read.csv(file.path(
    sensitivity_dir, "tcga_optimism_bootstrap_replicates.csv"
  ))
  optimism <- utils::read.csv(file.path(
    sensitivity_dir, "tcga_optimism_bootstrap_summary.csv"
  ))
  safe_write_csv(gsea, file.path(source_dir, "Figure6_hallmark_GSEA.csv"))
  safe_write_csv(immune, file.path(source_dir, "Figure6_immune_signatures.csv"))
  safe_write_csv(cox, file.path(source_dir, "Figure6_TCGA_cox.csv"))
  safe_write_csv(
    bootstrap,
    file.path(source_dir, "Figure6_TCGA_bootstrap.csv")
  )

  p1 <- ggplot2::ggplot(
    gsea,
    ggplot2::aes(
      x = NES,
      y = pathway,
      colour = disease,
      size = -log10(pmax(p.adjust, 1e-12))
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8BDC4") +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_size(range = c(1.4, 4.5), guide = "none") +
    ggplot2::labs(
      title = "Hallmark pathway activity",
      x = "Normalized enrichment score",
      y = NULL,
      colour = NULL,
      size = expression(-log[10]("FDR"))
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p2 <- ggplot2::ggplot(
    immune,
    ggplot2::aes(
      x = median_difference,
      y = stats::reorder(signature, median_difference),
      colour = disease,
      shape = significant
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8BDC4") +
    ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(0.35)) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Immune signature differences",
      x = "Disease − normal median score",
      y = NULL,
      colour = NULL,
      shape = "FDR < 0.05"
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p3 <- ggplot2::ggplot(
    cox,
    ggplot2::aes(
      x = hazard_ratio,
      y = stats::reorder(label, hazard_ratio)
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      width = 0.14,
      linewidth = 0.7,
      colour = submission_palette[["OC"]]
    ) +
    ggplot2::geom_point(size = 2.4, colour = submission_palette[["OC"]]) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "TCGA-OV exploratory survival associations",
      x = "Hazard ratio (95% CI; log scale)",
      y = NULL
    ) +
    submission_theme()
  p4 <- ggplot2::ggplot(
    bootstrap,
    ggplot2::aes(x = test_cindex)
  ) +
    ggplot2::geom_histogram(
      bins = 25,
      fill = "#BFD7EA",
      colour = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_vline(
      xintercept = optimism$optimism_corrected_cindex,
      colour = submission_palette[["OC"]],
      linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = optimism$optimism_corrected_cindex - 0.0005,
      y = Inf,
      label = sprintf(
        "Corrected C = %.3f",
        optimism$optimism_corrected_cindex
      ),
      hjust = 1,
      vjust = 1.4,
      size = 2.4,
      colour = submission_palette[["OC"]]
    ) +
    ggplot2::labs(
      title = "TCGA bootstrap C-index",
      x = "Test C-index",
      y = "Bootstrap models"
    ) +
    submission_theme()
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(figure, "Figure6_functional_and_tcga_context", figure_dir)
}

submission_build_supplementary1 <- function(
    project_root,
    figure_dir,
    source_dir
) {
  sensitivity_dir <- file.path(
    project_root, "results", "submission", "sensitivity"
  )
  deg <- utils::read.csv(file.path(
    sensitivity_dir, "deg_threshold_sensitivity_summary.csv"
  ))
  loso <- utils::read.csv(file.path(
    sensitivity_dir, "wgcna_module_trait_leave_one_out.csv"
  ))
  selection <- utils::read.csv(file.path(
    sensitivity_dir, "machine_learning_selection_frequency.csv"
  ))
  tcga <- utils::read.csv(file.path(
    sensitivity_dir, "tcga_lasso_selection_stability.csv"
  ))
  safe_write_csv(deg, file.path(source_dir, "SupplementaryFigure1_DEG.csv"))
  safe_write_csv(loso, file.path(source_dir, "SupplementaryFigure1_WGCNA.csv"))
  safe_write_csv(
    selection,
    file.path(source_dir, "SupplementaryFigure1_ML.csv")
  )
  safe_write_csv(tcga, file.path(source_dir, "SupplementaryFigure1_TCGA.csv"))

  deg$fdr <- factor(deg$fdr_threshold)
  p1 <- ggplot2::ggplot(
    deg,
    ggplot2::aes(
      x = absolute_log2fc_threshold,
      y = retained_hub_count,
      colour = fdr,
      group = fdr
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(
      `0.01` = submission_palette[["OA"]],
      `0.05` = submission_palette[["OC"]]
    )) +
    ggplot2::scale_y_continuous(breaks = 0:10, limits = c(0, 10.5)) +
    ggplot2::labs(
      title = "Hub retention across DEG rules",
      x = "Absolute log2 FC threshold",
      y = "Retained hubs",
      colour = "FDR"
    ) +
    submission_theme()
  p2 <- ggplot2::ggplot(
    loso,
    ggplot2::aes(x = dataset_id, y = correlation, fill = disease)
  ) +
    ggplot2::geom_boxplot(width = 0.55, outlier.size = 0.8) +
    ggplot2::scale_fill_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::labs(
      title = "LOSO module–trait correlations",
      x = NULL,
      y = "Correlation",
      fill = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "none")
  top <- do.call(rbind, lapply(
    split(selection, interaction(selection$dataset_id, selection$model)),
    function(table) {
      head(table[order(table$selection_frequency, decreasing = TRUE), ], 12)
    }
  ))
  top$label <- paste(top$disease, top$model, sep = " · ")
  p3 <- ggplot2::ggplot(
    top,
    ggplot2::aes(
      x = selection_frequency,
      y = stats::reorder(gene, selection_frequency),
      colour = disease
    )
  ) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~label, scales = "free_y", ncol = 2) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(
      title = "Most frequently selected discovery features",
      x = "Selection frequency",
      y = NULL,
      colour = NULL
    ) +
    submission_theme(7.5) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  tcga$gene <- factor(tcga$gene)
  p4 <- ggplot2::ggplot(
    tcga,
    ggplot2::aes(
      x = selection_frequency,
      y = stats::reorder(gene, selection_frequency),
      colour = lambda_rule
    )
  ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~lambda_rule, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::scale_colour_manual(values = c(
      "lambda.1se" = submission_palette[["OA"]],
      "lambda.min" = submission_palette[["OC"]]
    )) +
    ggplot2::labs(
      title = "TCGA LASSO bootstrap selection stability",
      x = "Selection frequency",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "SupplementaryFigure1_sensitivity_details",
    figure_dir,
    height_mm = 155
  )
}

submission_build_supplementary2 <- function(
    project_root,
    figure_dir,
    source_dir
) {
  mr_dir <- file.path(project_root, "results", "mr")
  estimates <- utils::read.csv(file.path(
    mr_dir, "MR_combined_estimates.csv"
  ), check.names = FALSE)
  estimates$direction <- ifelse(
    estimates$id.exposure == "ebi-a-GCST007092",
    "OA → OC",
    "OC → OA"
  )
  estimates$or <- exp(estimates$b)
  estimates$lower <- exp(estimates$b - 1.96 * estimates$se)
  estimates$upper <- exp(estimates$b + 1.96 * estimates$se)
  estimates$method <- factor(
    estimates$method,
    levels = rev(c(
      "Inverse variance weighted",
      "Weighted median",
      "MR Egger",
      "Weighted mode",
      "Simple mode"
    ))
  )
  heterogeneity <- rbind(
    utils::read.csv(file.path(
      mr_dir,
      "MR_ebi-a-GCST007092__ieu-a-1120_heterogeneity.csv"
    )),
    utils::read.csv(file.path(
      mr_dir,
      "MR_ieu-a-1120__ebi-a-GCST007092_heterogeneity.csv"
    ))
  )
  heterogeneity$direction <- ifelse(
    heterogeneity$id.exposure == "ebi-a-GCST007092",
    "OA → OC",
    "OC → OA"
  )
  heterogeneity <- heterogeneity[
    heterogeneity$method == "Inverse variance weighted",
  ]
  pleiotropy <- rbind(
    utils::read.csv(file.path(
      mr_dir,
      "MR_ebi-a-GCST007092__ieu-a-1120_pleiotropy.csv"
    )),
    utils::read.csv(file.path(
      mr_dir,
      "MR_ieu-a-1120__ebi-a-GCST007092_pleiotropy.csv"
    ))
  )
  pleiotropy$direction <- ifelse(
    pleiotropy$id.exposure == "ebi-a-GCST007092",
    "OA → OC",
    "OC → OA"
  )
  diagnostics <- merge(
    heterogeneity[, c("direction", "Q", "Q_df", "Q_pval")],
    pleiotropy[, c("direction", "egger_intercept", "pval")],
    by = "direction"
  )
  safe_write_csv(
    estimates,
    file.path(source_dir, "SupplementaryFigure2_MR_estimates.csv")
  )
  safe_write_csv(
    diagnostics,
    file.path(source_dir, "SupplementaryFigure2_MR_diagnostics.csv")
  )
  p1 <- ggplot2::ggplot(
    estimates,
    ggplot2::aes(
      x = or,
      y = method,
      colour = direction
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#9CA3AF"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0.13,
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::geom_point(
      size = 2.2,
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::scale_colour_manual(values = c(
      "OA → OC" = submission_palette[["OA"]],
      "OC → OA" = submission_palette[["OC"]]
    )) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Bidirectional MR estimates",
      x = "Odds ratio (95% CI; log scale)",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  diagnostic_text <- paste(
    sprintf(
      "%s\nIVW Q = %.2f (df %d), P = %s\nEgger intercept P = %.3f",
      diagnostics$direction,
      diagnostics$Q,
      diagnostics$Q_df,
      formatC(diagnostics$Q_pval, format = "g", digits = 3),
      diagnostics$pval
    ),
    collapse = "\n\n"
  )
  p2 <- ggplot2::ggplot() +
    ggplot2::annotate(
      "label",
      x = 0,
      y = 0.80,
      label = diagnostic_text,
      hjust = 0,
      vjust = 1,
      size = 2.4,
      family = "Arial",
      lineheight = 1.05,
      linewidth = 0.25,
      fill = "#F8FAFC"
    ) +
    ggplot2::annotate(
      "text",
      x = 0,
      y = -0.31,
      label = paste(
        "Both IVW CIs include OR = 1.",
        "Reverse MR: heterogeneity + 4 PRESSO outliers.",
        "No detected causal evidence;",
        "not proof of absence.",
        sep = "\n"
      ),
      hjust = 0,
      vjust = 0,
      size = 2.2,
      family = "Arial",
      colour = "#374151"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.38, 0.85)) +
    ggplot2::labs(title = "Sensitivity diagnostics and boundary") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9)
    )
  figure <- submission_panel_tag(p1 | p2)
  submission_save_plot(
    figure,
    "SupplementaryFigure2_negative_bidirectional_MR",
    figure_dir,
    height_mm = 90
  )
}

submission_read_umap <- function(path, dataset_id, max_cells, seed) {
  require_namespace("data.table", "single-cell UMAP montage")
  table <- data.table::fread(
    path,
    select = intersect(
      c(
        "UMAP1", "UMAP2", "cell_type", "celltype", "published_cell_type",
        "analysis_cluster"
      ),
      names(data.table::fread(path, nrows = 0L))
    ),
    data.table = FALSE
  )
  label_column <- intersect(
    c("cell_type", "celltype", "published_cell_type", "analysis_cluster"),
    names(table)
  )[[1L]]
  table$label <- as.character(table[[label_column]])
  table <- table[
    is.finite(table$UMAP1) & is.finite(table$UMAP2) & nzchar(table$label),
  ]
  if (nrow(table) > max_cells) {
    set.seed(seed)
    table <- table[sample(seq_len(nrow(table)), max_cells), ]
  }
  table$dataset_id <- dataset_id
  table[, c("dataset_id", "UMAP1", "UMAP2", "label")]
}

submission_build_supplementary3 <- function(
    project_root,
    figure_dir,
    source_dir
) {
  base <- file.path(project_root, "results", "single_cell_downstream")
  files <- c(
    GSE104782 = file.path(base, "GSE104782", "cell_annotations.tsv.gz"),
    GSE169454 = file.path(base, "GSE169454", "cell_annotations.tsv.gz"),
    GSE255460 = file.path(base, "GSE255460", "visualization_subsample.tsv.gz"),
    GSE154600 = file.path(base, "GSE154600", "cell_annotations.tsv.gz"),
    GSE180661 = file.path(base, "GSE180661", "cell_annotations.tsv.gz")
  )
  max_cells <- c(
    GSE104782 = 2000,
    GSE169454 = 30000,
    GSE255460 = 30000,
    GSE154600 = 30000,
    GSE180661 = 40000
  )
  umap <- do.call(rbind, lapply(seq_along(files), function(index) {
    submission_read_umap(
      files[[index]],
      names(files)[[index]],
      max_cells[[index]],
      20260726L + index
    )
  }))
  safe_write_csv(
    umap,
    file.path(source_dir, "SupplementaryFigure3_UMAP_subsamples.csv")
  )
  label_levels <- sort(unique(umap$label))
  colours <- grDevices::hcl.colors(
    length(label_levels),
    palette = "Dark 3"
  )
  names(colours) <- label_levels
  centroids <- aggregate(
    cbind(UMAP1, UMAP2) ~ dataset_id + label,
    data = umap,
    FUN = stats::median
  )
  plot <- ggplot2::ggplot(
    umap,
    ggplot2::aes(x = UMAP1, y = UMAP2, colour = label)
  ) +
    ggplot2::geom_point(size = 0.12, alpha = 0.45) +
    ggrepel::geom_text_repel(
      data = centroids,
      ggplot2::aes(label = label),
      size = 1.8,
      colour = "#111827",
      min.segment.length = 0,
      max.overlaps = 60,
      segment.size = 0.2
    ) +
    ggplot2::facet_wrap(~dataset_id, scales = "free", ncol = 3) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(
      title = "Dataset-specific single-cell embeddings and conservative labels",
      subtitle = "Large atlases are deterministically subsampled for display",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    submission_theme(7.5) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
  submission_save_plot(
    plot,
    "SupplementaryFigure3_single_cell_UMAPs",
    figure_dir,
    height_mm = 135
  )
}

submission_write_figure_legends <- function(figure_dir) {
  legends <- c(
    "# Figure legends",
    "",
    "## Figure 1. Study design and audited data resources",
    "",
    "**A,** Analysis flow from disease-specific bulk discovery through evidence-bounded synthesis. Bidirectional Mendelian randomization (MR) is retained as a null supplementary analysis. **B,** Bulk discovery and external-validation sample counts; opacity distinguishes training from validation cohorts. **C,** Total and QC-pass single-cell counts for five datasets on a logarithmic scale. Open points denote input cells and filled points denote QC-pass cells.",
    "",
    "## Figure 2. Cross-disease bulk transcriptomic discovery",
    "",
    "**A–B,** Differential-expression volcano plots for OA (GSE114007) and OC (GSE18520). Dashed lines indicate FDR 0.05 and absolute log2 fold change 1; final hub genes are labelled. **C,** OA and OC log2 fold changes for genes measured in both training datasets; purple points are the 286 primary shared DEGs. **D,** Shared-gene counts under the six predeclared FDR and absolute log2 fold-change rules. No rule was selected after viewing the results.",
    "",
    "## Figure 3. Network and machine-learning stability",
    "",
    "**A,** Absolute correlation between the most trait-associated WGCNA module and disease status under soft-power perturbation; diamonds mark the primary powers. **B,** Fraction of primary-module genes retained under each perturbation. **C,** Feature-selection frequency across 50 repeated five-fold cross-validation runs (250 outer training folds per disease/model). Candidate screening preceded resampling, so these frequencies and the perfect internal AUCs are exploratory rather than independent performance estimates. **D,** Disease-specific discovery log2 fold changes for the final ten genes.",
    "",
    "## Figure 4. Direction-fixed external validation",
    "",
    "**A,** Gene-level AUCs after fixing the expected expression direction from the corresponding training dataset; asterisks mark DeLong 95% intervals excluding 0.5. **B,** AUCs and 95% intervals for a transparent signed ten-gene score. **C,** Comparison of fixed-direction and legacy automatically oriented AUCs; points above the diagonal quantify optimism introduced by selecting ROC direction in each validation cohort. **D,** Minimum direction-fixed AUC across two external cohorts per disease; filled symbols indicate AUC at least 0.60 in both cohorts.",
    "",
    "## Figure 5. Single-cell localization of hub-gene evidence",
    "",
    "**A,** QC-pass fractions after count-level audit and mandatory doublet detection where supported. **B,** Cell-count-weighted fraction of cells expressing each hub gene in each dataset. **C,** Number of annotated cell types with at least 25% detection. **D,** Final hub genes reaching FDR below 0.05 in eligible pseudobulk contrasts; datasets or cell types without biological replication were not tested.",
    "",
    "## Figure 6. Functional and exploratory prognostic context",
    "",
    "**A,** Normalized enrichment scores for the twelve Hallmark sets with the largest absolute enrichment across OA and OC training data. **B,** Disease-minus-normal differences in rank-based immune-signature scores; filled points indicate FDR below 0.05. **C,** Exploratory TCGA-OV Cox estimates for the three selected genes and the age/stage-adjusted continuous risk score. **D,** C-index distribution when 200 bootstrap-fitted models were evaluated in the original cohort; the vertical line is the optimism-corrected C-index. Proportional-hazards diagnostics and a time-varying coefficient sensitivity model are reported separately.",
    "",
    "## Supplementary Figure 1. Prespecified non-MR sensitivity analyses",
    "",
    "**A,** Retention of final hub genes across DEG thresholds. **B,** Leave-one-sample-out correlations for fixed primary WGCNA module eigengenes. **C,** Most frequently selected genes in repeated outer cross-validation. **D,** TCGA LASSO gene-selection frequencies across 200 bootstrap samples under lambda.1se and lambda.min.",
    "",
    "## Supplementary Figure 2. Negative bidirectional MR results",
    "",
    "**A,** Odds-ratio estimates across five MR estimators for OA to OC and OC to OA. **B,** Heterogeneity and Egger-intercept diagnostics. Both IVW estimates included the null; reverse MR showed significant heterogeneity and MR-PRESSO outliers. The result is reported as no detected causal evidence, not proof of absence.",
    "",
    "## Supplementary Figure 3. Single-cell embeddings",
    "",
    "Dataset-specific UMAP or released coordinates with author-provided or conservative reference-transfer labels. Large datasets were deterministically subsampled only for display; all quantitative summaries used the complete QC-pass cells. OA and OC datasets were not forced into a shared latent space."
  )
  write_utf8(legends, file.path(figure_dir, "figure_legends.md"))
}

run_submission_figures <- function(project_root) {
  for (
    package in c(
      "ggplot2", "patchwork", "ggrepel", "scales",
      "ragg", "data.table"
    )
  ) {
    require_namespace(package, "submission figure assembly")
  }
  figure_dir <- ensure_dir(file.path(
    project_root,
    "results",
    "submission",
    "figures"
  ))
  source_dir <- ensure_dir(file.path(figure_dir, "source_data"))
  log_info("Building unified main figures.")
  submission_build_figure1(project_root, figure_dir, source_dir)
  submission_build_figure2(project_root, figure_dir, source_dir)
  submission_build_figure3(project_root, figure_dir, source_dir)
  submission_build_figure4(project_root, figure_dir, source_dir)
  submission_build_figure5(project_root, figure_dir, source_dir)
  submission_build_figure6(project_root, figure_dir, source_dir)
  log_info("Building unified supplementary figures.")
  submission_build_supplementary1(project_root, figure_dir, source_dir)
  submission_build_supplementary2(project_root, figure_dir, source_dir)
  submission_build_supplementary3(project_root, figure_dir, source_dir)
  submission_write_figure_legends(figure_dir)
  safe_write_csv(
    data.frame(
      setting = c(
        "target_journal",
        "width_mm",
        "png_dpi",
        "font_family",
        "palette",
        "vector_format"
      ),
      value = c(
        "not yet specified; provisional general biomedical layout",
        "180",
        "300",
        "Arial",
        "Okabe-Ito-derived disease colours plus redundant shapes/text",
        "PDF"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(figure_dir, "figure_style_manifest.csv")
  )
  log_info("Submission figures completed.")
  invisible(figure_dir)
}
