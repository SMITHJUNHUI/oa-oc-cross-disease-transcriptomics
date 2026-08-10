v30_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v30"))
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

v30_copy_pair <- function(source_dir, source_stem, target_dir, target_stem) {
  for (extension in c("png", "pdf")) {
    source <- file.path(source_dir, paste0(source_stem, ".", extension))
    target <- file.path(target_dir, paste0(target_stem, ".", extension))
    if (!file.exists(source)) {
      stop("Missing V3.0 figure source: ", source, call. = FALSE)
    }
    if (!file.copy(source, target, overwrite = TRUE, copy.date = TRUE)) {
      stop("Could not copy V3.0 figure: ", source, call. = FALSE)
    }
  }
}

v30_prepare_baseline <- function(project_root, paths) {
  source_root <- file.path(project_root, "results", "submission_v24")
  if (!dir.exists(source_root)) {
    stop("V2.4 baseline outputs are required for V3.0.", call. = FALSE)
  }
  v22_copy_tree(source_root, paths$root)

  v30_copy_pair(
    source_root,
    "figures/Figure3_network_and_ml_stability",
    paths$figures,
    "SupplementaryFigure8_candidate_evidence_stability"
  )

  existing_main <- list.files(
    paths$figures,
    pattern = "^Figure[1-6]_.*\\.(png|pdf)$",
    full.names = TRUE
  )
  if (length(existing_main) > 0L) {
    unlink(existing_main)
  }
  v30_copy_pair(
    source_root,
    "figures/Figure2_bulk_discovery",
    paths$figures,
    "Figure2_gene_direction"
  )
  invisible(paths)
}

v30_build_figure1 <- function(paths) {
  top_nodes <- data.frame(
    x = c(1.0, 3.0, 5.0),
    y = c(3.45, 3.45, 3.45),
    label = c(
      "OA cartilage / synovium\ntranscriptomes",
      "Partial shared\ntranscriptional landscape",
      "Ovarian tumor / reference\ntranscriptomes"
    ),
    fill = c("#DCEAF7", "#E7E2F3", "#FBE4D5"),
    stringsAsFactors = FALSE
  )
  questions <- data.frame(
    x = c(1, 3, 5, 1, 3, 5),
    y = c(2.15, 2.15, 2.15, 1.05, 1.05, 1.05),
    label = c(
      "Gene-level overlap\nand direction",
      "Pathway-level\nconcordance / divergence",
      "Cell-state-specific\nlocalization",
      "Cross-cohort\nmolecular separability",
      "Bidirectional genetic-\nliability assessment",
      "Focused regulatory\nhypotheses"
    ),
    domain = c(
      "Transcriptome", "Transcriptome", "Cell context",
      "Validation task", "Genetic boundary", "Exploratory"
    ),
    stringsAsFactors = FALSE
  )
  domain_colours <- c(
    Transcriptome = "#DCEAF7",
    `Cell context` = "#DDF0E3",
    `Validation task` = "#FBE4D5",
    `Genetic boundary` = "#E5E7EB",
    Exploratory = "#F4EDC9"
  )
  top_arrows <- data.frame(
    x = c(1.62, 4.38),
    xend = c(2.38, 3.62),
    y = 3.45,
    yend = 3.45
  )
  branch_lines <- rbind(
    data.frame(x = 3, xend = 3, y = 3.08, yend = 2.72),
    data.frame(x = 1, xend = 5, y = 2.72, yend = 2.72),
    data.frame(x = c(1, 3, 5), xend = c(1, 3, 5), y = 2.72, yend = 2.48),
    data.frame(x = 3, xend = 3, y = 1.82, yend = 1.60),
    data.frame(x = 1, xend = 5, y = 1.60, yend = 1.60),
    data.frame(x = c(1, 3, 5), xend = c(1, 3, 5), y = 1.60, yend = 1.38)
  )
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = top_arrows,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"),
      linewidth = 0.65,
      colour = "#667085"
    ) +
    ggplot2::geom_segment(
      data = branch_lines,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.48,
      colour = "#98A2B3"
    ) +
    ggplot2::geom_tile(
      data = top_nodes,
      ggplot2::aes(x = x, y = y),
      width = 1.15,
      height = 0.64,
      fill = top_nodes$fill,
      colour = "#475467",
      linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = top_nodes,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.75,
      lineheight = 0.95,
      fontface = c("plain", "bold", "plain"),
      colour = "#1F2937"
    ) +
    ggplot2::geom_tile(
      data = questions,
      ggplot2::aes(x = x, y = y, fill = domain),
      width = 1.35,
      height = 0.64,
      colour = "#667085",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = questions,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.55,
      lineheight = 0.94,
      colour = "#1F2937"
    ) +
    ggplot2::annotate(
      "label",
      x = 3,
      y = 0.22,
      label = paste0(
        "Shared membership does not imply shared direction, pathway state, ",
        "cellular meaning, clinical utility, or causality."
      ),
      size = 2.65,
      fontface = "bold",
      fill = "#FFF7ED",
      colour = "#9A3412",
      linewidth = 0.3,
      label.padding = grid::unit(1.4, "mm")
    ) +
    ggplot2::scale_fill_manual(values = domain_colours) +
    ggplot2::coord_cartesian(xlim = c(0.20, 5.80), ylim = c(-0.05, 3.90)) +
    ggplot2::labs(
      title = "Systems-level question",
      subtitle = paste0(
        "To what extent do OA and OC overlap transcriptionally, and how does ",
        "interpretation change across analytical contexts?"
      )
    ) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 10, colour = "#111827"
      ),
      plot.subtitle = ggplot2::element_text(size = 7.3, colour = "#475467"),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 2, 4)
    )

  resources <- data.frame(
    resource = factor(
      c(
        "Bulk discovery", "External bulk", "Single-cell",
        "HPA / TCGA-OV", "GWAS", "Regulatory resources"
      ),
      levels = rev(c(
        "Bulk discovery", "External bulk", "Single-cell",
        "HPA / TCGA-OV", "GWAS", "Regulatory resources"
      ))
    ),
    scale = c(
      "2 cohorts\n101 samples",
      "4 cohorts\n112 samples",
      "5 datasets\n1,025,361 QC-pass cells",
      "Normal reference\n307 TCGA-OV samples",
      "2 GWAS\n21 and 11 instruments",
      "KnockTF target sets\nmiRTarBase interactions"
    ),
    role = c(
      "Disease-specific discovery",
      "Secondary molecular separability",
      "Dataset-specific localization",
      "Tissue/composition\nand survival context",
      "Bidirectional genetic-\nliability assessment",
      "Focused upstream hypotheses"
    ),
    boundary = c(
      "Association",
      "Not diagnostic validation",
      "Not a shared cell state",
      "Not tissue specificity\nor prognosis validation",
      "Not shared heritability\nor genetic correlation",
      "Not TF/miRNA activity\nor causal regulation"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    top_nodes,
    file.path(paths$source, "Figure1_disease_tracks.csv")
  )
  safe_write_csv(
    questions,
    file.path(paths$source, "Figure1_analytic_questions.csv")
  )
  safe_write_csv(
    resources,
    file.path(paths$source, "Figure1_resource_roles.csv")
  )
  p2 <- ggplot2::ggplot(resources, ggplot2::aes(y = resource)) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 1),
      width = 0.97, height = 0.78,
      fill = "#F8FAFC", colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 2),
      width = 0.97, height = 0.78,
      fill = "#F9FAFB", colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 3),
      width = 0.97, height = 0.78,
      fill = "#FFF9F2", colour = "#D0D5DD"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 1, label = scale),
      size = 2.25, lineheight = 0.95, colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2, label = role),
      size = 2.25, lineheight = 0.95, colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 3, label = boundary),
      size = 2.20, lineheight = 0.95, colour = "#9A3412"
    ) +
    ggplot2::annotate(
      "text",
      x = c(1, 2, 3), y = 6.45,
      label = c("Scale", "Analytic role", "Inference boundary"),
      fontface = "bold", size = 3.0, colour = "#344054"
    ) +
    ggplot2::scale_x_continuous(limits = c(0.48, 3.52), breaks = NULL) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(0.25, 0.95))
    ) +
    ggplot2::labs(
      title = "Audited evidence layers and boundaries",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = "bold", colour = "#344054", size = 7.0
      ),
      plot.title = ggplot2::element_text(
        face = "bold", size = 10, colour = "#111827"
      ),
      plot.margin = ggplot2::margin(2, 4, 4, 4)
    )
  figure <- submission_panel_tag(
    p1 / p2 + patchwork::plot_layout(heights = c(1.35, 1))
  )
  submission_save_plot(
    figure,
    "Figure1_systems_framework",
    paths$figures,
    height_mm = 190
  )
}

v30_update_pathway_table <- function(paths) {
  path <- file.path(
    paths$tables,
    "Table_S18_Hallmark_pathway_direction_matrix.csv"
  )
  table <- utils::read.csv(path, check.names = FALSE)
  table$paired_direction_index_definition <- paste0(
    "OA NES multiplied by OC NES; positive indicates matching signs and ",
    "negative indicates opposite signs. Magnitude is descriptive."
  )
  table$quantitative_interpretation <- ifelse(
    table$paired_direction_index > 0,
    "matching NES signs",
    "opposite NES signs"
  )
  table$inference_boundary <- paste0(
    "The paired direction index summarizes independently estimated NES signs ",
    "and magnitudes; it is not a gene-level concordance proportion and does ",
    "not establish a shared pathway mechanism."
  )
  safe_write_csv(table, path)
  safe_write_csv(
    table,
    file.path(paths$source, "Figure3_pathway_direction.csv")
  )
  table
}

v30_build_figure3 <- function(paths, paired) {
  paired$display_class <- factor(
    ifelse(
      paired$both_significant,
      paste0("Both FDR <0.05: ", paired$direction_class),
      ifelse(
        paired$OA_significant | paired$OC_significant,
        "One disease FDR <0.05",
        "Neither disease FDR <0.05"
      )
    ),
    levels = c(
      "Both FDR <0.05: concordant",
      "Both FDR <0.05: discordant",
      "One disease FDR <0.05",
      "Neither disease FDR <0.05"
    )
  )
  colors <- c(
    "Both FDR <0.05: concordant" = "#009E73",
    "Both FDR <0.05: discordant" = "#D55E00",
    "One disease FDR <0.05" = "#6B7280",
    "Neither disease FDR <0.05" = "#D1D5DB"
  )
  shapes <- c(
    "Both FDR <0.05: concordant" = 16,
    "Both FDR <0.05: discordant" = 17,
    "One disease FDR <0.05" = 1,
    "Neither disease FDR <0.05" = 3
  )
  label_rows <- paired[paired$both_significant, , drop = FALSE]
  p1 <- ggplot2::ggplot(
    paired,
    ggplot2::aes(
      x = OA_NES, y = OC_NES,
      colour = display_class, shape = display_class
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.35) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      colour = "#D1D5DB", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_rows,
      ggplot2::aes(label = pathway),
      size = 2.0,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.23,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(
      title = "Hallmark direction across diseases",
      subtitle = "50 sets; 10 significant in both diseases",
      x = "OA normalized enrichment score",
      y = "OC normalized enrichment score",
      colour = NULL, shape = NULL
    ) +
    submission_theme(7.1) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 5.5)
    )

  selected <- paired[paired$both_significant, , drop = FALSE]
  selected <- selected[
    order(-abs(selected$paired_direction_index)),
    ,
    drop = FALSE
  ]
  selected$pathway_label <- factor(
    selected$pathway,
    levels = rev(selected$pathway)
  )
  long <- rbind(
    data.frame(
      pathway_label = selected$pathway_label,
      disease = "OA", NES = selected$OA_NES
    ),
    data.frame(
      pathway_label = selected$pathway_label,
      disease = "OC", NES = selected$OC_NES
    )
  )
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#9CA3AF",
      linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = selected,
      ggplot2::aes(
        x = OA_NES, xend = OC_NES,
        y = pathway_label, yend = pathway_label
      ),
      colour = "#CBD5E1", linewidth = 0.65
    ) +
    ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = NES, y = pathway_label, colour = disease, shape = disease),
      size = 2.35
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17)) +
    ggplot2::labs(
      title = "Jointly significant pathways",
      subtitle = "Independent disease-specific NES",
      x = "Normalized enrichment score",
      y = NULL, colour = NULL, shape = NULL
    ) +
    submission_theme(7.0) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 5.8)
    )

  p3 <- ggplot2::ggplot(
    selected,
    ggplot2::aes(
      x = paired_direction_index,
      y = pathway_label,
      fill = direction_class
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0, colour = "#6B7280",
      linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_col(width = 0.68, na.rm = TRUE) +
    ggplot2::scale_fill_manual(
      values = c(concordant = "#009E73", discordant = "#D55E00")
    ) +
    ggplot2::labs(
      title = "Paired direction index",
      subtitle = "OA NES * OC NES; signed descriptive summary",
      x = "Paired direction index",
      y = NULL, fill = NULL
    ) +
    submission_theme(7.0) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(margin = ggplot2::margin(l = 10))
    )
  safe_write_csv(
    selected,
    file.path(paths$source, "Figure3_jointly_significant_pathways.csv")
  )
  figure <- submission_panel_tag(p1 / (p2 | p3))
  submission_save_plot(
    figure,
    "Figure3_pathway_direction",
    paths$figures,
    height_mm = 180
  )

  ordered <- paired[
    order(paired$paired_direction_index, paired$pathway),
    ,
    drop = FALSE
  ]
  ordered$pathway_label <- factor(
    ordered$pathway,
    levels = ordered$pathway
  )
  full_long <- rbind(
    data.frame(
      pathway_label = ordered$pathway_label,
      disease = "OA", NES = ordered$OA_NES
    ),
    data.frame(
      pathway_label = ordered$pathway_label,
      disease = "OC", NES = ordered$OC_NES
    )
  )
  s1 <- ggplot2::ggplot(
    full_long,
    ggplot2::aes(x = disease, y = pathway_label, fill = NES)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = submission_palette[["negative"]],
      mid = "white",
      high = submission_palette[["positive"]],
      midpoint = 0
    ) +
    ggplot2::labs(
      title = "Complete Hallmark NES matrix",
      x = NULL, y = NULL, fill = "NES"
    ) +
    submission_theme(6.4) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 4.7),
      legend.position = "top"
    )
  ordered$direction_colour <- ifelse(
    ordered$direction_class == "concordant", "#009E73", "#D55E00"
  )
  s2 <- ggplot2::ggplot(
    ordered,
    ggplot2::aes(x = paired_direction_index, y = pathway_label)
  ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      colour = "#98A2B3", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = paired_direction_index, yend = pathway_label),
      colour = ordered$direction_colour,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(colour = ordered$direction_colour, size = 1.25) +
    ggplot2::labs(
      title = "Signed paired direction index",
      x = "OA NES * OC NES",
      y = NULL
    ) +
    submission_theme(6.4) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  supplementary <- submission_panel_tag(s1 | s2)
  submission_save_plot(
    supplementary,
    "SupplementaryFigure7_complete_pathway_direction",
    paths$figures,
    height_mm = 185
  )
}

v30_single_cell_detection_contrast <- function(paths) {
  path <- file.path(
    paths$tables,
    "Table_S10_single_cell_hub_gene_evidence.csv"
  )
  table <- utils::read.csv(path, check.names = FALSE)
  numeric_columns <- c("cells", "expression_value", "fraction_detected")
  table[numeric_columns] <- lapply(table[numeric_columns], as.numeric)
  table$detected_cells <- table$cells * table$fraction_detected
  table$other_cells <- NA_real_
  table$other_fraction_detected <- NA_real_
  table$within_atlas_detection_contrast <- NA_real_
  keys <- interaction(
    table$dataset_id,
    table$context,
    table$gene,
    drop = TRUE
  )
  for (key in levels(keys)) {
    index <- which(keys == key)
    total_cells <- sum(table$cells[index])
    total_detected <- sum(table$detected_cells[index])
    other_cells <- total_cells - table$cells[index]
    other_detected <- total_detected - table$detected_cells[index]
    other_fraction <- ifelse(
      other_cells > 0,
      other_detected / other_cells,
      NA_real_
    )
    table$other_cells[index] <- other_cells
    table$other_fraction_detected[index] <- other_fraction
    table$within_atlas_detection_contrast[index] <-
      table$fraction_detected[index] - other_fraction
  }
  table$specificity_status <- ifelse(
    table$within_atlas_detection_contrast >= 0.05,
    "higher detection than other annotated cells",
    ifelse(
      table$within_atlas_detection_contrast <= -0.05,
      "lower detection than other annotated cells",
      "similar detection to other annotated cells"
    )
  )
  table$metric_definition <- paste0(
    "Detection fraction in the annotated cell type minus the count-weighted ",
    "detection fraction across all other annotated cell types in the same ",
    "dataset and gene."
  )
  table$inference_boundary <- paste0(
    "Within-atlas descriptive detection contrast; values are not compared ",
    "across datasets or diseases and do not establish cell-type-specific ",
    "function, activity, or causality."
  )
  safe_write_csv(table, path)
  safe_write_csv(
    table,
    file.path(paths$source, "Figure4_within_atlas_detection_contrasts.csv")
  )
  table
}

v30_build_figure4 <- function(paths, table) {
  genes <- c(
    "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1",
    "CFI", "DDIT3", "DIRAS3", "EFEMP1", "HK2"
  )
  plot_panel <- function(dataset_id, disease, title) {
    part <- table[
      table$dataset_id == dataset_id & table$gene %in% genes,
      ,
      drop = FALSE
    ]
    part$gene <- factor(part$gene, levels = rev(genes))
    part$cell_type <- factor(
      part$cell_type,
      levels = unique(part$cell_type[order(part$cell_type)])
    )
    panel <- ggplot2::ggplot(
      part,
      ggplot2::aes(
        x = cell_type,
        y = gene,
        size = fraction_detected,
        fill = within_atlas_detection_contrast
      )
    ) +
      ggplot2::geom_point(
        shape = 21,
        colour = "#475467",
        stroke = 0.25,
        alpha = 0.95
      ) +
      ggplot2::scale_size_continuous(
        range = c(0.5, 5.1),
        limits = c(0, 1),
        labels = scales::percent
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-1, 1)
      ) +
      ggplot2::labs(
        title = title,
        subtitle = paste0(
          dataset_id,
          "; color = detection fraction minus other cells in the same atlas"
        ),
        x = NULL,
        y = NULL,
        size = "Detected",
        fill = "Within-atlas\ncontrast"
      ) +
      submission_theme(7.0) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45, hjust = 1, size = 5.9
        ),
        axis.text.y = ggplot2::element_text(size = 6.4),
        legend.position = "bottom"
      )
    if (length(unique(part$context)) > 1L) {
      panel <- panel +
        ggplot2::facet_grid(
          cols = ggplot2::vars(context),
          scales = "free_x",
          space = "free_x"
        )
    }
    panel
  }
  p1 <- plot_panel(
    "GSE255460",
    "OA",
    "OA cartilage candidate-detection context"
  )
  p2 <- plot_panel(
    "GSE154600",
    "OC",
    "OC tumor candidate-detection context"
  )

  pseudobulk <- utils::read.csv(file.path(
    paths$source,
    "Figure5_hub_pseudobulk_evidence.csv"
  ))
  pseudobulk$display <- paste0(
    pseudobulk$dataset_id, " | ", pseudobulk$cell_type,
    " | ", pseudobulk$contrast
  )
  pseudobulk$display <- factor(
    pseudobulk$display,
    levels = rev(unique(pseudobulk$display))
  )
  p3 <- ggplot2::ggplot(
    pseudobulk,
    ggplot2::aes(
      x = logFC,
      y = display,
      colour = gene
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "#98A2B3"
    ) +
    ggplot2::geom_point(size = 2.0, alpha = 0.9) +
    ggplot2::labs(
      title = "Eligible dataset-specific pseudobulk effects",
      subtitle = "Only contrasts with biological replication and FDR <0.05 are shown",
      x = "log2 fold change",
      y = NULL,
      colour = "Gene"
    ) +
    submission_theme(6.7) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.y = ggplot2::element_text(size = 5.3)
    )
  safe_write_csv(
    pseudobulk,
    file.path(paths$source, "Figure4_eligible_pseudobulk_effects.csv")
  )
  figure <- submission_panel_tag((p1 | p2) / p3)
  submission_save_plot(
    figure,
    "Figure4_single_cell_context",
    paths$figures,
    height_mm = 185
  )
}

v30_effect_size_panel <- function(effects) {
  effects$display <- factor(
    effects$dataset_id,
    levels = rev(unique(effects$dataset_id))
  )
  ggplot2::ggplot(
    effects,
    ggplot2::aes(
      x = Hedges_g,
      y = display,
      colour = disease,
      shape = estimate_type
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "#98A2B3"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      width = 0.18,
      linewidth = 0.55,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 2.35) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(
      cohort = 16,
      `random-effects summary` = 18
    )) +
    ggplot2::labs(
      title = "Direction-fixed signed-score contrasts",
      subtitle = "Hedges g; summaries across two cohorts are descriptive",
      x = "Hedges g (disease minus reference)",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.0) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 5.8)
    )
}

v30_build_figure5 <- function(paths) {
  pca <- utils::read.csv(file.path(
    paths$source,
    "Figure4_GSE54388_unsupervised_PCA.csv"
  ))
  permutations <- utils::read.csv(file.path(
    paths$source,
    "Figure4_permutation_AUC.csv"
  ))
  curves <- utils::read.csv(file.path(
    paths$source,
    "Figure4_direction_fixed_ROC_curves.csv"
  ))
  effects <- utils::read.csv(file.path(
    paths$tables,
    "Table_S22a_external_signed_score_effect_sizes.csv"
  ))
  variance <- utils::read.csv(file.path(
    paths$analysis,
    "GSE54388_unsupervised_PCA_variance.csv"
  ))
  pc1 <- 100 * variance$variance_fraction[[1L]]
  pc2 <- 100 * variance$variance_fraction[[2L]]
  p1 <- ggplot2::ggplot(
    pca,
    ggplot2::aes(x = PC1, y = PC2, colour = group, shape = group)
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = group),
      type = "norm",
      linewidth = 0.55,
      linetype = 2,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 2.35, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = c(
      Normal = "#7A8793",
      Disease = submission_palette[["OC"]]
    )) +
    ggplot2::scale_shape_manual(values = c(Normal = 17, Disease = 16)) +
    ggplot2::labs(
      title = "GSE54388 unsupervised PCA",
      subtitle = "Top 2,000 variable genes; labels used only for display",
      x = sprintf("PC1 (%.1f%%)", pc1),
      y = sprintf("PC2 (%.1f%%)", pc2),
      colour = NULL, shape = NULL
    ) +
    submission_theme(7.1) +
    ggplot2::theme(legend.position = "top")
  p2 <- v30_effect_size_panel(effects)

  levels_dataset <- c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
  permutations$dataset_id <- factor(
    permutations$dataset_id,
    levels = levels_dataset
  )
  p3 <- ggplot2::ggplot(
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
      ggplot2::aes(xintercept = observed_auc, colour = disease),
      linewidth = 0.75
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
    submission_theme(6.9) +
    ggplot2::theme(legend.position = "none")

  curves$legend_label <- factor(
    curves$legend_label,
    levels = c(
      "GSE117999 (AUC 0.520)",
      "GSE82107 (AUC 0.629)",
      "GSE54388 (AUC 1.000)",
      "GSE12470 (AUC 0.979)"
    )
  )
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
      title = "Secondary ROC display",
      subtitle = "Retrospective molecular separability",
      x = "1 - specificity",
      y = "Sensitivity",
      colour = NULL,
      linetype = NULL
    ) +
    submission_theme(6.9) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = "none"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 5.6)
    )
  safe_write_csv(
    effects,
    file.path(paths$source, "Figure5_external_effect_sizes.csv")
  )
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure5_molecular_separability",
    paths$figures,
    height_mm = 175
  )
}

v30_build_supplementary_figure10 <- function(paths) {
  context <- utils::read.csv(file.path(
    paths$tables,
    "Table_S21_cross_cohort_molecular_separability_context.csv"
  ))
  leave_one_out <- utils::read.csv(file.path(
    paths$source,
    "Figure4_leave_one_out_AUC.csv"
  ))
  calibration <- utils::read.csv(file.path(
    paths$source,
    "SupplementaryFigure10_calibration_bins.csv"
  ))
  schematic <- data.frame(
    disease = c("OA", "OC"),
    y = c(2, 1),
    left = c(
      "Reference cartilage / synovium",
      "Nonmalignant ovarian / peritoneal reference"
    ),
    right = c(
      "OA cartilage / synovium",
      "Serous ovarian carcinoma / HGSOC"
    ),
    contrast = c(
      "Within-tissue chronic degenerative / inflammatory contrast",
      "Malignant tumor-reference contrast"
    ),
    colour = c(submission_palette[["OA"]], submission_palette[["OC"]]),
    stringsAsFactors = FALSE
  )
  p1 <- ggplot2::ggplot(schematic) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 1.65, xend = 3.35, y = y, yend = y, colour = disease),
      arrow = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed"),
      linewidth = 0.8
    ) +
    ggplot2::geom_label(
      ggplot2::aes(x = 1, y = y, label = left),
      size = 2.25,
      lineheight = 0.95,
      fill = "#F8FAFC",
      colour = "#1F2937",
      linewidth = 0.3
    ) +
    ggplot2::geom_label(
      ggplot2::aes(x = 4, y = y, label = right),
      size = 2.25,
      lineheight = 0.95,
      fill = "#FFF9F2",
      colour = "#1F2937",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2.5, y = y - 0.28, label = contrast),
      size = 2.05,
      colour = "#475467"
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::coord_cartesian(xlim = c(0.25, 4.75), ylim = c(0.55, 2.45)) +
    ggplot2::labs(
      title = "Validation contrasts differed biologically",
      subtitle = "Schematic summarizes actual tissue/comparator classes; not a disease continuum"
    ) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 9),
      plot.subtitle = ggplot2::element_text(size = 6.8, colour = "#475467")
    )

  levels_dataset <- c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
  leave_one_out$dataset_id <- factor(
    leave_one_out$dataset_id,
    levels = levels_dataset
  )
  p2 <- ggplot2::ggplot(
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
    ggplot2::geom_jitter(
      height = 0.12, width = 0,
      alpha = 0.65, size = 1.35
    ) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 18,
      size = 3.0
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
    submission_theme(7.0) +
    ggplot2::theme(legend.position = "top")

  p3 <- ggplot2::ggplot(
    calibration,
    ggplot2::aes(
      x = mean_predicted,
      y = observed_fraction,
      colour = disease,
      group = dataset_id
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      colour = "#98A2B3"
    ) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(
      ggplot2::aes(size = n_in_bin),
      alpha = 0.88
    ) +
    ggplot2::facet_wrap(~dataset_id, ncol = 2) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Cross-fitted calibration sensitivity",
      subtitle = "No locked probability model was transported from discovery",
      x = "Mean cross-fitted probability",
      y = "Observed disease fraction",
      colour = NULL,
      size = "Bin n"
    ) +
    submission_theme(6.8) +
    ggplot2::theme(legend.position = "top")
  safe_write_csv(
    schematic,
    file.path(paths$source, "SupplementaryFigure10_contrast_schematic.csv")
  )
  safe_write_csv(
    context,
    file.path(paths$source, "SupplementaryFigure10_validation_context.csv")
  )
  figure <- submission_panel_tag(p1 / (p2 | p3))
  submission_save_plot(
    figure,
    "SupplementaryFigure10_design_and_reliability",
    paths$figures,
    height_mm = 175
  )
}

v30_build_figure6 <- function(paths) {
  mr <- utils::read.csv(file.path(
    paths$tables,
    "Table_S12a_MR_estimates_and_provenance.csv"
  ))
  methods <- c(
    "Inverse variance weighted",
    "Weighted median",
    "MR Egger",
    "Weighted mode",
    "Simple mode"
  )
  mr$method <- factor(mr$method, levels = rev(methods))
  mr$direction_short <- ifelse(
    mr$direction == "OA to ovarian cancer",
    "OA liability -> OC risk",
    "OC liability -> OA risk"
  )
  p1 <- ggplot2::ggplot(
    mr,
    ggplot2::aes(
      x = odds_ratio,
      y = method,
      colour = direction_short
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#667085"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_95, xmax = ci_upper_95),
      width = 0.18,
      linewidth = 0.55,
      orientation = "y",
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::geom_point(
      size = 2.1,
      position = ggplot2::position_dodge(width = 0.45)
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(
      "OA liability -> OC risk" = submission_palette[["OA"]],
      "OC liability -> OA risk" = submission_palette[["OC"]]
    )) +
    ggplot2::labs(
      title = "Bidirectional genetic-liability estimates",
      subtitle = "21 OA instruments and 11 OC instruments",
      x = "Odds ratio (95% CI; log scale)",
      y = NULL,
      colour = NULL
    ) +
    submission_theme(7.0) +
    ggplot2::theme(legend.position = "top")

  boundary <- data.frame(
    layer = factor(
      c(
        "Transcriptomic evidence",
        "Regulatory context",
        "Bidirectional MR"
      ),
      levels = rev(c(
        "Transcriptomic evidence",
        "Regulatory context",
        "Bidirectional MR"
      ))
    ),
    evidence = c(
      "Gene, pathway, cell-context,\nand external-cohort associations",
      "KnockTF target-set direction;\nmiRTarBase interaction coverage",
      "OA liability -> OC risk;\nOC liability -> OA risk"
    ),
    interpretation = c(
      "Partial convergence with\ndirectional and cellular heterogeneity",
      "Focused hypotheses;\nactivity was not measured",
      "No evidence under available\ninstruments and assumptions"
    ),
    boundary = c(
      "Not a shared mechanism",
      "Not causal regulation",
      "Not a test of shared heritability\nor genetic correlation"
    ),
    stringsAsFactors = FALSE
  )
  p2 <- ggplot2::ggplot(boundary, ggplot2::aes(y = layer)) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 1),
      width = 0.97, height = 0.80,
      fill = "#F8FAFC", colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 2),
      width = 0.97, height = 0.80,
      fill = "#F3F8F5", colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 3),
      width = 0.97, height = 0.80,
      fill = "#FFF9F2", colour = "#D0D5DD"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 1, label = evidence),
      size = 2.2, lineheight = 0.95, colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2, label = interpretation),
      size = 2.2, lineheight = 0.95, colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 3, label = boundary),
      size = 2.15, lineheight = 0.95, colour = "#9A3412"
    ) +
    ggplot2::annotate(
      "text",
      x = c(1, 2, 3),
      y = 3.45,
      label = c("Evidence", "Bounded interpretation", "What is not established"),
      fontface = "bold",
      size = 2.8,
      colour = "#344054"
    ) +
    ggplot2::scale_x_continuous(limits = c(0.48, 3.52), breaks = NULL) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(1.0, 0.9))
    ) +
    ggplot2::labs(
      title = "Integrative inference boundary",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = "bold", size = 6.5, colour = "#344054"
      ),
      plot.title = ggplot2::element_text(
        face = "bold", size = 9, colour = "#111827"
      )
    )
  safe_write_csv(
    boundary,
    file.path(paths$source, "Figure6_inference_boundary.csv")
  )
  figure <- submission_panel_tag(p1 / p2)
  submission_save_plot(
    figure,
    "Figure6_genetic_liability_boundary",
    paths$figures,
    height_mm = 175
  )
}

v30_write_figure_legends <- function(paths) {
  legends <- c(
    "## Figure 1. Systems-level framework for shared-but-nonidentical transcriptional features",
    "",
    paste0(
      "**A,** OA and OC were analyzed on disease-specific tracks that met at ",
      "a partial shared transcriptional landscape. Six analytical questions ",
      "then separated gene direction, pathway direction, cell-state localization, ",
      "cross-cohort separability, bidirectional genetic liability, and focused ",
      "regulatory hypotheses. **B,** Audited evidence layers, intended roles, ",
      "and inference boundaries. Shared membership was not treated as proof of ",
      "shared direction, pathway state, cellular meaning, clinical utility, or causality."
    ),
    "",
    "## Figure 2. Gene-level convergence with an approximately balanced directional structure",
    "",
    paste0(
      "**A-B,** OA and OC differential-expression volcano plots. **C,** All ",
      "commonly measured genes with the 286 primary shared DEGs colored by ",
      "direction; the inset isolates the ten-gene interpretable evidence ",
      "summary. The four quadrants contain 112, 34, 86, and 54 genes, yielding ",
      "146 concordant and 140 discordant genes. **D,** Shared-gene counts under ",
      "six prespecified thresholds."
    ),
    "",
    "## Figure 3. Quantitative pathway-direction heterogeneity",
    "",
    paste0(
      "**A,** OA and OC normalized enrichment scores (NES) for all 50 Hallmark ",
      "sets. **B,** Independent disease-specific NES for the ten sets significant ",
      "in both diseases; six had opposite signs and four had matching signs. ",
      "**C,** Signed paired direction index (OA NES multiplied by OC NES). The ",
      "index is a descriptive summary of independently estimated pathway ",
      "directions, not a gene-level concordance proportion or shared mechanism."
    ),
    "",
    "## Figure 4. Within-atlas candidate detection and cell-context localization",
    "",
    paste0(
      "**A-B,** Candidate detection fractions and within-atlas detection ",
      "contrasts in the largest count-level OA (GSE255460) and OC (GSE154600) ",
      "atlases used for the main display. Contrast equals the annotated cell-type ",
      "detection fraction minus the count-weighted detection fraction in all ",
      "other annotated cells from the same dataset and gene. Values were not ",
      "compared across datasets or diseases. **C,** Significant eligible ",
      "dataset-specific pseudobulk effects. OC datasets mainly contributed ",
      "localization because suitable reference replication was unavailable."
    ),
    "",
    "## Figure 5. Cross-cohort molecular separability across different contrast scales",
    "",
    paste0(
      "**A,** Unsupervised PCA of GSE54388. **B,** Direction-fixed signed-score ",
      "Hedges g estimates; two-cohort summaries are descriptive. **C,** Null AUC ",
      "distributions from 1,000 label permutations. **D,** ROC curves, shown last ",
      "as a secondary display. These panels describe retrospective molecular ",
      "separability in biologically different validation tasks, not a universal ",
      "diagnostic model."
    ),
    "",
    "## Figure 6. Bidirectional genetic-liability assessment and integrative boundary",
    "",
    paste0(
      "**A,** Five MR estimators for genetic liability to hip/knee OA on OC risk ",
      "and genetic liability to OC on OA risk. No estimator supported a ",
      "liability-to-outcome association under the selected datasets, instruments, ",
      "and assumptions. **B,** Evidence layers and their explicit boundaries. ",
      "Bidirectional MR is not a test of shared heritability, genetic correlation, ",
      "or common susceptibility loci; detailed diagnostics are in Figure S2 and ",
      "focused regulatory matrices are in Figure S11."
    ),
    "",
    "## Supplementary Figure 1. Prespecified non-MR sensitivity analyses",
    "",
    paste0(
      "DEG threshold retention, WGCNA leave-one-out stability, strict nested ",
      "feature frequency, and TCGA LASSO selection stability."
    ),
    "",
    "## Supplementary Figure 2. Detailed bidirectional MR estimates and diagnostics",
    "",
    paste0(
      "Five MR estimators and heterogeneity/pleiotropy diagnostics in both ",
      "directions. Null estimates mean no liability-to-outcome evidence was ",
      "detected under the available instruments and assumptions, not proof of ",
      "absence or a test of shared genetic architecture."
    ),
    "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings",
    "",
    paste0(
      "Released or recomputed embeddings with dataset-specific labels. ",
      "Quantitative summaries used all QC-pass cells; OA and OC were not ",
      "integrated into a shared latent space."
    ),
    "",
    "## Supplementary Figure 4. HPA normal-tissue context",
    "",
    paste0(
      "Normal-tissue specificity and distribution categories. Cartilage is ",
      "absent, so this audit cannot establish OA-OC specificity."
    ),
    "",
    "## Supplementary Figure 5. TCGA-OV relative stromal and immune context audit",
    "",
    paste0(
      "**A,** Spearman correlations between candidate expression and published ",
      "ESTIMATE rank-based stromal, immune, and combined scores. **B-C,** ",
      "Representative strongest absolute stromal and immune associations. Scores ",
      "are relative within-cohort transcriptomic proxies; absolute tumor purity ",
      "and histologic cell fractions were not inferred from RNA-seq."
    ),
    "",
    "## Supplementary Figure 6. Exploratory cell-type functional annotation",
    "",
    paste0(
      "**A,** Gene Ontology Biological Process over-representation among top ",
      "cluster markers from the GSE104782 OA cartilage atlas. **B,** Corresponding ",
      "annotation in GSE154600 ovarian tumors. Labels remain dataset specific, ",
      "so these panels localize functional themes without demonstrating conserved ",
      "function or mechanism."
    ),
    "",
    "## Supplementary Figure 7. Complete Hallmark direction and paired direction index",
    "",
    paste0(
      "**A,** Complete OA and OC Hallmark NES matrix. **B,** Signed paired ",
      "direction index for all 50 sets. Positive and negative values indicate ",
      "matching and opposite NES signs, respectively; magnitude is descriptive."
    ),
    "",
    "## Supplementary Figure 8. WGCNA and strict nested candidate-evidence stability",
    "",
    paste0(
      "**A,** WGCNA module-trait association under soft-power perturbation. ",
      "**B,** Primary-module gene retention. **C,** Strict nested feature-selection ",
      "frequency summary. **D,** Disease-specific candidate discovery effects. ",
      "These analyses support evidence stability, not a diagnostic panel."
    ),
    "",
    "## Supplementary Figure 9. Candidate-centered Hallmark contexts",
    "",
    paste0(
      "Disease-status-adjusted residual association rankings for SOX9, DDIT3, ",
      "BNC1, and AKAP12 against Hallmark sets. This transparent post hoc analysis ",
      "is exploratory and does not represent single-gene perturbation, mediation, ",
      "pathway activation, or mechanism."
    ),
    "",
    "## Supplementary Figure 10. External-cohort design, sample influence, and calibration sensitivity",
    "",
    paste0(
      "**A,** Factual schematic of the different OA and OC validation contrast ",
      "classes; it is not a disease continuum. **B,** Leave-one-sample-out AUCs. ",
      "**C,** Three-bin cross-fitted calibration sensitivity. Because no locked ",
      "probability model or clinical decision threshold was transported, a ",
      "nomogram and decision-curve analysis were not performed."
    ),
    "",
    "## Supplementary Figure 11. Focused upstream regulatory context",
    "",
    paste0(
      "**A,** OA and OC enrichment of KnockTF perturbational target sets for ",
      "transcription factors connected to four candidate exemplars; NES does not ",
      "estimate transcription-factor activity. **B,** Curated miRTarBase ",
      "interaction coverage for multi-candidate miRNAs. miRNA abundance and ",
      "activity were not measured, so disease-specific miRNA regulation was not inferred."
    ),
    ""
  )
  writeLines(
    legends,
    file.path(paths$figures, "figure_legends.md"),
    useBytes = TRUE
  )
}

v30_update_table_index <- function(paths) {
  path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(path, check.names = FALSE)
  row10 <- index$table_id == "Table S10"
  index$title[row10] <- "Dataset-specific candidate detection and within-atlas contrast"
  index$contents[row10] <- paste0(
    "Detection fractions, count-weighted other-cell detection fractions, ",
    "within-atlas detection contrasts, and eligible pseudobulk results."
  )
  row18 <- index$table_id == "Table S18"
  index$title[row18] <- "Hallmark pathway direction and paired direction index"
  index$contents[row18] <- paste0(
    "Complete paired OA/OC Hallmark NES, FDR, direction class, and signed ",
    "descriptive paired direction index with explicit definition and boundary."
  )
  safe_write_csv(index, path)
}

v30_update_claim_registry <- function(paths) {
  source <- file.path(paths$root, "claim_evidence_registry_v24.csv")
  registry <- utils::read.csv(source, check.names = FALSE)
  update <- function(id, claim = NULL, location = NULL, allowed = NULL, prohibited = NULL) {
    row <- registry$claim_id == id
    if (!is.null(claim)) registry$manuscript_claim[row] <<- claim
    if (!is.null(location)) registry$figure_or_table[row] <<- location
    if (!is.null(allowed)) registry$allowed_wording[row] <<- allowed
    if (!is.null(prohibited)) registry$prohibited_wording[row] <<- prohibited
  }
  update(
    "C02",
    paste0(
      "Of 286 shared DEGs, 146 (51.0%) were directionally concordant and ",
      "140 (49.0%) were discordant, indicating an approximately balanced ",
      "directional structure."
    ),
    "Figure 2; Tables S2-S3",
    "approximately balanced directional structure; partial overlap",
    "majority discordant; uniform shared program"
  )
  update("C03", location = "Figure 2; Figure S8; Tables S2 and S16")
  update("C04", location = "Figures S1 and S8; Table S4")
  update("C05", location = "Figures S1 and S8; Table S5")
  update("C06", location = "Figure 5; Figure S10; Tables S6 and S21-S22")
  update("C07", location = "Figure S5; Figure S1; Table S7")
  update(
    "C08",
    paste0(
      "Within-atlas detection contrasts and dataset-specific localization ",
      "place candidates in distinct annotated OA and OC cellular contexts; ",
      "eligible OA pseudobulk contrasts support cell-state associations."
    ),
    "Figure 4; Figures S3 and S6; Tables S9-S10 and S19",
    "within-atlas detection contrast; dataset-specific cellular context",
    "cross-dataset cell specificity; cell type causes disease"
  )
  update("C09", location = "Figures S5-S7; Tables S8 and S11")
  update(
    "C10",
    paste0(
      "Bidirectional MR found no evidence that genetic liability to OA affected ",
      "OC risk, or that genetic liability to OC affected OA risk, under the ",
      "selected GWAS datasets, instruments, and assumptions."
    ),
    "Figure 6; Figure S2; Table S12",
    "no detected liability-to-outcome evidence under available instruments",
    "no shared genetic architecture; MR proves absence"
  )
  update("C11", location = "Figure 6; Figure S2; Table S12")
  update(
    "C12",
    paste0(
      "The total evidence supports partial context-dependent transcriptomic ",
      "convergence, not a uniform shared program or a supported disease-to-disease ",
      "germline-liability relationship."
    ),
    "Figures 1-6",
    "partial context-dependent convergence with explicit inference boundaries",
    "shared mechanism; shared genetic causality"
  )
  update("C15", location = "Figure 5; Table S6")
  update(
    "C17",
    registry$manuscript_claim[registry$claim_id == "C17"],
    "Figure 6; Figure S2; Tables S12a-b",
    "verified MR provenance and no detected liability-to-outcome evidence",
    "shared heritability test; proof that causality is absent"
  )
  update("C19", location = "Figure 4; Figure S6; Table S17")
  update("C21", location = "Figure 3; Figure S7; Table S18")
  update("C22", location = "Figure 4; Table S19")
  update("C23", location = "Figure 5; Figure S10; Tables S21-S22")
  safe_write_csv(
    registry,
    file.path(paths$root, "claim_evidence_registry_v30.csv")
  )
  unlink(source)

  checklist_source <- file.path(
    paths$root,
    "reproducibility_checklist_v24.csv"
  )
  checklist <- utils::read.csv(checklist_source, check.names = FALSE)
  checklist <- rbind(
    checklist,
    data.frame(
      item_id = c("R30", "R31", "R32"),
      domain = c("narrative", "single-cell", "model boundary"),
      item = c(
        "Six-main-figure systems-level narrative is rebuilt without adding models.",
        "Within-atlas detection contrast is computed without cross-dataset integration.",
        "DCA and nomogram are excluded because no locked probability model or clinical decision threshold exists."
      ),
      status = "complete",
      evidence = c(
        "Figures 1-6 and figure_legends.md",
        "Table S10 and Figure 4 source data",
        "Methods, Figure S10 legend, and response matrix"
      ),
      stringsAsFactors = FALSE
    )
  )
  safe_write_csv(
    checklist,
    file.path(paths$root, "reproducibility_checklist_v30.csv")
  )
  unlink(checklist_source)
}

run_reviewer_v30 <- function(project_root) {
  paths <- v30_output_paths(project_root)
  v30_prepare_baseline(project_root, paths)
  log_info("Building V3.0 systems-level framework.")
  v30_build_figure1(paths)

  log_info("Retaining audited gene-direction quadrants without adding a duplicate analysis.")
  paired <- v30_update_pathway_table(paths)
  log_info("Promoting quantitative pathway direction to main Figure 3.")
  v30_build_figure3(paths, paired)

  log_info("Computing bounded within-atlas single-cell detection contrasts.")
  single_cell <- v30_single_cell_detection_contrast(paths)
  v30_build_figure4(paths, single_cell)

  log_info("Rebuilding external molecular separability with effect sizes before ROC.")
  v30_build_figure5(paths)
  v30_build_supplementary_figure10(paths)

  log_info("Promoting bidirectional genetic-liability assessment with explicit boundaries.")
  v30_build_figure6(paths)
  v30_write_figure_legends(paths)
  v30_update_table_index(paths)
  v30_update_claim_registry(paths)

  log_info(
    "V3.0 positioning upgrade completed: six main figures, ",
    "eleven supplementary figures, and no added prediction algorithm."
  )
  invisible(paths)
}
