v32_base_runner <- run_reviewer_v31

v32_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v32"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs"))
  )
}

v32_theme <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, colour = "#111111"),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.6, colour = "#4B5563"),
      axis.title = ggplot2::element_text(face = "plain", colour = "#222222"),
      axis.text = ggplot2::element_text(colour = "#333333"),
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "#222222"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "#111111"),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

v32_wrap <- function(x, width = 34L) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1L))
}

v32_copy_file <- function(source, target) {
  if (!file.exists(source)) stop("Missing V3.1 source: ", source, call. = FALSE)
  if (!file.copy(source, target, overwrite = TRUE)) {
    stop("Could not copy ", source, " to ", target, call. = FALSE)
  }
  invisible(target)
}

v32_prepare_scope <- function(project_root, paths) {
  source_root <- file.path(project_root, "results", "submission_v31")
  if (!dir.exists(source_root)) stop("V3.1 outputs are required for V3.2.", call. = FALSE)
  retained_tables <- c(
    "Table_S1_data_sources_and_cohorts.csv",
    "Table_S2_shared_differentially_expressed_genes.csv",
    "Table_S3a_DEG_threshold_summary.csv",
    "Table_S3b_DEG_threshold_membership.csv",
    "Table_S4_WGCNA_stability.csv",
    "Table_S5_machine_learning_resampling.csv",
    "Table_S6_external_validation.csv",
    "Table_S9_single_cell_QC_and_status.csv",
    "Table_S10_single_cell_hub_gene_evidence.csv",
    "Table_S11a_Hallmark_GSEA.csv",
    "Table_S11b_GO_shared_genes.csv",
    "Table_S11c_KEGG_shared_genes.csv",
    "Table_S16_candidate_prioritization_matrix.csv",
    "Table_S18_Hallmark_pathway_direction_matrix.csv",
    "Table_S19_gene_cell_function_context_matrix.csv",
    "Table_S21_cross_cohort_molecular_separability_context.csv",
    "Table_S22a_external_signed_score_effect_sizes.csv",
    "Table_S24a_dataset_context_CCSS.csv",
    "Table_S24b_disease_consensus_CCSS.csv",
    "Table_S24c_sample_aware_UCell_summary.csv",
    "Table_S25a_STRING_mapping_audit.csv",
    "Table_S25b_direction_aware_STRING_edges.csv",
    "Table_S25c_STRING_node_topology.csv",
    "Table_S25d_STRING_network_topology.csv",
    "Table_S25e_STRING_direction_label_permutation.csv",
    "Table_S29a_panel_size_composition.csv",
    "Table_S29b_panel_size_direction_sensitivity.csv",
    "Table_S30a_extended_gene_cell_detection.csv",
    "Table_S30b_panel_size_localization_sensitivity.csv",
    "Table_S30c_panel_size_localization_profile_comparison.csv",
    "Table_S31a_bulk_unsupervised_PCA.csv",
    "Table_S31b_bulk_sample_correlation_QC.csv"
  )
  for (name in retained_tables) {
    v32_copy_file(
      file.path(source_root, "supplementary_tables", name),
      file.path(paths$tables, name)
    )
  }
  scope <- data.frame(
    module = c(
      "Bulk DEG overlap", "Functional enrichment", "Direction analysis",
      "WGCNA", "STRING", "Candidate ranking", "Single-cell localization",
      "Exploratory classification", "Mendelian randomization",
      "CellChat/NicheNet", "TF-miRNA", "TCGA/immune/HPA context",
      "DCA/nomogram"
    ),
    v32_status = c(
      "main", "main", "main", "main, stability in supplement",
      "auxiliary main panel", "main", "main",
      "supplement only", "excluded from V3.2 scope", "excluded from V3.2 scope",
      "excluded from V3.2 scope", "excluded from V3.2 scope", "not performed"
    ),
    rationale = c(
      "directly answers whether transcriptomic alterations overlap",
      "summarizes significant biological themes without implying one conserved mechanism",
      "central evidence that shared membership does not imply shared state",
      "supports disease-associated co-expression context without dominating the narrative",
      "database association context only; topology did not select candidates",
      "transparent descriptive candidate set; not a diagnostic signature",
      "answers which exact source-defined cell labels carry candidate expression",
      "secondary cohort-dependent separation check; no clinical model claim",
      "does not directly answer the revised descriptive transcriptomic question",
      "no direct fixed-candidate ligand/receptor axis and no supported ligand-activity inference",
      "database-derived regulatory context did not yield a validated disease-specific axis",
      "not required for the discovery-direction-localization narrative",
      "no locked clinical probability model or decision threshold"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(scope, file.path(paths$analysis, "V32_scope_decisions.csv"))
  invisible(source_root)
}

v32_copy_supplementary_figures <- function(project_root, paths) {
  source_dir <- file.path(project_root, "results", "submission_v31", "figures")
  mapping <- c(
    "SupplementaryFigure15_panel_size_sensitivity" = "SupplementaryFigure2_candidate_set_sensitivity",
    "Figure5_molecular_separability" = "SupplementaryFigure3_exploratory_classification",
    "SupplementaryFigure7_complete_pathway_direction" = "SupplementaryFigure4_complete_Hallmark_direction",
    "SupplementaryFigure3_single_cell_UMAPs" = "SupplementaryFigure5_all_single_cell_UMAPs",
    "SupplementaryFigure16_bulk_PCA_and_QC" = "SupplementaryFigure6_bulk_PCA_QC",
    "SupplementaryFigure8_candidate_evidence_stability" = "SupplementaryFigure7_candidate_stability"
  )
  for (extension in c("png", "pdf")) {
    for (source_stem in names(mapping)) {
      v32_copy_file(
        file.path(source_dir, paste0(source_stem, ".", extension)),
        file.path(paths$figures, paste0(mapping[[source_stem]], ".", extension))
      )
    }
  }
  source_data_dir <- file.path(source_dir, "source_data")
  selected_source <- c(
    "SupplementaryFigure15_direction_composition.csv",
    "SupplementaryFigure15_Hallmark_profile.csv",
    "SupplementaryFigure15_cell_localization.csv",
    "Figure4_GSE54388_unsupervised_PCA.csv",
    "Figure4_permutation_AUC.csv",
    "Figure4_direction_fixed_ROC_curves.csv",
    "Figure5_external_effect_sizes.csv",
    "SupplementaryFigure7_pathway_direction.csv",
    "SupplementaryFigure3_UMAP_subsamples.csv",
    "SupplementaryFigure16_bulk_PCA.csv",
    "SupplementaryFigure16_bulk_sample_correlation.csv",
    "Figure3_WGCNA_bootstrap.csv",
    "Figure3_WGCNA_powers.csv",
    "Figure3_ML_selection_frequency.csv",
    "Figure3_feature_stability_summary.csv"
  )
  for (name in selected_source) {
    v32_copy_file(file.path(source_data_dir, name), file.path(paths$source, name))
  }
}

v32_volcano <- function(data, title) {
  data$status <- "Not primary"
  primary <- data$adj.P.Val < 0.05 & abs(data$logFC) >= 1
  data$status[primary & data$logFC > 0] <- "Higher"
  data$status[primary & data$logFC < 0] <- "Lower"
  data$status <- factor(data$status, levels = c("Not primary", "Lower", "Higher"))
  ggplot2::ggplot(data, ggplot2::aes(logFC, -log10(pmax(adj.P.Val, .Machine$double.xmin)))) +
    ggplot2::geom_point(ggplot2::aes(colour = status), size = 0.62, alpha = 0.62) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.3, colour = "#777777") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "#777777") +
    ggplot2::scale_colour_manual(values = c("Not primary" = "#C9CDD2", "Lower" = "#2C7FB8", "Higher" = "#D95F0E")) +
    ggplot2::labs(title = title, x = "log2 fold change", y = "-log10 FDR", colour = NULL) +
    v32_theme(7.1) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 5.8))
}

v32_build_figure1 <- function(paths) {
  nodes <- data.frame(
    x = c(1, 2, 3, 1, 2, 3), y = c(2.3, 2.3, 2.3, 1.0, 1.0, 1.0),
    stage = c("GEO discovery", "Disease-specific DEG", "286 shared genes", "Functional and\ndirection analysis", "Candidate-gene\ncharacterization", "Single-cell\nlocalization"),
    detail = c("OA and OC analyzed separately", "FDR <0.05; |log2FC| >=1", "membership overlap", "GO/KEGG/Hallmark", "WGCNA + secondary evidence", "exact source labels"),
    fill = c("#E8F1F8", "#E8F1F8", "#F2ECF7", "#FBE9E3", "#F3F3F3", "#E9F4EF"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(1.42, 2.42, 3.00, 2.58, 1.58), xend = c(1.58, 2.58, 3.00, 2.42, 1.42),
    y = c(2.3, 2.3, 1.92, 1.0, 1.0), yend = c(2.3, 2.3, 1.38, 1.0, 1.0)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
                          arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"), linewidth = 0.55, colour = "#333333") +
    ggplot2::geom_tile(data = nodes, ggplot2::aes(x, y, fill = fill), width = 0.82, height = 0.78, colour = "#333333", linewidth = 0.45) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y + 0.13, label = stage), fontface = "bold", size = 3.0, lineheight = 0.92) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y - 0.18, label = detail), size = 2.2, colour = "#4B5563", lineheight = 0.9) +
    ggplot2::annotate("text", x = 2, y = 0.35,
      label = "Primary question: how can partially shared transcriptomic alterations occupy different directions and cellular contexts?",
      size = 2.75, fontface = "italic", colour = "#374151") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.15, 2.8), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
  safe_write_csv(nodes, file.path(paths$source, "Figure1_workflow_nodes.csv"))
  submission_save_plot(p, "Figure1_study_workflow", paths$figures, 185, 112)
}

v32_build_figure2 <- function(project_root, paths) {
  source <- file.path(project_root, "results", "submission_v31", "figures", "source_data")
  oa <- utils::read.csv(file.path(source, "Figure2_OA_DEG.csv"))
  oc <- utils::read.csv(file.path(source, "Figure2_OC_DEG.csv"))
  shared <- utils::read.csv(file.path(source, "Figure2_common_gene_effects_quadrants.csv"))
  shared <- shared[as.logical(shared$primary_shared), , drop = FALSE]
  p1 <- v32_volcano(oa, "OA discovery")
  p2 <- v32_volcano(oc, "OC discovery")

  circle <- function(cx, cy, radius, set_name) {
    angle <- seq(0, 2 * pi, length.out = 240L)
    data.frame(x = cx + radius * cos(angle), y = cy + radius * sin(angle), set = set_name)
  }
  circles <- rbind(circle(1.15, 1, 0.86, "OA"), circle(1.95, 1, 0.86, "OC"))
  venn <- ggplot2::ggplot(circles, ggplot2::aes(x, y, group = set, fill = set)) +
    ggplot2::geom_polygon(alpha = 0.34, colour = "#333333", linewidth = 0.45) +
    ggplot2::annotate("text", x = 0.62, y = 1.0, label = "1,722", fontface = "bold", size = 4.0) +
    ggplot2::annotate("text", x = 1.55, y = 1.0, label = "286", fontface = "bold", size = 4.3) +
    ggplot2::annotate("text", x = 2.47, y = 1.0, label = "2,024", fontface = "bold", size = 4.0) +
    ggplot2::annotate("text", x = 0.70, y = 1.92, label = "OA DEGs\n2,008", size = 3.0, fontface = "bold") +
    ggplot2::annotate("text", x = 2.42, y = 1.92, label = "OC DEGs\n2,310", size = 3.0, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E"), guide = "none") +
    ggplot2::coord_equal(xlim = c(0.1, 3.0), ylim = c(0.05, 2.15), clip = "off") +
    ggplot2::labs(title = "Primary DEG overlap", subtitle = "286 shared genes") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 8), plot.subtitle = ggplot2::element_text(size = 6.6, colour = "#4B5563"))

  quadrant_levels <- c("OA higher / OC higher", "OA lower / OC lower", "OA higher / OC lower", "OA lower / OC higher")
  shared$direction_quadrant <- factor(shared$direction_quadrant, levels = quadrant_levels)
  shared <- shared[order(shared$direction_quadrant, -abs(shared$logFC_OA) - abs(shared$logFC_OC)), ]
  shared$row_index <- seq_len(nrow(shared))
  heat <- rbind(
    data.frame(gene = shared$gene, row_index = shared$row_index, disease = "OA", log2FC = shared$logFC_OA, quadrant = shared$direction_quadrant),
    data.frame(gene = shared$gene, row_index = shared$row_index, disease = "OC", log2FC = shared$logFC_OC, quadrant = shared$direction_quadrant)
  )
  separators <- cumsum(as.integer(table(shared$direction_quadrant)))
  heatmap <- ggplot2::ggplot(heat, ggplot2::aes(disease, row_index, fill = log2FC)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = separators[-length(separators)] + 0.5, linewidth = 0.45, colour = "#333333") +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0,
                                  limits = c(-5.5, 5.5), oob = scales::squish) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::labs(title = "Direction-ordered shared genes", subtitle = "286 genes; four quadrants separated", x = NULL, y = "Shared genes", fill = "log2FC") +
    v32_theme(7.1) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), legend.position = "right")

  safe_write_csv(shared, file.path(paths$source, "Figure2_shared_gene_order.csv"))
  safe_write_csv(heat, file.path(paths$source, "Figure2_shared_gene_heatmap.csv"))
  safe_write_csv(data.frame(set = c("OA only", "Shared", "OC only"), genes = c(1722L, 286L, 2024L)), file.path(paths$source, "Figure2_overlap_counts.csv"))
  figure <- submission_panel_tag((p1 | p2) / (venn | heatmap) + patchwork::plot_layout(heights = c(1, 0.95)))
  submission_save_plot(figure, "Figure2_shared_transcriptomic_alterations", paths$figures, 185, 170)
}

v32_build_figure3 <- function(paths) {
  go <- utils::read.csv(file.path(paths$tables, "Table_S11b_GO_shared_genes.csv"))
  go_ids <- c("GO:0098813", "GO:0044839", "GO:0045786", "GO:0031012", "GO:0005201", "GO:0002697", "GO:0002367", "GO:0140467", "GO:0001666", "GO:0070482")
  go_theme <- c(
    "Cell cycle", "Cell cycle", "Cell cycle", "Matrix", "Matrix",
    "Immune", "Immune", "Stress", "Stress", "Stress"
  )
  selected_go <- go[match(go_ids, go$ID), , drop = FALSE]
  if (anyNA(selected_go$ID)) stop("A representative V3.2 GO term is missing.", call. = FALSE)
  selected_go$theme <- go_theme
  selected_go$term <- factor(v32_wrap(selected_go$Description, 30L), levels = rev(v32_wrap(selected_go$Description, 30L)))
  p1 <- ggplot2::ggplot(selected_go, ggplot2::aes(FoldEnrichment, term, size = Count, fill = theme)) +
    ggplot2::geom_point(shape = 21, colour = "#333333", stroke = 0.35) +
    ggplot2::scale_fill_manual(values = c("Cell cycle" = "#7B6EA8", Matrix = "#2C7FB8", Immune = "#2A9D8F", Stress = "#D95F0E")) +
    ggplot2::labs(title = "Representative non-redundant GO terms", subtitle = "All shown terms FDR <0.05; full results in Table S11b", x = "Fold enrichment", y = NULL, size = "Genes", fill = NULL) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, nrow = 1),
      fill = ggplot2::guide_legend(order = 2, nrow = 1)
    ) +
    v32_theme(6.7) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = ggplot2::element_text(size = 5.2),
      legend.spacing.y = grid::unit(0.3, "mm"),
      axis.text.y = ggplot2::element_text(size = 5.8)
    )

  kegg <- utils::read.csv(file.path(paths$tables, "Table_S11c_KEGG_shared_genes.csv"))
  kegg$term <- v32_wrap(kegg$Description, 30L)
  p2 <- ggplot2::ggplot(kegg, ggplot2::aes(FoldEnrichment, term, size = Count, fill = -log10(p.adjust))) +
    ggplot2::geom_point(shape = 21, colour = "#333333", stroke = 0.4) +
    ggplot2::scale_fill_gradient(low = "#DDEBF7", high = "#2C7FB8") +
    ggplot2::scale_x_continuous(limits = c(0, max(kegg$FoldEnrichment) * 1.35)) +
    ggplot2::annotate("text", x = max(kegg$FoldEnrichment) * 0.70, y = 1.28, label = "One KEGG term passed FDR <0.05", size = 2.35, colour = "#4B5563") +
    ggplot2::labs(title = "KEGG enrichment", x = "Fold enrichment", y = NULL, size = "Genes", fill = "-log10 FDR") +
    v32_theme(6.9) + ggplot2::theme(legend.position = "bottom")

  pathways <- utils::read.csv(file.path(paths$tables, "Table_S18_Hallmark_pathway_direction_matrix.csv"))
  both <- pathways[as.logical(pathways$both_significant), ]
  hall <- rbind(
    data.frame(pathway = both$pathway, disease = "OA", NES = both$OA_NES),
    data.frame(pathway = both$pathway, disease = "OC", NES = both$OC_NES)
  )
  order_path <- both$pathway[order(both$paired_direction_index)]
  hall$pathway <- factor(v32_wrap(hall$pathway, 24L), levels = v32_wrap(order_path, 24L))
  p3 <- ggplot2::ggplot(hall, ggplot2::aes(NES, pathway, colour = disease, shape = disease)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(group = pathway), colour = "#C9CDD2", linewidth = 0.45) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E")) +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17)) +
    ggplot2::labs(title = "Hallmark states shared by significance", subtitle = "10 pathways significant in both diseases", x = "Normalized enrichment score", y = NULL, colour = NULL, shape = NULL) +
    v32_theme(6.5) + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 5.4))

  safe_write_csv(selected_go, file.path(paths$source, "Figure3_representative_GO_terms.csv"))
  safe_write_csv(kegg, file.path(paths$source, "Figure3_KEGG_terms.csv"))
  safe_write_csv(hall, file.path(paths$source, "Figure3_joint_Hallmark_states.csv"))
  figure <- submission_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(1, 1.05), widths = c(1.35, 0.65)))
  submission_save_plot(figure, "Figure3_functional_characterization", paths$figures, 185, 178)
}

v32_build_figure4 <- function(paths) {
  shared <- utils::read.csv(file.path(paths$tables, "Table_S2_shared_differentially_expressed_genes.csv"))
  if (!"direction_quadrant" %in% names(shared)) {
    shared$direction_quadrant <- ifelse(shared$logFC_OA > 0 & shared$logFC_OC > 0, "OA higher / OC higher",
      ifelse(shared$logFC_OA < 0 & shared$logFC_OC < 0, "OA lower / OC lower",
        ifelse(shared$logFC_OA > 0, "OA higher / OC lower", "OA lower / OC higher")))
  }
  levels_q <- c("OA higher / OC higher", "OA lower / OC lower", "OA higher / OC lower", "OA lower / OC higher")
  quadrants <- as.data.frame(table(factor(shared$direction_quadrant, levels = levels_q)), stringsAsFactors = FALSE)
  names(quadrants) <- c("quadrant", "genes")
  quadrants$direction <- c("Concordant", "Concordant", "Discordant", "Discordant")
  quadrants$label <- c("OA+ / OC+", "OA- / OC-", "OA+ / OC-", "OA- / OC+")
  p1 <- ggplot2::ggplot(quadrants, ggplot2::aes(label, genes, fill = direction)) +
    ggplot2::geom_col(width = 0.66, colour = "#333333", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = genes), vjust = -0.35, size = 2.8) +
    ggplot2::scale_fill_manual(values = c(Concordant = "#2C7FB8", Discordant = "#D95F0E")) +
    ggplot2::labs(title = "Shared-gene direction classes", subtitle = "146 concordant; 140 discordant", x = NULL, y = "Genes", fill = NULL) +
    v32_theme(7.0) + ggplot2::theme(legend.position = "bottom")

  candidates <- utils::read.csv(file.path(paths$tables, "Table_S16_candidate_prioritization_matrix.csv"))
  candidates <- candidates[order(candidates$prioritization_rank), ]
  effects <- rbind(
    data.frame(gene = candidates$gene, disease = "OA", log2FC = candidates$log2FC_OA),
    data.frame(gene = candidates$gene, disease = "OC", log2FC = candidates$log2FC_OC)
  )
  effects$gene <- factor(effects$gene, levels = rev(candidates$gene))
  p2 <- ggplot2::ggplot(effects, ggplot2::aes(disease, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", log2FC)), size = 2.2) +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-5.5, 5.5), oob = scales::squish) +
    ggplot2::labs(title = "Candidate-gene transcriptional divergence", x = NULL, y = NULL, fill = "log2FC") +
    v32_theme(7.0) + ggplot2::theme(panel.grid = ggplot2::element_blank())

  pathway <- utils::read.csv(file.path(paths$tables, "Table_S18_Hallmark_pathway_direction_matrix.csv"))
  pathway$status <- "Not jointly significant"
  pathway$status[as.logical(pathway$both_significant) & pathway$direction_class == "concordant"] <- "Concordant"
  pathway$status[as.logical(pathway$both_significant) & pathway$direction_class == "discordant"] <- "Discordant"
  label_data <- pathway[as.logical(pathway$both_significant), ]
  p3 <- ggplot2::ggplot(pathway, ggplot2::aes(OA_NES, OC_NES)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "#B6BBC2", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = status, shape = status), size = 1.8, alpha = 0.9) +
    ggrepel::geom_text_repel(data = label_data, ggplot2::aes(label = v32_wrap(pathway, 18L)), size = 2.0,
                             max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, box.padding = 0.15) +
    ggplot2::scale_colour_manual(values = c("Not jointly significant" = "#C9CDD2", Concordant = "#2C7FB8", Discordant = "#D95F0E")) +
    ggplot2::scale_shape_manual(values = c("Not jointly significant" = 1, Concordant = 16, Discordant = 17)) +
    ggplot2::labs(title = "Pathway-direction comparison", subtitle = "6 of 10 jointly significant Hallmarks had opposite signs", x = "OA NES", y = "OC NES", colour = NULL, shape = NULL) +
    v32_theme(7.0) + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 5.6))

  safe_write_csv(quadrants, file.path(paths$source, "Figure4_direction_counts.csv"))
  safe_write_csv(effects, file.path(paths$source, "Figure4_candidate_effects.csv"))
  safe_write_csv(pathway, file.path(paths$source, "Figure4_pathway_direction.csv"))
  figure <- submission_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(0.95, 1.15)))
  submission_save_plot(figure, "Figure4_transcriptional_divergence", paths$figures, 185, 175)
}

v32_build_figure5 <- function(paths) {
  wgcna <- utils::read.csv(file.path(paths$tables, "Table_S4_WGCNA_stability.csv"))
  bootstrap <- wgcna[wgcna$source_table == "bootstrap", ]
  bootstrap$disease <- factor(bootstrap$disease, levels = c("OA", "OC"))
  p1 <- ggplot2::ggplot(bootstrap, ggplot2::aes(disease, observed_correlation, colour = disease, shape = disease)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = bootstrap_lower_95, ymax = bootstrap_upper_95), width = 0.12, linewidth = 0.65) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_colour_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E"), guide = "none") +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(-1, 0)) +
    ggplot2::labs(title = "Disease-associated WGCNA modules", subtitle = "Observed correlation with 95% bootstrap interval", x = NULL, y = "Module-trait correlation") +
    v32_theme(7.0)

  edges <- utils::read.csv(file.path(paths$tables, "Table_S25b_direction_aware_STRING_edges.csv"))
  edges <- edges[edges$network_type == "high-confidence physical", ]
  nodes <- utils::read.csv(file.path(paths$tables, "Table_S25c_STRING_node_topology.csv"))
  nodes <- nodes[nodes$network_type == "high-confidence physical" & nodes$degree > 0, ]
  graph <- igraph::graph_from_data_frame(edges[, c("gene_A", "gene_B")], directed = FALSE, vertices = nodes$gene)
  set.seed(20260726L)
  layout <- igraph::layout_with_fr(graph)
  node_plot <- data.frame(gene = igraph::V(graph)$name, x = layout[, 1], y = layout[, 2], stringsAsFactors = FALSE)
  node_plot <- merge(node_plot, nodes[, c("gene", "degree", "candidate", "direction_binary")], by = "gene", all.x = TRUE)
  edge_plot <- data.frame(
    x = node_plot$x[match(igraph::as_edgelist(graph)[, 1], node_plot$gene)],
    y = node_plot$y[match(igraph::as_edgelist(graph)[, 1], node_plot$gene)],
    xend = node_plot$x[match(igraph::as_edgelist(graph)[, 2], node_plot$gene)],
    yend = node_plot$y[match(igraph::as_edgelist(graph)[, 2], node_plot$gene)]
  )
  label_nodes <- node_plot[node_plot$candidate | rank(-node_plot$degree, ties.method = "first") <= 8L, ]
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_plot, ggplot2::aes(x, y, xend = xend, yend = yend), colour = "#B8BDC5", linewidth = 0.35) +
    ggplot2::geom_point(data = node_plot, ggplot2::aes(x, y, fill = direction_binary, size = degree, shape = candidate), colour = "#333333", stroke = 0.3, alpha = 0.9) +
    ggrepel::geom_text_repel(data = label_nodes, ggplot2::aes(x, y, label = gene), size = 2.0, min.segment.length = 0, max.overlaps = Inf, segment.size = 0.2) +
    ggplot2::scale_fill_manual(values = c(concordant = "#2C7FB8", discordant = "#D95F0E")) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.4, 4.8)) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        title = "Direction",
        nrow = 1,
        override.aes = list(shape = 21, size = 3, colour = "#333333")
      ),
      size = ggplot2::guide_legend(order = 2, title = "Degree", nrow = 1)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Auxiliary STRING association context", subtitle = "46 connected products; only JUNB among ten candidates was connected", fill = "Direction", size = "Degree") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 7.8),
      plot.subtitle = ggplot2::element_text(size = 6.0, colour = "#4B5563"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.3, "mm"),
      legend.text = ggplot2::element_text(size = 5.2),
      legend.title = ggplot2::element_text(size = 5.2)
    )

  candidates <- utils::read.csv(file.path(paths$tables, "Table_S16_candidate_prioritization_matrix.csv"))
  candidates <- candidates[order(candidates$prioritization_rank), ]
  effects <- rbind(
    data.frame(gene = candidates$gene, disease = "OA", log2FC = candidates$log2FC_OA),
    data.frame(gene = candidates$gene, disease = "OC", log2FC = candidates$log2FC_OC)
  )
  effects$gene <- factor(effects$gene, levels = rev(candidates$gene))
  p3 <- ggplot2::ggplot(effects, ggplot2::aes(disease, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", log2FC)), size = 2.1) +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-5.5, 5.5), oob = scales::squish) +
    ggplot2::labs(title = "Candidate bulk-expression effects", x = NULL, y = NULL, fill = "log2FC") +
    v32_theme(6.8) + ggplot2::theme(panel.grid = ggplot2::element_blank())

  node_candidate <- nodes[nodes$gene %in% candidates$gene, c("gene", "degree")]
  candidates$ppi_connected <- candidates$gene %in% node_candidate$gene[node_candidate$degree > 0]
  evidence <- rbind(
    data.frame(gene = candidates$gene, layer = "Shared DEG", value = 1),
    data.frame(gene = candidates$gene, layer = "OA WGCNA", value = as.numeric(as.logical(candidates$OA_primary_WGCNA_module))),
    data.frame(gene = candidates$gene, layer = "OC WGCNA", value = as.numeric(as.logical(candidates$OC_primary_WGCNA_module))),
    data.frame(gene = candidates$gene, layer = "Model votes", value = pmin(candidates$original_model_vote_count / 4, 1)),
    data.frame(gene = candidates$gene, layer = "PPI connected", value = as.numeric(candidates$ppi_connected)),
    data.frame(gene = candidates$gene, layer = "OA scRNA", value = as.numeric(candidates$OA_top_context_detection_fraction > 0)),
    data.frame(gene = candidates$gene, layer = "OC scRNA", value = as.numeric(candidates$OC_top_context_detection_fraction > 0))
  )
  evidence$gene <- factor(evidence$gene, levels = rev(candidates$gene))
  evidence$layer <- factor(evidence$layer, levels = c("Shared DEG", "OA WGCNA", "OC WGCNA", "Model votes", "PPI connected", "OA scRNA", "OC scRNA"))
  p4 <- ggplot2::ggplot(evidence, ggplot2::aes(layer, gene, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_point(data = evidence[evidence$value > 0, ], shape = 21, fill = "white", colour = "#333333", size = 1.4, stroke = 0.3) +
    ggplot2::scale_fill_gradient(low = "#F1F3F5", high = "#2C7FB8", limits = c(0, 1)) +
    ggplot2::labs(title = "Transparent evidence map", subtitle = "Model votes are secondary ranking evidence", x = NULL, y = NULL, fill = "Support") +
    v32_theme(6.3) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), panel.grid = ggplot2::element_blank())

  safe_write_csv(bootstrap, file.path(paths$source, "Figure5_WGCNA_primary_modules.csv"))
  safe_write_csv(node_plot, file.path(paths$source, "Figure5_STRING_nodes.csv"))
  safe_write_csv(edge_plot, file.path(paths$source, "Figure5_STRING_edges.csv"))
  safe_write_csv(effects, file.path(paths$source, "Figure5_candidate_effects.csv"))
  safe_write_csv(evidence, file.path(paths$source, "Figure5_candidate_evidence_map.csv"))
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4) + patchwork::plot_layout(heights = c(1.1, 1)))
  submission_save_plot(figure, "Figure5_candidate_gene_characterization", paths$figures, 185, 178)
}

v32_build_figure6 <- function(project_root, paths) {
  umap <- utils::read.csv(file.path(project_root, "results", "submission_v31", "figures", "source_data", "SupplementaryFigure3_UMAP_subsamples.csv"))
  umap <- umap[umap$dataset_id %in% c("GSE255460", "GSE154600"), ]
  labels <- sort(unique(umap$label))
  colours <- grDevices::hcl.colors(length(labels), palette = "Dark 3")
  names(colours) <- labels
  make_umap <- function(dataset_id, title) {
    data <- umap[umap$dataset_id == dataset_id, ]
    centroids <- aggregate(cbind(UMAP1, UMAP2) ~ label, data, stats::median)
    ggplot2::ggplot(data, ggplot2::aes(UMAP1, UMAP2, colour = label)) +
      ggplot2::geom_point(size = 0.13, alpha = 0.48) +
      ggrepel::geom_text_repel(data = centroids, ggplot2::aes(label = label), size = 1.9, colour = "#111111", min.segment.length = 0, max.overlaps = 40, segment.size = 0.18) +
      ggplot2::scale_colour_manual(values = colours, guide = "none") +
      ggplot2::labs(title = title, subtitle = "Exact source labels; deterministic display subsample", x = "UMAP 1", y = "UMAP 2") +
      v32_theme(7.0) + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  }
  p1 <- make_umap("GSE255460", "OA atlas: GSE255460")
  p2 <- make_umap("GSE154600", "OC atlas: GSE154600")

  detection <- utils::read.csv(file.path(paths$tables, "Table_S30a_extended_gene_cell_detection.csv"))
  composition <- utils::read.csv(file.path(paths$tables, "Table_S29a_panel_size_composition.csv"))
  detection <- merge(detection, composition[, c("gene", "original_rank")], by = "gene")
  detection <- detection[detection$original_rank <= 10L, ]
  type_score <- aggregate(fraction_detected ~ dataset_id + disease + cell_type, detection, mean)
  selected_types <- do.call(rbind, lapply(split(type_score, type_score$dataset_id), function(data) head(data[order(-data$fraction_detected), ], 6L)))
  selected_key <- paste(selected_types$dataset_id, selected_types$cell_type)
  detection <- detection[paste(detection$dataset_id, detection$cell_type) %in% selected_key, ]
  detection$gene <- factor(detection$gene, levels = rev(composition$gene[composition$original_rank <= 10L]))
  detection$atlas <- ifelse(detection$disease == "OA", "OA: GSE255460", "OC: GSE154600")
  p3 <- ggplot2::ggplot(detection, ggplot2::aes(cell_type, gene, size = fraction_detected, fill = fraction_detected)) +
    ggplot2::geom_point(shape = 21, colour = "#333333", stroke = 0.3) +
    ggplot2::facet_grid(. ~ atlas, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient(low = "#F5F5F5", high = "#D95F0E") +
    ggplot2::scale_size_continuous(range = c(0.6, 5.0)) +
    ggplot2::labs(title = "Candidate detection across leading source-defined cell labels", subtitle = "Within-atlas detection fractions; no numerical OA-OC comparison", x = NULL, y = NULL, size = "Detected", fill = "Detected") +
    v32_theme(6.6) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "bottom", panel.grid = ggplot2::element_blank())

  safe_write_csv(umap, file.path(paths$source, "Figure6_representative_UMAPs.csv"))
  safe_write_csv(detection, file.path(paths$source, "Figure6_candidate_cell_detection.csv"))
  figure <- submission_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(1.05, 1)))
  submission_save_plot(figure, "Figure6_single_cell_localization", paths$figures, 185, 176)
}

v32_build_figure7 <- function(paths) {
  boxes <- data.frame(
    x = c(2, 2, 1, 3), y = c(2.75, 1.75, 0.65, 0.65),
    width = c(1.65, 1.85, 1.45, 1.45), height = c(0.62, 0.72, 0.72, 0.72),
    title = c("Partially shared transcriptomic alterations", "Directionally divergent states\nin distinct cellular contexts", "OA phenotype", "OC phenotype"),
    detail = c("286 shared DEGs; ECM, immune, stress, and cell-cycle themes", "~49% discordant genes; 6/10 joint Hallmarks opposite", "cartilage degeneration context", "malignant tissue context"),
    fill = c("#F2ECF7", "#FBE9E3", "#E8F1F8", "#FBE9E3"), stringsAsFactors = FALSE
  )
  arrows <- data.frame(x = c(2, 1.78, 2.22), y = c(2.42, 1.40, 1.40), xend = c(2, 1, 3), yend = c(2.12, 1.02, 1.02))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend), arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed"), linewidth = 0.6, colour = "#333333") +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill), colour = "#333333", linewidth = 0.5) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y + 0.10, label = title), fontface = "bold", size = 3.05, lineheight = 0.92) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y - 0.16, label = detail), size = 2.25, colour = "#4B5563", lineheight = 0.92) +
    ggplot2::annotate("text", x = 2, y = 0.05,
      label = "Interpretation: shared membership can reflect recurrent biological themes without implying one conserved disease mechanism.",
      fontface = "italic", size = 2.65, colour = "#374151") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.1, 3.9), ylim = c(-0.1, 3.2), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
  safe_write_csv(boxes, file.path(paths$source, "Figure7_summary_model.csv"))
  submission_save_plot(p, "Figure7_summary_model", paths$figures, 185, 118)
}

v32_build_supplementary1 <- function(paths) {
  deg <- utils::read.csv(file.path(paths$tables, "Table_S3a_DEG_threshold_summary.csv"))
  deg$threshold <- factor(
    paste0("|log2FC| >= ", deg$absolute_log2fc_threshold),
    levels = paste0("|log2FC| >= ", sort(unique(deg$absolute_log2fc_threshold)))
  )
  p1 <- ggplot2::ggplot(
    deg,
    ggplot2::aes(
      absolute_log2fc_threshold,
      shared_count,
      colour = factor(fdr_threshold),
      shape = factor(fdr_threshold)
    )
  ) +
    ggplot2::geom_line(linewidth = 0.6) + ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = c("0.01" = "#2C7FB8", "0.05" = "#D95F0E")) +
    ggplot2::scale_shape_manual(values = c("0.01" = 16, "0.05" = 17)) +
    ggplot2::labs(title = "Shared-DEG threshold sensitivity", x = "Absolute log2FC threshold", y = "Shared genes", colour = "FDR", shape = "FDR") +
    v32_theme(7.0) + ggplot2::theme(legend.position = "bottom")
  wgcna <- utils::read.csv(file.path(paths$tables, "Table_S4_WGCNA_stability.csv"))
  perturb <- wgcna[wgcna$source_table == "soft_power_perturbation", ]
  p2 <- ggplot2::ggplot(perturb, ggplot2::aes(soft_power, absolute_module_trait_correlation, colour = disease, shape = is_primary)) +
    ggplot2::geom_line(ggplot2::aes(group = disease), linewidth = 0.55) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E")) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16)) +
    ggplot2::labs(title = "WGCNA soft-power stability", x = "Soft power", y = "|Module-trait correlation|", colour = NULL, shape = "Primary") +
    v32_theme(7.0) + ggplot2::theme(legend.position = "bottom")
  safe_write_csv(deg, file.path(paths$source, "SupplementaryFigure1_DEG_sensitivity.csv"))
  safe_write_csv(perturb, file.path(paths$source, "SupplementaryFigure1_WGCNA_sensitivity.csv"))
  figure <- submission_panel_tag(p1 | p2)
  submission_save_plot(figure, "SupplementaryFigure1_core_sensitivity", paths$figures, 185, 92)
}

v32_write_documentation <- function(paths) {
  legends <- c(
    "## Figure 1. Study workflow focused on discovery, divergence, and cellular localization", "",
    "OA and OC were analyzed separately through disease-specific differential expression. The shared-gene set was then characterized by functional enrichment, transcriptional direction, transparent candidate evidence, and exact-label single-cell localization. Genetic-causal and communication-network analyses were outside the V3.2 submission scope.", "",
    "## Figure 2. Shared transcriptomic alterations between OA and OC", "",
    "**A-B,** Discovery-cohort volcano plots using FDR <0.05 and absolute log2 fold change >=1. **C,** Primary DEG overlap. **D,** Direction-ordered heatmap of the 286 shared DEGs. Heatmap rows are genes sorted by direction class; no cross-disease normalization is implied.", "",
    "## Figure 3. Functional characterization of the shared-gene set", "",
    "**A,** Ten representative, non-redundant significant GO terms selected to span cell-cycle, matrix, immune, and stress categories; the complete enrichment table remains in Table S11b. **B,** The single KEGG term that passed FDR <0.05. **C,** OA and OC normalized enrichment scores for the ten Hallmark pathways significant in both diseases. These displays summarize recurring biological themes rather than one conserved mechanism.", "",
    "## Figure 4. Transcriptional divergence within the shared-gene set", "",
    "**A,** Four directional quadrants among the 286 shared DEGs. **B,** OA and OC log2 fold changes for the ten transparently ranked candidate genes. **C,** OA versus OC Hallmark normalized enrichment scores; jointly significant concordant and discordant pathways are highlighted. Directional differences are descriptive and do not establish regulatory decoupling or mechanism.", "",
    "## Figure 5. Candidate-gene characterization across co-expression and association contexts", "",
    "**A,** Primary disease-associated WGCNA module-trait correlations with bootstrap 95% intervals. **B,** Connected high-confidence physical STRING association graph among shared-DEG products; only JUNB among the ten candidates was connected, and topology did not determine candidate selection. **C,** Candidate bulk-expression effects. **D,** Transparent evidence matrix. The ten genes are descriptive candidates, not a diagnostic signature or therapeutic target set.", "",
    "## Figure 6. Representative single-cell atlases reveal distinct cellular localization", "",
    "**A-B,** Representative OA GSE255460 and OC GSE154600 embeddings with exact source labels and deterministic display subsampling. **C,** Candidate detection fractions across the six leading labels within each atlas. Values support within-atlas localization only and are not compared numerically across OA and OC.", "",
    "## Figure 7. Summary model of partial overlap and context-dependent divergence", "",
    "Partially shared transcriptomic alterations can recur across diseases while taking different directions and occupying different cell contexts. The summary is an observational interpretation and does not imply a common cause, conserved disease mechanism, or shared treatment target.", "",
    "## Supplementary Figure 1. Core DEG and WGCNA sensitivity analyses", "",
    "Shared-DEG retention across prespecified thresholds and disease-specific WGCNA soft-power perturbation.", "",
    "## Supplementary Figure 2. Candidate-set size sensitivity", "",
    "Direction composition, descriptive Hallmark profiles, and within-atlas localization for deterministic top-5, top-10, and top-15 candidate sets.", "",
    "## Supplementary Figure 3. Exploratory cohort-dependent classification", "",
    "Unsupervised structure, standardized score contrasts, label-permutation distributions, and ROC curves. This is a supplementary reproducibility check, not a clinical prediction model.", "",
    "## Supplementary Figure 4. Complete Hallmark direction matrix", "",
    "Complete paired Hallmark normalized enrichment scores and direction index. Matching significance does not imply matching pathway state.", "",
    "## Supplementary Figure 5. Dataset-specific single-cell embeddings", "",
    "All five single-cell atlases with exact source labels; OA and OC were not integrated into one latent space.", "",
    "## Supplementary Figure 6. Discovery-cohort bulk PCA and sample-correlation QC", "",
    "Separate unsupervised PCA and sample-correlation audits for OA and OC discovery cohorts; no outcome-informed sample exclusion was performed.", "",
    "## Supplementary Figure 7. Candidate-ranking and nested-stability details", "",
    "WGCNA stability, original penalized-regression/random-forest evidence, and strict outer-fold selection frequencies. These are ranking and sensitivity layers, not clinical model validation."
  )
  writeLines(legends, file.path(paths$figures, "figure_legends.md"), useBytes = TRUE)

  table_index <- data.frame(
    table_id = c("S1", "S2", "S3a-b", "S4", "S5", "S6", "S9", "S10", "S11a-c", "S16", "S18", "S19", "S21-S22", "S24a-c", "S25a-e", "S29-S30", "S31a-b"),
    scope = c(
      "datasets and cohorts", "shared DEGs", "DEG threshold sensitivity", "WGCNA stability",
      "secondary candidate-ranking stability", "external exploratory classification", "single-cell QC",
      "candidate detection and eligible pseudobulk evidence", "Hallmark, GO, and KEGG enrichment",
      "transparent candidate evidence matrix", "complete Hallmark direction matrix", "gene-cell-function context",
      "external-cohort separability and effect sizes", "CCSS and sample-aware UCell summaries",
      "STRING mapping, edges, topology, and permutation", "candidate-size and localization sensitivity",
      "bulk discovery PCA and sample-correlation QC"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(table_index, file.path(paths$tables, "supplementary_table_index_v32.csv"))

  registry <- data.frame(
    claim_id = sprintf("V32-C%02d", 1:10),
    claim = c(
      "OA and OC share 286 primary DEGs", "Shared-gene direction is approximately balanced",
      "Ten Hallmarks are significant in both diseases", "Six jointly significant Hallmarks have opposite signs",
      "Representative significant GO terms include matrix, immune, stress, and cell-cycle categories",
      "The ten candidates are a transparent descriptive set", "STRING topology did not select the candidates",
      "Candidate localization is atlas and label specific", "Exploratory classification is cohort dependent",
      "Shared membership does not establish one conserved disease mechanism"
    ),
    evidence = c(
      "Figure 2; Table S2", "Figure 4A; Table S2", "Figures 3C/4C; Table S18", "Figures 3C/4C; Table S18",
      "Figure 3A; Table S11b", "Figure 5; Table S16", "Figure 5B; Tables S25a-e", "Figure 6; Tables S24/S30",
      "Figure S3; Tables S6/S21-S22", "Figures 4/6/7"
    ),
    prohibited_overstatement = c(
      "shared cause", "uniform shared program", "identical pathway activation", "pathway decoupling mechanism",
      "aging enrichment when not observed", "diagnostic signature", "PPI-selected hub genes", "homologous cell states",
      "clinical diagnostic performance", "common pathogenic mechanism"
    ),
    status = "audited",
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v32.csv"))
  checklist <- data.frame(
    item = c("Seven main figures", "Seven supplementary figures", "MR removed from submission scope", "CellChat/NicheNet removed", "TF-miRNA removed", "ROC supplementary only", "Exact source labels retained", "Common journal-style figure system"),
    status = "complete",
    evidence = c("Figures 1-7", "Figures S1-S7", "V32_scope_decisions.csv", "V32_scope_decisions.csv", "V32_scope_decisions.csv", "Figure S3", "Figure 6", "white background; Arial; consistent line weights; blue/orange palette; panel labels"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(checklist, file.path(paths$root, "reproducibility_checklist_v32.csv"))
}

run_reviewer_v32 <- function(project_root) {
  v32_base_runner(project_root)
  paths <- v32_output_paths(project_root)
  log_info("V3.2: applying strategic submission-scope simplification.")
  v32_prepare_scope(project_root, paths)
  v32_copy_supplementary_figures(project_root, paths)
  log_info("V3.2: rebuilding seven discovery-direction-localization main figures.")
  v32_build_figure1(paths)
  v32_build_figure2(project_root, paths)
  v32_build_figure3(paths)
  v32_build_figure4(paths)
  v32_build_figure5(paths)
  v32_build_figure6(project_root, paths)
  v32_build_figure7(paths)
  v32_build_supplementary1(paths)
  v32_write_documentation(paths)
  log_info("V3.2 strategic reconstruction completed: seven main and seven supplementary figures; deleted modules retained only in historical V3.1 outputs.")
  invisible(paths)
}
