v33_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v33"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs"))
  )
}

v33_copy <- function(source, target) {
  if (!file.exists(source)) stop("Missing V3.3 source: ", source, call. = FALSE)
  if (!file.copy(source, target, overwrite = TRUE)) stop("Could not copy: ", source, call. = FALSE)
  invisible(target)
}

v33_prepare_scope <- function(project_root, paths) {
  v32 <- file.path(project_root, "results", "submission_v32")
  blood <- file.path(project_root, "results", "blood_v33_internal")
  decision <- utils::read.csv(file.path(blood, "blood_module_decision.csv"))
  if (!identical(decision$manuscript_status[[1L]], "ELIGIBLE") || decision$positive_gene_count[[1L]] < 1L) {
    stop("The prespecified blood module is not eligible for V3.3.", call. = FALSE)
  }

  retained <- c(
    "Table_S1_data_sources_and_cohorts.csv",
    "Table_S2_shared_differentially_expressed_genes.csv",
    "Table_S3a_DEG_threshold_summary.csv",
    "Table_S3b_DEG_threshold_membership.csv",
    "Table_S4_WGCNA_stability.csv",
    "Table_S9_single_cell_QC_and_status.csv",
    "Table_S10_single_cell_hub_gene_evidence.csv",
    "Table_S11a_Hallmark_GSEA.csv",
    "Table_S11b_GO_shared_genes.csv",
    "Table_S11c_KEGG_shared_genes.csv",
    "Table_S18_Hallmark_pathway_direction_matrix.csv",
    "Table_S19_gene_cell_function_context_matrix.csv",
    "Table_S24a_dataset_context_CCSS.csv",
    "Table_S24b_disease_consensus_CCSS.csv",
    "Table_S24c_sample_aware_UCell_summary.csv",
    "Table_S25a_STRING_mapping_audit.csv",
    "Table_S25b_direction_aware_STRING_edges.csv",
    "Table_S25c_STRING_node_topology.csv",
    "Table_S25d_STRING_network_topology.csv",
    "Table_S25e_STRING_direction_label_permutation.csv",
    "Table_S30a_extended_gene_cell_detection.csv",
    "Table_S31a_bulk_unsupervised_PCA.csv",
    "Table_S31b_bulk_sample_correlation_QC.csv"
  )
  for (name in retained) {
    v33_copy(file.path(v32, "supplementary_tables", name), file.path(paths$tables, name))
  }

  representative <- utils::read.csv(file.path(v32, "supplementary_tables", "Table_S16_candidate_prioritization_matrix.csv"))
  representative <- representative[, c(
    "gene", "log2FC_OA", "log2FC_OC", "direction_class",
    "OA_primary_WGCNA_module", "OC_primary_WGCNA_module",
    "OA_top_cell_context", "OA_top_context_detection_fraction",
    "OC_top_cell_context", "OC_top_context_detection_fraction"
  )]
  names(representative)[1L] <- "representative_shared_gene"
  safe_write_csv(representative, file.path(paths$tables, "Table_S12_representative_shared_gene_context.csv"))

  dataset_audit <- utils::read.csv(file.path(blood, "blood_dataset_audit.csv"))
  positive <- utils::read.csv(file.path(blood, "FDR_supported_systemic_component_genes.csv"))
  attrition <- utils::read.csv(file.path(blood, "blood_screen_attrition.csv"))
  safe_write_csv(dataset_audit, file.path(paths$tables, "Table_S13_blood_dataset_audit.csv"))
  safe_write_csv(positive, file.path(paths$tables, "Table_S14_FDR_supported_systemic_component.csv"))
  safe_write_csv(attrition, file.path(paths$tables, "Table_S15_blood_screen_attrition.csv"))

  scope <- data.frame(
    module = c(
      "Bulk tissue DEG", "Functional and direction analysis", "WGCNA", "STRING",
      "Representative shared-gene localization", "Peripheral-blood systemic-component screen",
      "AUC/complex machine learning", "MR", "CellChat/NicheNet", "TF-miRNA", "DCA/nomogram"
    ),
    status = c(
      "main", "main", "supplementary stability", "auxiliary",
      "main", "main only because dual-cohort FDR rule passed",
      "excluded", "excluded", "excluded", "excluded", "not performed"
    ),
    rationale = c(
      "primary discovery", "biological characterization with explicit directions",
      "within-cohort co-expression support", "database association context; not candidate selection",
      "cellular localization only", "systemic-component validation, not a blood biomarker study",
      "not a prediction study", "outside the focused transcriptomic question",
      "no direct validated communication axis", "prediction-only regulatory evidence",
      "no locked clinical probability model"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(scope, file.path(paths$analysis, "V33_scope_decisions.csv"))
  invisible(list(v32 = v32, blood = blood, positive = positive, attrition = attrition, dataset_audit = dataset_audit))
}

v33_copy_retained_figures <- function(project_root, paths) {
  v32_figures <- file.path(project_root, "results", "submission_v32", "figures")
  mapping <- c(
    "Figure2_shared_transcriptomic_alterations" = "Figure2_shared_transcriptomic_alterations",
    "Figure3_functional_characterization" = "Figure3_functional_characterization",
    "Figure4_transcriptional_divergence" = "Figure4_transcriptional_divergence",
    "Figure6_single_cell_localization" = "Figure6_single_cell_localization",
    "SupplementaryFigure1_core_sensitivity" = "SupplementaryFigure1_core_sensitivity",
    "SupplementaryFigure4_complete_Hallmark_direction" = "SupplementaryFigure2_complete_Hallmark_direction",
    "SupplementaryFigure5_all_single_cell_UMAPs" = "SupplementaryFigure3_all_single_cell_UMAPs",
    "SupplementaryFigure6_bulk_PCA_QC" = "SupplementaryFigure4_bulk_PCA_QC"
  )
  for (extension in c("png", "pdf")) {
    for (source_stem in names(mapping)) {
      v33_copy(
        file.path(v32_figures, paste0(source_stem, ".", extension)),
        file.path(paths$figures, paste0(mapping[[source_stem]], ".", extension))
      )
    }
    v33_copy(
      file.path(project_root, "results", "submission_v31", "figures", paste0("SupplementaryFigure12_direction_aware_STRING_network.", extension)),
      file.path(paths$figures, paste0("SupplementaryFigure5_direction_aware_STRING_network.", extension))
    )
  }
}

v33_rebuild_hallmark_supplement <- function(paths) {
  hallmark <- utils::read.csv(
    file.path(paths$tables, "Table_S18_Hallmark_pathway_direction_matrix.csv"),
    check.names = FALSE
  )
  hallmark <- hallmark[order(hallmark$paired_direction_index, decreasing = TRUE), , drop = FALSE]
  hallmark$pathway_label <- factor(hallmark$pathway, levels = rev(hallmark$pathway))
  long <- rbind(
    data.frame(pathway_label = hallmark$pathway_label, disease = "OA", NES = hallmark$OA_NES),
    data.frame(pathway_label = hallmark$pathway_label, disease = "OC", NES = hallmark$OC_NES)
  )
  p1 <- ggplot2::ggplot(long, ggplot2::aes(disease, pathway_label, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.22) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#1B9E77", midpoint = 0) +
    ggplot2::labs(title = "A   Complete Hallmark NES matrix", x = NULL, y = NULL, fill = "NES") +
    submission_theme(6.3) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 4.6),
      legend.position = "top",
      plot.title = ggplot2::element_text(face = "bold", size = 7.2, hjust = 0, margin = ggplot2::margin(b = 4))
    )
  hallmark$direction_colour <- ifelse(hallmark$direction_class == "concordant", "#009E73", "#D55E00")
  p2 <- ggplot2::ggplot(hallmark, ggplot2::aes(paired_direction_index, pathway_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#98A2B3", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = paired_direction_index, yend = pathway_label),
      colour = hallmark$direction_colour, linewidth = 0.45
    ) +
    ggplot2::geom_point(colour = hallmark$direction_colour, size = 1.2) +
    ggplot2::labs(title = "B   Signed paired direction index", x = "OA NES x OC NES", y = NULL) +
    submission_theme(6.3) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 7.2, hjust = 0, margin = ggplot2::margin(b = 4))
    )
  combined <- p1 | p2
  submission_save_plot(
    combined,
    "SupplementaryFigure2_complete_Hallmark_direction",
    paths$figures,
    width_mm = 185,
    height_mm = 185
  )
}

v33_build_figure1 <- function(paths) {
  boxes <- data.frame(
    x = c(1, 3, 2, 1, 3, 2, 2),
    y = c(3.1, 3.1, 2.35, 1.55, 1.55, 0.82, 0.1),
    width = c(1.45, 1.45, 1.65, 1.45, 1.45, 1.65, 1.85),
    height = c(0.56, 0.56, 0.60, 0.58, 0.58, 0.60, 0.52),
    title = c(
      "OA tissue", "OC tissue", "286 shared DEGs",
      "Direction + function", "Single-cell localization",
      "Prespecified blood screen", "One FDR-supported systemic candidate"
    ),
    detail = c(
      "GSE114007", "GSE18520", "membership overlap",
      "GO / KEGG / Hallmark", "exact source labels",
      "OA GSE48556 + OC GSE31682", "G0S2"
    ),
    fill = c("#E8F1F8", "#FBE9E3", "#F2ECF7", "#F7F7F7", "#E8F4F1", "#FFF3E8", "#FBE9E3"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(1.45, 2.55, 2, 1.45, 2.55, 2, 2),
    y = c(2.92, 2.92, 2.05, 2.05, 2.05, 1.25, 0.52),
    xend = c(1.85, 2.15, 1.45, 1, 2.55, 2, 2),
    yend = c(2.57, 2.57, 1.85, 1.85, 1.85, 1.12, 0.36)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"),
      linewidth = 0.55, colour = "#333333"
    ) +
    ggplot2::geom_rect(
      data = boxes,
      ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill),
      colour = "#333333", linewidth = 0.45
    ) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y + 0.08, label = title), fontface = "bold", size = 2.7) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y - 0.13, label = detail), colour = "#4B5563", size = 2.0) +
    ggplot2::annotate(
      "text", x = 2, y = -0.32,
      label = "Blood was a conditional validation filter, not an independent biomarker study.",
      fontface = "italic", size = 2.3, colour = "#374151"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.05, 3.95), ylim = c(-0.42, 3.55), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
  safe_write_csv(boxes, file.path(paths$source, "Figure1_workflow_nodes.csv"))
  submission_save_plot(p, "Figure1_study_workflow_with_blood", paths$figures, 185, 135)
}

v33_build_figure5 <- function(paths, blood_data) {
  positive <- blood_data$positive
  attrition <- blood_data$attrition
  audit <- blood_data$dataset_audit
  if (nrow(positive) != 1L || positive$gene[[1L]] != "G0S2") {
    stop("V3.3 blood figure expects the audited single positive gene G0S2.", call. = FALSE)
  }

  dataset_plot <- data.frame(
    x = c(1, 2), y = 1,
    dataset = audit$dataset_id,
    disease = audit$disease,
    detail = paste0(audit$biospecimen, "\n", audit$disease_n, " disease / ", audit$control_n, " control"),
    fill = c("#E8F1F8", "#FBE9E3")
  )
  p1 <- ggplot2::ggplot(dataset_plot) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.58, ymax = 1.42, fill = fill), colour = "#333333", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(x, 1.18, label = paste0(disease, ": ", dataset)), fontface = "bold", size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x, 0.88, label = detail), size = 2.15, colour = "#4B5563", lineheight = 0.95) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.45, 2.55), ylim = c(0.45, 1.55), clip = "off") +
    ggplot2::labs(title = "Eligible peripheral-blood cohorts", subtitle = "Independent microarrays; disease versus healthy control") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 8), plot.subtitle = ggplot2::element_text(size = 6.2, colour = "#4B5563"))

  attrition$stage <- factor(v32_wrap(attrition$stage, 28L), levels = rev(v32_wrap(attrition$stage, 28L)))
  attrition$highlight <- ifelse(attrition$genes == min(attrition$genes), "Retained", "Screened")
  p2 <- ggplot2::ggplot(attrition, ggplot2::aes(genes, stage, fill = highlight)) +
    ggplot2::geom_col(width = 0.64, colour = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = genes), hjust = -0.25, size = 2.5) +
    ggplot2::scale_fill_manual(values = c(Screened = "#B8C4D0", Retained = "#D95F0E"), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, max(attrition$genes) * 1.15), clip = "off") +
    ggplot2::labs(title = "Prespecified four-layer attrition", x = "Genes", y = NULL) +
    v32_theme(6.7) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5))

  heat <- data.frame(
    gene = "G0S2",
    context = factor(c("OA tissue", "OC tissue", "OA blood", "OC blood"), levels = c("OA tissue", "OC tissue", "OA blood", "OC blood")),
    log2FC = c(positive$logFC_OA, positive$logFC_OC, positive$OA_blood_logFC, positive$OC_blood_logFC)
  )
  p3 <- ggplot2::ggplot(heat, ggplot2::aes(context, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", log2FC)), size = 2.7) +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-2.2, 2.2), oob = scales::squish) +
    ggplot2::labs(title = "Four-context direction", subtitle = "All effects were lower in disease", x = NULL, y = NULL, fill = "log2FC") +
    v32_theme(7.0) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), panel.grid = ggplot2::element_blank())

  forest <- data.frame(
    cohort = factor(c("OA blood: GSE48556", "OC blood: GSE31682"), levels = rev(c("OA blood: GSE48556", "OC blood: GSE31682"))),
    effect = c(positive$OA_blood_hedges_g, positive$OC_blood_hedges_g),
    lower = c(positive$OA_blood_hedges_g_lower, positive$OC_blood_hedges_g_lower),
    upper = c(positive$OA_blood_hedges_g_upper, positive$OC_blood_hedges_g_upper),
    disease = c("OA", "OC")
  )
  p4 <- ggplot2::ggplot(forest, ggplot2::aes(effect, cohort, colour = disease, shape = disease)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = lower, xmax = upper),
      orientation = "y", width = 0.13, linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E"), guide = "none") +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17), guide = "none") +
    ggplot2::labs(title = "Standardized blood effects", subtitle = "Hedges g with 95% confidence interval; no cross-disease pooling", x = "Hedges g", y = NULL) +
    v32_theme(6.7)

  safe_write_csv(dataset_plot, file.path(paths$source, "Figure5_blood_datasets.csv"))
  safe_write_csv(attrition, file.path(paths$source, "Figure5_blood_attrition.csv"))
  safe_write_csv(heat, file.path(paths$source, "Figure5_G0S2_four_context_effects.csv"))
  safe_write_csv(forest, file.path(paths$source, "Figure5_G0S2_blood_standardized_effects.csv"))
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4) + patchwork::plot_layout(heights = c(1, 0.9)))
  submission_save_plot(figure, "Figure5_peripheral_blood_systemic_component", paths$figures, 185, 145)
}

v33_build_figure7 <- function(paths) {
  boxes <- data.frame(
    x = c(2, 2, 1, 3, 2), y = c(3.0, 2.0, 1.0, 1.0, 0.05),
    width = c(1.8, 1.8, 1.45, 1.45, 1.8), height = c(0.58, 0.68, 0.68, 0.68, 0.48),
    title = c(
      "Partially shared tissue alterations", "Directionally heterogeneous states",
      "Distinct cellular contexts", "One blood-replicated alteration", "Bounded interpretation"
    ),
    detail = c(
      "286 shared DEGs", "140/286 discordant; 6/10 Hallmarks opposite",
      "exact OA and OC source labels", "G0S2 lower in all four contexts",
      "limited systemic component, not a biomarker or common mechanism"
    ),
    fill = c("#F2ECF7", "#FBE9E3", "#E8F1F8", "#FFF3E8", "#F7F7F7"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(2, 1.78, 2.22, 1, 3), y = c(2.70, 1.64, 1.64, 0.64, 0.64),
    xend = c(2, 1, 3, 1.75, 2.25), yend = c(2.35, 1.36, 1.36, 0.30, 0.30)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend), arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"), linewidth = 0.55, colour = "#333333") +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill), colour = "#333333", linewidth = 0.45) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y + 0.09, label = title), fontface = "bold", size = 2.75) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x, y = y - 0.13, label = v32_wrap(detail, 42L)), colour = "#4B5563", size = 2.0, lineheight = 0.95) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.05, 3.95), ylim = c(-0.35, 3.45), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
  safe_write_csv(boxes, file.path(paths$source, "Figure7_integrated_model.csv"))
  submission_save_plot(p, "Figure7_integrated_summary_model", paths$figures, 185, 125)
}

v33_write_documentation <- function(paths) {
  legends <- c(
    "## Figure 1. Study workflow with a conditional peripheral-blood validation filter", "",
    "OA and OC tissue cohorts were analyzed separately to identify shared DEGs, functional themes, transcriptional direction, and cell localization. Peripheral blood was evaluated only as a prespecified validation filter. It entered the submitted narrative because one gene passed the dual-cohort FDR and four-context direction rule; no blood biomarker model was built.", "",
    "## Figure 2. Shared transcriptomic alterations between OA and OC", "",
    "**A-B,** Discovery-cohort volcano plots using FDR <0.05 and absolute log2 fold change >=1. **C,** Primary DEG overlap. **D,** Direction-ordered heatmap of the 286 shared DEGs. Rows are genes sorted by direction class; no cross-disease normalization is implied.", "",
    "## Figure 3. Functional characterization of the shared-gene set", "",
    "**A,** Ten representative, non-redundant significant GO terms spanning cell-cycle, matrix, immune, and stress categories. **B,** The only KEGG term passing FDR <0.05. **C,** OA and OC NES values for the ten Hallmark pathways significant in both diseases. These displays summarize recurring themes rather than one conserved mechanism.", "",
    "## Figure 4. Transcriptional divergence within the shared-gene set", "",
    "**A,** Four directional quadrants among 286 shared DEGs. **B,** OA and OC log2 fold changes for ten representative shared genes used for cellular localization. **C,** OA versus OC Hallmark NES values; jointly significant concordant and discordant pathways are highlighted. Directional differences are descriptive and do not establish regulatory mechanism.", "",
    "## Figure 5. Peripheral-blood validation identifies one limited systemic component", "",
    "**A,** Eligible OA PBMC and OC blood-cell-fraction case-control microarrays. **B,** Prespecified attrition from 286 tissue-shared DEGs to one gene meeting tissue concordance, measurement in both blood cohorts, four-context sign agreement, and FDR <0.05 in each blood cohort. **C,** G0S2 log2 fold changes in OA tissue, OC tissue, OA blood, and OC blood. **D,** Within-cohort standardized blood effects. Effects were not pooled across diseases. G0S2 is a blood-replicated systemic-component candidate, not a blood biomarker or shared causal mechanism.", "",
    "## Figure 6. Representative single-cell atlases reveal distinct cellular localization", "",
    "**A-B,** Representative OA GSE255460 and OC GSE154600 embeddings with exact source labels and deterministic display subsampling. **C,** Detection fractions for ten representative shared genes across six leading labels within each atlas. Values support within-atlas localization only and are not compared numerically across OA and OC.", "",
    "## Figure 7. Integrated interpretation of tissue overlap, cellular divergence, and limited blood replication", "",
    "Partially shared tissue alterations coexist with frequent direction differences and distinct source-defined cell contexts. G0S2 was the sole gene satisfying the prespecified dual-blood-cohort FDR rule. The model does not imply a conserved disease mechanism, blood diagnostic signature, or common therapeutic target.", "",
    "## Supplementary Figure 1. Core DEG and WGCNA sensitivity analyses", "",
    "Shared-DEG retention across prespecified thresholds and disease-specific WGCNA soft-power perturbation.", "",
    "## Supplementary Figure 2. Complete Hallmark direction matrix", "",
    "Complete paired Hallmark NES values and direction index. Matching significance does not imply matching pathway state.", "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings", "",
    "All five single-cell atlases with exact source labels; OA and OC were not integrated into one latent space.", "",
    "## Supplementary Figure 4. Discovery-cohort bulk PCA and sample-correlation QC", "",
    "Separate unsupervised PCA and sample-correlation audits for OA and OC discovery cohorts; no outcome-informed sample exclusion was performed.", "",
    "## Supplementary Figure 5. Direction-aware STRING association landscape", "",
    "High-confidence physical STRING association graph, descriptive topology, direction-aware subgraphs, and fixed-size label permutations. STRING topology did not select the representative genes or G0S2."
  )
  writeLines(legends, file.path(paths$figures, "figure_legends.md"), useBytes = TRUE)

  registry <- data.frame(
    claim = c(
      "286 tissue-shared DEGs", "146 concordant and 140 discordant genes",
      "10 jointly significant Hallmarks, 6 with opposite signs",
      "G0S2 was the sole dual-blood-FDR positive gene", "Blood is a limited systemic component only",
      "Single-cell localization is atlas specific", "No prediction or causal claim"
    ),
    evidence = c(
      "Figure 2; Table S2", "Figure 4; Table S2", "Figures 3-4; Table S18",
      "Figure 5; Tables S13-S15", "Figure 5; scope decision", "Figure 6; Tables S9/S12/S24/S30",
      "scope decision; Discussion"
    ),
    boundary = c(
      "not shared mechanism", "membership does not imply direction", "NES is transcriptional state",
      "one candidate, not a signature", "blood cohorts use different platforms and blood fractions",
      "no cross-disease numerical comparison", "observational and retrospective"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v33.csv"))
}

run_reviewer_v33 <- function(project_root) {
  paths <- v33_output_paths(project_root)
  log_info("V3.3: adding a conditional, prespecified peripheral-blood systemic-component validation.")
  blood_data <- v33_prepare_scope(project_root, paths)
  v33_copy_retained_figures(project_root, paths)
  v33_rebuild_hallmark_supplement(paths)
  v33_build_figure1(paths)
  v33_build_figure5(paths, blood_data)
  v33_build_figure7(paths)
  v33_write_documentation(paths)
  log_info("V3.3 completed: seven main figures, five supplementary figures, one dual-blood-FDR systemic candidate, and no AUC/ML submission module.")
  invisible(paths)
}
