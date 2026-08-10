v42_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v42"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs")),
    reference_audit = ensure_dir(file.path(root, "reference_audit"))
  )
}

v42_build_supplementary_figure1 <- function(paths) {
  table <- utils::read.csv(file.path(paths$tables, "Table_S3a_DEG_threshold_summary.csv"), check.names = FALSE)
  table$fdr_label <- factor(
    sprintf("FDR %.2f", table$fdr_threshold),
    levels = c("FDR 0.01", "FDR 0.05")
  )
  primary <- table[as.logical(table$is_primary), , drop = FALSE]
  p <- ggplot2::ggplot(
    table,
    ggplot2::aes(
      x = absolute_log2fc_threshold,
      y = shared_count,
      colour = fdr_label,
      group = fdr_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_point(
      data = primary,
      shape = 21,
      size = 4.2,
      stroke = 0.8,
      fill = "white",
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = primary,
      ggplot2::aes(label = paste0("Primary: ", shared_count)),
      nudge_x = 0.09,
      nudge_y = 65,
      colour = "black",
      size = 2.35,
      hjust = 0
    ) +
    ggplot2::scale_colour_manual(values = c("FDR 0.01" = "#3C8DBC", "FDR 0.05" = "#E66101")) +
    ggplot2::scale_x_continuous(breaks = c(0.5, 1.0, 1.5), limits = c(0.45, 1.62)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(
      x = "Absolute log2-fold-change threshold",
      y = "Shared genes",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 7.2, base_family = "Arial") +
    ggplot2::theme(
      legend.position = "top",
      axis.title = ggplot2::element_text(size = 6.8),
      axis.text = ggplot2::element_text(size = 6.2),
      legend.text = ggplot2::element_text(size = 6.2),
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
  safe_write_csv(table, file.path(paths$source, "SupplementaryFigure1_DEG_threshold_sensitivity.csv"))
  v34_save_plot(p, "SupplementaryFigure1_DEG_threshold_sensitivity", paths$figures, 170, 102)
}

v42_build_figure7 <- function(paths) {
  nodes <- data.frame(
    x = rep(2, 6),
    y = c(5.5, 4.5, 3.5, 2.5, 1.5, 0.5),
    width = rep(2.25, 6),
    height = rep(0.64, 6),
    title = c(
      "Tissue overlap", "Directional heterogeneity", "External replication",
      "Source-defined cell localization", "Blood persistence", "G0S2 association"
    ),
    detail = c(
      "286 shared DEGs",
      "146 concordant; 140 discordant",
      "stronger in OC; variable in OA",
      "separate OA and OC atlases",
      "one gene met both blood-cohort FDR thresholds",
      "association; source and function unresolved"
    ),
    fill = c("#ECE7F2", "#E8F2EC", "#DDEAF4", "#E7F1EF", "#F8E1D8", "#FBE8D6"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = rep(2, 5), y = c(5.16, 4.16, 3.16, 2.16, 1.16),
    xend = rep(2, 5), yend = c(4.84, 3.84, 2.84, 1.84, 0.84)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
      linewidth = 0.45, colour = "#4A4A4A"
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill),
      colour = "#3C3C3C", linewidth = 0.4
    ) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y + 0.09, label = title), size = 2.45, fontface = "bold") +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y - 0.14, label = v32_wrap(detail, 42L)), size = 1.9, colour = "#59636E", lineheight = 0.92) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.05, 5.95), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 8, 6, 8))
  safe_write_csv(nodes, file.path(paths$source, "Figure7_integrated_model.csv"))
  v34_save_plot(p, "Figure7_integrated_interpretation", paths$figures, 170, 104)
}

v42_write_documentation <- function(project_root, paths) {
  legends_source <- file.path(project_root, "manuscript", "figure_legends_v42.md")
  if (!file.copy(legends_source, file.path(paths$figures, "figure_legends.md"), overwrite = TRUE)) {
    stop("Could not install V4.2 figure legends.", call. = FALSE)
  }
  terminology <- data.frame(
    concept = c("cross-disease overlap", "direction result", "five-gene display", "blood result", "single-cell populations", "overall interpretation"),
    canonical_term = c("shared transcriptomic overlap", "directional heterogeneity", "illustrative genes", "blood-persistent G0S2 association", "exact source-defined cell populations", "partial molecular convergence with context dependence"),
    avoided_term = c("shared disease mechanism", "common disease state", "predictive panel", "systemic biomarker", "shared immune-cell origin", "causal relationship"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(terminology, file.path(paths$analysis, "V42_terminology_ledger.csv"))
  registry <- data.frame(
    claim = c("286 shared tissue genes", "overlap does not define a uniform disease state", "external replication differs between OA and OC", "illustrative genes have source-defined localization", "G0S2 is the sole dual-blood-FDR result"),
    evidence = c("Figure 2; Table S2", "Figure 2; Table S2", "Figure 4; Tables S4-S6", "Figure 5; Tables S8a-S8b", "Figure 6; Tables S9-S11"),
    boundary = c("overlap is threshold dependent", "direction is descriptive, not a pathogenesis assay", "OA replication is weaker and cohort dependent", "localization does not establish function or blood-cell source", "source, function and prospective reproducibility remain unresolved"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v42.csv"))
  scope <- data.frame(
    module = c("Tissue DEG", "Functional enrichment", "External tissue replication", "Illustrative genes", "Single-cell localization", "Peripheral blood", "WGCNA", "STRING", "Machine learning/ROC", "MR", "CellChat/NicheNet", "TF-miRNA", "DCA/nomogram"),
    manuscript_status = c("main", "main but compact", "main", "main evidence summary", "main", "main validation layer", "excluded", "supplementary", "excluded", "excluded", "excluded", "excluded", "excluded"),
    rationale = c("primary discovery", "bounded biological interpretation", "cross-cohort reproducibility", "transparent illustration", "cellular localization only", "prespecified blood persistence screen", "did not advance the central evidence chain", "auxiliary database context", "not a prediction study", "outside the focused transcriptomic question", "no validated communication axis", "prediction-only layer", "no clinical probability model"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(scope, file.path(paths$analysis, "V42_scope_decisions.csv"))
}

run_reviewer_v42 <- function(project_root) {
  source_root <- file.path(project_root, "results", "submission_v41")
  paths <- v42_output_paths(project_root)
  v40_copy_tree(source_root, paths$root)

  unlink(file.path(paths$tables, "Table_S15_WGCNA_stability.csv"), force = TRUE)
  for (suffix in c("a_STRING_mapping_audit.csv", "b_STRING_edges.csv", "c_STRING_node_topology.csv")) {
    source <- file.path(paths$tables, paste0("Table_S16", suffix))
    target <- file.path(paths$tables, paste0("Table_S15", suffix))
    unlink(target, force = TRUE)
    if (!file.rename(source, target)) stop("Could not renumber STRING table: ", source, call. = FALSE)
  }
  unlink(list.files(paths$figures, pattern = "^SupplementaryFigure1_core_sensitivity\\.", full.names = TRUE), force = TRUE)
  v42_build_supplementary_figure1(paths)
  v42_build_figure7(paths)
  v42_write_documentation(project_root, paths)

  log_info("V4.2 figure package completed: WGCNA was removed from the submission scope and the batch/G0S2 interpretation was tightened.")
  invisible(paths)
}
