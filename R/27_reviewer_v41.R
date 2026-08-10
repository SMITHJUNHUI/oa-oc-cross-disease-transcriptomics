v41_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v41"))
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

v41_reword_function <- function(fun, replacements) {
  source_text <- paste(deparse(fun), collapse = "\n")
  for (old in names(replacements)) {
    source_text <- gsub(old, replacements[[old]], source_text, fixed = TRUE)
  }
  eval(parse(text = source_text), envir = .GlobalEnv)
}

v41_build_figure7 <- function(paths) {
  nodes <- data.frame(
    x = rep(2, 6),
    y = c(5.5, 4.5, 3.5, 2.5, 1.5, 0.5),
    width = c(2.25, 2.25, 2.25, 2.25, 2.25, 2.25),
    height = rep(0.64, 6),
    title = c(
      "Tissue overlap", "Directional heterogeneity", "External replication",
      "Source-defined cell localization", "Blood persistence", "G0S2 candidate"
    ),
    detail = c(
      "286 shared DEGs",
      "146 concordant; 140 discordant",
      "stronger in OC; variable in OA",
      "separate OA and OC atlases",
      "one gene passed both blood FDR filters",
      "limited systemic component requiring confirmation"
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

v41_write_documentation <- function(project_root, paths) {
  legends_source <- file.path(project_root, "manuscript", "figure_legends_v41.md")
  if (!file.copy(legends_source, file.path(paths$figures, "figure_legends.md"), overwrite = TRUE)) {
    stop("Could not install V4.1 figure legends.", call. = FALSE)
  }
  terminology <- data.frame(
    concept = c("cross-disease overlap", "functional recurrence", "five-gene display", "blood result", "single-cell populations", "overall interpretation"),
    canonical_term = c("shared molecular features", "recurring biological themes", "illustrative genes", "candidate systemic molecular signal", "exact source-defined cell populations", "partial molecular convergence with context dependence"),
    avoided_term = c("shared disease mechanism", "activated mechanism", "predictive panel", "blood biomarker", "forced homologous cell labels", "causal relationship"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(terminology, file.path(paths$analysis, "V41_terminology_ledger.csv"))
  registry <- data.frame(
    claim = c("286 shared tissue genes", "146 concordant and 140 discordant genes", "external replication differs between OA and OC", "illustrative genes have source-defined localization", "G0S2 is the sole dual-blood-FDR result"),
    evidence = c("Figure 2; Table S2", "Figure 2; Table S2", "Figure 4; Tables S4-S6", "Figure 5; Tables S8a-S8b", "Figure 6; Tables S9-S11"),
    boundary = c("overlap is threshold dependent", "direction is not a mechanistic assay", "OA replication is weaker and cohort dependent", "localization does not establish different function", "prospective and protein-level confirmation is required"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v41.csv"))
}

run_reviewer_v41 <- function(project_root) {
  source_root <- file.path(project_root, "results", "submission_v40")
  paths <- v41_output_paths(project_root)
  v40_copy_tree(source_root, paths$root)

  external <- v34_external_validation(project_root, paths)
  genes <- c("G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3")
  unlink(file.path(paths$tables, "Table_S8b_representative_gene_single_cell_detection.csv"), force = TRUE)
  unlink(list.files(paths$figures, pattern = "^Figure4_external_tissue_and_representative_genes\\.", full.names = TRUE), force = TRUE)
  unlink(file.path(paths$source, c("Figure4_representative_gene_effects.csv", "Figure5_representative_UMAPs.csv")), force = TRUE)
  single_builder <- v41_reword_function(
    v34_single_cell_candidates,
    c("Table_S8b_representative_gene_single_cell_detection.csv" = "Table_S8b_illustrative_gene_single_cell_detection.csv")
  )
  single_cell <- single_builder(project_root, paths, genes)
  candidate_builder <- v41_reword_function(
    v34_candidate_evidence,
    c("Representative evidence summary; not an optimized signature" = "Illustrative evidence summary; not an optimized panel")
  )
  candidates <- candidate_builder(paths, external, single_cell, genes)

  v41_reword_function(v34_build_figure1, c("Representative genes" = "Illustrative genes"))(paths)
  v41_reword_function(v34_build_figure3, c("Representative significant GO terms" = "Illustrative significant GO terms"))(project_root, paths)
  v41_reword_function(
    v34_build_figure4,
    c(
      "Representative gene effects" = "Illustrative gene effects",
      "Five genes selected to summarize distinct evidence roles" = "Five genes display distinct evidence roles",
      "Representative genes, not a predictive signature" = "Illustrative genes across evidence layers",
      "Figure4_representative_gene_effects.csv" = "Figure4_illustrative_gene_effects.csv",
      "Figure4_external_tissue_and_representative_genes" = "Figure4_external_tissue_and_illustrative_genes"
    )
  )(paths, external, candidates)
  v41_reword_function(
    v34_build_figure5,
    c(
      "Representative-gene detection by source-defined cell label" = "Illustrative-gene detection by source-defined cell label",
      "Figure5_representative_UMAPs.csv" = "Figure5_illustrative_UMAPs.csv"
    )
  )(project_root, paths, single_cell, genes)
  v41_build_figure7(paths)
  v41_write_documentation(project_root, paths)

  log_info("V4.1 figure package completed: data are unchanged; terminology and the integrated evidence sequence were strategically compressed.")
  invisible(paths)
}
