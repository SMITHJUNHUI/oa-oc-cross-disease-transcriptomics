submission_workflow_plot <- function() {
  nodes <- data.frame(
    x = 1:6,
    label = c(
      "Bulk discovery\nOA + OC",
      "Direction-aware\nshared DEGs",
      "WGCNA + strict\nnested ML",
      "Cross-cohort\nreproducibility",
      "TCGA +\nsingle-cell context",
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
      size = 2.55,
      fontface = "bold",
      linewidth = 0,
      label.padding = grid::unit(2.2, "mm"),
      lineheight = 0.95
    ) +
    ggplot2::annotate(
      "label",
      x = 3.5,
      y = 1.43,
      label = paste(
        "Hypothesis: aging-associated, mesenchymal,",
        "and immune-remodeling programs may partially converge"
      ),
      fill = "#F5F0FA",
      colour = submission_palette[["shared"]],
      size = 2.35,
      fontface = "bold",
      linewidth = 0.25,
      label.r = grid::unit(1.5, "mm")
    ) +
    ggplot2::annotate(
      "label",
      x = 3.5,
      y = 0.48,
      label = "Bidirectional MR: null supplementary evidence",
      fill = "white",
      colour = submission_palette[["neutral"]],
      size = 2.3,
      linewidth = 0.25,
      label.r = grid::unit(1.5, "mm")
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.45, 6.55), ylim = c(0.25, 1.65)) +
    ggplot2::theme_void(base_family = "Arial")
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
  permutations <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_permutation_auc.csv"
  ))
  leave_one_out <- utils::read.csv(file.path(
    sensitivity_dir, "external_validation_leave_one_out_auc.csv"
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
  composite$dataset_id <- factor(
    composite$dataset_id,
    levels = levels(fixed$dataset_id)
  )
  permutations$dataset_id <- factor(
    permutations$dataset_id,
    levels = levels(fixed$dataset_id)
  )
  leave_one_out$dataset_id <- factor(
    leave_one_out$dataset_id,
    levels = levels(fixed$dataset_id)
  )
  safe_write_csv(
    fixed,
    file.path(source_dir, "Figure4_direction_fixed_external_AUC.csv")
  )
  safe_write_csv(
    composite,
    file.path(source_dir, "Figure4_signed_composite_AUC.csv")
  )
  safe_write_csv(
    permutations,
    file.path(source_dir, "Figure4_permutation_AUC.csv")
  )
  safe_write_csv(
    leave_one_out,
    file.path(source_dir, "Figure4_leave_one_out_AUC.csv")
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
      title = "Direction-fixed gene-level reproducibility",
      subtitle = "* 95% CI excludes AUC = 0.5",
      x = NULL,
      y = NULL,
      fill = "AUC"
    ) +
    submission_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  p2 <- ggplot2::ggplot(
    composite,
    ggplot2::aes(x = auc, y = dataset_id, colour = disease)
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
      title = "Fixed signed ten-gene score",
      subtitle = "Molecular separation only",
      x = "AUC (95% CI)",
      y = NULL,
      colour = NULL
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "top")
  p3 <- ggplot2::ggplot(
    permutations,
    ggplot2::aes(x = auc, fill = disease)
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      alpha = 0.75,
      position = "identity",
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
      subtitle = "Vertical lines are observed AUCs; 1,000 permutations/cohort",
      x = "Permuted-label AUC",
      y = "Frequency",
      fill = NULL,
      colour = NULL
    ) +
    submission_theme(7.5) +
    ggplot2::theme(legend.position = "none")
  p4 <- ggplot2::ggplot(
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
    ggplot2::geom_jitter(height = 0.12, width = 0, alpha = 0.65, size = 1.5) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 18,
      size = 3.2
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
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure4_external_validation",
    figure_dir,
    height_mm = 155
  )
}

submission_build_supplementary4 <- function(
    project_root,
    figure_dir,
    source_dir
) {
  context <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "sensitivity",
    "hpa_normal_tissue_context.csv"
  ))
  context$rna_tissue_specificity <- factor(
    context$rna_tissue_specificity,
    levels = c(
      "Low tissue specificity",
      "Tissue enhanced",
      "Group enriched",
      "Tissue enriched"
    )
  )
  context$rna_cell_type_specificity <- factor(
    context$rna_cell_type_specificity,
    levels = c(
      "Low cell type specificity",
      "Cell type enhanced",
      "Group enriched",
      "Cell type enriched"
    )
  )
  context$gene <- factor(
    context$gene,
    levels = rev(submission_load_cache(
      project_root, "07_machine_learning.rds"
    )$final_genes)
  )
  safe_write_csv(
    context,
    file.path(source_dir, "SupplementaryFigure4_HPA_context.csv")
  )
  p1 <- ggplot2::ggplot(
    context,
    ggplot2::aes(
      x = rna_tissue_specificity,
      y = gene,
      fill = ovary_listed_as_specific
    )
  ) +
    ggplot2::geom_point(shape = 21, size = 3, colour = "#374151") +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = submission_palette[["OC"]], `FALSE` = "#D1D5DB"),
      labels = c(`TRUE` = "Ovary listed", `FALSE` = "Ovary not listed")
    ) +
    ggplot2::labs(
      title = "Normal-tissue RNA specificity",
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    submission_theme() +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )
  p2 <- ggplot2::ggplot(
    context,
    ggplot2::aes(
      x = rna_cell_type_specificity,
      y = gene,
      colour = rna_tissue_distribution
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = c(
        "Detected in all" = submission_palette[["positive"]],
        "Detected in many" = submission_palette[["OC"]],
        "Detected in some" = "#7570B3"
      ),
      labels = c(
        "Detected in all" = "All",
        "Detected in many" = "Many",
        "Detected in some" = "Some"
      )
    ) +
    ggplot2::labs(
      title = "Normal single-cell type specificity",
      subtitle = "HPA categories describe normal reference data",
      x = NULL,
      y = NULL,
      colour = "Detected in tissues"
    ) +
    submission_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )
  note <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = -0.95,
      y = 0,
      label = paste0(
        "Interpretation boundary: HPA 25.1 integrates normal HPA/GTEx data.\n",
        "Cartilage is absent; only DIRAS3 listed ovary among its ",
        "group-enriched tissues.\n",
        "These annotations audit background expression and do not ",
        "establish OA–OC disease specificity."
      ),
      hjust = 0,
      size = 2.8,
      lineheight = 1.15,
      colour = "#374151"
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-0.5, 0.5)) +
    ggplot2::theme_void(base_family = "Arial")
  figure <- submission_panel_tag(
    ((p1 | p2) / note) +
      patchwork::plot_layout(heights = c(4, 1))
  )
  submission_save_plot(
    figure,
    "SupplementaryFigure4_HPA_normal_tissue_context",
    figure_dir,
    height_mm = 150
  )
}

submission_write_figure_legends <- function(figure_dir) {
  legends <- c(
    "# Figure legends",
    "",
    "## Figure 1. Study design, biological hypothesis, and audited data resources",
    "",
    "**A,** The study tests whether aging-associated, mesenchymal-remodeling, and immune-remodeling programs partially converge between osteoarthritis (OA) and ovarian cancer (OC), while keeping disease-specific analyses separate. Bidirectional Mendelian randomization (MR) is retained as supplementary null evidence. **B,** Bulk discovery and external-cohort sample counts. **C,** Total and quality-control (QC)-pass single-cell counts.",
    "",
    "## Figure 2. Direction-aware cross-disease transcriptomic discovery",
    "",
    "**A–B,** Differential-expression volcano plots for OA (GSE114007) and OC (GSE18520). **C,** Disease-specific log2 fold changes for commonly measured genes; purple points are the 286 shared differentially expressed genes. **D,** Shared-gene counts under six prespecified thresholds. Only 146/286 primary shared genes changed concordantly.",
    "",
    "## Figure 3. Robust identification of molecular candidates",
    "",
    "**A,** Absolute association between the most disease-associated weighted gene co-expression network analysis (WGCNA) module and phenotype under soft-power perturbation. **B,** Primary-module gene retention. **C,** Selection frequencies from 50 repeated five-fold outer resamples in which univariate screening was repeated using all measured genes inside each outer training fold; LASSO lambda and random-forest mtry were tuned without access to the outer test fold. **D,** Disease-specific discovery effects for the ten prioritized genes.",
    "",
    "## Figure 4. Direction-fixed cross-cohort molecular reproducibility",
    "",
    "**A,** Gene-level AUCs with expression direction fixed from the corresponding discovery cohort. **B,** AUCs and DeLong 95% intervals for the fixed signed ten-gene score; these summarize molecular separation, not clinical diagnostic accuracy. **C,** Null AUC distributions from 1,000 label permutations per cohort; vertical lines indicate observed values. **D,** AUCs after sequentially omitting each sample and recomputing the score; diamonds mark medians.",
    "",
    "## Figure 5. Candidate-gene localization across distinct cellular contexts",
    "",
    "**A,** QC-pass fractions after count-level audit and supported doublet detection. **B,** Cell-count-weighted detection fractions. **C,** Number of annotated cell types with at least 25% detection. **D,** Candidate genes meeting FDR below 0.05 in eligible pseudobulk contrasts. OA and OC cells were not forced into a common latent space; the panel shows that bulk overlap can arise in different cellular contexts.",
    "",
    "## Figure 6. Functional and exploratory prognostic context",
    "",
    "**A,** Hallmark normalized enrichment scores. **B,** Disease-minus-reference differences in rank-based immune-signature scores. **C,** Exploratory TCGA-OV Cox estimates. **D,** Original-cohort C-index distribution from 200 bootstrap-fitted models; the vertical line is the optimism-corrected estimate.",
    "",
    "## Supplementary Figure 1. Prespecified non-MR sensitivity analyses",
    "",
    "**A,** Candidate retention across DEG thresholds. **B,** Leave-one-sample-out WGCNA correlations. **C,** Most frequently selected genes under strict nested resampling. **D,** TCGA LASSO selection stability.",
    "",
    "## Supplementary Figure 2. Negative bidirectional MR results",
    "",
    "**A,** Estimates from five MR methods in both directions. **B,** Heterogeneity and Egger-intercept diagnostics. Null estimates mean no causal effect was detected under the available instruments and assumptions, not proof of absence.",
    "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings",
    "",
    "Released or recomputed dataset-specific embeddings with author-provided or conservative labels. All quantitative summaries used complete QC-pass cells. OA and OC datasets were not integrated into a shared latent space.",
    "",
    "## Supplementary Figure 4. Human Protein Atlas normal-tissue context",
    "",
    "**A,** RNA tissue-specificity categories for the ten candidates; fill indicates whether ovary was listed among the specific tissues. **B,** RNA single-cell-type specificity and tissue-distribution categories. Human Protein Atlas version 25.1 integrates normal HPA and GTEx expression. Cartilage is absent; therefore, these results audit background expression but cannot establish OA–OC disease specificity."
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
  log_info("Building revised main figures.")
  submission_build_figure1(project_root, figure_dir, source_dir)
  submission_build_figure2(project_root, figure_dir, source_dir)
  submission_build_figure3(project_root, figure_dir, source_dir)
  submission_build_figure4(project_root, figure_dir, source_dir)
  submission_build_figure5(project_root, figure_dir, source_dir)
  submission_build_figure6(project_root, figure_dir, source_dir)
  log_info("Building revised supplementary figures.")
  submission_build_supplementary1(project_root, figure_dir, source_dir)
  submission_build_supplementary2(project_root, figure_dir, source_dir)
  submission_build_supplementary3(project_root, figure_dir, source_dir)
  submission_build_supplementary4(project_root, figure_dir, source_dir)
  submission_write_figure_legends(figure_dir)
  safe_write_csv(
    data.frame(
      setting = c(
        "target_journal", "width_mm", "png_dpi",
        "font_family", "palette", "vector_format"
      ),
      value = c(
        "not yet specified; provisional general biomedical layout",
        "180", "300", "Arial",
        "Okabe-Ito-derived disease colours plus shapes/text", "PDF"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(figure_dir, "figure_style_manifest.csv")
  )
  log_info("Revised submission figures completed.")
  invisible(figure_dir)
}
