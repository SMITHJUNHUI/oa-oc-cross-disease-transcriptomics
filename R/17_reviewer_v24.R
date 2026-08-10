v24_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v24"))
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

v24_build_figure1 <- function(project_root, paths) {
  stages <- data.frame(
    stage = 1:9,
    x = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
    y = c(3, 3, 3, 2, 2, 2, 1, 1, 1),
    label = c(
      "Disease-specific\nbulk discovery",
      "Direction-aware\ngene overlap",
      "Paired pathway\ndirection",
      "WGCNA and\ncandidate evidence",
      "Strict nested\nfeature stability",
      "Cross-cohort\nmolecular separability",
      "Single-cell and\ngene-cell-function\ncontext",
      "HPA and TCGA\ncontext audits",
      "Bidirectional MR and\nbounded synthesis"
    ),
    domain = c(
      "Discovery", "Discovery", "Discovery",
      "Prioritization", "Prioritization", "Reproducibility",
      "Cell context", "Context", "Causality boundary"
    ),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(1.35, 2.35, 3, 1.35, 2.35, 3, 1.35, 2.35),
    xend = c(1.65, 2.65, 1, 1.65, 2.65, 1, 1.65, 2.65),
    y = c(3, 3, 2.82, 2, 2, 1.82, 1, 1),
    yend = c(3, 3, 2.18, 2, 2, 1.18, 1, 1),
    stringsAsFactors = FALSE
  )
  domain_colours <- c(
    Discovery = "#DCEAF7",
    Prioritization = "#E7E2F3",
    Reproducibility = "#FBE4D5",
    `Cell context` = "#DDF0E3",
    Context = "#F4EDC9",
    `Causality boundary` = "#E5E7EB"
  )
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"),
      linewidth = 0.55,
      colour = "#667085"
    ) +
    ggplot2::geom_tile(
      data = stages,
      ggplot2::aes(x = x, y = y, fill = domain),
      width = 0.70,
      height = 0.60,
      colour = "#475467",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = stages,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.65,
      lineheight = 0.95,
      colour = "#1F2937"
    ) +
    ggplot2::geom_label(
      data = stages,
      ggplot2::aes(x = x - 0.29, y = y + 0.23, label = stage),
      size = 2.5,
      label.size = 0,
      label.padding = grid::unit(0.9, "mm"),
      fill = "#344054",
      colour = "white",
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = domain_colours) +
    ggplot2::coord_cartesian(xlim = c(0.55, 3.45), ylim = c(0.58, 3.42)) +
    ggplot2::labs(
      title = "Linear evidence chain",
      subtitle = paste0(
        "Each layer answers a distinct question; no layer is treated as ",
        "proof of a shared mechanism."
      )
    ) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 10,
        colour = "#111827"
      ),
      plot.subtitle = ggplot2::element_text(size = 7.5, colour = "#475467"),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 2, 4)
    )

  resources <- data.frame(
    resource = factor(
      c(
        "Bulk discovery",
        "External bulk",
        "Single-cell",
        "HPA / TCGA-OV",
        "GWAS"
      ),
      levels = rev(c(
        "Bulk discovery",
        "External bulk",
        "Single-cell",
        "HPA / TCGA-OV",
        "GWAS"
      ))
    ),
    scale = c(
      "2 cohorts\n101 samples",
      "4 cohorts\n112 samples",
      "5 datasets\n1,025,361 QC-pass cells",
      "Normal reference\n307 TCGA-OV samples",
      "2 datasets\n21 and 11 instruments"
    ),
    role = c(
      "Disease-specific discovery",
      "Secondary molecular separability",
      "Cellular localization",
      "Tissue/composition\nand survival context",
      "Inherited-causality assessment"
    ),
    boundary = c(
      "Association",
      "Not diagnostic validation",
      "Not a shared cell state",
      "Not tissue specificity\nor prognosis validation",
      "Null under available instruments"
    ),
    stringsAsFactors = FALSE
  )
  resources$y <- as.numeric(resources$resource)
  safe_write_csv(stages, file.path(paths$source, "Figure1_linear_workflow.csv"))
  safe_write_csv(resources, file.path(paths$source, "Figure1_resource_roles.csv"))

  p2 <- ggplot2::ggplot(resources, ggplot2::aes(y = resource)) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 1),
      width = 0.97,
      height = 0.78,
      fill = "#F8FAFC",
      colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 2),
      width = 0.97,
      height = 0.78,
      fill = "#F9FAFB",
      colour = "#D0D5DD"
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(x = 3),
      width = 0.97,
      height = 0.78,
      fill = "#FFF9F2",
      colour = "#D0D5DD"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 1, label = scale),
      size = 2.35,
      lineheight = 0.95,
      colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2, label = role),
      size = 2.35,
      lineheight = 0.95,
      colour = "#1F2937"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = 3, label = boundary),
      size = 2.3,
      lineheight = 0.95,
      colour = "#9A3412"
    ) +
    ggplot2::annotate(
      "text",
      x = c(1, 2, 3),
      y = 5.45,
      label = c("Scale", "Analytic role", "Inference boundary"),
      fontface = "bold",
      size = 3.1,
      colour = "#344054"
    ) +
    ggplot2::scale_x_continuous(limits = c(0.48, 3.52), breaks = NULL) +
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(0.25, 0.95))
    ) +
    ggplot2::labs(
      title = "Audited data resources and intended roles",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_family = "Arial", base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = "bold",
        colour = "#344054",
        size = 7.2
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 10,
        colour = "#111827"
      ),
      plot.margin = ggplot2::margin(2, 4, 4, 4)
    )

  figure <- submission_panel_tag(
    p1 / p2 + patchwork::plot_layout(heights = c(1.5, 1))
  )
  submission_save_plot(
    figure,
    "Figure1_study_design",
    paths$figures,
    height_mm = 188
  )
}

v24_candidate_selection <- function(paths) {
  candidates <- utils::read.csv(file.path(
    paths$tables,
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  selected <- candidates[candidates$gene %in% c(
    "SOX9", "DDIT3", "BNC1", "AKAP12"
  ), , drop = FALSE]
  reasons <- c(
    SOX9 = paste0(
      "Only cross-disease model consensus; discordant bulk direction and ",
      "distinct OA/OC single-cell contexts."
    ),
    DDIT3 = paste0(
      "Highest strict nested OA selection frequency among the ten candidates; ",
      "discordant bulk direction."
    ),
    BNC1 = paste0(
      "Highest strict nested OC selection frequency among the ten candidates; ",
      "OC-skewed model evidence and stromal context association."
    ),
    AKAP12 = paste0(
      "Concordant bulk direction with non-zero OA nested stability and ",
      "explicitly different OA/OC gene-cell-function contexts."
    )
  )
  selected$exploratory_selection_rationale <- unname(reasons[selected$gene])
  selected$selection_scope <- paste0(
    "Transparent post hoc exemplars for candidate-centered pathway context; ",
    "not a revised signature or new candidate-ranking algorithm."
  )
  selected <- selected[match(
    c("SOX9", "DDIT3", "BNC1", "AKAP12"),
    selected$gene
  ), , drop = FALSE]
  safe_write_csv(
    selected,
    file.path(paths$analysis, "candidate_centered_selection.csv")
  )
  selected
}

v24_partial_association_rank <- function(dataset, anchor_gene) {
  expression <- dataset$expression
  if (!anchor_gene %in% rownames(expression)) {
    stop("Candidate not present in ", dataset$id, ": ", anchor_gene, call. = FALSE)
  }
  group <- stats::relevel(factor(dataset$group), ref = "Normal")
  design <- stats::model.matrix(~group)
  residuals <- qr.resid(qr(design), t(expression))
  anchor <- residuals[, match(anchor_gene, rownames(expression))]
  numerator <- as.numeric(crossprod(anchor, residuals))
  denominator <- sqrt(sum(anchor^2) * colSums(residuals^2))
  association <- numerator / denominator
  association[!is.finite(association)] <- 0
  association <- pmax(pmin(association, 0.999999), -0.999999)
  degrees_freedom <- nrow(residuals) - qr(design)$rank - 1
  statistic <- association * sqrt(
    degrees_freedom / pmax(1 - association^2, .Machine$double.eps)
  )
  names(statistic) <- rownames(expression)
  statistic <- statistic[
    is.finite(statistic) & names(statistic) != anchor_gene
  ]
  sort(statistic, decreasing = TRUE)
}

v24_run_candidate_gsea <- function(
  ranked_statistic,
  gene_sets,
  seed
) {
  require_namespace("BiocParallel", "deterministic serial enrichment")
  BiocParallel::register(BiocParallel::SerialParam())
  term_to_gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unname(unlist(gene_sets, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  result <- suppressMessages(clusterProfiler::GSEA(
    geneList = ranked_statistic,
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

v24_candidate_centered_pathways <- function(project_root, paths, selected) {
  bulk <- submission_load_cache(project_root, "02_bulk_training.rds")
  enrichment <- submission_load_cache(project_root, "05_enrichment.rds")
  gene_sets <- enrichment$gsea$oa_train$hallmark@geneSets
  datasets <- list(OA = bulk$oa_train, OC = bulk$oc_train)
  output <- list()
  ranked_output <- list()
  index <- 1L
  ranked_index <- 1L
  for (candidate_index in seq_len(nrow(selected))) {
    candidate <- selected$gene[[candidate_index]]
    for (disease in names(datasets)) {
      dataset <- datasets[[disease]]
      ranking <- v24_partial_association_rank(dataset, candidate)
      ranked_output[[ranked_index]] <- data.frame(
        candidate = candidate,
        disease = disease,
        dataset_id = dataset$id,
        gene = names(ranking),
        partial_association_t = as.numeric(ranking),
        stringsAsFactors = FALSE
      )
      ranked_index <- ranked_index + 1L
      result <- v24_run_candidate_gsea(
        ranking,
        gene_sets,
        20260730L + candidate_index * 10L + match(disease, names(datasets))
      )
      result <- result[, c(
        "ID", "Description", "setSize", "NES", "pvalue", "p.adjust",
        "qvalue", "rank", "leading_edge", "core_enrichment"
      )]
      names(result) <- c(
        "pathway_id", "pathway", "set_size", "NES", "P", "FDR",
        "q_value", "peak_rank", "leading_edge_summary", "leading_edge_genes"
      )
      complete_pathways <- data.frame(
        pathway_id = names(gene_sets),
        stringsAsFactors = FALSE
      )
      result <- merge(
        complete_pathways,
        result,
        by = "pathway_id",
        all.x = TRUE,
        sort = FALSE
      )
      result$calculation_status <- ifelse(
        is.finite(result$NES) & is.finite(result$FDR),
        "estimated",
        paste0(
          "not estimated by fgsea multilevel because the residual ranking ",
          "was insufficiently balanced for this gene set"
        )
      )
      result$candidate <- candidate
      result$disease <- disease
      result$dataset_id <- dataset$id
      result$adjustment <- "Residual association after adjustment for disease/reference group"
      result$association_direction <- ifelse(
        result$NES > 0,
        "positive residual association",
        "negative residual association"
      )
      result$exploratory_selection_rationale <-
        selected$exploratory_selection_rationale[[candidate_index]]
      result$inference_boundary <- paste0(
        "Exploratory candidate-centered association enrichment; not ",
        "single-gene perturbation, pathway activation, mediation, or mechanism."
      )
      output[[index]] <- result
      index <- index + 1L
    }
  }
  combined <- do.call(rbind, output)
  combined <- combined[, c(
    "candidate", "disease", "dataset_id", "pathway_id", "pathway",
    "set_size", "NES", "P", "FDR", "q_value", "peak_rank",
    "association_direction", "calculation_status", "adjustment",
    "exploratory_selection_rationale", "leading_edge_summary",
    "leading_edge_genes", "inference_boundary"
  )]
  combined$pathway <- v23_title_case_pathway(combined$pathway_id)
  combined <- combined[order(
    match(combined$candidate, selected$gene),
    combined$disease,
    combined$FDR,
    -abs(combined$NES)
  ), , drop = FALSE]
  safe_write_csv(
    combined,
    file.path(
      paths$tables,
      "Table_S20_candidate_centered_Hallmark_context.csv"
    )
  )
  ranked <- do.call(rbind, ranked_output)
  safe_write_csv(
    ranked,
    file.path(paths$analysis, "candidate_centered_partial_association_ranks.csv")
  )

  plot_data <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
    candidate <- selected$gene[[i]]
    part <- combined[combined$candidate == candidate, , drop = FALSE]
    part <- part[
      is.finite(part$FDR) & is.finite(part$NES),
      ,
      drop = FALSE
    ]
    summary_fdr <- stats::aggregate(
      FDR ~ pathway_id,
      data = part,
      FUN = min
    )
    part$abs_NES <- abs(part$NES)
    summary_abs <- stats::aggregate(
      abs_NES ~ pathway_id,
      data = part,
      FUN = max
    )
    summary <- merge(summary_fdr, summary_abs, by = "pathway_id")
    summary <- summary[order(summary$FDR, -summary$abs_NES), , drop = FALSE]
    keep <- head(summary$pathway_id, 6L)
    part[part$pathway_id %in% keep, , drop = FALSE]
  }))
  plot_data$significant <- ifelse(
    plot_data$FDR < 0.05,
    "FDR < 0.05",
    "FDR >= 0.05"
  )
  plot_data$pathway_label <- factor(
    plot_data$pathway,
    levels = rev(unique(plot_data$pathway))
  )
  safe_write_csv(
    plot_data,
    file.path(
      paths$source,
      "SupplementaryFigure9_candidate_centered_Hallmark_context.csv"
    )
  )
  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = NES, y = pathway_label, colour = disease)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "#98A2B3"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = pathway_id),
      colour = "#C7CDD4",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = significant),
      size = 2.3,
      stroke = 0.7
    ) +
    ggplot2::facet_wrap(
      ~candidate,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(
      "FDR < 0.05" = 16,
      "FDR >= 0.05" = 1
    )) +
    ggplot2::labs(
      title = "Exploratory candidate-centered Hallmark contexts",
      subtitle = paste0(
        "Residual gene association after disease/reference adjustment; ",
        "six strongest pathways per candidate"
      ),
      x = "Normalized enrichment score",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.1) +
    ggplot2::theme(
      legend.position = "top",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 5.9)
    )
  submission_save_plot(
    figure,
    "SupplementaryFigure9_candidate_centered_Hallmark_context",
    paths$figures,
    height_mm = 175
  )
  combined
}

v24_cross_cohort_context_table <- function(paths) {
  validation <- utils::read.csv(file.path(
    paths$tables,
    "Table_S6_external_validation.csv"
  ))
  summary <- validation[
    validation$source_table == "signed_score",
    ,
    drop = FALSE
  ]
  manifest <- utils::read.csv(file.path(
    paths$tables,
    "Table_S1_data_sources_and_cohorts.csv"
  ))
  manifest <- manifest[
    match(summary$dataset_id, manifest$source_id),
    ,
    drop = FALSE
  ]
  contexts <- data.frame(
    dataset_id = summary$dataset_id,
    disease = summary$disease,
    tissue_and_comparator = c(
      "Human cartilage; OA versus reference cartilage",
      "Human synovium; OA versus reference synovium",
      "Ovarian/peritoneal tissue; tumor versus reference tissue",
      "Ovarian/peritoneal tissue; serous carcinoma versus normal peritoneum"
    ),
    validation_task_scale = c(
      "Within-tissue chronic degenerative contrast",
      "Within-tissue chronic inflammatory/degenerative contrast",
      "Malignant transformation and tumor-microenvironment contrast",
      "Malignant tissue versus non-malignant comparator-tissue contrast"
    ),
    n_reference = summary$n_normal,
    n_disease = summary$n_disease,
    fixed_direction_AUC = summary$auc,
    CI_lower = summary$ci_lower,
    CI_upper = summary$ci_upper,
    permutation_empirical_P = summary$permutation_empirical_p,
    leave_one_out_AUC_min = summary$leave_one_out_auc_minimum,
    leave_one_out_AUC_median = summary$leave_one_out_auc_median,
    leave_one_out_AUC_max = summary$leave_one_out_auc_maximum,
    source_design = manifest$source_design,
    task_interpretation = c(
      "Modest molecular separability; interval includes chance",
      "Modest molecular separability; interval includes chance",
      "Near-complete retrospective molecular separability",
      "Near-complete retrospective molecular separability"
    ),
    inference_boundary = paste0(
      "The same fixed shared molecular summary was evaluated in biologically ",
      "different validation tasks. AUC describes cohort-specific separability, ",
      "not clinical utility, disease severity, or a shared mechanism."
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    contexts,
    file.path(
      paths$tables,
      "Table_S21_cross_cohort_molecular_separability_context.csv"
    )
  )
  safe_write_csv(
    contexts,
    file.path(paths$source, "Figure4_validation_task_context.csv")
  )
  contexts
}

v24_hedges_g <- function(score, status) {
  disease <- score[status == 1L]
  reference <- score[status == 0L]
  n1 <- length(disease)
  n0 <- length(reference)
  pooled_sd <- sqrt(
    ((n1 - 1) * stats::var(disease) + (n0 - 1) * stats::var(reference)) /
      (n1 + n0 - 2)
  )
  cohen_d <- (mean(disease) - mean(reference)) / pooled_sd
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  hedges_g <- correction * cohen_d
  variance <- (n1 + n0) / (n1 * n0) +
    hedges_g^2 / (2 * (n1 + n0 - 2))
  c(
    hedges_g = hedges_g,
    standard_error = sqrt(variance),
    variance = variance
  )
}

v24_external_effect_sizes <- function(paths) {
  validation <- utils::read.csv(file.path(
    paths$tables,
    "Table_S6_external_validation.csv"
  ))
  scores <- validation[
    validation$source_table == "sample_scores",
    c("dataset_id", "disease", "sample_id", "group", "signed_score"),
    drop = FALSE
  ]
  datasets <- unique(scores$dataset_id)
  estimates <- do.call(rbind, lapply(datasets, function(dataset_id) {
    part <- scores[scores$dataset_id == dataset_id, , drop = FALSE]
    status <- as.integer(part$group == "Disease")
    estimate <- v24_hedges_g(part$signed_score, status)
    data.frame(
      disease = part$disease[[1L]],
      dataset_id = dataset_id,
      estimate_type = "cohort",
      Hedges_g = estimate[["hedges_g"]],
      standard_error = estimate[["standard_error"]],
      variance = estimate[["variance"]],
      CI_lower = estimate[["hedges_g"]] - 1.96 * estimate[["standard_error"]],
      CI_upper = estimate[["hedges_g"]] + 1.96 * estimate[["standard_error"]],
      n_reference = sum(status == 0L),
      n_disease = sum(status == 1L),
      tau_squared = NA_real_,
      I_squared = NA_real_,
      interpretation = paste0(
        "Direction-fixed signed-score standardized mean difference; ",
        "positive values indicate higher disease-group score."
      ),
      inference_boundary = paste0(
        "Cohort-level molecular contrast; not diagnostic accuracy or a ",
        "shared biological effect size across different tissues."
      ),
      stringsAsFactors = FALSE
    )
  }))
  pooled <- do.call(rbind, lapply(c("OA", "OC"), function(disease) {
    part <- estimates[estimates$disease == disease, , drop = FALSE]
    weight_fixed <- 1 / part$variance
    pooled_fixed <- sum(weight_fixed * part$Hedges_g) / sum(weight_fixed)
    Q <- sum(weight_fixed * (part$Hedges_g - pooled_fixed)^2)
    k <- nrow(part)
    denominator <- sum(weight_fixed) -
      sum(weight_fixed^2) / sum(weight_fixed)
    tau_squared <- max(0, (Q - (k - 1)) / denominator)
    weight_random <- 1 / (part$variance + tau_squared)
    estimate <- sum(weight_random * part$Hedges_g) / sum(weight_random)
    standard_error <- sqrt(1 / sum(weight_random))
    I_squared <- if (Q > 0) max(0, 100 * (Q - (k - 1)) / Q) else 0
    data.frame(
      disease = disease,
      dataset_id = paste0(disease, " random-effects summary"),
      estimate_type = "random-effects summary",
      Hedges_g = estimate,
      standard_error = standard_error,
      variance = standard_error^2,
      CI_lower = estimate - 1.96 * standard_error,
      CI_upper = estimate + 1.96 * standard_error,
      n_reference = sum(part$n_reference),
      n_disease = sum(part$n_disease),
      tau_squared = tau_squared,
      I_squared = I_squared,
      interpretation = paste0(
        "Illustrative within-disease random-effects summary across two ",
        "heterogeneous validation cohorts."
      ),
      inference_boundary = paste0(
        "With only two cohorts, heterogeneity and pooled estimates are ",
        "imprecise and remain descriptive."
      ),
      stringsAsFactors = FALSE
    )
  }))
  output <- rbind(estimates, pooled)
  safe_write_csv(
    output,
    file.path(
      paths$tables,
      "Table_S22a_external_signed_score_effect_sizes.csv"
    )
  )
  output
}

v24_ridge_logit_predict <- function(x_train, y_train, x_test, penalty = 0.1) {
  center <- mean(x_train)
  scale_value <- stats::sd(x_train)
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }
  x_scaled <- (x_train - center) / scale_value
  objective <- function(parameters) {
    eta <- parameters[[1L]] + parameters[[2L]] * x_scaled
    log_partition <- ifelse(
      eta > 0,
      eta + log1p(exp(-eta)),
      log1p(exp(eta))
    )
    sum(log_partition - y_train * eta) +
      0.5 * penalty * parameters[[2L]]^2
  }
  fit <- stats::optim(
    par = c(stats::qlogis(mean(y_train)), 0),
    fn = objective,
    method = "BFGS"
  )
  eta_test <- fit$par[[1L]] +
    fit$par[[2L]] * ((x_test - center) / scale_value)
  stats::plogis(eta_test)
}

v24_cross_fitted_calibration <- function(paths) {
  validation <- utils::read.csv(file.path(
    paths$tables,
    "Table_S6_external_validation.csv"
  ))
  scores <- validation[
    validation$source_table == "sample_scores",
    c("dataset_id", "disease", "sample_id", "group", "signed_score"),
    drop = FALSE
  ]
  predictions <- do.call(rbind, lapply(unique(scores$dataset_id), function(id) {
    part <- scores[scores$dataset_id == id, , drop = FALSE]
    status <- as.integer(part$group == "Disease")
    probability <- rep(NA_real_, nrow(part))
    for (omitted in seq_len(nrow(part))) {
      retained <- setdiff(seq_len(nrow(part)), omitted)
      probability[[omitted]] <- v24_ridge_logit_predict(
        x_train = part$signed_score[retained],
        y_train = status[retained],
        x_test = part$signed_score[[omitted]],
        penalty = 0.1
      )
    }
    data.frame(
      dataset_id = id,
      disease = part$disease,
      sample_id = part$sample_id,
      group = part$group,
      observed = status,
      signed_score = part$signed_score,
      cross_fitted_probability = probability,
      model = paste0(
        "Leave-one-sample-out ridge logistic recalibration of the fixed ",
        "signed score; slope L2 penalty=0.1"
      ),
      stringsAsFactors = FALSE
    )
  }))
  metrics <- do.call(rbind, lapply(unique(predictions$dataset_id), function(id) {
    part <- predictions[predictions$dataset_id == id, , drop = FALSE]
    probability <- pmax(
      pmin(part$cross_fitted_probability, 1 - 1e-6),
      1e-6
    )
    logit <- stats::qlogis(probability)
    calibration <- suppressWarnings(stats::glm(
      part$observed ~ logit,
      family = stats::binomial()
    ))
    coefficients <- stats::coef(calibration)
    raw_intercept <- unname(coefficients[[1L]])
    raw_slope <- unname(coefficients[[2L]])
    logit_sd <- stats::sd(logit)
    calibration_status <- if (
      !is.finite(raw_slope) || logit_sd < 0.1 || abs(raw_slope) > 10
    ) {
      "numerically unstable; do not interpret intercept or slope"
    } else if (raw_slope <= 0) {
      "directionally unstable; do not interpret as calibration"
    } else {
      "estimable descriptive sensitivity"
    }
    interpretable <- identical(
      calibration_status,
      "estimable descriptive sensitivity"
    )
    data.frame(
      dataset_id = id,
      disease = part$disease[[1L]],
      n = nrow(part),
      Brier_score = mean((probability - part$observed)^2),
      calibration_intercept = if (interpretable) raw_intercept else NA_real_,
      calibration_slope = if (interpretable) raw_slope else NA_real_,
      raw_calibration_intercept = raw_intercept,
      raw_calibration_slope = raw_slope,
      cross_fitted_logit_SD = logit_sd,
      calibration_status = calibration_status,
      assessment = paste0(
        "Cross-fitted cohort-specific score recalibration; descriptive ",
        "sensitivity only."
      ),
      inference_boundary = paste0(
        "This is not external probability calibration because no locked ",
        "probability model was transported from discovery."
      ),
      stringsAsFactors = FALSE
    )
  }))
  predictions$calibration_bin <- ave(
    predictions$cross_fitted_probability,
    predictions$dataset_id,
    FUN = function(value) {
      pmin(3L, ceiling(rank(value, ties.method = "first") / (length(value) / 3)))
    }
  )
  binned <- stats::aggregate(
    cbind(
      mean_predicted = predictions$cross_fitted_probability,
      observed_fraction = predictions$observed
    ),
    by = list(
      dataset_id = predictions$dataset_id,
      disease = predictions$disease,
      calibration_bin = predictions$calibration_bin
    ),
    FUN = mean
  )
  binned$n_in_bin <- as.integer(table(
    interaction(
      predictions$dataset_id,
      predictions$disease,
      predictions$calibration_bin,
      drop = TRUE
    )
  )[interaction(
    binned$dataset_id,
    binned$disease,
    binned$calibration_bin,
    drop = TRUE
  )])
  safe_write_csv(
    metrics,
    file.path(
      paths$tables,
      "Table_S22b_cross_fitted_calibration_metrics.csv"
    )
  )
  safe_write_csv(
    predictions,
    file.path(paths$analysis, "cross_fitted_calibration_predictions.csv")
  )
  safe_write_csv(
    binned,
    file.path(paths$source, "SupplementaryFigure10_calibration_bins.csv")
  )
  list(metrics = metrics, predictions = predictions, binned = binned)
}

v24_build_external_evaluation_figure <- function(paths, effects, calibration) {
  effects$display <- effects$dataset_id
  effects$display <- factor(effects$display, levels = rev(effects$display))
  p1 <- ggplot2::ggplot(
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
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_lower, xmax = CI_upper),
      height = 0.16,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_shape_manual(values = c(
      cohort = 16,
      `random-effects summary` = 18
    )) +
    ggplot2::labs(
      title = "Direction-fixed signed-score effect sizes",
      subtitle = "Diamonds are descriptive two-cohort random-effects summaries",
      x = "Hedges g (disease minus reference)",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    submission_theme(7.4) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 6.4)
    )
  p2 <- ggplot2::ggplot(
    calibration$binned,
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
      subtitle = paste0(
        "LOO ridge recalibration of the fixed score; not transport of a ",
        "locked probability model"
      ),
      x = "Mean cross-fitted probability",
      y = "Observed disease fraction",
      colour = NULL,
      size = "Bin n"
    ) +
    submission_theme(7.1) +
    ggplot2::theme(legend.position = "top")
  safe_write_csv(
    effects,
    file.path(paths$source, "SupplementaryFigure10_effect_sizes.csv")
  )
  figure <- submission_panel_tag(p1 / p2)
  submission_save_plot(
    figure,
    "SupplementaryFigure10_external_evaluation_context",
    paths$figures,
    height_mm = 180
  )
}

v24_upstream_regulatory_context <- function(project_root, paths, selected) {
  gmt_path <- file.path(
    project_root,
    "data",
    "raw",
    "regulatory",
    "human_500diff.gmt"
  )
  if (!file.exists(gmt_path)) {
    stop(
      "Project-local KnockTF GMT input is missing: ",
      gmt_path,
      call. = FALSE
    )
  }
  gmt <- clusterProfiler::read.gmt(gmt_path)
  names(gmt) <- c("tf", "gene")
  gmt <- unique(gmt)
  enrichment <- submission_load_cache(project_root, "05_enrichment.rds")
  rankings <- list(
    OA = enrichment$gsea$oa_train$hallmark@geneList,
    OC = enrichment$gsea$oc_train$hallmark@geneList
  )
  results <- lapply(seq_along(rankings), function(index) {
    disease <- names(rankings)[[index]]
    result <- v24_run_candidate_gsea(
      rankings[[index]],
      split(gmt$gene, gmt$tf),
      20260820L + index
    )
    result$disease <- disease
    result
  })
  results <- do.call(rbind, results)
  results <- results[, c("ID", "NES", "pvalue", "p.adjust", "disease")]
  names(results) <- c("tf", "NES", "P", "FDR", "disease")
  links <- utils::read.csv(file.path(
    paths$tables,
    "Table_S13_supplementary_KnockTF_GMT_hub_interactions.csv"
  ))
  links <- unique(links[
    links$target_normalized %in% selected$gene,
    c("tf", "target_normalized"),
    drop = FALSE
  ])
  linked <- stats::aggregate(
    links$target_normalized,
    by = list(tf = links$tf),
    FUN = function(value) paste(sort(unique(value)), collapse = ";")
  )
  names(linked)[[2L]] <- "linked_candidates"
  linked$candidate_count <- lengths(strsplit(linked$linked_candidates, ";", fixed = TRUE))
  paired <- merge(
    results[results$disease == "OA", c("tf", "NES", "P", "FDR")],
    results[results$disease == "OC", c("tf", "NES", "P", "FDR")],
    by = "tf",
    all = TRUE,
    suffixes = c("_OA", "_OC")
  )
  paired <- merge(paired, linked, by = "tf", all = FALSE)
  paired$direction_class <- ifelse(
    sign(paired$NES_OA) == sign(paired$NES_OC),
    "concordant target-set direction",
    "discordant target-set direction"
  )
  paired$resource <- paste0(
    "Project-local KnockTF perturbational top-differential target sets; ",
    "data/raw/regulatory/human_500diff.gmt"
  )
  paired$inference_boundary <- paste0(
    "Target-set enrichment does not establish TF activity, direct regulation ",
    "in OA or OC, or causal explanation of candidate expression."
  )
  paired <- paired[order(
    pmin(paired$FDR_OA, paired$FDR_OC, na.rm = TRUE),
    -paired$candidate_count,
    -pmax(abs(paired$NES_OA), abs(paired$NES_OC), na.rm = TRUE)
  ), , drop = FALSE]
  safe_write_csv(
    paired,
    file.path(
      paths$tables,
      "Table_S23a_KnockTF_candidate_regulatory_context.csv"
    )
  )

  mirna <- utils::read.csv(file.path(
    paths$tables,
    "Table_S13_supplementary_miRTarBase_hub_interactions.csv"
  ))
  mirna <- unique(mirna[
    mirna$target_normalized %in% selected$gene,
    c("miRNA", "target_normalized"),
    drop = FALSE
  ])
  coverage <- stats::aggregate(
    mirna$target_normalized,
    by = list(miRNA = mirna$miRNA),
    FUN = function(value) paste(sort(unique(value)), collapse = ";")
  )
  names(coverage)[[2L]] <- "linked_candidates"
  coverage$candidate_count <- lengths(strsplit(
    coverage$linked_candidates,
    ";",
    fixed = TRUE
  ))
  coverage$resource <- "miRTarBase locally released interaction catalogue"
  coverage$disease_specific_evidence <- "not available in this study"
  coverage$inference_boundary <- paste0(
    "Interaction presence is hypothesis-generating; miRNA abundance, ",
    "direction, and activity were not measured in the analyzed cohorts."
  )
  coverage <- coverage[order(
    -coverage$candidate_count,
    coverage$miRNA
  ), , drop = FALSE]
  safe_write_csv(
    coverage,
    file.path(
      paths$tables,
      "Table_S23b_miRTarBase_candidate_regulatory_context.csv"
    )
  )

  top_tf <- head(paired$tf, 10L)
  tf_plot <- results[results$tf %in% top_tf, , drop = FALSE]
  tf_plot <- merge(
    tf_plot,
    paired[, c("tf", "linked_candidates", "candidate_count")],
    by = "tf",
    all.x = TRUE
  )
  tf_plot$significant <- tf_plot$FDR < 0.05
  tf_order <- paired$tf[paired$tf %in% top_tf]
  tf_plot$tf <- factor(tf_plot$tf, levels = rev(tf_order))
  p1 <- ggplot2::ggplot(
    tf_plot,
    ggplot2::aes(x = disease, y = tf, fill = NES)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(significant, "*", "")),
      size = 4,
      colour = "#111827"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    ggplot2::labs(
      title = "KnockTF perturbational target-set context",
      subtitle = "Asterisks denote FDR <0.05; NES is not TF activity",
      x = NULL,
      y = NULL,
      fill = "NES"
    ) +
    submission_theme(7.4) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.4))

  top_mirna <- head(
    coverage$miRNA[coverage$candidate_count >= 2L],
    12L
  )
  if (length(top_mirna) < 6L) {
    top_mirna <- head(coverage$miRNA, 12L)
  }
  mirna_plot <- mirna[mirna$miRNA %in% top_mirna, , drop = FALSE]
  mirna_plot$miRNA <- factor(mirna_plot$miRNA, levels = rev(top_mirna))
  mirna_plot$target_normalized <- factor(
    mirna_plot$target_normalized,
    levels = selected$gene
  )
  p2 <- ggplot2::ggplot(
    mirna_plot,
    ggplot2::aes(x = target_normalized, y = miRNA)
  ) +
    ggplot2::geom_tile(
      fill = "#DDEAF3",
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      shape = 21,
      fill = "#0072B2",
      colour = "white",
      size = 2.2,
      stroke = 0.4
    ) +
    ggplot2::labs(
      title = "Curated miRNA-candidate interaction coverage",
      subtitle = "Interaction presence only; activity was not measured",
      x = "Candidate gene",
      y = NULL
    ) +
    submission_theme(7.4) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.2))
  safe_write_csv(
    tf_plot,
    file.path(paths$source, "SupplementaryFigure11_KnockTF_context.csv")
  )
  safe_write_csv(
    mirna_plot,
    file.path(paths$source, "SupplementaryFigure11_miRNA_context.csv")
  )
  figure <- submission_panel_tag(p1 | p2)
  submission_save_plot(
    figure,
    "SupplementaryFigure11_upstream_regulatory_context",
    paths$figures,
    height_mm = 145
  )
  list(tf = paired, mirna = coverage)
}

v24_update_figure_materials <- function(paths) {
  legend_path <- file.path(paths$figures, "figure_legends.md")
  legends <- readLines(legend_path, warn = FALSE, encoding = "UTF-8")
  legends <- v23_replace_legend(
    legends,
    "Figure 1. Study design, biological hypothesis, and audited data resources",
    paste0(
      "**A,** Linear evidence chain from disease-specific discovery through ",
      "direction-aware overlap, pathway comparison, network and nested-model ",
      "stability, external molecular separability, single-cell localization, ",
      "context audits, and bidirectional MR. **B,** Audited resource scale, ",
      "analytic role, and inference boundary. Each layer answers a distinct ",
      "question and no layer is treated as proof of a shared mechanism."
    )
  )
  legends[legends == "## Figure 1. Study design, biological hypothesis, and audited data resources"] <-
    "## Figure 1. Linear study design and audited evidence boundaries"
  legends <- v23_replace_legend(
    legends,
    "Figure 4. Direction-fixed cross-cohort molecular reproducibility",
    paste0(
      "**A,** Unsupervised PCA of GSE54388. **B,** Null AUC distributions ",
      "from 1,000 label permutations. **C,** Leave-one-sample-out AUCs. ",
      "**D,** ROC curves for the fixed signed ten-gene molecular summary, ",
      "shown last as a secondary display. The same summary was evaluated in ",
      "validation tasks with different tissues, comparator groups, and ",
      "biological scales (Table S21); these panels show retrospective ",
      "cross-cohort molecular separability, not clinical performance."
    )
  )
  legends[legends == "## Figure 4. Direction-fixed cross-cohort molecular reproducibility"] <-
    "## Figure 4. Direction-fixed cross-cohort molecular separability"
  additions <- c(
    "## Supplementary Figure 9. Candidate-centered Hallmark contexts",
    "",
    paste0(
      "Disease-status-adjusted residual association rankings for SOX9, DDIT3, ",
      "BNC1, and AKAP12 were analyzed against Hallmark sets. Six pathways per ",
      "candidate are shown; filled points denote FDR <0.05. This transparent ",
      "post hoc analysis is exploratory and does not represent single-gene ",
      "perturbation, mediation, pathway activation, or mechanism."
    ),
    "",
    "## Supplementary Figure 10. External evaluation context",
    "",
    paste0(
      "**A,** Direction-fixed signed-score standardized mean differences by ",
      "cohort with descriptive two-cohort random-effects summaries by disease. ",
      "**B,** Three-bin calibration sensitivity from leave-one-sample-out ridge ",
      "recalibration of the fixed score. Because a locked probability model was ",
      "not transported from discovery, this is a cross-fitted sensitivity ",
      "analysis rather than external clinical calibration."
    ),
    "",
    "## Supplementary Figure 11. Focused upstream regulatory context",
    "",
    paste0(
      "**A,** OA and OC enrichment of KnockTF perturbational target sets for ",
      "transcription factors connected to the four candidate exemplars; NES ",
      "does not estimate TF activity. **B,** Curated miRTarBase interaction ",
      "coverage for multi-candidate miRNAs. miRNA abundance and activity were ",
      "not measured, so these matrices define testable regulatory context rather ",
      "than an inferred disease mechanism."
    )
  )
  if (!"## Supplementary Figure 11. Focused upstream regulatory context" %in%
      legends) {
    legends <- c(legends, "", additions, "")
  }
  writeLines(enc2utf8(legends), legend_path, useBytes = TRUE)

  style_path <- file.path(paths$figures, "figure_style_manifest.csv")
  style <- utils::read.csv(style_path, stringsAsFactors = FALSE)
  style$value[style$setting == "revision"] <- "V2.4"
  safe_write_csv(style, style_path)
}

v24_update_table_index <- function(paths) {
  index_path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  additions <- data.frame(
    table_id = c(
      "Table S20",
      "Table S21",
      "Table S22a",
      "Table S22b",
      "Table S23a",
      "Table S23b"
    ),
    filename = c(
      "Table_S20_candidate_centered_Hallmark_context.csv",
      "Table_S21_cross_cohort_molecular_separability_context.csv",
      "Table_S22a_external_signed_score_effect_sizes.csv",
      "Table_S22b_cross_fitted_calibration_metrics.csv",
      "Table_S23a_KnockTF_candidate_regulatory_context.csv",
      "Table_S23b_miRTarBase_candidate_regulatory_context.csv"
    ),
    title = c(
      "Exploratory candidate-centered Hallmark context",
      "Cross-cohort molecular separability context",
      "External signed-score standardized effect sizes",
      "Cross-fitted calibration sensitivity metrics",
      "KnockTF candidate regulatory context",
      "miRTarBase candidate regulatory context"
    ),
    contents = c(
      paste0(
        "Disease-status-adjusted candidate-gene residual association GSEA ",
        "for SOX9, DDIT3, BNC1, and AKAP12."
      ),
      paste0(
        "Cohort tissue/comparator scale, fixed AUC with confidence interval, ",
        "permutation P value, leave-one-out range, and inference boundary."
      ),
      paste0(
        "Hedges g by cohort and descriptive two-cohort random-effects ",
        "summaries within disease."
      ),
      paste0(
        "Brier score, calibration intercept, and slope from leave-one-sample-",
        "out ridge recalibration of the fixed signed score."
      ),
      paste0(
        "Paired OA/OC enrichment of KnockTF perturbational target sets linked ",
        "to four transparent candidate exemplars."
      ),
      paste0(
        "Curated miRNA interaction coverage for the four candidate exemplars; ",
        "disease-specific activity is explicitly unavailable."
      )
    ),
    source = c(
      "GSE114007 and GSE18520 normalized expression; Hallmark gene sets",
      "Table S1 metadata and Table S6 direction-fixed validation results",
      "Table S6 external signed-score sample values",
      "Table S6 external signed-score sample values",
      "KnockTF GMT release and full OA/OC discovery rankings",
      "miRTarBase local release"
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
      "All tables are UTF-8 CSV files. Candidate-centered enrichment and ",
      "interaction tables remain hypothesis-generating and are not treatment ",
      "recommendations."
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
  writeLines(enc2utf8(readme), file.path(paths$tables, "README.md"), useBytes = TRUE)
}

v24_update_registries <- function(project_root, paths) {
  baseline <- file.path(project_root, "results", "submission_v23")
  claims <- utils::read.csv(file.path(
    baseline,
    "claim_evidence_registry_v23.csv"
  ))
  claims <- lapply(claims, function(column) {
    if (is.character(column)) {
      gsub("submission_v23", "submission_v24", column, fixed = TRUE)
    } else {
      column
    }
  })
  claims <- as.data.frame(claims, stringsAsFactors = FALSE)
  additions <- data.frame(
    claim_id = c("C23", "C24"),
    manuscript_claim = c(
      paste0(
        "The same molecular summary showed different retrospective ",
        "separability across OA and OC validation tasks."
      ),
      paste0(
        "Four transparent candidate exemplars showed disease-specific ",
        "candidate-centered Hallmark association contexts."
      )
    ),
    primary_data = c(
      "results/submission_v24/supplementary_tables/Table_S21_cross_cohort_molecular_separability_context.csv",
      "results/submission_v24/supplementary_tables/Table_S20_candidate_centered_Hallmark_context.csv"
    ),
    figure_or_table = c(
      "Figure 4; Table S21",
      "Figure S9; Table S20"
    ),
    allowed_wording = c(
      "cross-cohort molecular separability; validation-task heterogeneity",
      "exploratory residual association enrichment; candidate context"
    ),
    prohibited_wording = c(
      "clinical prediction; OA model failure; OC diagnostic signature",
      "single-gene mechanism; perturbation; pathway activation"
    ),
    status = c(
      "verified with task-scale boundary",
      "verified as exploratory post hoc analysis"
    ),
    stringsAsFactors = FALSE
  )
  claims <- claims[!claims$claim_id %in% additions$claim_id, , drop = FALSE]
  claims <- rbind(claims, additions)
  safe_write_csv(
    claims,
    file.path(paths$root, "claim_evidence_registry_v24.csv")
  )

  checklist <- utils::read.csv(file.path(
    baseline,
    "reproducibility_checklist_v23.csv"
  ))
  checklist$evidence <- gsub(
    "submission_v23",
    "submission_v24",
    checklist$evidence,
    fixed = TRUE
  )
  checklist$item[checklist$item_id == "R19"] <-
    "V2.4 submission strengthening is one-command reproducible."
  checklist$evidence[checklist$item_id == "R19"] <- "run_submission_v24.ps1"
  checklist$item[checklist$item_id == "R21"] <-
    "Six main and eleven supplementary figures have paired PDF/PNG outputs."
  checklist$evidence[checklist$item_id == "R21"] <-
    "results/submission_v24/figures/"
  additions <- data.frame(
    item_id = c("R28", "R29", "R30", "R31", "R32"),
    domain = c(
      "validation interpretation",
      "candidate context",
      "analysis scope",
      "calibration",
      "regulatory context"
    ),
    item = c(
      "External AUCs are interpreted by tissue, comparator, and validation-task scale.",
      "Candidate-centered GSEA adjusts for disease/reference status and is explicitly exploratory.",
      "No DCA, nomogram, PPI, new model family, or LASSO-RF intersection was added.",
      "Calibration is cross-fitted and explicitly not presented as transported external probability calibration.",
      "TF/miRNA evidence is a focused matrix with perturbational/curated provenance and explicit activity boundaries."
    ),
    status = c("complete", "complete", "complete", "complete", "complete"),
    evidence = c(
      "Figure 4; Table S21",
      "Figure S9; Table S20",
      "V2.4 response matrix and source tree",
      "Figure S10; Table S22b",
      "Figure S11; Tables S23a-S23b"
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
    file.path(paths$root, "reproducibility_checklist_v24.csv")
  )
}

run_reviewer_v24 <- function(project_root) {
  for (package in c(
    "ggplot2", "patchwork", "ragg", "clusterProfiler", "BiocParallel"
  )) {
    require_namespace(package, "V2.4 submission strengthening")
  }
  baseline <- file.path(project_root, "results", "submission_v23")
  if (!dir.exists(baseline)) {
    stop("V2.3 baseline outputs are required before V2.4.", call. = FALSE)
  }
  paths <- v24_output_paths(project_root)
  v22_copy_tree(
    baseline,
    paths$root,
    skip = c(
      "claim_evidence_registry_v23.csv",
      "reproducibility_checklist_v23.csv",
      "submission_audit_v23.json"
    )
  )
  log_info("Rebuilding Figure 1 as a linear evidence chain.")
  v24_build_figure1(project_root, paths)
  log_info("Selecting transparent candidate-centered exemplars.")
  selected <- v24_candidate_selection(paths)
  log_info("Running disease-status-adjusted candidate-centered Hallmark analysis.")
  candidate_pathways <- v24_candidate_centered_pathways(
    project_root,
    paths,
    selected
  )
  log_info("Building cross-cohort molecular-separability context table.")
  validation_context <- v24_cross_cohort_context_table(paths)
  log_info("Estimating external signed-score effect sizes and calibration sensitivity.")
  external_effects <- v24_external_effect_sizes(paths)
  calibration <- v24_cross_fitted_calibration(paths)
  v24_build_external_evaluation_figure(
    paths,
    external_effects,
    calibration
  )
  log_info("Building focused TF/miRNA upstream regulatory context.")
  regulatory <- v24_upstream_regulatory_context(
    project_root,
    paths,
    selected
  )
  v24_update_figure_materials(paths)
  v24_update_table_index(paths)
  v24_update_registries(project_root, paths)
  log_info(
    "V2.4 strengthening completed: ",
    nrow(candidate_pathways),
    " candidate-pathway rows and ",
    nrow(validation_context),
    " validation-task rows, ",
    nrow(regulatory$tf),
    " TF rows, and ",
    nrow(regulatory$mirna),
    " miRNA rows."
  )
  invisible(paths)
}
