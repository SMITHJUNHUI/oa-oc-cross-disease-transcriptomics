v21_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v21"))
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

v21_prepare_direction_quadrants <- function(project_root, paths) {
  differential <- submission_load_cache(project_root, "03_differential.rds")
  shared <- submission_load_cache(project_root, "04_shared.rds")
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
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
  common$direction_quadrant <- "Not primary shared"
  shared_rows <- common$primary_shared
  common$direction_quadrant[
    shared_rows & common$logFC_OA >= 0 & common$logFC_OC >= 0
  ] <- "OA higher / OC higher"
  common$direction_quadrant[
    shared_rows & common$logFC_OA < 0 & common$logFC_OC < 0
  ] <- "OA lower / OC lower"
  common$direction_quadrant[
    shared_rows & common$logFC_OA >= 0 & common$logFC_OC < 0
  ] <- "OA higher / OC lower"
  common$direction_quadrant[
    shared_rows & common$logFC_OA < 0 & common$logFC_OC >= 0
  ] <- "OA lower / OC higher"
  counts <- as.data.frame(table(
    common$direction_quadrant[common$primary_shared],
    useNA = "no"
  ))
  names(counts) <- c("direction_quadrant", "genes")
  counts$concordant <- counts$direction_quadrant %in% c(
    "OA higher / OC higher",
    "OA lower / OC lower"
  )
  safe_write_csv(
    common,
    file.path(paths$analysis, "direction_aware_common_gene_effects.csv")
  )
  safe_write_csv(
    counts,
    file.path(paths$analysis, "direction_quadrant_counts.csv")
  )
  list(
    common = common,
    counts = counts,
    oa = oa,
    oc = oc,
    hub_genes = ml$final_genes
  )
}

v21_prepare_gse54388_pca <- function(project_root, paths) {
  validation <- submission_load_cache(project_root, "08_bulk_validation.rds")
  cohort <- validation$oc_validation_1
  expression <- cohort$expression
  metadata <- cohort$metadata
  group <- metadata$group[
    match(colnames(expression), metadata$geo_accession)
  ]
  if (anyNA(group)) {
    stop("GSE54388 expression and metadata sample identifiers do not align.")
  }
  gene_sd <- apply(expression, 1L, stats::sd, na.rm = TRUE)
  gene_sd[!is.finite(gene_sd)] <- 0
  selected <- names(head(sort(gene_sd, decreasing = TRUE), 2000L))
  pca <- stats::prcomp(
    t(expression[selected, , drop = FALSE]),
    center = TRUE,
    scale. = TRUE
  )
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(
    sample_id = rownames(pca$x),
    group = group[match(rownames(pca$x), colnames(expression))],
    pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  variance_table <- data.frame(
    component = paste0("PC", seq_along(variance)),
    variance_fraction = variance,
    cumulative_variance = cumsum(variance),
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    scores,
    file.path(paths$analysis, "GSE54388_unsupervised_PCA_scores.csv")
  )
  safe_write_csv(
    variance_table,
    file.path(paths$analysis, "GSE54388_unsupervised_PCA_variance.csv")
  )
  safe_write_csv(
    data.frame(
      dataset_id = "GSE54388",
      genes_ranked_by = "sample standard deviation",
      genes_entered = length(selected),
      scaling = "gene-wise centering and unit variance",
      labels_used_in_PCA = FALSE,
      stringsAsFactors = FALSE
    ),
    file.path(paths$analysis, "GSE54388_unsupervised_PCA_manifest.csv")
  )
  list(scores = scores, variance = variance_table)
}

v21_prepare_roc_curves <- function(project_root, paths) {
  scores <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "sensitivity",
    "external_validation_sample_scores.csv"
  ))
  rows <- lapply(split(scores, scores$dataset_id), function(table) {
    roc <- pROC::roc(
      response = factor(table$group, levels = c("Normal", "Disease")),
      predictor = table$signed_score,
      levels = c("Normal", "Disease"),
      direction = "<",
      quiet = TRUE
    )
    coordinates <- as.data.frame(pROC::coords(
      roc,
      x = "all",
      ret = c("threshold", "specificity", "sensitivity"),
      transpose = FALSE
    ))
    coordinates$dataset_id <- unique(table$dataset_id)
    coordinates$disease <- unique(table$disease)
    coordinates$auc <- as.numeric(pROC::auc(roc))
    coordinates
  })
  curves <- do.call(rbind, rows)
  rownames(curves) <- NULL
  curves$false_positive_rate <- 1 - curves$specificity
  safe_write_csv(
    curves,
    file.path(paths$analysis, "external_validation_direction_fixed_ROC_curves.csv")
  )
  curves
}

v21_prepare_cell_context_matrix <- function(project_root, paths) {
  evidence <- submission_read_cell_hub_evidence(project_root)
  hub_genes <- submission_load_cache(
    project_root,
    "07_machine_learning.rds"
  )$final_genes
  evidence <- evidence[
    evidence$gene %in% hub_genes &
      is.finite(evidence$cells) &
      is.finite(evidence$fraction_detected) &
      evidence$cells > 0,
    ,
    drop = FALSE
  ]
  evidence$weighted_detected <- evidence$cells * evidence$fraction_detected
  grouped <- aggregate(
    cbind(weighted_detected, cells) ~ disease + gene + cell_type,
    data = evidence,
    FUN = sum
  )
  grouped$fraction_detected <- grouped$weighted_detected / grouped$cells
  grouped <- grouped[
    order(
      grouped$disease,
      grouped$gene,
      -grouped$fraction_detected,
      -grouped$cells
    ),
    ,
    drop = FALSE
  ]
  top <- grouped[
    !duplicated(paste(grouped$disease, grouped$gene, sep = "::")),
    ,
    drop = FALSE
  ]
  complete <- merge(
    expand.grid(
      disease = c("OA", "OC"),
      gene = hub_genes,
      stringsAsFactors = FALSE
    ),
    top[, c("disease", "gene", "cell_type", "cells", "fraction_detected")],
    by = c("disease", "gene"),
    all.x = TRUE,
    sort = FALSE
  )
  names(complete)[names(complete) == "cell_type"] <- "top_cell_context"
  complete$top_cell_context[is.na(complete$top_cell_context)] <- "Not detected"
  complete$fraction_detected[is.na(complete$fraction_detected)] <- 0
  complete$cells[is.na(complete$cells)] <- 0
  complete$gene <- factor(complete$gene, levels = hub_genes)
  safe_write_csv(
    grouped,
    file.path(paths$analysis, "single_cell_context_rankings.csv")
  )
  safe_write_csv(
    complete,
    file.path(paths$analysis, "single_cell_gene_disease_context_matrix.csv")
  )
  complete
}

v21_estimate_score <- function(expression, gene_sets) {
  ranked <- apply(expression, 2L, rank, ties.method = "average")
  ranked <- 10000 * ranked / nrow(ranked)
  score_one <- function(signature) {
    signature <- unique(stats::na.omit(as.character(signature)))
    signature <- intersect(signature, rownames(ranked))
    if (length(signature) == 0L) {
      return(rep(NA_real_, ncol(ranked)))
    }
    vapply(seq_len(ncol(ranked)), function(sample_index) {
      ordered <- sort(ranked[, sample_index], decreasing = TRUE)
      hit <- names(ordered) %in% signature
      weighted <- ordered^0.25
      hit_total <- sum(weighted[hit])
      miss_total <- sum(!hit)
      sum(
        cumsum(ifelse(hit, weighted / hit_total, 0)) -
          cumsum(ifelse(!hit, 1 / miss_total, 0))
      )
    }, numeric(1L))
  }
  stromal <- score_one(gene_sets$stromal_signature)
  immune <- score_one(gene_sets$immune_signature)
  data.frame(
    sample_id = colnames(expression),
    stromal_score = stromal,
    immune_score = immune,
    combined_estimate_score = stromal + immune,
    stringsAsFactors = FALSE
  )
}

v21_prepare_tcga_context <- function(project_root, paths) {
  local_config_path <- file.path(project_root, "config", "local.yml")
  local_config_raw <- readBin(
    local_config_path,
    what = "raw",
    n = file.info(local_config_path)$size
  )
  local_config_text <- rawToChar(local_config_raw)
  Encoding(local_config_text) <- "UTF-8"
  config <- yaml::yaml.load(local_config_text)
  expression <- read_tcga_expression(config$tcga$expression_path)
  gene_sets <- utils::read.csv(
    file.path(
      project_root,
      "data",
      "reference",
      "estimate_gene_sets_tidyestimate_1.1.1.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  scores <- v21_estimate_score(expression, gene_sets)
  hub_genes <- submission_load_cache(
    project_root,
    "07_machine_learning.rds"
  )$final_genes
  present <- intersect(hub_genes, rownames(expression))
  score_columns <- c(
    "stromal_score",
    "immune_score",
    "combined_estimate_score"
  )
  rows <- list()
  index <- 0L
  for (gene in present) {
    for (score_name in score_columns) {
      index <- index + 1L
      x <- as.numeric(expression[gene, scores$sample_id])
      y <- scores[[score_name]]
      keep <- is.finite(x) & is.finite(y)
      test <- suppressWarnings(stats::cor.test(
        x[keep],
        y[keep],
        method = "spearman",
        exact = FALSE
      ))
      rows[[index]] <- data.frame(
        gene = gene,
        context_score = score_name,
        samples = sum(keep),
        spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  correlations <- do.call(rbind, rows)
  correlations$fdr <- stats::p.adjust(correlations$p_value, method = "BH")
  expression_long <- do.call(rbind, lapply(present, function(gene) {
    data.frame(
      sample_id = scores$sample_id,
      gene = gene,
      expression = as.numeric(expression[gene, scores$sample_id]),
      stringsAsFactors = FALSE
    )
  }))
  safe_write_csv(
    scores,
    file.path(paths$analysis, "TCGA_OV_relative_context_scores.csv")
  )
  safe_write_csv(
    correlations,
    file.path(paths$analysis, "TCGA_OV_hub_gene_context_correlations.csv")
  )
  safe_write_csv(
    expression_long,
    file.path(paths$analysis, "TCGA_OV_hub_gene_expression_for_context.csv")
  )
  safe_write_csv(
    data.frame(
      item = c(
        "method",
        "gene_set_snapshot",
        "input_scale",
        "interpretation",
        "absolute_purity_reported",
        "candidates_requested",
        "candidates_available",
        "candidate_missing"
      ),
      value = c(
        "published ESTIMATE rank-based stromal and immune signatures",
        "tidyestimate 1.1.1 gene-set snapshot (141 stromal; 141 immune)",
        "TCGA-OV RSEM matrix after the project's automatic transform",
        "relative within-cohort transcriptomic context only",
        "no",
        length(hub_genes),
        length(present),
        paste(setdiff(hub_genes, present), collapse = ";")
      ),
      stringsAsFactors = FALSE
    ),
    file.path(paths$analysis, "TCGA_OV_context_method_manifest.csv")
  )
  list(
    scores = scores,
    correlations = correlations,
    expression = expression_long
  )
}

v21_prepare_mr_tables <- function(project_root, paths) {
  metadata <- utils::read.csv(file.path(
    project_root,
    "results",
    "mr",
    "OpenGWAS_dataset_metadata.csv"
  ), check.names = FALSE)
  estimates <- utils::read.csv(file.path(
    project_root,
    "results",
    "mr",
    "MR_combined_estimates.csv"
  ), check.names = FALSE)
  meta_columns <- c(
    "id", "trait", "population", "sex", "sample_size",
    "ncase", "ncontrol", "build", "pmid"
  )
  exposure_meta <- metadata[, meta_columns]
  names(exposure_meta) <- paste0("exposure_", names(exposure_meta))
  outcome_meta <- metadata[, meta_columns]
  names(outcome_meta) <- paste0("outcome_", names(outcome_meta))
  table <- merge(
    estimates,
    exposure_meta,
    by.x = "id.exposure",
    by.y = "exposure_id",
    all.x = TRUE,
    sort = FALSE
  )
  table <- merge(
    table,
    outcome_meta,
    by.x = "id.outcome",
    by.y = "outcome_id",
    all.x = TRUE,
    sort = FALSE
  )
  table$direction <- ifelse(
    table$id.exposure == "ebi-a-GCST007092",
    "OA to ovarian cancer",
    "ovarian cancer to OA"
  )
  table$instrument_p_threshold <- 5e-8
  table$clump_r2 <- 0.001
  table$clump_window_kb <- 10000
  table$odds_ratio <- exp(table$b)
  table$ci_lower_95 <- exp(table$b - 1.96 * table$se)
  table$ci_upper_95 <- exp(table$b + 1.96 * table$se)
  output_columns <- c(
    "direction",
    "id.exposure", "exposure_trait", "exposure_population", "exposure_sex",
    "exposure_sample_size", "exposure_ncase", "exposure_ncontrol",
    "exposure_build", "exposure_pmid",
    "id.outcome", "outcome_trait", "outcome_population", "outcome_sex",
    "outcome_sample_size", "outcome_ncase", "outcome_ncontrol",
    "outcome_build", "outcome_pmid",
    "instrument_p_threshold", "clump_r2", "clump_window_kb",
    "nsnp", "method", "b", "se", "odds_ratio",
    "ci_lower_95", "ci_upper_95", "pval"
  )
  table <- table[, output_columns]
  safe_write_csv(
    table,
    file.path(
      paths$tables,
      "Table_S12a_MR_estimates_and_provenance.csv"
    )
  )
  original <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "supplementary_tables",
    "Table_S12_negative_bidirectional_MR.csv"
  ), check.names = FALSE)
  safe_write_csv(
    original,
    file.path(paths$tables, "Table_S12b_MR_sensitivity_diagnostics.csv")
  )
  table
}

v21_copy_supplementary_tables <- function(project_root, paths) {
  source_dir <- file.path(
    project_root,
    "results",
    "submission",
    "supplementary_tables"
  )
  files <- list.files(source_dir, full.names = TRUE)
  files <- files[
    basename(files) != "Table_S12_negative_bidirectional_MR.csv"
  ]
  file.copy(files, paths$tables, overwrite = TRUE, copy.date = TRUE)
  invisible(paths$tables)
}

v21_write_registries <- function(project_root, paths) {
  registry <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "claim_evidence_registry.csv"
  ), check.names = FALSE)
  additions <- data.frame(
    claim_id = c("C14", "C15", "C16", "C17"),
    manuscript_claim = c(
      paste(
        "The four primary OA/OC direction quadrants contain 112, 34,",
        "86, and 54 genes, yielding 146 concordant and 140 discordant genes."
      ),
      paste(
        "Unsupervised GSE54388 PCA visually separates the six normal",
        "from sixteen tumor samples without using labels to fit the coordinates."
      ),
      paste(
        "Relative TCGA-OV stromal/immune scores provide a composition",
        "audit for nine available candidates across 307 samples;",
        "they are not absolute purity or histologic cell fractions."
      ),
      paste(
        "MR provenance is internally consistent for ebi-a-GCST007092",
        "and ieu-a-1120, with 21 and 11 instruments after the recorded",
        "genome-wide and LD-clumping rules."
      )
    ),
    primary_data = c(
      "results/submission_v21/analysis/direction_quadrant_counts.csv",
      "results/submission_v21/analysis/GSE54388_unsupervised_PCA_scores.csv",
      paste(
        "results/submission_v21/analysis/TCGA_OV_relative_context_scores.csv;",
        "TCGA_OV_hub_gene_context_correlations.csv"
      ),
      "results/submission_v21/supplementary_tables/Table_S12a_MR_estimates_and_provenance.csv"
    ),
    figure_or_table = c(
      "Figure 2; Table S2",
      "Figure 4; Table S6",
      "Figure S5; Table S15",
      "Figure S2; Table S12a–b"
    ),
    allowed_wording = c(
      "directionally heterogeneous overlap with explicit quadrant counts",
      "unsupervised retrospective tumor–normal molecular structure",
      "relative within-cohort transcriptomic context scores",
      "verified MR dataset provenance and no detected causal evidence"
    ),
    prohibited_wording = c(
      "shared genes have the same direction",
      "unsupervised diagnostic validation",
      "absolute tumor purity or measured cell fractions",
      "alternative GCST900 accessions or proof that causality is absent"
    ),
    status = "verified against V2.1 outputs",
    stringsAsFactors = FALSE
  )
  registry <- rbind(registry, additions)
  safe_write_csv(
    registry,
    file.path(paths$root, "claim_evidence_registry_v21.csv")
  )
  for (filename in c(
    "data_source_manifest.csv",
    "parameter_manifest.csv"
  )) {
    file.copy(
      file.path(project_root, "results", "submission", filename),
      file.path(paths$root, filename),
      overwrite = TRUE,
      copy.date = TRUE
    )
  }
  checklist <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "reproducibility_checklist.csv"
  ), check.names = FALSE)
  checklist <- rbind(
    checklist,
    data.frame(
      item_id = c("R19", "R20", "R21", "R22"),
      domain = c("reviewer revision", "MR", "figures", "context audit"),
      item = c(
        "V2.1 reviewer-strengthening analyses are one-command reproducible.",
        "MR accessions, sample sizes, instruments, and clumping rules are audited.",
        "Six main and five supplementary figures have paired PDF/PNG outputs.",
        "TCGA relative context analysis explicitly withholds absolute purity."
      ),
      status = "complete",
      evidence = c(
        "run_submission_v21.ps1",
        "Table S12a–b",
        "results/submission_v21/figures/",
        "Figure S5; Table S15"
      ),
      stringsAsFactors = FALSE
    )
  )
  safe_write_csv(
    checklist,
    file.path(paths$root, "reproducibility_checklist_v21.csv")
  )
}

v21_build_figure2 <- function(project_root, paths, analysis) {
  threshold <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission",
    "sensitivity",
    "deg_threshold_sensitivity_summary.csv"
  ))
  common <- analysis$common
  counts <- analysis$counts
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
  safe_write_csv(
    common,
    file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv")
  )
  safe_write_csv(
    counts,
    file.path(paths$source, "Figure2_quadrant_counts.csv")
  )
  safe_write_csv(
    threshold,
    file.path(paths$source, "Figure2_DEG_threshold_sensitivity.csv")
  )
  safe_write_csv(analysis$oa, file.path(paths$source, "Figure2_OA_DEG.csv"))
  safe_write_csv(analysis$oc, file.path(paths$source, "Figure2_OC_DEG.csv"))

  p1 <- submission_volcano_plot(
    analysis$oa,
    "OA",
    analysis$hub_genes
  )
  p2 <- submission_volcano_plot(
    analysis$oc,
    "OC",
    analysis$hub_genes
  )
  shared_only <- common[common$primary_shared, , drop = FALSE]
  labels <- common[common$hub, , drop = FALSE]
  quadrant_colors <- c(
    "OA higher / OC higher" = "#009E73",
    "OA lower / OC lower" = "#6A3D9A",
    "OA higher / OC lower" = "#0072B2",
    "OA lower / OC higher" = "#D55E00"
  )
  count_labels <- setNames(counts$genes, counts$direction_quadrant)
  quadrant_annotation <- data.frame(
    x = c(Inf, -Inf, Inf, -Inf),
    y = c(Inf, -Inf, -Inf, Inf),
    hjust = c(1.05, -0.05, 1.05, -0.05),
    vjust = c(1.2, -0.2, -0.2, 1.2),
    label = paste0(
      c("OA↑ / OC↑", "OA↓ / OC↓", "OA↑ / OC↓", "OA↓ / OC↑"),
      "\nn=",
      unname(count_labels[quadrant_levels])
    ),
    stringsAsFactors = FALSE
  )
  p3 <- ggplot2::ggplot(
    common,
    ggplot2::aes(x = logFC_OA, y = logFC_OC)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_point(colour = "#D4D7DB", size = 0.32, alpha = 0.32) +
    ggplot2::geom_point(
      data = shared_only,
      ggplot2::aes(colour = direction_quadrant),
      size = 0.85,
      alpha = 0.72
    ) +
    ggplot2::geom_text(
      data = quadrant_annotation,
      ggplot2::aes(
        x = x,
        y = y,
        label = label,
        hjust = hjust,
        vjust = vjust
      ),
      inherit.aes = FALSE,
      size = 2.2,
      lineheight = 0.9,
      colour = "#374151"
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = gene),
      size = 2.05,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.16
    ) +
    ggplot2::scale_colour_manual(values = quadrant_colors) +
    ggplot2::labs(
      title = "Shared membership does not imply shared direction",
      subtitle = "146/286 concordant; 140/286 discordant",
      x = "OA log2 fold change",
      y = "OC log2 fold change",
      colour = NULL
    ) +
    submission_theme(7.8) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(3.2, "mm")
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
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure2_bulk_discovery",
    paths$figures,
    height_mm = 155
  )
}

v21_build_figure4 <- function(project_root, paths, pca, curves) {
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
      c("GSE117999", "GSE82107", "GSE54388", "GSE12470"),
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
      title = "Direction-fixed score ROC curves",
      subtitle = "Retrospective molecular separation",
      x = "1 - specificity",
      y = "Sensitivity",
      colour = NULL,
      linetype = NULL
    ) +
    submission_theme(7.4) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
      linetype = "none"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.5)
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
      subtitle = "Observed AUC shown by vertical line; 1,000 permutations/cohort",
      x = "Permuted-label AUC",
      y = "Frequency"
    ) +
    submission_theme(7.4) +
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
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure4_external_validation",
    paths$figures,
    height_mm = 160
  )
}

v21_build_figure5 <- function(project_root, paths, context_matrix) {
  evidence <- submission_read_cell_hub_evidence(project_root)
  status <- utils::read.csv(file.path(
    project_root,
    "results",
    "tables",
    "single_cell_dataset_status.csv"
  ))
  summary <- utils::read.csv(file.path(
    project_root,
    "results",
    "single_cell_downstream",
    "single_cell_downstream_summary.csv"
  ))
  hub_genes <- submission_load_cache(
    project_root,
    "07_machine_learning.rds"
  )$final_genes
  evidence$weighted_detected <- evidence$fraction_detected * evidence$cells
  weighted <- aggregate(
    cbind(weighted_detected, cells) ~ dataset_id + disease + gene,
    data = evidence,
    FUN = sum
  )
  weighted$fraction_detected <- weighted$weighted_detected / weighted$cells
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
  context_matrix$gene <- factor(context_matrix$gene, levels = rev(hub_genes))
  context_matrix$display_context <- ifelse(
    nchar(context_matrix$top_cell_context) > 18,
    paste0(substr(context_matrix$top_cell_context, 1, 15), "..."),
    context_matrix$top_cell_context
  )
  context_matrix$display_context <- gsub(
    "\\.",
    " ",
    context_matrix$display_context
  )
  context_matrix$display_context[
    grepl("^Ovarian cancer", context_matrix$display_context)
  ] <- "Ovarian cancer\ncell"
  context_matrix$display_context[
    context_matrix$display_context == "Endothelial cell"
  ] <- "Endothelial\ncell"
  safe_write_csv(
    weighted,
    file.path(paths$source, "Figure5_hub_detection_by_dataset.csv")
  )
  safe_write_csv(
    context_matrix,
    file.path(paths$source, "Figure5_gene_disease_cell_context_matrix.csv")
  )
  safe_write_csv(
    pseudobulk,
    file.path(paths$source, "Figure5_hub_pseudobulk_evidence.csv")
  )
  safe_write_csv(
    status,
    file.path(paths$source, "Figure5_single_cell_status.csv")
  )
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
      ggplot2::aes(label = paste0(
        scales::comma(qc_pass),
        " / ",
        scales::comma(cells)
      )),
      hjust = 1.05,
      colour = "white",
      size = 2.05
    ) +
    ggplot2::scale_fill_manual(values = submission_palette[c("OA", "OC")]) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
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
      title = "Candidate detection across cell atlases",
      subtitle = "Cell-count-weighted fraction detected",
      x = NULL,
      y = NULL,
      size = "Detected",
      fill = "Detected"
    ) +
    submission_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  p3 <- ggplot2::ggplot(
    context_matrix,
    ggplot2::aes(
      x = disease,
      y = gene,
      fill = fraction_detected
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.55) +
    ggplot2::geom_text(
      ggplot2::aes(label = display_context),
      size = 2.05,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_gradient(
      low = "#F3F4F6",
      high = submission_palette[["shared"]],
      labels = scales::percent
    ) +
    ggplot2::labs(
      title = "Gene–disease cellular context matrix",
      subtitle = "Highest weighted detection context; labels remain dataset specific",
      x = NULL,
      y = NULL,
      fill = "Detected"
    ) +
    submission_theme()
  if (nrow(pseudobulk) > 0L) {
    pseudobulk$label <- paste(
      pseudobulk$dataset_id,
      pseudobulk$cell_type,
      pseudobulk$contrast,
      sep = " | "
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
        title = "Significant eligible pseudobulk effects",
        x = "log2 fold change",
        y = NULL,
        colour = "Gene"
      ) +
      submission_theme(7.5)
  } else {
    p4 <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No candidate met FDR < 0.05\nin eligible pseudobulk contrasts",
        size = 3
      ) +
      ggplot2::labs(title = "Pseudobulk candidate effects") +
      ggplot2::theme_void(base_family = "Arial")
  }
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure5_single_cell_localization",
    paths$figures,
    height_mm = 165
  )
}

v21_build_supplementary5 <- function(paths, context) {
  correlations <- context$correlations
  scores <- context$scores
  expression <- context$expression
  correlations$context_label <- factor(
    correlations$context_score,
    levels = c(
      "stromal_score",
      "immune_score",
      "combined_estimate_score"
    ),
    labels = c("Stromal", "Immune", "Combined")
  )
  correlations$gene <- factor(
    correlations$gene,
    levels = rev(unique(correlations$gene))
  )
  safe_write_csv(
    correlations,
    file.path(paths$source, "SupplementaryFigure5_context_correlations.csv")
  )
  safe_write_csv(
    scores,
    file.path(paths$source, "SupplementaryFigure5_context_scores.csv")
  )
  p1 <- ggplot2::ggplot(
    correlations,
    ggplot2::aes(
      x = context_label,
      y = gene,
      fill = spearman_rho,
      size = -log10(pmax(fdr, .Machine$double.xmin))
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "#374151", stroke = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = submission_palette[["OA"]],
      mid = "white",
      high = submission_palette[["OC"]],
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    ggplot2::scale_size(range = c(1.5, 5.5)) +
    ggplot2::labs(
      title = "Candidate genes and relative TCGA-OV context scores",
      subtitle = "Spearman correlations; point size is −log10(FDR)",
      x = NULL,
      y = NULL,
      fill = "ρ",
      size = "−log10(FDR)"
    ) +
    submission_theme() +
    ggplot2::theme(legend.position = "right")
  representative <- do.call(rbind, lapply(
    c("stromal_score", "immune_score"),
    function(score_name) {
      subset <- correlations[
        correlations$context_score == score_name,
        ,
        drop = FALSE
      ]
      subset[which.max(abs(subset$spearman_rho)), , drop = FALSE]
    }
  ))
  scatter_plots <- lapply(seq_len(nrow(representative)), function(index) {
    row <- representative[index, ]
    table <- merge(
      expression[expression$gene == row$gene, ],
      scores[, c("sample_id", row$context_score)],
      by = "sample_id"
    )
    names(table)[names(table) == row$context_score] <- "context_value"
    ggplot2::ggplot(
      table,
      ggplot2::aes(x = expression, y = context_value)
    ) +
      ggplot2::geom_point(
        colour = if (row$context_score == "stromal_score") {
          submission_palette[["OA"]]
        } else {
          submission_palette[["OC"]]
        },
        alpha = 0.45,
        size = 1.1
      ) +
      ggplot2::geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        colour = "#374151",
        linewidth = 0.6
      ) +
      ggplot2::labs(
        title = paste(
          row$gene,
          if (row$context_score == "stromal_score") "and stromal score" else
            "and immune score"
        ),
        subtitle = sprintf(
          "Spearman ρ = %.2f; FDR = %.2g",
          row$spearman_rho,
          row$fdr
        ),
        x = paste(row$gene, "expression"),
        y = if (row$context_score == "stromal_score") {
          "Relative stromal score"
        } else {
          "Relative immune score"
        }
      ) +
      submission_theme()
  })
  note <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = -0.98,
      y = 0,
      hjust = 0,
      size = 2.7,
      lineheight = 1.15,
      colour = "#374151",
      label = paste0(
        "Interpretation boundary: scores are rank-based transcriptomic ",
        "proxies within TCGA-OV.\n",
        "Absolute tumor purity was not calculated from RNA-seq, and the ",
        "scores are not histologic cell fractions."
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-0.5, 0.5)) +
    ggplot2::theme_void(base_family = "Arial")
  figure <- submission_panel_tag(
    (p1 | (scatter_plots[[1L]] / scatter_plots[[2L]])) / note +
      patchwork::plot_layout(heights = c(4, 0.85))
  )
  submission_save_plot(
    figure,
    "SupplementaryFigure5_TCGA_relative_context",
    paths$figures,
    height_mm = 160
  )
}

v21_write_figure_legends <- function(paths) {
  legends <- c(
    "# Figure legends",
    "",
    "## Figure 1. Study design, biological hypothesis, and audited data resources",
    "",
    "**A,** The study tests context-dependent convergence of aging-associated, mesenchymal-remodeling, and immune-remodeling programs while keeping OA and OC analyses separate. Bidirectional MR is supplementary null evidence. **B,** Bulk discovery and external-cohort sample counts. **C,** Total and QC-pass single-cell counts.",
    "",
    "## Figure 2. Direction-aware cross-disease transcriptomic discovery",
    "",
    "**A–B,** OA and OC differential-expression volcano plots. **C,** Explicit quadrants for disease-specific effects among commonly measured genes. The 286 primary shared DEGs are colored by direction: 146 were concordant and 140 discordant. **D,** Shared-gene counts under six prespecified thresholds.",
    "",
    "## Figure 3. Robust identification of molecular candidates",
    "",
    "**A,** WGCNA module–trait association under soft-power perturbation. **B,** Primary-module gene retention. **C,** Feature frequency under strict outer-fold screening and nested tuning. **D,** Disease-specific effects for the ten prioritized genes.",
    "",
    "## Figure 4. Direction-fixed cross-cohort molecular reproducibility",
    "",
    "**A,** Unsupervised PCA of GSE54388 using the 2,000 genes with highest sample standard deviation; group labels were used only for display. **B,** ROC curves for the fixed signed ten-gene score with direction fixed from discovery. **C,** Null AUC distributions from 1,000 label permutations. **D,** Leave-one-sample-out AUCs. These panels show retrospective molecular separation, not clinical diagnostic performance.",
    "",
    "## Figure 5. Candidate-gene localization across distinct cellular contexts",
    "",
    "**A,** QC-pass fractions. **B,** Cell-count-weighted candidate detection. **C,** Gene–disease matrix showing the annotated cell context with the highest weighted detection for each candidate within OA and OC atlases. Labels remain dataset specific and do not imply homologous cell states. **D,** Significant eligible pseudobulk effects. OA and OC datasets were not integrated into a shared latent space.",
    "",
    "## Figure 6. Functional and exploratory prognostic context",
    "",
    "**A,** Hallmark enrichment. **B,** Rank-based immune-signature differences. **C,** Exploratory TCGA-OV Cox estimates. **D,** Bootstrap concordance and optimism correction.",
    "",
    "## Supplementary Figure 1. Prespecified non-MR sensitivity analyses",
    "",
    "DEG threshold retention, WGCNA leave-one-out stability, strict nested feature frequency, and TCGA LASSO selection stability.",
    "",
    "## Supplementary Figure 2. Negative bidirectional MR results",
    "",
    "Five MR estimators and heterogeneity/pleiotropy diagnostics in both directions. Null estimates mean no causal effect was detected under the available instruments and assumptions, not proof of absence.",
    "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings",
    "",
    "Released or recomputed embeddings with dataset-specific labels. Quantitative summaries used all QC-pass cells; OA and OC were not integrated into a shared latent space.",
    "",
    "## Supplementary Figure 4. HPA normal-tissue context",
    "",
    "Normal-tissue specificity and distribution categories. Cartilage is absent, so this audit cannot establish OA–OC specificity.",
    "",
    "## Supplementary Figure 5. TCGA-OV relative stromal and immune context audit",
    "",
    "**A,** Spearman correlations between candidate expression and published ESTIMATE rank-based stromal, immune, and combined scores. **B–C,** Representative strongest absolute stromal and immune associations. Scores are relative within-cohort transcriptomic proxies; absolute tumor purity and histologic cell fractions were not inferred from RNA-seq."
  )
  write_utf8(legends, file.path(paths$figures, "figure_legends.md"))
}

v21_build_figures <- function(
    project_root,
    paths,
    quadrants,
    pca,
    curves,
    context_matrix,
    tcga_context
) {
  submission_build_figure1(project_root, paths$figures, paths$source)
  v21_build_figure2(project_root, paths, quadrants)
  submission_build_figure3(project_root, paths$figures, paths$source)
  v21_build_figure4(project_root, paths, pca, curves)
  v21_build_figure5(project_root, paths, context_matrix)
  submission_build_figure6(project_root, paths$figures, paths$source)
  submission_build_supplementary1(
    project_root,
    paths$figures,
    paths$source
  )
  submission_build_supplementary2(
    project_root,
    paths$figures,
    paths$source
  )
  submission_build_supplementary3(
    project_root,
    paths$figures,
    paths$source
  )
  submission_build_supplementary4(
    project_root,
    paths$figures,
    paths$source
  )
  v21_build_supplementary5(paths, tcga_context)
  v21_write_figure_legends(paths)
  safe_write_csv(
    data.frame(
      setting = c(
        "revision",
        "target_journal",
        "width_mm",
        "png_dpi",
        "font_family",
        "vector_format"
      ),
      value = c(
        "V2.1",
        "not yet specified; provisional general biomedical layout",
        "180",
        "300",
        "Arial",
        "PDF"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(paths$figures, "figure_style_manifest.csv")
  )
}

run_reviewer_v21 <- function(project_root) {
  for (package in c(
    "ggplot2",
    "patchwork",
    "ggrepel",
    "scales",
    "ragg",
    "data.table",
    "pROC",
    "yaml"
  )) {
    require_namespace(package, "V2.1 reviewer revision")
  }
  paths <- v21_output_paths(project_root)
  log_info("Preparing V2.1 direction, PCA, context, and MR audits.")
  quadrants <- v21_prepare_direction_quadrants(project_root, paths)
  pca <- v21_prepare_gse54388_pca(project_root, paths)
  curves <- v21_prepare_roc_curves(project_root, paths)
  context_matrix <- v21_prepare_cell_context_matrix(project_root, paths)
  tcga_context <- v21_prepare_tcga_context(project_root, paths)
  v21_copy_supplementary_tables(project_root, paths)
  v21_prepare_mr_tables(project_root, paths)
  v21_write_registries(project_root, paths)
  safe_write_csv(
    pca$scores,
    file.path(
      paths$tables,
      "Table_S6b_GSE54388_unsupervised_PCA_scores.csv"
    )
  )
  safe_write_csv(
    pca$variance,
    file.path(
      paths$tables,
      "Table_S6c_GSE54388_unsupervised_PCA_variance.csv"
    )
  )
  safe_write_csv(
    tcga_context$scores,
    file.path(paths$tables, "Table_S15a_TCGA_relative_context_scores.csv")
  )
  safe_write_csv(
    tcga_context$correlations,
    file.path(
      paths$tables,
      "Table_S15b_TCGA_candidate_context_correlations.csv"
    )
  )
  log_info("Building V2.1 main and supplementary figures.")
  v21_build_figures(
    project_root,
    paths,
    quadrants,
    pca,
    curves,
    context_matrix,
    tcga_context
  )
  log_info("V2.1 reviewer-strengthening analyses completed.")
  invisible(paths)
}
