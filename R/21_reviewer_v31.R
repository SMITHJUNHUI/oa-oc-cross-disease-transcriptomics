v31_base_runner <- run_reviewer_v30

v31_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v31"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs"))
  )
}

v31_prepare_baseline <- function(project_root, paths) {
  source_root <- file.path(project_root, "results", "submission_v30")
  if (!dir.exists(source_root)) {
    stop("V3.0 outputs are required for V3.1.", call. = FALSE)
  }
  v22_copy_tree(source_root, paths$root)
  obsolete <- file.path(paths$figures, c(
    "Figure1_systems_framework.png", "Figure1_systems_framework.pdf",
    "Figure2_gene_direction.png", "Figure2_gene_direction.pdf",
    "Figure6_genetic_liability_boundary.png", "Figure6_genetic_liability_boundary.pdf"
  ))
  unlink(obsolete[file.exists(obsolete)])
  invisible(paths)
}

v31_panel_ranking <- function(project_root, paths) {
  model_files <- list(
    OA_lasso = file.path(project_root, "results", "tables", "ML_OA_GSE114007_lasso_coefficients.csv"),
    OA_rf = file.path(project_root, "results", "tables", "ML_OA_GSE114007_random_forest_importance.csv"),
    OC_lasso = file.path(project_root, "results", "tables", "ML_OC_GSE18520_lasso_coefficients.csv"),
    OC_rf = file.path(project_root, "results", "tables", "ML_OC_GSE18520_random_forest_importance.csv")
  )
  if (any(!file.exists(unlist(model_files)))) {
    stop("Original model evidence files are missing.", call. = FALSE)
  }
  oa_lasso <- utils::read.csv(model_files$OA_lasso)
  oa_rf <- utils::read.csv(model_files$OA_rf)
  oc_lasso <- utils::read.csv(model_files$OC_lasso)
  oc_rf <- utils::read.csv(model_files$OC_rf)
  support <- list(
    OA_lasso = toupper(oa_lasso$gene[oa_lasso$coefficient != 0]),
    OA_rf = toupper(head(oa_rf$gene[order(-oa_rf$importance)], 20L)),
    OC_lasso = toupper(oc_lasso$gene[oc_lasso$coefficient != 0]),
    OC_rf = toupper(head(oc_rf$gene[order(-oc_rf$importance)], 20L))
  )
  votes <- sort(table(unlist(support, use.names = FALSE)), decreasing = TRUE)
  ranked <- unique(c("SOX9", names(votes)))
  expected <- c(
    "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1", "CFI", "DDIT3", "DIRAS3",
    "EFEMP1", "HK2", "KIT", "MYZAP", "NOD2", "OGN", "RTN1"
  )
  top15 <- head(ranked, 15L)
  if (!identical(top15, expected)) {
    stop(
      "The reconstructed original model-vote ranking changed: ",
      paste(top15, collapse = ", "),
      call. = FALSE
    )
  }
  shared <- utils::read.csv(file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv"))
  shared <- shared[as.logical(shared$primary_shared), , drop = FALSE]
  rows <- lapply(seq_along(expected), function(index) {
    gene <- expected[[index]]
    record <- shared[match(gene, toupper(shared$gene)), , drop = FALSE]
    data.frame(
      gene = gene,
      original_rank = index,
      original_model_vote_count = as.integer(votes[[gene]] %||% 0L),
      cross_disease_model_consensus = gene == "SOX9",
      top5_member = index <= 5L,
      top10_member = index <= 10L,
      top15_member = TRUE,
      log2FC_OA = record$logFC_OA,
      log2FC_OC = record$logFC_OC,
      direction_class = ifelse(sign(record$logFC_OA) == sign(record$logFC_OC), "concordant", "discordant"),
      ranking_provenance = ifelse(
        gene == "SOX9",
        "original cross-disease model consensus",
        "deterministic extension of the original four-model vote ranking"
      ),
      interpretation = "panel-size sensitivity only; not an optimized predictive signature",
      stringsAsFactors = FALSE
    )
  })
  composition <- do.call(rbind, rows)
  safe_write_csv(composition, file.path(paths$tables, "Table_S29a_panel_size_composition.csv"))
  composition
}

v31_panel_direction_pathway <- function(project_root, paths, composition) {
  panel_sizes <- c(5L, 10L, 15L)
  direction_rows <- lapply(panel_sizes, function(size) {
    selected <- composition[composition$original_rank <= size, , drop = FALSE]
    counts <- table(factor(selected$direction_class, levels = c("concordant", "discordant")))
    data.frame(
      panel = paste0("Top ", size),
      panel_size = size,
      direction_class = names(counts),
      genes = as.integer(counts),
      proportion = as.integer(counts) / size,
      median_absolute_log2FC_OA = stats::median(abs(selected$log2FC_OA)),
      median_absolute_log2FC_OC = stats::median(abs(selected$log2FC_OC)),
      stringsAsFactors = FALSE
    )
  })
  direction <- do.call(rbind, direction_rows)

  cache <- readRDS(file.path(project_root, "results", "cache", "05_enrichment.rds"))
  if (is.list(cache) && "value" %in% names(cache)) cache <- cache$value
  gene_sets <- cache$gsea$oa_train$hallmark@geneSets
  shared <- utils::read.csv(file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv"))
  background <- unique(toupper(shared$gene[as.logical(shared$primary_shared)]))
  universe_size <- length(background)
  ora_rows <- list()
  row_index <- 0L
  for (size in panel_sizes) {
    selected <- toupper(composition$gene[composition$original_rank <= size])
    for (pathway_id in names(gene_sets)) {
      members <- intersect(toupper(gene_sets[[pathway_id]]), background)
      overlap <- intersect(selected, members)
      row_index <- row_index + 1L
      ora_rows[[row_index]] <- data.frame(
        panel = paste0("Top ", size),
        panel_size = size,
        pathway_id = pathway_id,
        pathway = tools::toTitleCase(tolower(gsub("^HALLMARK_|_", " ", pathway_id))),
        overlap_genes = paste(overlap, collapse = ";"),
        overlap_count = length(overlap),
        pathway_genes_in_shared_background = length(members),
        shared_background_size = universe_size,
        fold_enrichment = if (length(members) > 0L) {
          (length(overlap) / size) / (length(members) / universe_size)
        } else {
          NA_real_
        },
        p_value = if (length(members) > 0L) {
          stats::phyper(length(overlap) - 1L, length(members), universe_size - length(members), size, lower.tail = FALSE)
        } else {
          1
        },
        stringsAsFactors = FALSE
      )
    }
  }
  ora <- do.call(rbind, ora_rows)
  ora$FDR <- ave(ora$p_value, ora$panel, FUN = stats::p.adjust, method = "BH")
  ora$analysis_scope <- "descriptive Hallmark over-representation against the 286 shared-DEG background"
  ora$inference_boundary <- "small-panel sensitivity; not pathway activity or a new candidate-selection rule"
  safe_write_csv(direction, file.path(paths$tables, "Table_S29b_panel_size_direction_sensitivity.csv"))
  safe_write_csv(ora, file.path(paths$tables, "Table_S29c_panel_size_Hallmark_sensitivity.csv"))

  pair_rows <- list()
  pair_index <- 0L
  for (pair in list(c("Top 5", "Top 10"), c("Top 5", "Top 15"), c("Top 10", "Top 15"))) {
    left <- ora[ora$panel == pair[[1L]], c("pathway_id", "fold_enrichment", "p_value")]
    right <- ora[ora$panel == pair[[2L]], c("pathway_id", "fold_enrichment", "p_value")]
    joined <- merge(left, right, by = "pathway_id", suffixes = c("_left", "_right"))
    left_top <- head(joined$pathway_id[order(joined$p_value_left, -joined$fold_enrichment_left)], 10L)
    right_top <- head(joined$pathway_id[order(joined$p_value_right, -joined$fold_enrichment_right)], 10L)
    pair_index <- pair_index + 1L
    pair_rows[[pair_index]] <- data.frame(
      panel_left = pair[[1L]],
      panel_right = pair[[2L]],
      spearman_fold_enrichment = suppressWarnings(stats::cor(
        joined$fold_enrichment_left, joined$fold_enrichment_right,
        method = "spearman", use = "pairwise.complete.obs"
      )),
      top10_jaccard = length(intersect(left_top, right_top)) / length(union(left_top, right_top)),
      stringsAsFactors = FALSE
    )
  }
  comparisons <- do.call(rbind, pair_rows)
  safe_write_csv(comparisons, file.path(paths$tables, "Table_S29d_panel_size_pathway_profile_comparison.csv"))
  list(direction = direction, ora = ora, comparisons = comparisons)
}

v31_panel_localization <- function(project_root, paths, composition) {
  require_namespace("data.table", "panel-size single-cell localization sensitivity")
  cache_dir <- ensure_dir(file.path(project_root, "results", "cache", "submission_v31"))
  cache_path <- file.path(cache_dir, "panel_size_detection_localization.rds")
  genes <- composition$gene
  if (file.exists(cache_path)) {
    evidence <- readRDS(cache_path)
  } else {
    evidence_parts <- list()
    annotation_path <- file.path(
      project_root, "results", "single_cell_downstream", "GSE154600", "cell_annotations.tsv.gz"
    )
    annotation <- data.table::fread(annotation_path, showProgress = FALSE)
    sce_paths <- sort(list.files(
      file.path(project_root, "results", "cache", "single_cell", "GSE154600"),
      pattern = "_qc_sce\\.rds$", full.names = TRUE
    ))
    for (index in seq_along(sce_paths)) {
      log_info("V3.1 localization sensitivity: GSE154600 sample ", index, "/", length(sce_paths), ".")
      sce <- readRDS(sce_paths[[index]])
      pass <- as.logical(SummarizedExperiment::colData(sce)$passes_QC)
      pass[is.na(pass)] <- FALSE
      raw_counts <- SummarizedExperiment::assay(sce, "counts")[, pass, drop = FALSE]
      counts <- .scd_collapse_gene_symbols(
        raw_counts,
        SummarizedExperiment::rowData(sce)$gene_symbol,
        SummarizedExperiment::rowData(sce)$gene_id
      )
      ann <- annotation[match(colnames(counts), cell_id)]
      metadata <- data.frame(
        cell_id = colnames(counts), sample = ann$batch, context = "HGSOC_all",
        cell_type = ann$cell_type, stringsAsFactors = FALSE
      )
      evidence_parts[[length(evidence_parts) + 1L]] <- v30_candidate_count_summary(
        counts, metadata, genes, "GSE154600", "OC", "all QC-pass tumor cells"
      )
      rm(sce, raw_counts, counts)
      invisible(gc())
    }

    local_config <- file.path(project_root, "config", "local.yml")
    config <- read_project_config(project_root, if (file.exists(local_config)) "config/local.yml" else "config/config.yml")
    dataset <- .scd_dataset_config(config, "GSE255460")
    source_metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
    bundle <- .sc_gse255460_ensure_sparse_bundle(dataset, config)
    validated <- .sc_gse255460_validate_manifest(bundle, source_metadata)
    annotation <- data.table::fread(file.path(
      project_root, "results", "single_cell_downstream", "GSE255460",
      "cell_annotations_all_QC_pass.tsv.gz"
    ), showProgress = FALSE)
    annotation <- annotation[trait == "OA"]
    partitions <- sort(unique(as.character(annotation$ID)))
    for (index in seq_along(partitions)) {
      partition_id <- partitions[[index]]
      log_info("V3.1 localization sensitivity: GSE255460 partition ", index, "/", length(partitions), ".")
      imported <- .sc_read_gse255460_partition(bundle, validated, partition_id, source_metadata)
      ann <- annotation[ID == partition_id]
      cell_index <- match(ann$cell_id, colnames(imported$counts))
      if (anyNA(cell_index)) stop("GSE255460 panel sensitivity annotation mismatch.", call. = FALSE)
      counts <- imported$counts[, cell_index, drop = FALSE]
      metadata <- data.frame(
        cell_id = colnames(counts), sample = ann$donor,
        context = paste0("OA:", ann$group), cell_type = ann$celltype,
        stringsAsFactors = FALSE
      )
      evidence_parts[[length(evidence_parts) + 1L]] <- v30_candidate_count_summary(
        counts, metadata, genes, "GSE255460", "OA", "all QC-pass OA cells"
      )
      rm(imported, counts)
      invisible(gc())
    }
    evidence <- data.table::rbindlist(evidence_parts, use.names = TRUE, fill = TRUE)
    saveRDS(evidence, cache_path, compress = "xz")
  }
  evidence <- data.table::as.data.table(evidence)
  evidence[, detected_cells := fraction_detected * cells]
  aggregated <- evidence[, .(
    cells = sum(cells),
    detected_cells = sum(detected_cells),
    samples = sum(samples)
  ), by = .(dataset_id, disease, cell_type, gene)]
  aggregated[, fraction_detected := detected_cells / cells]
  aggregated <- aggregated[
    cells >= 100L & !cell_type %in% c("Unassigned", "Ambiguous", "Other")
  ]
  panel_membership <- data.table::rbindlist(lapply(c(5L, 10L, 15L), function(size) {
    data.table::data.table(
      panel = paste0("Top ", size), panel_size = size,
      gene = composition$gene[composition$original_rank <= size]
    )
  }))
  scored <- merge(aggregated, panel_membership, by = "gene", allow.cartesian = TRUE)
  localization <- scored[, .(
    panel_mean_detection = mean(fraction_detected),
    panel_median_detection = stats::median(fraction_detected),
    cells = max(cells),
    measured_genes = .N
  ), by = .(dataset_id, disease, cell_type, panel, panel_size)]
  localization[, within_dataset_panel_rank := rank(-panel_mean_detection, ties.method = "min"),
               by = .(dataset_id, panel)]
  localization[, top_cell_type := within_dataset_panel_rank == 1L]
  safe_write_csv(as.data.frame(aggregated), file.path(paths$tables, "Table_S30a_extended_gene_cell_detection.csv"))
  safe_write_csv(as.data.frame(localization), file.path(paths$tables, "Table_S30b_panel_size_localization_sensitivity.csv"))

  correlation_rows <- list()
  counter <- 0L
  for (current_dataset in unique(localization$dataset_id)) {
    subset <- localization[dataset_id == current_dataset]
    for (pair in list(c("Top 5", "Top 10"), c("Top 5", "Top 15"), c("Top 10", "Top 15"))) {
      left <- subset[panel == pair[[1L]], .(cell_type, left = panel_mean_detection)]
      right <- subset[panel == pair[[2L]], .(cell_type, right = panel_mean_detection)]
      joined <- merge(left, right, by = "cell_type")
      counter <- counter + 1L
      correlation_rows[[counter]] <- data.frame(
        dataset_id = current_dataset,
        panel_left = pair[[1L]], panel_right = pair[[2L]],
        cell_types = nrow(joined),
        spearman_localization_profile = suppressWarnings(stats::cor(joined$left, joined$right, method = "spearman")),
        top_cell_type_left = left$cell_type[which.max(left$left)],
        top_cell_type_right = right$cell_type[which.max(right$right)],
        stringsAsFactors = FALSE
      )
    }
  }
  correlations <- do.call(rbind, correlation_rows)
  safe_write_csv(correlations, file.path(paths$tables, "Table_S30c_panel_size_localization_profile_comparison.csv"))
  list(evidence = aggregated, localization = localization, correlations = correlations)
}

v31_reference_theme <- function(base_size = 8.5) {
  submission_theme(base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#111111", size = base_size + 0.8),
      plot.subtitle = ggplot2::element_text(colour = "#4B5563", size = base_size - 0.6),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "#111111"),
      axis.line = ggplot2::element_line(colour = "#222222", linewidth = 0.35),
      panel.grid.major.y = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.24)
    )
}

v31_build_figure1 <- function(paths) {
  top <- data.frame(
    xmin = c(0.25, 4.15), xmax = c(2.05, 5.95), ymin = 3.25, ymax = 3.95,
    label = c("OA tissue\ntranscriptomes", "OC tissue\ntranscriptomes"),
    fill = c("#DDEBF7", "#F8D9D2"), stringsAsFactors = FALSE
  )
  shared <- data.frame(
    xmin = 2.25, xmax = 3.95, ymin = 3.25, ymax = 3.95,
    label = "Partial shared\ntranscriptome", fill = "#E9E2F3"
  )
  questions <- data.frame(
    xmin = c(0.10, 1.62, 3.14, 4.66), xmax = c(1.44, 2.96, 4.48, 6.00),
    ymin = 1.15, ymax = 2.55,
    question = c(
      "Q1\nAre alterations\nshared?", "Q2\nAre directions and\npathways concordant?",
      "Q3\nDo signals occupy\nsimilar cell contexts?", "Q4\nIs convergence\ngenetically mediated?"
    ),
    answer = c(
      "286 shared DEGs", "~49% discordant genes;\n6/10 Hallmarks opposite",
      "Dataset-specific\ncellular localization", "No liability-to-outcome\nevidence detected"
    ),
    fill = c("#E7F0F8", "#FCE8E3", "#EFE8F6", "#EEEEEE"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(2.05, 4.15, 3.10), xend = c(2.25, 3.95, 3.10),
    y = c(3.60, 3.60, 3.25), yend = c(3.60, 3.60, 2.70)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = top, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "#222222", linewidth = 0.5) +
    ggplot2::geom_rect(data = shared, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "#222222", linewidth = 0.55) +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                          arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"), colour = "#333333", linewidth = 0.55) +
    ggplot2::geom_rect(data = questions, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "#333333", linewidth = 0.45) +
    ggplot2::geom_text(data = top, ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), size = 3.5, fontface = "bold", lineheight = 0.95) +
    ggplot2::geom_text(data = shared, ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label), size = 3.5, fontface = "bold", lineheight = 0.95) +
    ggplot2::geom_text(data = questions, ggplot2::aes(x = (xmin + xmax) / 2, y = ymax - 0.20, label = question), vjust = 1, size = 2.85, fontface = "bold", lineheight = 0.92) +
    ggplot2::geom_text(data = questions, ggplot2::aes(x = (xmin + xmax) / 2, y = ymin + 0.24, label = answer), vjust = 0, size = 2.35, lineheight = 0.93, colour = "#333333") +
    ggplot2::annotate("text", x = 3.05, y = 0.72,
      label = "Interpretive boundary: shared expression membership is not shared mechanism, clinical utility, or causality",
      fontface = "italic", size = 2.85, colour = "#444444") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0, 6.1), ylim = c(0.45, 4.15), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 7, 7, 7))
  safe_write_csv(rbind(
    data.frame(element = "disease", label = top$label, result = NA_character_),
    data.frame(element = "shared", label = shared$label, result = NA_character_),
    data.frame(element = "question", label = questions$question, result = questions$answer)
  ), file.path(paths$source, "Figure1_question_driven_framework.csv"))
  submission_save_plot(p, "Figure1_question_driven_framework", paths$figures, 185, 118)
}

v31_volcano_plot <- function(data, disease, colour) {
  data$primary <- data$adj.P.Val < 0.05 & abs(data$logFC) >= 1
  data$display <- ifelse(data$primary, ifelse(data$logFC > 0, "Higher", "Lower"), "Not primary")
  ggplot2::ggplot(data, ggplot2::aes(logFC, -log10(pmax(adj.P.Val, .Machine$double.xmin)))) +
    ggplot2::geom_point(ggplot2::aes(alpha = primary), colour = colour, size = 0.75) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.16, `TRUE` = 0.78), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "#777777", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "#777777", linewidth = 0.3) +
    ggplot2::labs(title = paste0(disease, " discovery"), x = "log2 fold change", y = "-log10 FDR") +
    v31_reference_theme(7.5)
}

v31_build_figure2 <- function(paths, composition) {
  oa <- utils::read.csv(file.path(paths$source, "Figure2_OA_DEG.csv"))
  oc <- utils::read.csv(file.path(paths$source, "Figure2_OC_DEG.csv"))
  shared <- utils::read.csv(file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv"))
  shared <- shared[as.logical(shared$primary_shared), , drop = FALSE]
  shared$evidence_panel <- toupper(shared$gene) %in% composition$gene[composition$original_rank <= 10L]
  shared$direction <- ifelse(sign(shared$logFC_OA) == sign(shared$logFC_OC), "Concordant", "Discordant")
  p1 <- v31_volcano_plot(oa, "OA", "#0072B2")
  p2 <- v31_volcano_plot(oc, "OC", "#D55E00")
  p3 <- ggplot2::ggplot(shared, ggplot2::aes(logFC_OA, logFC_OC)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(fill = direction), shape = 21, colour = "white", stroke = 0.15, size = 2.1, alpha = 0.82) +
    ggplot2::geom_point(data = shared[shared$evidence_panel, ], shape = 21, fill = NA, colour = "#111111", stroke = 0.75, size = 3.0) +
    ggplot2::scale_fill_manual(values = c(Concordant = "#0072B2", Discordant = "#D55E00")) +
    ggplot2::labs(title = "Direction of 286 shared DEGs", subtitle = "Black rings: ten-gene evidence panel", x = "OA log2 fold change", y = "OC log2 fold change", fill = NULL) +
    v31_reference_theme(7.5) + ggplot2::theme(legend.position = "bottom")
  quadrants <- as.data.frame(table(shared$direction_quadrant), stringsAsFactors = FALSE)
  names(quadrants) <- c("direction_quadrant", "genes")
  quadrants$direction_quadrant <- factor(
    quadrants$direction_quadrant,
    levels = c(
      "OA higher / OC higher", "OA lower / OC lower",
      "OA higher / OC lower", "OA lower / OC higher"
    )
  )
  quadrants$direction <- ifelse(
    grepl("higher / OC higher|lower / OC lower", quadrants$direction_quadrant),
    "Concordant", "Discordant"
  )
  p4 <- ggplot2::ggplot(quadrants, ggplot2::aes(direction_quadrant, genes, fill = direction)) +
    ggplot2::geom_col(width = 0.68, colour = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = genes), vjust = -0.35, size = 2.8) +
    ggplot2::scale_fill_manual(values = c(Concordant = "#0072B2", Discordant = "#D55E00"), guide = "none") +
    ggplot2::labs(title = "Balanced directional structure", subtitle = "146 concordant; 140 discordant", x = NULL, y = "Shared DEGs") +
    v31_reference_theme(7.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4) + patchwork::plot_layout(heights = c(1, 1.08)))
  safe_write_csv(shared, file.path(paths$source, "Figure2_shared_DEG_direction.csv"))
  safe_write_csv(quadrants, file.path(paths$source, "Figure2_direction_quadrant_counts.csv"))
  submission_save_plot(figure, "Figure2_gene_and_direction_heterogeneity", paths$figures, 185, 177)
}

v31_build_panel_sensitivity_figure <- function(paths, panel, localization) {
  direction <- panel$direction
  direction$panel <- factor(direction$panel, levels = c("Top 5", "Top 10", "Top 15"))
  direction$direction_class <- factor(direction$direction_class, levels = c("concordant", "discordant"))
  p1 <- ggplot2::ggplot(direction, ggplot2::aes(panel, proportion, fill = direction_class)) +
    ggplot2::geom_col(width = 0.65, colour = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(genes, "/", panel_size)), position = ggplot2::position_stack(vjust = 0.5), size = 2.7, colour = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(concordant = "#0072B2", discordant = "#D55E00"), labels = c("Concordant", "Discordant")) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = "Direction composition", x = NULL, y = "Proportion", fill = NULL) +
    v31_reference_theme(7.2) + ggplot2::theme(legend.position = "bottom")

  ora <- panel$ora
  ora <- ora[ora$overlap_count > 0L, , drop = FALSE]
  top_ids <- unique(unlist(lapply(split(ora, ora$panel), function(x) {
    head(x$pathway_id[order(x$p_value, -x$fold_enrichment)], 6L)
  })))
  ora <- ora[ora$pathway_id %in% top_ids, , drop = FALSE]
  ora$pathway <- factor(ora$pathway, levels = rev(unique(ora$pathway[order(ora$p_value)])))
  ora$panel <- factor(ora$panel, levels = c("Top 5", "Top 10", "Top 15"))
  p2 <- ggplot2::ggplot(ora, ggplot2::aes(panel, pathway, size = overlap_count, fill = -log10(p_value))) +
    ggplot2::geom_point(shape = 21, colour = "#333333") +
    ggplot2::scale_fill_gradient(low = "#DDEBF7", high = "#0072B2") +
    ggplot2::labs(title = "Descriptive Hallmark profile", subtitle = "Shared-DEG background; no pathway-activity claim", x = NULL, y = NULL, size = "Genes", fill = "-log10 P") +
    v31_reference_theme(6.7) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.7))

  loc <- as.data.frame(localization$localization)
  loc$panel <- factor(loc$panel, levels = c("Top 5", "Top 10", "Top 15"))
  top_types <- unique(unlist(lapply(split(loc, loc$dataset_id), function(x) {
    score <- aggregate(panel_mean_detection ~ cell_type, x, max)
    head(score$cell_type[order(-score$panel_mean_detection)], 7L)
  })))
  loc <- loc[loc$cell_type %in% top_types, , drop = FALSE]
  loc$dataset_display <- ifelse(loc$dataset_id == "GSE255460", "OA atlas: GSE255460", "OC atlas: GSE154600")
  p3 <- ggplot2::ggplot(loc, ggplot2::aes(panel, cell_type, size = panel_mean_detection, fill = panel_mean_detection)) +
    ggplot2::geom_point(shape = 21, colour = "#333333") +
    ggplot2::facet_wrap(~dataset_display, scales = "free_y", nrow = 1) +
    ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#D55E00") +
    ggplot2::labs(title = "Detection-based cell-context localization", subtitle = "Within-atlas only; exact source labels retained", x = NULL, y = NULL, size = "Mean detection", fill = "Mean detection") +
    v31_reference_theme(6.8) + ggplot2::theme(legend.position = "bottom")
  figure <- submission_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(1.05, 1.0)))
  safe_write_csv(panel$direction, file.path(paths$source, "SupplementaryFigure15_direction_composition.csv"))
  safe_write_csv(panel$ora, file.path(paths$source, "SupplementaryFigure15_Hallmark_profile.csv"))
  safe_write_csv(as.data.frame(localization$localization), file.path(paths$source, "SupplementaryFigure15_cell_localization.csv"))
  submission_save_plot(figure, "SupplementaryFigure15_panel_size_sensitivity", paths$figures, 185, 180)
}

v31_bulk_qc_dataset <- function(dataset, dataset_label) {
  expression <- dataset$expression
  variances <- apply(expression, 1L, stats::var, na.rm = TRUE)
  top <- names(sort(variances[is.finite(variances)], decreasing = TRUE))[seq_len(min(1000L, sum(is.finite(variances))))]
  matrix <- expression[top, , drop = FALSE]
  pca <- stats::prcomp(t(matrix), center = TRUE, scale. = TRUE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  group <- as.character(dataset$group[colnames(expression)])
  scores <- data.frame(
    dataset = dataset_label, sample = colnames(expression), group = group,
    PC1 = pca$x[, 1L], PC2 = pca$x[, 2L],
    PC1_variance = variance[[1L]], PC2_variance = variance[[2L]],
    genes_used = length(top), stringsAsFactors = FALSE
  )
  correlation <- stats::cor(matrix, use = "pairwise.complete.obs", method = "pearson")
  diag(correlation) <- NA_real_
  sample_qc <- data.frame(
    dataset = dataset_label, sample = colnames(expression), group = group,
    median_sample_correlation = apply(correlation, 2L, stats::median, na.rm = TRUE),
    minimum_sample_correlation = apply(correlation, 2L, min, na.rm = TRUE),
    genes_used = length(top), stringsAsFactors = FALSE
  )
  sample_qc$lowest_three <- rank(sample_qc$median_sample_correlation, ties.method = "first") <= 3L
  list(scores = scores, sample_qc = sample_qc, correlation = correlation)
}

v31_build_bulk_qc <- function(project_root, paths) {
  cache <- readRDS(file.path(project_root, "results", "cache", "02_bulk_training.rds"))
  if (is.list(cache) && "value" %in% names(cache)) cache <- cache$value
  oa <- v31_bulk_qc_dataset(cache$oa_train, "OA discovery")
  oc <- v31_bulk_qc_dataset(cache$oc_train, "OC discovery")
  scores <- rbind(oa$scores, oc$scores)
  qc <- rbind(oa$sample_qc, oc$sample_qc)
  safe_write_csv(scores, file.path(paths$tables, "Table_S31a_bulk_unsupervised_PCA.csv"))
  safe_write_csv(qc, file.path(paths$tables, "Table_S31b_bulk_sample_correlation_QC.csv"))
  safe_write_csv(scores, file.path(paths$source, "SupplementaryFigure16_bulk_PCA.csv"))
  safe_write_csv(qc, file.path(paths$source, "SupplementaryFigure16_bulk_sample_correlation.csv"))
  make_pca <- function(data, title, colour) {
    xlab <- paste0("PC1 (", sprintf("%.1f", unique(data$PC1_variance)), "%)")
    ylab <- paste0("PC2 (", sprintf("%.1f", unique(data$PC2_variance)), "%)")
    ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, fill = group, shape = group)) +
      ggplot2::geom_point(size = 2.5, colour = "#222222", stroke = 0.35) +
      ggplot2::stat_ellipse(ggplot2::aes(colour = group), type = "norm", linewidth = 0.45, linetype = 2, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c(Normal = "white", Disease = colour)) +
      ggplot2::scale_colour_manual(values = c(Normal = "#777777", Disease = colour)) +
      ggplot2::scale_shape_manual(values = c(Normal = 24, Disease = 21)) +
      ggplot2::labs(title = title, subtitle = "Top 1,000 variable genes; labels used only for display", x = xlab, y = ylab, fill = NULL, shape = NULL) +
      v31_reference_theme(7.2) + ggplot2::theme(legend.position = "bottom")
  }
  make_qc <- function(data, title, colour) {
    data$sample_order <- seq_len(nrow(data))
    ggplot2::ggplot(data, ggplot2::aes(sample_order, median_sample_correlation, fill = group, shape = group)) +
      ggplot2::geom_point(size = 2.25, colour = "#222222", stroke = 0.3) +
      ggplot2::geom_text(data = data[data$lowest_three, ], ggplot2::aes(label = sample), vjust = -0.65, size = 2.0, check_overlap = TRUE) +
      ggplot2::scale_fill_manual(values = c(Normal = "white", Disease = colour)) +
      ggplot2::scale_shape_manual(values = c(Normal = 24, Disease = 21)) +
      ggplot2::labs(title = title, subtitle = "Median correlation to all other samples; three lowest labeled", x = "Sample index", y = "Median Pearson correlation", fill = NULL, shape = NULL) +
      v31_reference_theme(7.2) + ggplot2::theme(legend.position = "bottom")
  }
  p1 <- make_pca(oa$scores, "OA unsupervised PCA", "#0072B2")
  p2 <- make_pca(oc$scores, "OC unsupervised PCA", "#D55E00")
  p3 <- make_qc(oa$sample_qc, "OA sample-correlation QC", "#0072B2")
  p4 <- make_qc(oc$sample_qc, "OC sample-correlation QC", "#D55E00")
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(figure, "SupplementaryFigure16_bulk_PCA_and_QC", paths$figures, 185, 168)
  list(scores = scores, qc = qc)
}

v31_build_figure6 <- function(paths) {
  nodes <- data.frame(
    x = 1:4, y = 2.55,
    title = c("Gene layer", "Pathway layer", "Cellular layer", "Genetic-liability\nlayer"),
    result = c(
      "286 shared DEGs\n51% concordant",
      "10 shared Hallmarks\n6 opposite",
      "Five atlases\ncontext-specific localization",
      "Bidirectional MR\nno detected effect"
    ),
    fill = c("#DDEBF7", "#F8D9D2", "#E9E2F3", "#E8E8E8"), stringsAsFactors = FALSE
  )
  arrows <- data.frame(x = c(1.42, 2.42, 3.42), xend = c(1.58, 2.58, 3.58), y = 2.55, yend = 2.55)
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x = x, xend = xend, y = y, yend = yend), arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"), colour = "#333333", linewidth = 0.55) +
    ggplot2::geom_tile(data = nodes, ggplot2::aes(x = x, y = y, fill = fill), width = 0.84, height = 1.20, colour = "#222222", linewidth = 0.45) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y + 0.30, label = title), fontface = "bold", size = 2.75, lineheight = 0.92) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y - 0.20, label = result), size = 2.15, lineheight = 0.93, colour = "#333333") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.48, 4.52), ylim = c(1.75, 3.28), clip = "off") +
    ggplot2::labs(title = "Context-dependent convergence across evidence layers") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), plot.margin = ggplot2::margin(6, 7, 6, 7))
  boundary <- data.frame(
    evidence = factor(c("Shared genes", "Pathway overlap", "Cellular localization", "Null MR"), levels = rev(c("Shared genes", "Pathway overlap", "Cellular localization", "Null MR"))),
    supports = c("Partial molecular convergence", "Context-dependent pathway state", "Cellular contingency", "No detected disease-to-disease inherited effect"),
    does_not_support = c("Shared program", "Conserved mechanism", "Homologous cell state", "Absence of all genetic sharing"),
    stringsAsFactors = FALSE
  )
  boundary$index <- seq_len(nrow(boundary))
  p2 <- ggplot2::ggplot(boundary, ggplot2::aes(y = evidence)) +
    ggplot2::geom_tile(ggplot2::aes(x = 1, fill = "Supports"), width = 0.92, height = 0.74, colour = "white") +
    ggplot2::geom_tile(ggplot2::aes(x = 2, fill = "Does not support"), width = 0.92, height = 0.74, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = supports), size = 2.45, lineheight = 0.9) +
    ggplot2::geom_text(ggplot2::aes(x = 2, label = does_not_support), size = 2.45, lineheight = 0.9) +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("Supported interpretation", "Explicit boundary"), position = "top") +
    ggplot2::scale_fill_manual(values = c(Supports = "#DDEBF7", `Does not support` = "#F3F3F3"), guide = "none") +
    ggplot2::labs(title = "Evidence-to-claim boundary", x = NULL, y = NULL) +
    v31_reference_theme(7.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"), panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  figure <- submission_panel_tag(p1 / p2 + patchwork::plot_layout(heights = c(0.9, 1.15)))
  safe_write_csv(nodes, file.path(paths$source, "Figure6_integrated_evidence_layers.csv"))
  safe_write_csv(boundary, file.path(paths$source, "Figure6_evidence_claim_boundary.csv"))
  submission_save_plot(figure, "Figure6_integrated_context_model", paths$figures, 185, 155)
}

v31_update_documentation <- function(paths, panel, localization, bulk_qc) {
  legends <- c(
    "## Figure 1. Question-driven framework for context-dependent transcriptomic convergence",
    "",
    "OA and OC were analyzed on separate disease-specific tracks. Four prespecified questions distinguish shared transcriptional membership, direction and pathway state, cellular context, and disease-to-disease genetic liability. Each answer is paired with an explicit inference boundary.",
    "",
    "## Figure 2. Shared genes display an approximately balanced directional structure",
    "",
    "**A-B,** OA and OC discovery-cohort volcano plots. **C,** OA and OC log2 fold changes for 286 primary shared DEGs; black rings identify the fixed ten-gene evidence panel. **D,** Directional quadrant counts. The panel is an interpretable evidence summary, not an optimized predictive signature.",
    "",
    "## Figure 3. Quantitative pathway-direction heterogeneity",
    "",
    "**A,** OA and OC Hallmark normalized enrichment scores (NES). **B,** Paired-direction index for Hallmarks significant in both diseases. **C,** Leading-edge overlap. **D,** Complete direction classification. Pathway direction is descriptive and does not establish a conserved mechanism.",
    "",
    "## Figure 4. Within-atlas candidate detection and cell-context localization",
    "",
    "**A,** Candidate Cell Context Specificity Score (CCSS) with exact source labels. **B,** Sample-aware UCell ranks for the unsigned ten-gene evidence panel. **C,** Eligible sample-level pseudobulk effects. Scores are summarized within atlas and are not compared numerically across OA and OC.",
    "",
    "## Figure 5. Disease-context-dependent separability of the evidence panel",
    "",
    "**A,** Unsupervised PCA of GSE54388. **B,** Direction-fixed Hedges g estimates. **C,** Null AUC distributions from 1,000 label permutations. **D,** ROC curves shown last as a secondary display. These panels describe retrospective molecular separability across different comparator scales, not a universal diagnostic model.",
    "",
    "## Figure 6. Integrated context model and evidence-to-claim boundary",
    "",
    "**A,** Gene, pathway, cellular, and genetic-liability layers converge on a shared-but-nonidentical interpretation. **B,** Supported interpretations and their explicit boundaries. Detailed MR estimates and diagnostics remain in Figure S2.",
    "",
    "## Supplementary Figure 1. Prespecified non-MR sensitivity analyses",
    "",
    "DEG thresholds, WGCNA stability, machine-learning resampling, TCGA model sensitivity, and external-validation checks.",
    "",
    "## Supplementary Figure 2. Detailed bidirectional MR estimates and diagnostics",
    "",
    "Five MR estimators and heterogeneity/pleiotropy diagnostics in both directions. Null estimates mean no liability-to-outcome evidence was detected under the available instruments and assumptions, not proof of absence or a test of shared genetic architecture.",
    "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings",
    "",
    "Dataset-specific embeddings and exact source labels; OA and OC were not integrated into one latent space.",
    "",
    "## Supplementary Figure 4. HPA normal-tissue context",
    "",
    "Normal-tissue and normal-cell expression context; cartilage is not represented by the HPA/GTEx normal-tissue reference.",
    "",
    "## Supplementary Figure 5. TCGA-OV relative stromal and immune context audit",
    "",
    "Relative within-cohort transcriptomic context scores and candidate correlations; not absolute purity or histologic cell fractions.",
    "",
    "## Supplementary Figure 6. Exploratory cell-type functional annotation",
    "",
    "Marker-derived functional enrichment within eligible source-defined cell labels; descriptive only.",
    "",
    "## Supplementary Figure 7. Complete Hallmark direction and paired direction index",
    "",
    "Complete paired Hallmark NES matrix with matching and opposite directions identified explicitly.",
    "",
    "## Supplementary Figure 8. WGCNA and strict nested candidate-evidence stability",
    "",
    "Module stability, original LASSO/random-forest evidence, and strict outer-fold selection frequency. These layers prioritize evidence and assess stability; they do not define a clinical model.",
    "",
    "## Supplementary Figure 9. Candidate-centered Hallmark contexts",
    "",
    "Disease-status-adjusted residual association rankings for prespecified exemplars; not perturbation, mediation, or pathway activation.",
    "",
    "## Supplementary Figure 10. External-cohort design, sample influence, and calibration sensitivity",
    "",
    "Contrast schematic, leave-one-sample-out AUCs, and cross-fitted calibration. No nomogram or decision-curve analysis was performed because no locked clinical probability model or decision threshold exists.",
    "",
    "## Supplementary Figure 11. Focused upstream regulatory context",
    "",
    "KnockTF and miRTarBase evidence for prespecified exemplars; hypothesis generation only.",
    "",
    "## Supplementary Figure 12. Direction-aware STRING protein-association landscape",
    "",
    "High-confidence physical STRING graph, descriptive topology, direction-aware induced subgraphs, and 10,000 fixed-size label permutations. STRING is database-derived association context, not a tissue-specific interaction assay or mechanism.",
    "",
    "## Supplementary Figure 13. Sample-consensus CellChat communication context",
    "",
    "CellChat was fitted separately to biological samples and aggregated by context-specific consensus. OA and OC communication probabilities were not compared numerically.",
    "",
    "## Supplementary Figure 14. NicheNet prior-consistency overlay",
    "",
    "Consensus ligands projected onto the official NicheNet v2 human ligand-target prior for the fixed ten-gene panel. This is not ligand-activity inference, differential regulation, mediation, or causal signaling.",
    "",
    "## Supplementary Figure 15. Candidate-panel size sensitivity",
    "",
    "**A,** Direction composition of deterministic top-5, top-10, and top-15 extensions of the original model-vote ranking. **B,** Descriptive Hallmark over-representation profiles against the 286 shared-DEG background. **C,** Detection-based cell-context localization in the count-level GSE255460 OA and GSE154600 OC atlases. This analysis tests dependence on panel size; it does not optimize a signature or establish pathway activity.",
    "",
    "## Supplementary Figure 16. Discovery-cohort bulk PCA and sample-correlation QC",
    "",
    "**A-B,** Unsupervised PCA from the 1,000 most variable genes in each discovery cohort, with phenotype labels used only for display. **C-D,** Median correlation of each sample to all other samples, with the three lowest samples labeled for audit. OA and OC were analyzed separately."
  )
  writeLines(legends, file.path(paths$figures, "figure_legends.md"), useBytes = TRUE)

  index_path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  additions <- data.frame(
    table_id = c("Table S29a", "Table S29b", "Table S29c", "Table S29d", "Table S30a", "Table S30b", "Table S30c", "Table S31a", "Table S31b"),
    filename = c(
      "Table_S29a_panel_size_composition.csv", "Table_S29b_panel_size_direction_sensitivity.csv",
      "Table_S29c_panel_size_Hallmark_sensitivity.csv", "Table_S29d_panel_size_pathway_profile_comparison.csv",
      "Table_S30a_extended_gene_cell_detection.csv", "Table_S30b_panel_size_localization_sensitivity.csv",
      "Table_S30c_panel_size_localization_profile_comparison.csv", "Table_S31a_bulk_unsupervised_PCA.csv",
      "Table_S31b_bulk_sample_correlation_QC.csv"
    ),
    title = c(
      "Candidate-panel size composition", "Panel-size direction sensitivity", "Panel-size Hallmark sensitivity",
      "Panel-size Hallmark profile comparison", "Extended gene-cell detection", "Panel-size localization sensitivity",
      "Panel-size localization profile comparison", "Bulk unsupervised PCA", "Bulk sample-correlation QC"
    ),
    contents = c(
      "Deterministic top-5, top-10, and top-15 membership from the original model-vote ranking.",
      "Concordant/discordant composition and median absolute effects by panel size.",
      "Hallmark over-representation against the 286 shared-DEG background by panel size.",
      "Spearman and top-term Jaccard comparisons across panel sizes.",
      "Gene-level detection summaries in GSE154600 and GSE255460.",
      "Detection-based cell-context profiles for each panel size.",
      "Within-atlas cell-context profile correlations and top-label agreement.",
      "Unsupervised PCA source data for OA and OC discovery cohorts.",
      "Per-sample median and minimum correlation audit for each discovery cohort."
    ),
    source = c(
      "original four-model ranking", "Table S29a", "Hallmark sets and shared DEG background",
      "Table S29c", "count-level eligible single-cell atlases", "Table S30a", "Table S30b",
      "bulk discovery expression matrices", "top-1,000-variable-gene correlation matrices"
    ),
    stringsAsFactors = FALSE
  )
  index <- index[!sub("^Table ", "", index$table_id) %in% sub("^Table ", "", additions$table_id), , drop = FALSE]
  safe_write_csv(rbind(index, additions), index_path)

  registry_source <- file.path(paths$root, "claim_evidence_registry_v30.csv")
  registry_target <- file.path(paths$root, "claim_evidence_registry_v31.csv")
  registry <- utils::read.csv(registry_source, check.names = FALSE)
  new_claims <- data.frame(
    claim_id = c("C26", "C27"),
    manuscript_claim = c(
      "Direction and cell-context conclusions were audited under deterministic top-5, top-10, and top-15 evidence-panel definitions.",
      "Unsupervised PCA and sample-correlation summaries audited discovery-cohort structure without cross-disease integration."
    ),
    primary_data = c("original model ranking plus bulk and count-level single-cell data", "OA and OC bulk discovery matrices"),
    figure_or_table = c("Figure S15; Tables S29-S30", "Figure S16; Table S31"),
    allowed_wording = c("panel-size sensitivity; descriptive stability", "unsupervised structure; sample-correlation QC"),
    prohibited_wording = c("optimized signature; universal panel stability", "batch correction proves absence of confounding"),
    status = c("audited", "audited"),
    stringsAsFactors = FALSE
  )
  common <- intersect(names(registry), names(new_claims))
  registry <- registry[!registry$claim_id %in% new_claims$claim_id, , drop = FALSE]
  template <- registry[rep(1L, nrow(new_claims)), , drop = FALSE]
  template[,] <- NA
  template[, common] <- new_claims[, common]
  safe_write_csv(rbind(registry, template), registry_target)

  checklist <- data.frame(
    item = c(
      "Six question-driven main figures", "PPI retained as supplementary association landscape",
      "Candidate-panel size sensitivity", "Bulk discovery PCA/QC", "CellChat/NicheNet kept exploratory",
      "Clinical prediction template excluded", "Exact source labels retained"
    ),
    status = "complete",
    evidence = c(
      "Figures 1-6", "Figure S12", "Figure S15 and Tables S29-S30", "Figure S16 and Table S31",
      "Figures S13-S14 and Table S28", "No DCA or nomogram", "Figure 4 and Tables S24/S30"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(checklist, file.path(paths$root, "reproducibility_checklist_v31.csv"))
}

run_reviewer_v31 <- function(project_root) {
  v31_base_runner(project_root)
  paths <- v31_output_paths(project_root)
  v31_prepare_baseline(project_root, paths)
  log_info("V3.1: reconstructing deterministic candidate-panel ranks.")
  composition <- v31_panel_ranking(project_root, paths)
  panel <- v31_panel_direction_pathway(project_root, paths, composition)
  log_info("V3.1: auditing panel-size localization in eligible count-level atlases.")
  localization <- v31_panel_localization(project_root, paths, composition)
  log_info("V3.1: rebuilding question-driven main figures.")
  v31_build_figure1(paths)
  v31_build_figure2(paths, composition)
  v31_build_panel_sensitivity_figure(paths, panel, localization)
  log_info("V3.1: adding unsupervised bulk PCA and sample-correlation QC.")
  bulk_qc <- v31_build_bulk_qc(project_root, paths)
  v31_build_figure6(paths)
  v31_update_documentation(paths, panel, localization, bulk_qc)
  log_info("V3.1 refinement completed without adding a prediction algorithm or causal mechanism claim.")
  invisible(paths)
}
