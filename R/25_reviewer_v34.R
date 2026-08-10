v34_output_paths <- function(project_root) {
  root <- ensure_dir(file.path(project_root, "results", "submission_v34"))
  list(
    root = root,
    figures = ensure_dir(file.path(root, "figures")),
    source = ensure_dir(file.path(root, "figures", "source_data")),
    tables = ensure_dir(file.path(root, "supplementary_tables")),
    analysis = ensure_dir(file.path(root, "analysis")),
    logs = ensure_dir(file.path(root, "logs"))
  )
}

v34_copy <- function(source, target) {
  if (!file.exists(source)) stop("Missing V3.4 source: ", source, call. = FALSE)
  if (!file.copy(source, target, overwrite = TRUE)) stop("Could not copy: ", source, call. = FALSE)
  invisible(target)
}

v34_theme <- function(base_size = 6.6) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.4, colour = "#444444"),
      legend.title = ggplot2::element_text(size = base_size - 0.2, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 0.6),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

v34_panel_tag <- function(plot) {
  plot + patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(size = 8.2, face = "bold"))
}

v34_save_plot <- function(plot, stem, output_dir, width_mm = 170, height_mm = 116, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  tiff_path <- file.path(output_dir, paste0(stem, ".tiff"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "Arial")
  print(plot)
  grDevices::dev.off()
  svglite::svglite(svg_path, width = width, height = height)
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(tiff_path, width = width, height = height, units = "in", res = dpi, background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, svg = svg_path, png = png_path, tiff = tiff_path))
}

v34_prepare_tables <- function(project_root, paths) {
  stale_tables <- list.files(paths$tables, pattern = "\\.csv$", full.names = TRUE)
  if (length(stale_tables)) unlink(stale_tables)
  v33 <- file.path(project_root, "results", "submission_v33", "supplementary_tables")
  mapping <- c(
    "Table_S1_data_sources_and_cohorts.csv" = "Table_S1_tissue_data_sources_and_cohorts.csv",
    "Table_S2_shared_differentially_expressed_genes.csv" = "Table_S2_shared_tissue_DEGs.csv",
    "Table_S3a_DEG_threshold_summary.csv" = "Table_S3a_DEG_threshold_summary.csv",
    "Table_S3b_DEG_threshold_membership.csv" = "Table_S3b_DEG_threshold_membership.csv",
    "Table_S9_single_cell_QC_and_status.csv" = "Table_S8a_single_cell_QC_and_status.csv",
    "Table_S13_blood_dataset_audit.csv" = "Table_S9_blood_dataset_audit.csv",
    "Table_S14_FDR_supported_systemic_component.csv" = "Table_S10_FDR_supported_blood_component.csv",
    "Table_S15_blood_screen_attrition.csv" = "Table_S11_blood_screen_attrition.csv",
    "Table_S11b_GO_shared_genes.csv" = "Table_S12_GO_shared_genes.csv",
    "Table_S11c_KEGG_shared_genes.csv" = "Table_S13_KEGG_shared_genes.csv",
    "Table_S18_Hallmark_pathway_direction_matrix.csv" = "Table_S14_discovery_Hallmark_direction.csv",
    "Table_S4_WGCNA_stability.csv" = "Table_S15_WGCNA_stability.csv",
    "Table_S25a_STRING_mapping_audit.csv" = "Table_S16a_STRING_mapping_audit.csv",
    "Table_S25b_direction_aware_STRING_edges.csv" = "Table_S16b_STRING_edges.csv",
    "Table_S25c_STRING_node_topology.csv" = "Table_S16c_STRING_node_topology.csv"
  )
  for (source_name in names(mapping)) {
    v34_copy(file.path(v33, source_name), file.path(paths$tables, mapping[[source_name]]))
  }
  invisible(TRUE)
}

v34_fit_external_dataset <- function(dataset) {
  expression <- dataset$expression
  group <- droplevels(dataset$group)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(expression, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(Disease - Normal, levels = design))
  fit <- limma::eBayes(fit, trend = TRUE)
  table <- limma::topTable(fit, number = Inf, sort.by = "none")
  table$gene <- rownames(table)
  table$dataset_id <- dataset$id
  table$disease <- dataset$disease
  table
}

v34_external_validation <- function(project_root, paths) {
  set.seed(20260726L)
  validation <- readRDS(file.path(project_root, "results", "cache", "08_bulk_validation.rds"))$value
  shared <- readRDS(file.path(project_root, "results", "cache", "04_shared.rds"))$value$table
  result_list <- lapply(validation, v34_fit_external_dataset)
  names(result_list) <- vapply(validation, `[[`, character(1), "id")

  shared_results <- list()
  summary_rows <- list()
  for (id in names(result_list)) {
    full <- result_list[[id]]
    disease <- unique(full$disease)
    expected <- if (identical(disease, "OA")) shared$logFC_OA else shared$logFC_OC
    names(expected) <- shared$gene
    genes <- intersect(shared$gene, full$gene)
    table <- full[match(genes, full$gene), , drop = FALSE]
    table$discovery_logFC <- as.numeric(expected[table$gene])
    table$direction_match <- sign(table$logFC) == sign(table$discovery_logFC)
    table$external_FDR_lt_0_05 <- table$adj.P.Val < 0.05
    table$external_FDR_and_direction <- table$external_FDR_lt_0_05 & table$direction_match
    shared_results[[id]] <- table
    binomial <- stats::binom.test(sum(table$direction_match), nrow(table), p = 0.5)
    correlation <- suppressWarnings(stats::cor.test(table$logFC, table$discovery_logFC, method = "spearman", exact = FALSE))
    summary_rows[[id]] <- data.frame(
      dataset_id = id,
      disease = disease,
      measured_shared_genes = nrow(table),
      direction_matches = sum(table$direction_match),
      direction_agreement = mean(table$direction_match),
      binomial_P = binomial$p.value,
      external_FDR_genes = sum(table$external_FDR_lt_0_05),
      external_FDR_direction_matches = sum(table$external_FDR_and_direction),
      spearman_rho = unname(correlation$estimate),
      spearman_P = correlation$p.value,
      stringsAsFactors = FALSE
    )
  }
  gene_table <- do.call(rbind, shared_results)
  rownames(gene_table) <- NULL
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  safe_write_csv(gene_table, file.path(paths$tables, "Table_S4_external_tissue_gene_effects.csv"))
  safe_write_csv(summary_table, file.path(paths$tables, "Table_S5_external_tissue_direction_summary.csv"))

  local_config <- file.path(project_root, "config", "local.yml")
  project_config <- read_project_config(
    project_root,
    if (file.exists(local_config)) "config/local.yml" else "config/config.yml"
  )
  term_to_gene <- clusterProfiler::read.gmt(project_config$gene_sets$hallmark)
  gsea_rows <- list()
  for (id in names(result_list)) {
    ranked <- result_list[[id]]$logFC
    names(ranked) <- result_list[[id]]$gene
    ranked <- sort(ranked[is.finite(ranked) & !duplicated(names(ranked))], decreasing = TRUE)
    enrichment <- suppressWarnings(clusterProfiler::GSEA(
      geneList = ranked,
      TERM2GENE = term_to_gene,
      pvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    ))
    table <- as.data.frame(enrichment)
    table$dataset_id <- id
    table$disease <- unique(result_list[[id]]$disease)
    gsea_rows[[id]] <- table
  }
  external_gsea <- do.call(rbind, gsea_rows)
  rownames(external_gsea) <- NULL
  safe_write_csv(external_gsea, file.path(paths$tables, "Table_S6_external_tissue_Hallmark_GSEA.csv"))
  invisible(list(results = result_list, shared = shared_results, summary = summary_table, gsea = external_gsea, discovery = shared))
}

v34_sc_aggregate <- function(counts, symbols, cell_types, genes, dataset_id, disease) {
  output <- list()
  for (gene in genes) {
    rows <- which(toupper(symbols) == toupper(gene))
    if (!length(rows)) next
    gene_counts <- if (length(rows) == 1L) {
      as.numeric(counts[rows, ])
    } else {
      as.numeric(Matrix::colSums(counts[rows, , drop = FALSE]))
    }
    for (cell_type in unique(cell_types)) {
      index <- which(cell_types == cell_type)
      output[[length(output) + 1L]] <- data.frame(
        dataset_id = dataset_id,
        disease = disease,
        cell_type = cell_type,
        gene = gene,
        cells = length(index),
        total_counts = sum(gene_counts[index]),
        detected_cells = sum(gene_counts[index] > 0),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, output)
}

v34_single_cell_candidates <- function(project_root, paths, genes) {
  oa_annotations <- utils::read.delim(
    gzfile(file.path(project_root, "results", "single_cell_downstream", "GSE255460", "cell_annotations_all_QC_pass.tsv.gz")),
    check.names = FALSE
  )
  local_config <- file.path(project_root, "config", "local.yml")
  project_config <- read_project_config(
    project_root,
    if (file.exists(local_config)) "config/local.yml" else "config/config.yml"
  )
  dataset_ids <- vapply(project_config$single_cell$datasets, function(x) x$id, character(1))
  oa_dataset <- project_config$single_cell$datasets[[match("GSE255460", dataset_ids)]]
  oa_metadata_path <- oa_dataset$metadata_path
  oa_metadata <- .oa_sc_read_gse255460_metadata(oa_metadata_path)
  backed <- readRDS(file.path(project_root, "results", "single_cell", "GSE255460", "backed_count_manifest.rds"))
  bundle <- list(directory = backed$bundle_directory, manifest_path = backed$manifest_path)
  validated <- .sc_gse255460_validate_manifest(bundle, oa_metadata)
  oa_parts <- list()
  for (partition_id in validated$partitions$partition_id) {
    imported <- .sc_read_gse255460_partition(bundle, validated, partition_id, oa_metadata)
    annotation <- oa_annotations[oa_annotations$batch == partition_id, , drop = FALSE]
    index <- match(annotation$cell_id, colnames(imported$counts))
    if (anyNA(index)) stop("OA single-cell annotation mapping failed: ", partition_id, call. = FALSE)
    oa_parts[[partition_id]] <- v34_sc_aggregate(
      imported$counts[, index, drop = FALSE], imported$gene_id,
      annotation$celltype, genes, "GSE255460", "OA"
    )
    rm(imported)
    gc(FALSE)
  }
  oa <- do.call(rbind, oa_parts)
  oa <- stats::aggregate(cbind(cells, total_counts, detected_cells) ~ dataset_id + disease + cell_type + gene, data = oa, sum)

  oc_annotations <- utils::read.delim(
    gzfile(file.path(project_root, "results", "single_cell_downstream", "GSE154600", "cell_annotations.tsv.gz")),
    check.names = FALSE
  )
  sce_files <- list.files(
    file.path(project_root, "results", "cache", "single_cell", "GSE154600"),
    pattern = "qc_sce[.]rds$", full.names = TRUE
  )
  oc_parts <- list()
  for (sce_path in sce_files) {
    sce <- readRDS(sce_path)
    batch <- sub("_qc_sce[.]rds$", "", basename(sce_path))
    annotation <- oc_annotations[oc_annotations$batch == batch, , drop = FALSE]
    index <- match(annotation$cell_id, colnames(sce))
    if (anyNA(index)) stop("OC single-cell annotation mapping failed: ", batch, call. = FALSE)
    oc_parts[[batch]] <- v34_sc_aggregate(
      SummarizedExperiment::assay(sce, "counts")[, index, drop = FALSE],
      as.character(SummarizedExperiment::rowData(sce)$gene_symbol),
      annotation$cell_type, genes, "GSE154600", "OC"
    )
    rm(sce)
    gc(FALSE)
  }
  oc <- do.call(rbind, oc_parts)
  oc <- stats::aggregate(cbind(cells, total_counts, detected_cells) ~ dataset_id + disease + cell_type + gene, data = oc, sum)
  result <- rbind(oa, oc)
  result$fraction_detected <- result$detected_cells / result$cells
  result$mean_umi_per_cell <- result$total_counts / result$cells
  safe_write_csv(result, file.path(paths$tables, "Table_S8b_representative_gene_single_cell_detection.csv"))
  invisible(result)
}

v34_candidate_evidence <- function(paths, external, single_cell, genes) {
  shared <- external$discovery
  rows <- list()
  theme <- c(
    G0S2 = "Dual-blood-FDR candidate",
    EFEMP1 = "Extracellular matrix",
    AKAP12 = "Cell-context localization",
    SOX9 = "Tissue differentiation",
    DDIT3 = "Integrated stress response"
  )
  for (gene in genes) {
    direction_matches <- vapply(external$shared, function(x) {
      row <- x[x$gene == gene, , drop = FALSE]
      if (!nrow(row)) return(NA)
      isTRUE(row$direction_match[[1L]])
    }, logical(1))
    fdr_hits <- vapply(external$shared, function(x) {
      row <- x[x$gene == gene, , drop = FALSE]
      if (!nrow(row)) return(NA)
      isTRUE(row$external_FDR_lt_0_05[[1L]])
    }, logical(1))
    top_cell <- function(disease) {
      table <- single_cell[single_cell$gene == gene & single_cell$disease == disease, , drop = FALSE]
      if (!nrow(table)) return(NA_character_)
      table$cell_type[[which.max(table$fraction_detected)]]
    }
    shared_row <- shared[shared$gene == gene, , drop = FALSE]
    rows[[gene]] <- data.frame(
      gene = gene,
      OA_discovery_logFC = shared_row$logFC_OA,
      OC_discovery_logFC = shared_row$logFC_OC,
      direction_class = ifelse(shared_row$directionally_concordant, "Concordant", "Discordant"),
      external_direction_matches = sum(direction_matches, na.rm = TRUE),
      external_measured_cohorts = sum(!is.na(direction_matches)),
      external_FDR_hits = sum(fdr_hits, na.rm = TRUE),
      OA_top_cell_label = top_cell("OA"),
      OC_top_cell_label = top_cell("OC"),
      evidence_role = unname(theme[[gene]]),
      dual_blood_FDR = gene == "G0S2",
      interpretation = "Representative evidence summary; not an optimized signature",
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  safe_write_csv(result, file.path(paths$tables, "Table_S7_candidate_evidence_summary.csv"))
  result
}

v34_build_figure1 <- function(paths) {
  nodes <- data.frame(
    x = c(1, 3, 2, 2, 2, 2, 2),
    y = c(5.6, 5.6, 4.6, 3.6, 2.6, 1.6, 0.6),
    width = c(1.45, 1.45, 1.85, 1.85, 1.85, 1.85, 1.85),
    title = c("OA tissue", "OC tissue", "Shared DEGs", "Biological themes", "Representative genes", "Separate single-cell atlases", "Peripheral-blood filter"),
    detail = c("GSE114007", "GSE18520", "286 genes", "GO, KEGG, Hallmark", "transparent evidence summary", "exact source labels", "G0S2 retained"),
    fill = c("#DDEAF4", "#F8E1D8", "#ECE7F2", "#E8F2EC", "#F2F2F2", "#E7F1EF", "#FBE8D6"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(1.42, 2.58, 2, 2, 2, 2), y = c(5.35, 5.35, 4.26, 3.26, 2.26, 1.26),
    xend = c(1.82, 2.18, 2, 2, 2, 2), yend = c(4.92, 4.92, 3.93, 2.93, 1.93, 0.93)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
      linewidth = 0.45, colour = "#4A4A4A"
    ) +
    ggplot2::geom_rect(
      data = nodes,
      ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - 0.32, ymax = y + 0.32, fill = fill),
      linewidth = 0.4, colour = "#3C3C3C"
    ) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y + 0.08, label = title), size = 2.55, fontface = "bold") +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y - 0.13, label = detail), size = 1.95, colour = "#59636E") +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.05, 3.95), ylim = c(0.15, 6.05), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 8, 5, 8))
  safe_write_csv(nodes, file.path(paths$source, "Figure1_workflow_nodes.csv"))
  v34_save_plot(p, "Figure1_study_design", paths$figures, 170, 110)
}

v34_volcano <- function(table, disease) {
  table$class <- "Not primary"
  table$class[table$adj.P.Val < 0.05 & table$logFC <= -1] <- "Lower"
  table$class[table$adj.P.Val < 0.05 & table$logFC >= 1] <- "Higher"
  table$minus_log10_FDR <- -log10(pmax(table$adj.P.Val, 1e-300))
  ggplot2::ggplot(table, ggplot2::aes(logFC, minus_log10_FDR, colour = class)) +
    ggplot2::geom_point(size = 0.52, alpha = 0.68) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = "#777777") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "#777777") +
    ggplot2::scale_colour_manual(values = c("Not primary" = "#D8DDE2", Lower = "#2C7FB8", Higher = "#D95F0E")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 FDR", colour = NULL) +
    v34_theme(6.3) +
    ggplot2::theme(legend.position = "bottom")
}

v34_build_figure2 <- function(project_root, paths) {
  oa <- utils::read.csv(file.path(project_root, "results", "tables", "DEG_OA_GSE114007_all.csv"))
  oc <- utils::read.csv(file.path(project_root, "results", "tables", "DEG_OC_GSE18520_all.csv"))
  shared <- utils::read.csv(file.path(paths$tables, "Table_S2_shared_tissue_DEGs.csv"))
  p1 <- v34_volcano(oa, "OA")
  p2 <- v34_volcano(oc, "OC")

  angle <- seq(0, 2 * pi, length.out = 300)
  circles <- rbind(
    data.frame(x = 0.90 + 0.78 * cos(angle), y = 1 + 0.68 * sin(angle), disease = "OA"),
    data.frame(x = 1.60 + 0.78 * cos(angle), y = 1 + 0.68 * sin(angle), disease = "OC")
  )
  p3 <- ggplot2::ggplot(circles, ggplot2::aes(x, y, group = disease, fill = disease)) +
    ggplot2::geom_polygon(alpha = 0.42, colour = "#333333", linewidth = 0.4) +
    ggplot2::annotate("text", x = c(0.50, 1.25, 2.00), y = 1, label = c("1,722", "286", "2,024"), size = 3.0, fontface = "bold") +
    ggplot2::annotate("text", x = c(0.60, 1.90), y = 1.77, label = c("OA DEGs\n2,008", "OC DEGs\n2,310"), size = 2.5, fontface = "bold") +
    ggplot2::annotate("text", x = 1.25, y = 0.18, label = "286 shared genes", size = 2.35, colour = "#59636E") +
    ggplot2::scale_fill_manual(values = c(OA = "#79ADD2", OC = "#E8A179"), guide = "none") +
    ggplot2::coord_equal(xlim = c(0, 2.5), ylim = c(0.05, 1.95), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))

  shared$direction_class <- ifelse(shared$directionally_concordant %in% c(TRUE, "TRUE"), "Concordant", "Discordant")
  shared$quadrant <- paste0(ifelse(shared$logFC_OA > 0, "OA+", "OA-"), " / ", ifelse(shared$logFC_OC > 0, "OC+", "OC-"))
  shared <- shared[order(shared$quadrant, shared$logFC_OA), ]
  shared$row <- seq_len(nrow(shared))
  heat <- rbind(
    data.frame(row = shared$row, disease = "OA", log2FC = shared$logFC_OA),
    data.frame(row = shared$row, disease = "OC", log2FC = shared$logFC_OC)
  )
  p4 <- ggplot2::ggplot(heat, ggplot2::aes(disease, row, fill = log2FC)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = c(34.5, 88.5, 174.5), linewidth = 0.45, colour = "#333333") +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-5.5, 5.5), oob = scales::squish) +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = NULL, y = NULL, fill = "log2FC") +
    v34_theme(6.3) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right")

  safe_write_csv(shared, file.path(paths$source, "Figure2_shared_genes_ordered.csv"))
  safe_write_csv(oa, file.path(paths$source, "Figure2_OA_differential.csv"))
  safe_write_csv(oc, file.path(paths$source, "Figure2_OC_differential.csv"))
  top_row <- (p1 | p2) + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  figure <- v34_panel_tag(top_row / (p3 | p4) + patchwork::plot_layout(heights = c(1.0, 0.88)))
  v34_save_plot(figure, "Figure2_shared_tissue_transcriptomics", paths$figures, 170, 147)
}

v34_build_figure3 <- function(project_root, paths) {
  go <- utils::read.csv(file.path(project_root, "results", "submission_v32", "figures", "source_data", "Figure3_representative_GO_terms.csv"))
  go$theme <- c("Cell cycle", "Cell cycle", "Cell cycle", "Matrix", "Matrix", "Immune", "Immune", "Stress", "Stress", "Stress")
  go$label <- v32_wrap(go$Description, 32L)
  go$label <- factor(go$label, levels = rev(go$label))
  p1 <- ggplot2::ggplot(go, ggplot2::aes(FoldEnrichment, label, size = Count, colour = theme)) +
    ggplot2::geom_point(alpha = 0.95, stroke = 0.35) +
    ggplot2::scale_colour_manual(values = c("Cell cycle" = "#6A51A3", Matrix = "#2C7FB8", Immune = "#238B8E", Stress = "#D95F0E")) +
    ggplot2::scale_size(range = c(1.6, 5.0), breaks = c(10, 20, 30)) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, colour = NULL, size = "Genes") +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, nrow = 1, title.position = "left"),
      colour = ggplot2::guide_legend(order = 2, nrow = 1, title.position = "left", override.aes = list(size = 2.2))
    ) +
    v34_theme(6.4) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0.3, "mm"),
      axis.text.y = ggplot2::element_text(size = 5.5)
    )

  hallmark <- utils::read.csv(file.path(paths$tables, "Table_S14_discovery_Hallmark_direction.csv"), check.names = FALSE)
  hallmark <- hallmark[hallmark$both_significant %in% c(TRUE, "TRUE"), , drop = FALSE]
  hallmark$pathway_label <- factor(v32_wrap(hallmark$pathway, 25L), levels = rev(v32_wrap(hallmark$pathway, 25L)))
  hallmark$direction_label <- ifelse(hallmark$direction_class == "concordant", "Same sign", "Opposite sign")
  p2 <- ggplot2::ggplot(hallmark, ggplot2::aes(y = pathway_label)) +
    ggplot2::geom_segment(ggplot2::aes(x = OA_NES, xend = OC_NES, yend = pathway_label, colour = direction_label), linewidth = 0.65, alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = OA_NES), shape = 16, size = 2.0, colour = "#2C7FB8") +
    ggplot2::geom_point(ggplot2::aes(x = OC_NES), shape = 17, size = 2.2, colour = "#D95F0E") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
    ggplot2::scale_colour_manual(values = c("Same sign" = "#7A9E9F", "Opposite sign" = "#A98274")) +
    ggplot2::labs(x = "Normalized enrichment score", y = NULL, colour = NULL) +
    v34_theme(6.4) +
    ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 5.3))

  safe_write_csv(go, file.path(paths$source, "Figure3_GO_themes.csv"))
  safe_write_csv(hallmark, file.path(paths$source, "Figure3_joint_Hallmark.csv"))
  figure <- v34_panel_tag((p1 | p2) + patchwork::plot_layout(widths = c(1.03, 0.97)))
  v34_save_plot(figure, "Figure3_shared_biological_themes", paths$figures, 170, 104)
}

v34_build_figure4 <- function(paths, external, candidates) {
  summary <- external$summary
  summary$label <- paste0(summary$dataset_id, "\n", summary$direction_matches, "/", summary$measured_shared_genes)
  summary$disease <- factor(summary$disease, levels = c("OA", "OC"))
  summary$dataset_id <- factor(summary$dataset_id, levels = c("GSE117999", "GSE82107", "GSE54388", "GSE12470"))
  p1 <- ggplot2::ggplot(summary, ggplot2::aes(dataset_id, 100 * direction_agreement, fill = disease)) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
    ggplot2::geom_col(width = 0.64, colour = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", 100 * direction_agreement)), vjust = -0.35, size = 2.25) +
    ggplot2::scale_fill_manual(values = c(OA = "#79ADD2", OC = "#E8A179"), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 101), clip = "off") +
    ggplot2::labs(x = NULL, y = "Shared genes with matching direction (%)", fill = NULL) +
    v34_theme(6.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), legend.position = "none")

  contexts <- c("OA discovery", "GSE117999", "GSE82107", "OC discovery", "GSE54388", "GSE12470")
  heat_rows <- list()
  for (gene in candidates$gene) {
    discovery <- external$discovery[external$discovery$gene == gene, , drop = FALSE]
    heat_rows[[gene]] <- data.frame(
      gene = gene,
      context = contexts,
      log2FC = c(
        discovery$logFC_OA,
        external$results$GSE117999$logFC[match(gene, external$results$GSE117999$gene)],
        external$results$GSE82107$logFC[match(gene, external$results$GSE82107$gene)],
        discovery$logFC_OC,
        external$results$GSE54388$logFC[match(gene, external$results$GSE54388$gene)],
        external$results$GSE12470$logFC[match(gene, external$results$GSE12470$gene)]
      ),
      stringsAsFactors = FALSE
    )
  }
  heat <- do.call(rbind, heat_rows)
  heat$gene <- factor(heat$gene, levels = rev(candidates$gene))
  heat$context <- factor(heat$context, levels = contexts)
  p2 <- ggplot2::ggplot(heat, ggplot2::aes(context, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.finite(log2FC), sprintf("%.1f", log2FC), "NA")), size = 1.85) +
    ggplot2::geom_vline(xintercept = 3.5, linewidth = 0.55, colour = "#333333") +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-4.5, 4.5), oob = scales::squish, na.value = "#EEEEEE") +
    ggplot2::labs(x = NULL, y = NULL, fill = "log2FC") +
    v34_theme(6.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), panel.grid = ggplot2::element_blank(), legend.position = "right")

  evidence <- rbind(
    data.frame(gene = candidates$gene, evidence = "External direction", value = candidates$external_direction_matches / candidates$external_measured_cohorts, label = paste0(candidates$external_direction_matches, "/", candidates$external_measured_cohorts)),
    data.frame(gene = candidates$gene, evidence = "Single-cell localization", value = 1, label = "OA + OC"),
    data.frame(gene = candidates$gene, evidence = "Dual-blood FDR", value = ifelse(candidates$dual_blood_FDR, 1, 0), label = ifelse(candidates$dual_blood_FDR, "Yes", "No")),
    data.frame(gene = candidates$gene, evidence = "Interpretive role", value = 0.65, label = candidates$evidence_role)
  )
  evidence$gene <- factor(evidence$gene, levels = rev(candidates$gene))
  evidence$evidence <- factor(evidence$evidence, levels = c("External direction", "Single-cell localization", "Dual-blood FDR", "Interpretive role"))
  p3 <- ggplot2::ggplot(evidence, ggplot2::aes(evidence, gene, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = v32_wrap(label, 21L)), size = 1.85, lineheight = 0.9) +
    ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#5B8FA8", limits = c(0, 1), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    v34_theme(6.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), panel.grid = ggplot2::element_blank())

  safe_write_csv(summary, file.path(paths$source, "Figure4_external_direction_summary.csv"))
  safe_write_csv(heat, file.path(paths$source, "Figure4_representative_gene_effects.csv"))
  safe_write_csv(evidence, file.path(paths$source, "Figure4_candidate_evidence_matrix.csv"))
  figure <- v34_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(1, 0.88)))
  v34_save_plot(figure, "Figure4_external_tissue_and_representative_genes", paths$figures, 170, 138)
}

v34_build_figure5 <- function(project_root, paths, single_cell, genes) {
  umap <- utils::read.csv(file.path(project_root, "results", "submission_v32", "figures", "source_data", "Figure6_representative_UMAPs.csv"))
  selected <- list(
    GSE255460 = c("EC", "HomC", "HTC", "preHTC", "ProC", "RepC"),
    GSE154600 = c("B.cell", "Endothelial.cell", "Fibroblast", "Myeloid.cell", "Ovarian.cancer.cell", "T.cell")
  )
  palette <- c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377")
  build_umap <- function(dataset_id) {
    table <- umap[umap$dataset_id == dataset_id, , drop = FALSE]
    table$display_label <- ifelse(table$label %in% selected[[dataset_id]], table$label, "Other labels")
    levels <- c(selected[[dataset_id]], "Other labels")
    table$display_label <- factor(table$display_label, levels = levels)
    colours <- c(stats::setNames(palette, selected[[dataset_id]]), "Other labels" = "#D9D9D9")
    centres <- stats::aggregate(cbind(UMAP1, UMAP2) ~ display_label, data = table[table$display_label != "Other labels", ], median)
    base <- ggplot2::ggplot(table, ggplot2::aes(UMAP1, UMAP2, colour = display_label)) +
      ggplot2::geom_point(size = 0.15, alpha = 0.65) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE, guide = "none") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
      ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
      ggplot2::coord_equal(clip = "off") +
      v34_theme(6.2) +
      ggrepel::geom_text_repel(
        data = centres,
        ggplot2::aes(label = display_label),
        size = 2.0,
        colour = "#303030",
        seed = 20260810,
        box.padding = 0.4,
        point.padding = 0.15,
        force = 2,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.size = 0.18
      )
    base
  }
  p1 <- build_umap("GSE255460")
  p2 <- build_umap("GSE154600")

  selected_rows <- rbind(
    single_cell[single_cell$dataset_id == "GSE255460" & single_cell$cell_type %in% selected$GSE255460, ],
    single_cell[single_cell$dataset_id == "GSE154600" & single_cell$cell_type %in% selected$GSE154600, ]
  )
  selected_rows$gene <- factor(selected_rows$gene, levels = rev(genes))
  selected_rows$cell_label <- ifelse(selected_rows$dataset_id == "GSE255460", paste0("OA: ", selected_rows$cell_type), paste0("OC: ", selected_rows$cell_type))
  selected_rows$cell_label <- factor(selected_rows$cell_label, levels = c(paste0("OA: ", selected$GSE255460), paste0("OC: ", selected$GSE154600)))
  p3 <- ggplot2::ggplot(selected_rows, ggplot2::aes(cell_label, gene, size = fraction_detected, fill = fraction_detected)) +
    ggplot2::geom_point(shape = 21, colour = "#555555", stroke = 0.3) +
    ggplot2::geom_vline(xintercept = 6.5, linewidth = 0.45, colour = "#777777") +
    ggplot2::scale_size(
      name = "Detected fraction",
      range = c(0.8, 5.3),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00)
    ) +
    ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#D95F0E", limits = c(0, 1), guide = "none") +
    ggplot2::scale_x_discrete(labels = function(x) sub(": ", ":\n", x, fixed = TRUE)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(size = ggplot2::guide_legend(nrow = 1, title.position = "left", override.aes = list(fill = "#F4B183"))) +
    v34_theme(6.2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 32, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "italic"),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid = ggplot2::element_blank()
    )
  safe_write_csv(selected_rows, file.path(paths$source, "Figure5_single_cell_candidate_detection.csv"))
  safe_write_csv(umap, file.path(paths$source, "Figure5_representative_UMAPs.csv"))
  figure <- v34_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(1, 0.80)))
  v34_save_plot(figure, "Figure5_single_cell_localization", paths$figures, 170, 134)
}

v34_build_figure6 <- function(project_root, paths) {
  blood <- file.path(project_root, "results", "blood_v33_internal")
  attrition <- utils::read.csv(file.path(blood, "blood_screen_attrition.csv"))
  positive <- utils::read.csv(file.path(blood, "FDR_supported_systemic_component_genes.csv"))
  stage_labels <- c(
    "Shared tissue DEGs" = "Shared tissue DEGs",
    "Tissue-concordant DEGs" = "Tissue-concordant DEGs",
    "Measured in both blood cohorts" = "Measured in both\nblood cohorts",
    "All four effects same direction" = "All four effects in\nthe same direction",
    "Nominal P<0.05 in both blood cohorts" = "Nominal P < 0.05 in\nboth blood cohorts",
    "FDR<0.05 in both blood cohorts" = "FDR < 0.05 in\nboth blood cohorts"
  )
  attrition$stage <- factor(
    unname(stage_labels[as.character(attrition$stage)]),
    levels = rev(unname(stage_labels))
  )
  attrition$highlight <- ifelse(attrition$genes == min(attrition$genes), "Retained", "Screened")
  p1 <- ggplot2::ggplot(attrition, ggplot2::aes(genes, stage, fill = highlight)) +
    ggplot2::geom_col(width = 0.62, colour = "#333333", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = genes), hjust = -0.22, size = 2.2) +
    ggplot2::scale_fill_manual(values = c(Screened = "#B8C4D0", Retained = "#D95F0E"), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, max(attrition$genes) * 1.16), clip = "off") +
    ggplot2::labs(x = "Genes", y = NULL) +
    v34_theme(6.3) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.3))

  heat <- data.frame(
    context = factor(c("OA tissue", "OC tissue", "OA blood", "OC blood"), levels = c("OA tissue", "OC tissue", "OA blood", "OC blood")),
    gene = "G0S2",
    log2FC = c(positive$logFC_OA, positive$logFC_OC, positive$OA_blood_logFC, positive$OC_blood_logFC)
  )
  p2 <- ggplot2::ggplot(heat, ggplot2::aes(context, gene, fill = log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", log2FC)), size = 2.4) +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-2.2, 2.2), oob = scales::squish) +
    ggplot2::labs(x = NULL, y = NULL, fill = "log2FC") +
    v34_theme(6.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1), panel.grid = ggplot2::element_blank())

  forest <- data.frame(
    cohort = factor(c("OA blood: GSE48556", "OC blood: GSE31682"), levels = rev(c("OA blood: GSE48556", "OC blood: GSE31682"))),
    effect = c(positive$OA_blood_hedges_g, positive$OC_blood_hedges_g),
    lower = c(positive$OA_blood_hedges_g_lower, positive$OC_blood_hedges_g_lower),
    upper = c(positive$OA_blood_hedges_g_upper, positive$OC_blood_hedges_g_upper),
    disease = c("OA", "OC")
  )
  p3 <- ggplot2::ggplot(forest, ggplot2::aes(effect, cohort, colour = disease, shape = disease)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper), orientation = "y", width = 0.12, linewidth = 0.55) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(values = c(OA = "#2C7FB8", OC = "#D95F0E"), guide = "none") +
    ggplot2::scale_shape_manual(values = c(OA = 16, OC = 17), guide = "none") +
    ggplot2::labs(x = "Hedges g", y = NULL) +
    v34_theme(6.3)

  safe_write_csv(attrition, file.path(paths$source, "Figure6_blood_attrition.csv"))
  safe_write_csv(heat, file.path(paths$source, "Figure6_G0S2_context_effects.csv"))
  safe_write_csv(forest, file.path(paths$source, "Figure6_G0S2_blood_effects.csv"))
  right_column <- (p2 / p3) + patchwork::plot_layout(heights = c(0.80, 1.00))
  figure <- v34_panel_tag((p1 | right_column) + patchwork::plot_layout(widths = c(1.08, 0.92)))
  v34_save_plot(figure, "Figure6_peripheral_blood_validation", paths$figures, 170, 110)
}

v34_build_figure7 <- function(paths) {
  nodes <- data.frame(
    x = c(2, 2, 1, 3, 2), y = c(3.0, 2.0, 1.0, 1.0, 0.0),
    width = c(2.05, 2.05, 1.55, 1.55, 2.15), height = c(0.62, 0.68, 0.72, 0.72, 0.54),
    title = c("Shared tissue transcription", "Recurring biological themes", "Tissue-specific cellular context", "Limited systemic component", "Integrated interpretation"),
    detail = c("286 shared DEGs", "matrix, immune, stress and cell-cycle programs", "separate OA and OC atlases", "G0S2 replicated in both blood cohorts", "partial convergence with substantial context dependence"),
    fill = c("#ECE7F2", "#E8F2EC", "#DDEAF4", "#F8E1D8", "#F2F2F2"),
    stringsAsFactors = FALSE
  )
  arrows <- data.frame(
    x = c(2, 1.78, 2.22, 1, 3), y = c(2.66, 1.62, 1.62, 0.60, 0.60),
    xend = c(2, 1, 3, 1.72, 2.28), yend = c(2.36, 1.37, 1.37, 0.29, 0.29)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x, y, xend = xend, yend = yend), arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"), linewidth = 0.45, colour = "#4A4A4A") +
    ggplot2::geom_rect(data = nodes, ggplot2::aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill), colour = "#3C3C3C", linewidth = 0.4) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y + 0.09, label = title), size = 2.45, fontface = "bold") +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y = y - 0.14, label = v32_wrap(detail, 38L)), size = 1.9, colour = "#59636E", lineheight = 0.92) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(0.05, 3.95), ylim = c(-0.38, 3.45), clip = "off") +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 8, 6, 8))
  safe_write_csv(nodes, file.path(paths$source, "Figure7_integrated_model.csv"))
  v34_save_plot(p, "Figure7_integrated_interpretation", paths$figures, 170, 104)
}

v34_build_external_hallmark_supplement <- function(paths, external) {
  discovery <- utils::read.csv(file.path(paths$tables, "Table_S14_discovery_Hallmark_direction.csv"), check.names = FALSE)
  discovery <- discovery[discovery$both_significant %in% c(TRUE, "TRUE"), ]
  discovery_long <- rbind(
    data.frame(pathway_id = discovery$pathway_id, pathway = discovery$pathway, dataset_id = "OA discovery", disease = "OA", NES = discovery$OA_NES, FDR = discovery$OA_FDR),
    data.frame(pathway_id = discovery$pathway_id, pathway = discovery$pathway, dataset_id = "OC discovery", disease = "OC", NES = discovery$OC_NES, FDR = discovery$OC_FDR)
  )
  external_rows <- external$gsea[external$gsea$ID %in% discovery$pathway_id, , drop = FALSE]
  external_long <- data.frame(
    pathway_id = external_rows$ID,
    pathway = discovery$pathway[match(external_rows$ID, discovery$pathway_id)],
    dataset_id = external_rows$dataset_id,
    disease = external_rows$disease,
    NES = external_rows$NES,
    FDR = external_rows$p.adjust,
    stringsAsFactors = FALSE
  )
  table <- rbind(discovery_long, external_long)
  contexts <- c("OA discovery", "GSE117999", "GSE82107", "OC discovery", "GSE54388", "GSE12470")
  complete <- expand.grid(pathway_id = discovery$pathway_id, dataset_id = contexts, stringsAsFactors = FALSE)
  table <- merge(complete, table, by = c("pathway_id", "dataset_id"), all.x = TRUE)
  table$pathway <- discovery$pathway[match(table$pathway_id, discovery$pathway_id)]
  table$pathway_label <- factor(v32_wrap(table$pathway, 28L), levels = rev(v32_wrap(discovery$pathway, 28L)))
  table$dataset_id <- factor(table$dataset_id, levels = contexts)
  table$significance <- ifelse(is.finite(table$FDR) & table$FDR < 0.05, "FDR < 0.05", "Not FDR significant / unavailable")
  p <- ggplot2::ggplot(table, ggplot2::aes(dataset_id, pathway_label, fill = NES)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_point(data = table[table$significance != "FDR < 0.05", ], shape = 4, size = 1.2, stroke = 0.45, colour = "#555555") +
    ggplot2::geom_vline(xintercept = 3.5, linewidth = 0.55, colour = "#333333") +
    ggplot2::scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, limits = c(-3, 3), oob = scales::squish, na.value = "#EFEFEF") +
    ggplot2::labs(x = NULL, y = NULL, fill = "NES") +
    v34_theme(6.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), axis.text.y = ggplot2::element_text(size = 5.2), panel.grid = ggplot2::element_blank())
  safe_write_csv(table, file.path(paths$source, "SupplementaryFigure2_external_Hallmark_matrix.csv"))
  v34_save_plot(p, "SupplementaryFigure2_external_tissue_Hallmark", paths$figures, 170, 115)
}

v34_copy_supplementary_figures <- function(project_root, paths) {
  mapping <- list(
    SupplementaryFigure1_core_sensitivity = file.path(project_root, "results", "submission_v33", "figures", "SupplementaryFigure1_core_sensitivity"),
    SupplementaryFigure3_all_single_cell_UMAPs = file.path(project_root, "results", "submission_v33", "figures", "SupplementaryFigure3_all_single_cell_UMAPs"),
    SupplementaryFigure4_bulk_PCA_QC = file.path(project_root, "results", "submission_v33", "figures", "SupplementaryFigure4_bulk_PCA_QC"),
    SupplementaryFigure5_direction_aware_STRING_network = file.path(project_root, "results", "submission_v33", "figures", "SupplementaryFigure5_direction_aware_STRING_network")
  )
  for (target_stem in names(mapping)) {
    for (extension in c("pdf", "png")) {
      v34_copy(paste0(mapping[[target_stem]], ".", extension), file.path(paths$figures, paste0(target_stem, ".", extension)))
    }
  }
}

v34_write_documentation <- function(paths, external, candidates) {
  legends <- c(
    "## Figure 1. Study design and convergent evidence chain", "",
    "OA and OC tissue cohorts were analyzed separately. Shared tissue DEGs were characterized by functional enrichment, summarized with five representative genes, localized in separate single-cell atlases, and tested in independent blood cohorts. The blood layer was a validation filter rather than a biomarker-discovery study.", "",
    "## Figure 2. Shared tissue transcriptomic alterations in OA and OC", "",
    "**A-B,** Discovery-cohort volcano plots using FDR <0.05 and absolute log2 fold change >=1. **C,** Overlap of primary OA and OC DEGs. **D,** Direction-ordered heatmap of the 286 shared genes. Rows are ordered by the four OA/OC sign combinations; no cross-disease normalization is implied.", "",
    "## Figure 3. Shared genes converge on recurring biological themes", "",
    "**A,** Representative FDR-significant GO terms spanning cell-cycle, extracellular-matrix, immune, and stress categories. Point size denotes the number of shared genes. **B,** Hallmark pathways significant in both discovery cohorts. Blue circles show OA NES and orange triangles show OC NES; connector color indicates matching or opposite signs.", "",
    "## Figure 4. External tissue replication and representative-gene evidence", "",
    "**A,** Proportion of measured shared genes with the same direction as the corresponding OA or OC discovery contrast. Labels above bars report percentages; cohort-specific denominators are provided in Table S5. **B,** Discovery and external-cohort log2 fold changes for five representative genes. **C,** Transparent evidence matrix. The five genes summarize distinct evidence roles and do not constitute an optimized predictive signature.", "",
    "## Figure 5. Single-cell localization of representative genes", "",
    "**A-B,** Representative OA GSE255460 and OC GSE154600 UMAPs. Six exact source labels are highlighted in each atlas; remaining source labels are grey. **C,** Fraction of cells with detected G0S2, EFEMP1, AKAP12, SOX9, or DDIT3 expression within selected source labels. Fractions are interpreted within atlas and are not direct OA-OC quantitative comparisons. G0S2 was most frequently detected in ProC in the OA atlas and Myeloid.cell in the OC atlas.", "",
    "## Figure 6. Peripheral-blood validation identifies G0S2 as a limited systemic component", "",
    "**A,** Prespecified attrition from 286 shared tissue DEGs to one gene meeting the separate dual-cohort blood FDR rule. **B,** G0S2 log2 fold changes in OA tissue, OC tissue, OA PBMC, and OC blood-cell fraction. **C,** Within-cohort standardized blood effects shown as Hedges g with 95% confidence intervals. Effects were not pooled across diseases.", "",
    "## Figure 7. Integrated interpretation", "",
    "Shared tissue alterations converge on recurring biological themes but retain tissue, cohort, and cellular context. G0S2 was the sole dual-blood-FDR result. The integrated model supports partial molecular convergence with substantial context dependence, not a shared disease mechanism or circulating diagnostic signature.", "",
    "## Supplementary Figure 1. DEG-threshold and WGCNA sensitivity", "",
    "Shared-DEG retention across prespecified thresholds and disease-specific WGCNA soft-power sensitivity. Bootstrap and leave-one-out estimates are reported in Table S15.", "",
    "## Supplementary Figure 2. Hallmark states in discovery and external tissue cohorts", "",
    "NES values for the ten Hallmark pathways significant in both discovery cohorts. Crosses mark FDR >=0.05 or an unavailable enrichment estimate. The matrix shows strong OC replication and greater OA cohort dependence.", "",
    "## Supplementary Figure 3. Dataset-specific single-cell embeddings", "",
    "All five single-cell atlases with exact source labels. OA and OC datasets were not integrated into one latent space.", "",
    "## Supplementary Figure 4. Discovery-cohort bulk PCA and sample-correlation QC", "",
    "Separate unsupervised audits for OA and OC discovery cohorts. No outcome-informed sample exclusion was performed.", "",
    "## Supplementary Figure 5. Direction-aware STRING association landscape", "",
    "High-confidence physical STRING association graph and direction-aware topology. STRING was retained as auxiliary database context and did not select the representative genes."
  )
  writeLines(legends, file.path(paths$figures, "figure_legends.md"), useBytes = TRUE)

  registry <- data.frame(
    claim = c(
      "286 shared tissue DEGs", "Recurring matrix, immune, stress and cell-cycle themes",
      "External tissue direction agreement is disease and cohort dependent",
      "Five genes form an interpretable evidence summary", "G0S2 is the sole dual-blood-FDR result",
      "Representative signals occupy source-defined cell contexts", "No causal or diagnostic claim"
    ),
    evidence = c(
      "Figure 2; Table S2", "Figure 3; Tables S12-S14", "Figure 4; Tables S4-S6",
      "Figure 4; Table S7", "Figure 6; Tables S9-S11", "Figure 5; Tables S7-S8", "Figures 4-7; Discussion"
    ),
    boundary = c(
      "membership does not equal mechanism", "enrichment is not protein-level activity",
      "OA replication is weaker than OC replication", "not a signature or model",
      "one systemic-component candidate, not a biomarker", "no forced cell-type equivalence", "observational association only"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(registry, file.path(paths$root, "claim_evidence_registry_v34.csv"))

  scope <- data.frame(
    module = c("Tissue DEG", "Functional enrichment", "External tissue replication", "Representative genes", "Single-cell localization", "Peripheral blood", "WGCNA", "STRING", "Machine learning/ROC", "MR", "CellChat/NicheNet", "TF-miRNA", "DCA/nomogram"),
    manuscript_status = c("main", "main", "main", "main evidence summary", "main", "main validation layer", "supplementary", "supplementary", "excluded", "excluded", "excluded", "excluded", "excluded"),
    rationale = c("primary discovery", "biological interpretation", "cross-cohort reproducibility", "transparent illustration", "cellular source", "systemic-component validation", "within-cohort co-expression sensitivity", "database association context", "not a prediction study", "outside the focused transcriptomic question", "no validated communication axis", "prediction-only layer", "no clinical probability model"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(scope, file.path(paths$analysis, "V34_scope_decisions.csv"))

  terminology <- data.frame(
    concept = c("shared genes", "functional result", "five-gene set", "blood result", "single-cell labels", "overall inference"),
    canonical_term = c("shared tissue DEGs", "shared biological themes", "representative genes", "blood-replicated systemic component", "exact source-defined cell labels", "partial molecular convergence with context dependence"),
    prohibited_or_avoided = c("shared mechanism", "activated mechanism", "diagnostic signature", "blood biomarker", "forced CAF/tumor epithelial relabeling", "causal relationship"),
    stringsAsFactors = FALSE
  )
  safe_write_csv(terminology, file.path(paths$analysis, "V34_terminology_ledger.csv"))
}

run_reviewer_v34 <- function(project_root) {
  require_namespace("ggplot2", "V3.4 figures")
  require_namespace("patchwork", "V3.4 figure assembly")
  require_namespace("limma", "V3.4 external tissue models")
  require_namespace("clusterProfiler", "V3.4 external Hallmark analysis")
  require_namespace("svglite", "V3.4 vector export")
  require_namespace("ragg", "V3.4 raster preview export")
  require_namespace("ggrepel", "V3.4 UMAP labels")
  require_namespace("Matrix", "V3.4 sparse single-cell summaries")
  require_namespace("SummarizedExperiment", "V3.4 single-cell assay access")

  paths <- v34_output_paths(project_root)
  log_info("V3.4: converging the manuscript on tissue overlap, biological themes, cellular localization, external tissue replication, and conditional blood validation.")
  v34_prepare_tables(project_root, paths)
  external <- v34_external_validation(project_root, paths)
  genes <- c("G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3")
  single_cell <- v34_single_cell_candidates(project_root, paths, genes)
  candidates <- v34_candidate_evidence(paths, external, single_cell, genes)
  v34_build_figure1(paths)
  v34_build_figure2(project_root, paths)
  v34_build_figure3(project_root, paths)
  v34_build_figure4(paths, external, candidates)
  v34_build_figure5(project_root, paths, single_cell, genes)
  v34_build_figure6(project_root, paths)
  v34_build_figure7(paths)
  v34_build_external_hallmark_supplement(paths, external)
  v34_copy_supplementary_figures(project_root, paths)
  v34_write_documentation(paths, external, candidates)
  log_info("V3.4 completed: seven main figures, five supplementary figures, external tissue replication, five representative genes, and one dual-blood-FDR systemic-component candidate.")
  invisible(paths)
}
