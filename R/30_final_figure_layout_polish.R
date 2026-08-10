#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/30_final_figure_layout_polish.R <project_root>", call. = FALSE)
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
})

figure_dir <- file.path(project_root, "results", "submission_v42", "figures")
source_dir <- file.path(figure_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_submission <- function(base_size = 6.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#303030"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "#404040"),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.6),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", colour = "#111111"),
      plot.tag.position = c(0.01, 0.99),
      panel.grid = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(5, 7, 5, 6)
    )
}

save_bundle <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(figure_dir, paste0(stem, ".svg"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  tiff_path <- file.path(figure_dir, paste0(stem, ".tiff"))

  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()

  svglite::svglite(svg_path, width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = dpi, background = "white")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in", res = dpi,
                 background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()

  invisible(c(pdf = pdf_path, svg = svg_path, png = png_path, tiff = tiff_path))
}

wrap_lines <- function(x, width) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

# Supplementary Figure S1: threshold sensitivity. The highlighted open point
# is the prespecified primary threshold, not a third statistical series.
s1_data <- read.csv(
  file.path(source_dir, "SupplementaryFigure1_DEG_threshold_sensitivity.csv"),
  check.names = FALSE
)
s1_data$fdr_label <- factor(s1_data$fdr_label, levels = c("FDR 0.01", "FDR 0.05"))
s1_primary <- s1_data[s1_data$is_primary %in% c(TRUE, "TRUE"), , drop = FALSE]

s1 <- ggplot(s1_data, aes(absolute_log2fc_threshold, shared_count, colour = fdr_label, group = fdr_label)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.35) +
  geom_point(
    data = s1_primary,
    shape = 21,
    size = 4.1,
    stroke = 0.9,
    colour = "#111111",
    fill = "white",
    inherit.aes = FALSE,
    aes(x = absolute_log2fc_threshold, y = shared_count)
  ) +
  geom_text(
    data = s1_primary,
    aes(x = absolute_log2fc_threshold, y = shared_count, label = paste0("Primary: ", shared_count)),
    colour = "#222222",
    size = 2.15,
    hjust = -0.28,
    vjust = -0.65,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = c("FDR 0.01" = "#2C7FB8", "FDR 0.05" = "#D95F02"), name = NULL) +
  scale_x_continuous(breaks = c(0.5, 1.0, 1.5), expand = expansion(mult = c(0.04, 0.14))) +
  scale_y_continuous(breaks = c(300, 600, 900), expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = "Absolute log2-fold-change threshold", y = "Shared genes") +
  theme_submission(7.0) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.key.width = grid::unit(7.5, "mm")
  )
save_bundle(s1, "SupplementaryFigure1_DEG_threshold_sensitivity", 170, 102)

# Supplementary Figure S2: external Hallmark matrix. Tiles retain the available
# NES estimate; crosses flag estimates that are not FDR-significant or unavailable.
s2_data <- read.csv(
  file.path(source_dir, "SupplementaryFigure2_external_Hallmark_matrix.csv"),
  check.names = FALSE
)
s2_contexts <- c("OA discovery", "GSE117999", "GSE82107", "OC discovery", "GSE54388", "GSE12470")
s2_pathways <- c(
  "E2f Targets", "G2m Checkpoint", "Epithelial Mesenchymal\nTransition",
  "Coagulation", "Kras Signaling Up", "Mtorc1 Signaling", "Apoptosis",
  "Complement", "Glycolysis", "Mitotic Spindle"
)
s2_data$dataset_id <- factor(s2_data$dataset_id, levels = s2_contexts)
s2_data$pathway_label <- gsub("Epithelial Mesenchymal Transition", "Epithelial Mesenchymal\nTransition", s2_data$pathway_label, fixed = TRUE)
s2_data$pathway_label <- factor(s2_data$pathway_label, levels = rev(s2_pathways))
s2_data$flag <- is.na(s2_data$FDR) | s2_data$FDR >= 0.05

s2 <- ggplot(s2_data, aes(dataset_id, pathway_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(data = s2_data[s2_data$flag, , drop = FALSE], label = "×", colour = "#555555", size = 2.65) +
  geom_vline(xintercept = 3.5, colour = "#303030", linewidth = 0.65) +
  scale_fill_gradient2(
    low = "#2C7FB8", mid = "white", high = "#ED7D31", midpoint = 0,
    limits = c(-3, 3), oob = scales::squish, na.value = "#EEEEEE", name = "NES"
  ) +
  labs(x = NULL, y = NULL) +
  theme_submission(6.7) +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 5.8),
    axis.text.y = element_text(size = 5.5),
    legend.position = "right",
    plot.margin = margin(6, 8, 5, 7)
  )
save_bundle(s2, "SupplementaryFigure2_external_tissue_Hallmark", 170, 115)

# Supplementary Figure S3: five data sets remain in separate analytical spaces.
# Panel letters, rather than embedded titles, map each data set in the legend.
s3_data <- read.csv(
  file.path(project_root, "results", "submission_v32", "figures", "source_data", "SupplementaryFigure3_UMAP_subsamples.csv"),
  check.names = FALSE
)
s3_ids <- c("GSE104782", "GSE154600", "GSE169454", "GSE180661", "GSE255460")
s3_labels <- sort(unique(s3_data$label))
s3_colours <- setNames(scales::hue_pal(l = 60, c = 90)(length(s3_labels)), s3_labels)
s3_manual_labels <- rbind(
  data.frame(
    dataset_id = "GSE154600",
    label = c("B.cell", "Endothelial.cell", "Fibroblast", "Myeloid.cell", "Ovarian.cancer.cell", "Plasma.cell", "T.cell"),
    label_x = c(-9.5, 7.0, 11.5, -9.5, -3.5, -2.5, -3.5),
    label_y = c(-4.0, 8.3, -3.0, 10.5, 0.0, 8.4, -10.0),
    hjust = c(0, 0, 1, 0, 0.5, 1, 0.5),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset_id = "GSE180661",
    label = c("B.cell", "Dendritic.cell", "Endothelial.cell", "Fibroblast", "Mast.cell", "Myeloid.cell", "Ovarian.cancer.cell", "Plasma.cell", "T.cell"),
    label_x = c(1.0, -10.5, 1.5, 12.0, -11.5, -4.5, 8.5, 1.5, -10.0),
    label_y = c(16.0, -4.0, -3.8, -12.0, 2.0, -13.0, 8.5, 1.3, 10.5),
    hjust = c(0.5, 1, 0, 1, 0, 0.5, 1, 0.5, 0),
    stringsAsFactors = FALSE
  )
)

build_umap_panel <- function(dataset_id, tag) {
  dat <- s3_data[s3_data$dataset_id == dataset_id, , drop = FALSE]
  centres <- aggregate(cbind(UMAP1, UMAP2) ~ label, data = dat, FUN = median)
  p <- ggplot(dat, aes(UMAP1, UMAP2, colour = label)) +
    geom_point(size = 0.10, alpha = 0.70, stroke = 0)
  manual <- s3_manual_labels[s3_manual_labels$dataset_id == dataset_id, , drop = FALSE]
  if (nrow(manual)) {
    manual <- merge(centres, manual, by = "label", all.x = TRUE, sort = FALSE)
    p <- p +
      geom_segment(
        data = manual,
        aes(x = UMAP1, y = UMAP2, xend = label_x, yend = label_y),
        inherit.aes = FALSE,
        colour = "#777777",
        linewidth = 0.18
      ) +
      geom_text(
        data = manual,
        aes(x = label_x, y = label_y, label = label, hjust = hjust),
        inherit.aes = FALSE,
        size = 1.80,
        colour = "#333333"
      )
  } else {
    p <- p + geom_text_repel(
        data = centres,
        aes(label = label),
        seed = 20260810,
        size = 1.80,
        colour = "#333333",
        box.padding = 0.42,
        point.padding = 0.14,
        min.segment.length = 0,
        max.overlaps = Inf,
        force = 3,
        force_pull = 0.08,
        max.time = 5,
        max.iter = 100000,
        segment.size = 0.18,
        segment.color = "#777777",
        show.legend = FALSE
      )
  }
  p +
    scale_colour_manual(values = s3_colours, guide = "none") +
    coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = "UMAP 1", y = "UMAP 2", tag = tag) +
    theme_submission(6.2) +
    theme(
      axis.title = element_text(size = 5.8),
      axis.text = element_text(size = 5.0),
      plot.margin = margin(4, 5, 4, 5)
    )
}

s3_plots <- Map(build_umap_panel, s3_ids, LETTERS[1:5])
s3_top <- cowplot::plot_grid(plotlist = s3_plots[1:3], ncol = 3, align = "hv", axis = "tblr")
s3_bottom <- cowplot::plot_grid(
  NULL, s3_plots[[4]], s3_plots[[5]], NULL,
  ncol = 4,
  rel_widths = c(0.5, 1, 1, 0.5),
  align = "hv",
  axis = "tblr"
)
s3 <- cowplot::plot_grid(s3_top, s3_bottom, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
save_bundle(s3, "SupplementaryFigure3_all_single_cell_UMAPs", 170, 113)

# Supplementary Figure S4: unsupervised discovery-cohort QC. Labels are limited
# to the three lowest-correlation samples per cohort and repelled from plot edges.
s4_pca <- read.csv(
  file.path(project_root, "results", "submission_v33", "supplementary_tables", "Table_S31a_bulk_unsupervised_PCA.csv"),
  check.names = FALSE
)
s4_cor <- read.csv(
  file.path(project_root, "results", "submission_v33", "supplementary_tables", "Table_S31b_bulk_sample_correlation_QC.csv"),
  check.names = FALSE
)

build_pca_panel <- function(dataset, colour, tag) {
  dat <- s4_pca[s4_pca$dataset == dataset, , drop = FALSE]
  disease <- dat[dat$group == "Disease", , drop = FALSE]
  normal <- dat[dat$group != "Disease", , drop = FALSE]
  pc1 <- unique(dat$PC1_variance)[1]
  pc2 <- unique(dat$PC2_variance)[1]
  ggplot() +
    stat_ellipse(data = disease, aes(PC1, PC2), colour = colour, linetype = "dashed", linewidth = 0.55, level = 0.95) +
    stat_ellipse(data = normal, aes(PC1, PC2), colour = "#7A7A7A", linetype = "dashed", linewidth = 0.55, level = 0.95) +
    geom_point(data = disease, aes(PC1, PC2), shape = 21, fill = colour, colour = "#3B3B3B", size = 2.2, stroke = 0.45) +
    geom_point(data = normal, aes(PC1, PC2), shape = 24, fill = "white", colour = "#666666", size = 2.25, stroke = 0.55) +
    labs(x = sprintf("PC1 (%.1f%%)", pc1), y = sprintf("PC2 (%.1f%%)", pc2), tag = tag) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    theme_submission(6.6) +
    theme(plot.margin = margin(7, 8, 6, 8))
}

build_correlation_panel <- function(dataset, colour, tag) {
  dat <- s4_cor[s4_cor$dataset == dataset, , drop = FALSE]
  dat$sample_index <- seq_len(nrow(dat))
  disease <- dat[dat$group == "Disease", , drop = FALSE]
  normal <- dat[dat$group != "Disease", , drop = FALSE]
  low <- dat[dat$lowest_three %in% c(TRUE, "TRUE"), , drop = FALSE]
  low <- low[order(low$median_sample_correlation), , drop = FALSE]
  x_span <- diff(range(dat$sample_index))
  labels_on_left <- stats::median(low$sample_index) > mean(range(dat$sample_index))
  low$label_x <- if (labels_on_left) {
    low$sample_index - max(4, 0.10 * x_span)
  } else {
    low$sample_index + max(2, 0.06 * x_span)
  }
  low$label_y <- seq(
    max(0, min(low$median_sample_correlation) - 0.025),
    max(low$median_sample_correlation) + 0.075,
    length.out = nrow(low)
  )
  ggplot() +
    geom_point(data = disease, aes(sample_index, median_sample_correlation), shape = 21, fill = colour,
               colour = "#3B3B3B", size = 2.05, stroke = 0.42) +
    geom_point(data = normal, aes(sample_index, median_sample_correlation), shape = 24, fill = "white",
               colour = "#666666", size = 2.15, stroke = 0.52) +
    geom_segment(
      data = low,
      aes(x = sample_index, y = median_sample_correlation, xend = label_x, yend = label_y),
      linewidth = 0.20,
      colour = "#777777"
    ) +
    geom_text(
      data = low,
      aes(label_x, label_y, label = sample),
      hjust = if (labels_on_left) 1 else 0,
      size = 1.78,
      colour = "#454545",
      lineheight = 0.95
    ) +
    labs(x = "Sample index", y = "Median Pearson correlation", tag = tag) +
    scale_x_continuous(expand = expansion(mult = c(0.13, 0.16))) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.13))) +
    coord_cartesian(clip = "off") +
    theme_submission(6.6) +
    theme(plot.margin = margin(7, 11, 8, 10))
}

s4 <- (
  build_pca_panel("OA discovery", "#2C7FB8", "A") |
    build_pca_panel("OC discovery", "#D95F02", "B")
) / (
  build_correlation_panel("OA discovery", "#2C7FB8", "C") |
    build_correlation_panel("OC discovery", "#D95F02", "D")
) + plot_layout(heights = c(0.95, 1.05))
save_bundle(s4, "SupplementaryFigure4_bulk_PCA_QC", 170, 138)

# Supplementary Figure S5: transparent cross-layer roles for the illustrative
# genes. Column headings are axis labels, not subplot titles.
s5_data <- read.csv(file.path(source_dir, "Figure4_candidate_evidence_matrix.csv"), check.names = FALSE)
s5_genes <- c("G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3")
s5_data$gene <- factor(s5_data$gene, levels = rev(s5_genes))
s5_data$evidence <- factor(
  s5_data$evidence,
  levels = c("External direction", "Single-cell localization", "Dual-blood FDR", "Interpretive role")
)
s5 <- ggplot(s5_data, aes(evidence, gene, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = wrap_lines(label, 18L)), size = 2.35, lineheight = 0.92) +
  scale_fill_gradient(low = "#F2F2F2", high = "#5B8FA8", limits = c(0, 1), guide = "none") +
  scale_x_discrete(labels = c(
    "External direction" = "External\ndirection",
    "Single-cell localization" = "Single-cell\nlocalization",
    "Dual-blood FDR" = "Dual-blood\nFDR",
    "Interpretive role" = "Interpretive\nrole"
  )) +
  labs(x = NULL, y = NULL) +
  theme_submission(7.0) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6.2),
    axis.text.y = element_text(face = "italic", size = 6.4),
    axis.ticks = element_blank(),
    plot.margin = margin(6, 7, 7, 7)
  )
save_bundle(s5, "SupplementaryFigure5_illustrative_gene_evidence", 170, 94)

# Assemble the five R-generated figures into an A4, one-figure-per-page PDF.
# The page heading and caption remain outside the figure panels.
supplement_pdf <- file.path(project_root, "results", "submission_v42", "Additional_file_2_supplementary_figures.pdf")
page_specs <- list(
  list(
    heading = "Figure S1. Differential-expression threshold sensitivity",
    plot = s1,
    caption = paste0(
      "Shared-gene retention across six prespecified combinations of FDR and absolute log2-fold-change thresholds. ",
      "Blue and orange lines denote FDR thresholds of 0.01 and 0.05, respectively. The open point denotes the primary ",
      "threshold (FDR < 0.05 and absolute log2 fold change >= 1). Complete membership is provided in Additional file 1: Tables S3a and S3b."
    ),
    plot_y = 0.56, plot_h = 0.68
  ),
  list(
    heading = "Figure S2. Hallmark states in discovery and external tissue cohorts",
    plot = s2,
    caption = paste0(
      "Normalized enrichment scores for the ten Hallmark pathways significant in both discovery cohorts. Columns are ",
      "ordered as OA discovery, GSE117999, GSE82107, OC discovery, GSE54388 and GSE12470. Orange and blue denote positive ",
      "and negative enrichment, respectively; crosses mark FDR >= 0.05 or an unavailable estimate."
    ),
    plot_y = 0.57, plot_h = 0.70
  ),
  list(
    heading = "Figure S3. Dataset-specific single-cell embeddings",
    plot = s3,
    caption = paste0(
      "Dataset-specific single-cell embeddings with exact source labels: (A) GSE104782, (B) GSE154600, (C) GSE169454, ",
      "(D) GSE180661 and (E) GSE255460. Large atlases were deterministically subsampled for display. OA and OC data sets ",
      "were retained as separate analytical spaces; distances and cluster positions are therefore interpreted only within each panel."
    ),
    plot_y = 0.61, plot_h = 0.52
  ),
  list(
    heading = "Figure S4. Discovery-cohort PCA and sample-correlation quality control",
    plot = s4,
    caption = paste0(
      "Unsupervised quality-control displays for the OA and OC discovery cohorts. (A) OA and (B) OC PCA using the top ",
      "1,000 variable genes. Filled circles denote disease samples, open triangles denote reference samples and dashed ",
      "ellipses summarize group dispersion. (C) OA and (D) OC median sample correlations; only the three lowest samples ",
      "per cohort are labelled. These audits were not used for outcome-informed sample removal or batch correction."
    ),
    plot_y = 0.57, plot_h = 0.71
  ),
  list(
    heading = "Figure S5. Cross-layer evidence summary for the illustrative genes",
    plot = s5,
    caption = paste0(
      "Descriptive summary of external-direction agreement, source-defined single-cell localization, dual-blood FDR status ",
      "and interpretive role for G0S2, EFEMP1, AKAP12, SOX9 and DDIT3. These genes illustrate distinct evidence patterns ",
      "and do not constitute an optimized predictive signature."
    ),
    plot_y = 0.56, plot_h = 0.66
  )
)

grDevices::cairo_pdf(supplement_pdf, width = 210 / 25.4, height = 297 / 25.4, family = "Arial", onefile = TRUE)
for (spec in page_specs) {
  grid::grid.newpage()
  grid::grid.text(
    spec$heading,
    x = grid::unit(15, "mm"),
    y = grid::unit(282, "mm"),
    just = c("left", "top"),
    gp = grid::gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "#2B5C90")
  )
  print(
    spec$plot,
    newpage = FALSE,
    vp = grid::viewport(x = 0.5, y = spec$plot_y, width = 0.88, height = spec$plot_h)
  )
  caption <- wrap_lines(spec$caption, 118L)
  grid::grid.text(
    caption,
    x = grid::unit(15, "mm"),
    y = grid::unit(30, "mm"),
    just = c("left", "bottom"),
    gp = grid::gpar(fontfamily = "Arial", fontsize = 7.5, col = "#111111", lineheight = 1.18)
  )
}
grDevices::dev.off()

message("Final supplementary figures and the A4 supplementary PDF were rebuilt in R.")
