#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript R/29_final_discovery_iteration.R <project_root>", call. = FALSE)
}

project_root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
figure_dir <- file.path(project_root, "results", "submission_v42", "figures")
source_dir <- file.path(figure_dir, "source_data")

required <- c(
  "Figure4_external_direction_summary.csv",
  "Figure4_illustrative_gene_effects.csv",
  "Figure4_candidate_evidence_matrix.csv"
)
missing <- required[!file.exists(file.path(source_dir, required))]
if (length(missing)) stop("Missing source data: ", paste(missing, collapse = ", "), call. = FALSE)

palette <- c(OA = "#79ADD2", OC = "#E8A179")

theme_submission <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#303030"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#303030"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "#303030"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(size = 8.2, face = "bold", colour = "#202020"),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

wrap_text <- function(x, width = 20L) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

save_plot_bundle <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(file.path(figure_dir, paste0(stem, ".svg")), width = width_in, height = height_in)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    file.path(figure_dir, paste0(stem, ".pdf")),
    width = width_in,
    height = height_in,
    family = "Arial"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    file.path(figure_dir, paste0(stem, ".tiff")),
    width = width_in,
    height = height_in,
    units = "in",
    res = dpi,
    compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(figure_dir, paste0(stem, ".png")),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300
  )
  print(plot)
  dev.off()
}

# Figure 4 contract:
# Core conclusion: external reproducibility is high in OC and more variable in OA,
# while the illustrative genes show distinct cohort-level effect patterns.
# Archetype: two-panel quantitative validation grid. Panel A is the hero summary;
# panel B provides gene-level resolution. The interpretive matrix is moved to S5.
summary <- read.csv(file.path(source_dir, required[[1L]]), check.names = FALSE)
summary$disease <- factor(summary$disease, levels = c("OA", "OC"))
summary$dataset_id <- factor(
  summary$dataset_id,
  levels = c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
)

p4a <- ggplot(summary, aes(dataset_id, 100 * direction_agreement, fill = disease)) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_col(width = 0.64, colour = "#333333", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.1f%%", 100 * direction_agreement)),
    vjust = -0.45,
    size = 2.45
  ) +
  scale_fill_manual(values = palette, guide = "none") +
  coord_cartesian(ylim = c(0, 102), clip = "off") +
  labs(
    x = NULL,
    y = "Shared genes with matching direction (%)",
    fill = NULL
  ) +
  theme_submission(7.0) +
  theme(
    axis.text.x = element_text(angle = 24, hjust = 1),
    legend.position = "none"
  )

effects <- read.csv(file.path(source_dir, required[[2L]]), check.names = FALSE)
contexts <- c("OA discovery", "GSE117999", "GSE82107", "OC discovery", "GSE54388", "GSE12470")
genes <- c("G0S2", "EFEMP1", "AKAP12", "SOX9", "DDIT3")
effects$context <- factor(effects$context, levels = contexts)
effects$gene <- factor(effects$gene, levels = rev(genes))

p4b <- ggplot(effects, aes(context, gene, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(is.finite(log2FC), sprintf("%.1f", ifelse(abs(log2FC) < 0.05, 0, log2FC)), "NA")),
    size = 2.15
  ) +
  geom_vline(xintercept = 3.5, linewidth = 0.55, colour = "#333333") +
  scale_fill_gradient2(
    low = "#2C7FB8",
    mid = "white",
    high = "#D95F0E",
    midpoint = 0,
    limits = c(-4.5, 4.5),
    oob = scales::squish,
    na.value = "#EEEEEE"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "log2 fold change"
  ) +
  theme_submission(7.0) +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  )

figure4 <- (p4a | p4b) +
  plot_layout(widths = c(0.9, 1.25)) +
  plot_annotation(tag_levels = "A")
save_plot_bundle(figure4, "Figure4_external_tissue_and_illustrative_genes", 170, 84)

# Supplementary Figure S5 contract:
# Core conclusion: the five genes were selected to illustrate different evidence
# roles, not to define a predictive signature or shared mechanism.
evidence <- read.csv(file.path(source_dir, required[[3L]]), check.names = FALSE)
evidence$gene <- factor(evidence$gene, levels = rev(genes))
evidence$evidence <- factor(
  evidence$evidence,
  levels = c("External direction", "Single-cell localization", "Dual-blood FDR", "Interpretive role")
)

s5 <- ggplot(evidence, aes(evidence, gene, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = wrap_text(label, 19L)), size = 2.35, lineheight = 0.92) +
  scale_fill_gradient(low = "#F2F2F2", high = "#5B8FA8", limits = c(0, 1), guide = "none") +
  scale_x_discrete(labels = c(
    "External direction" = "External\ndirection",
    "Single-cell localization" = "Single-cell\nlocalization",
    "Dual-blood FDR" = "Dual-blood\nFDR",
    "Interpretive role" = "Interpretive\nrole"
  )) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_submission(7.2) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(face = "italic")
  )
save_plot_bundle(s5, "SupplementaryFigure5_illustrative_gene_evidence", 170, 94)

# Figure 7 contract:
# Core conclusion: the shared molecular landscape is progressively resolved by
# direction, cohort, cell context and blood persistence without implying a driver.
nodes <- data.frame(
  x = rep(2, 6),
  y = c(5.5, 4.5, 3.5, 2.5, 1.5, 0.5),
  width = rep(2.6, 6),
  height = rep(0.64, 6),
  title = c(
    "Shared tissue landscape", "Directional stratification", "External reproducibility",
    "Source-defined cell localization", "Blood persistence", "G0S2 candidate"
  ),
  detail = c(
    "286 shared differentially expressed genes",
    "conserved and disease-context-dependent patterns",
    "high agreement in OC; greater variability in OA",
    "separate OA and OC atlases",
    "one gene met both blood-cohort FDR thresholds",
    "candidate systemic feature requiring validation"
  ),
  fill = c("#ECE7F2", "#E8F2EC", "#DDEAF4", "#E7F1EF", "#F8E1D8", "#FBE8D6"),
  stringsAsFactors = FALSE
)
arrows <- data.frame(
  x = rep(2, 5),
  y = c(5.16, 4.16, 3.16, 2.16, 1.16),
  xend = rep(2, 5),
  yend = c(4.84, 3.84, 2.84, 1.84, 0.84)
)

figure7 <- ggplot() +
  geom_segment(
    data = arrows,
    aes(x, y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
    linewidth = 0.45,
    colour = "#4A4A4A"
  ) +
  geom_rect(
    data = nodes,
    aes(xmin = x - width / 2, xmax = x + width / 2, ymin = y - height / 2, ymax = y + height / 2, fill = fill),
    colour = "#3C3C3C",
    linewidth = 0.4
  ) +
  geom_text(data = nodes, aes(x, y = y + 0.09, label = title), size = 2.6, fontface = "bold") +
  geom_text(
    data = nodes,
    aes(x, y = y - 0.14, label = wrap_text(detail, 56L)),
    size = 2.05,
    colour = "#59636E",
    lineheight = 0.92
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.20, 3.80), ylim = c(0.05, 5.95), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(6, 8, 6, 8))
save_plot_bundle(figure7, "Figure7_integrated_interpretation", 170, 104)

write.csv(nodes, file.path(source_dir, "Figure7_integrated_model.csv"), row.names = FALSE)
message("Final discovery-led figure iteration completed.")
