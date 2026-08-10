read_geo_platform_annotation <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, warn = FALSE)
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  if (is.na(begin) || is.na(end) || end <= begin + 1L) {
    stop("Invalid GEO annotation table: ", path, call. = FALSE)
  }
  table_lines <- lines[(begin + 1L):(end - 1L)]
  utils::read.delim(
    textConnection(table_lines),
    sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    stringsAsFactors = FALSE, fill = TRUE
  )
}

collapse_geo_expression_to_gene <- function(expression, annotation) {
  annotation <- annotation[, c("ID", "Gene symbol"), drop = FALSE]
  names(annotation) <- c("probe", "gene")
  annotation$probe <- as.character(annotation$probe)
  annotation$gene <- trimws(as.character(annotation$gene))
  annotation <- annotation[
    nzchar(annotation$gene) &
      !grepl("///|//|;|,", annotation$gene) &
      grepl("^[A-Za-z0-9._-]+$", annotation$gene) &
      annotation$probe %in% rownames(expression),
    , drop = FALSE
  ]
  annotation <- annotation[!duplicated(annotation$probe), , drop = FALSE]
  probe_iqr <- apply(expression[annotation$probe, , drop = FALSE], 1L, stats::IQR, na.rm = TRUE)
  annotation$probe_iqr <- probe_iqr[annotation$probe]
  annotation <- annotation[order(annotation$gene, -annotation$probe_iqr, annotation$probe), , drop = FALSE]
  annotation <- annotation[!duplicated(annotation$gene), , drop = FALSE]
  gene_expression <- expression[annotation$probe, , drop = FALSE]
  rownames(gene_expression) <- annotation$gene
  list(expression = gene_expression, mapping = annotation)
}

fit_blood_limma <- function(expression, group, disease_label, control_label) {
  group <- factor(group, levels = c(control_label, disease_label))
  if (anyNA(group)) stop("Unrecognized blood cohort group.", call. = FALSE)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- c("Control", "Disease")
  fit <- limma::lmFit(expression, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(Disease - Control, levels = design))
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  result <- limma::topTable(fit, number = Inf, sort.by = "none")
  result$gene <- rownames(result)
  result <- result[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

  disease <- group == disease_label
  control <- group == control_label
  n_disease <- sum(disease)
  n_control <- sum(control)
  disease_mean <- rowMeans(expression[, disease, drop = FALSE], na.rm = TRUE)
  control_mean <- rowMeans(expression[, control, drop = FALSE], na.rm = TRUE)
  disease_sd <- apply(expression[, disease, drop = FALSE], 1L, stats::sd, na.rm = TRUE)
  control_sd <- apply(expression[, control, drop = FALSE], 1L, stats::sd, na.rm = TRUE)
  pooled_sd <- sqrt(((n_disease - 1) * disease_sd^2 + (n_control - 1) * control_sd^2) /
    (n_disease + n_control - 2))
  cohen_d <- (disease_mean - control_mean) / pooled_sd
  correction <- 1 - 3 / (4 * (n_disease + n_control) - 9)
  hedges_g <- correction * cohen_d
  variance_g <- (n_disease + n_control) / (n_disease * n_control) +
    hedges_g^2 / (2 * (n_disease + n_control - 2))
  effects <- data.frame(
    gene = rownames(expression),
    n_disease = n_disease,
    n_control = n_control,
    disease_mean = disease_mean,
    control_mean = control_mean,
    hedges_g = hedges_g,
    hedges_g_se = sqrt(variance_g),
    hedges_g_lower = hedges_g - 1.96 * sqrt(variance_g),
    hedges_g_upper = hedges_g + 1.96 * sqrt(variance_g),
    stringsAsFactors = FALSE
  )
  merge(result, effects, by = "gene", all.x = TRUE, sort = FALSE)
}

run_blood_systemic_v33 <- function(project_root) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery is required.", call. = FALSE)
  if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required.", call. = FALSE)
  cache <- file.path(project_root, "results", "cache", "blood_v33")
  output <- ensure_dir(file.path(project_root, "results", "blood_v33_internal"))
  required <- c(
    "GSE48556_series_matrix.txt.gz", "GSE31682_series_matrix.txt.gz",
    "GPL6947.annot.gz", "GPL2986.annot.gz"
  )
  missing <- required[!file.exists(file.path(cache, required))]
  if (length(missing)) stop("Missing blood cache files: ", paste(missing, collapse = ", "), call. = FALSE)

  oa_eset <- GEOquery::getGEO(
    filename = file.path(cache, "GSE48556_series_matrix.txt.gz"),
    getGPL = FALSE
  )
  oc_eset <- GEOquery::getGEO(
    filename = file.path(cache, "GSE31682_series_matrix.txt.gz"),
    getGPL = FALSE
  )
  oa_group <- as.character(Biobase::pData(oa_eset)[["disease state:ch1"]])
  oc_group <- as.character(Biobase::pData(oc_eset)[["disease:ch1"]])
  if (!identical(as.integer(table(oa_group)[c("case", "control")]), c(106L, 33L))) {
    stop("GSE48556 group audit failed.", call. = FALSE)
  }
  if (!identical(as.integer(table(oc_group)[c("epithelial ovarian cancer", "none")]), c(48L, 20L))) {
    stop("GSE31682 group audit failed.", call. = FALSE)
  }

  oa_annotation <- read_geo_platform_annotation(file.path(cache, "GPL6947.annot.gz"))
  oc_annotation <- read_geo_platform_annotation(file.path(cache, "GPL2986.annot.gz"))
  oa_gene <- collapse_geo_expression_to_gene(Biobase::exprs(oa_eset), oa_annotation)
  oc_gene <- collapse_geo_expression_to_gene(Biobase::exprs(oc_eset), oc_annotation)
  oa_de <- fit_blood_limma(oa_gene$expression, oa_group, "case", "control")
  oc_de <- fit_blood_limma(oc_gene$expression, oc_group, "epithelial ovarian cancer", "none")
  names(oa_de)[names(oa_de) != "gene"] <- paste0("OA_blood_", names(oa_de)[names(oa_de) != "gene"])
  names(oc_de)[names(oc_de) != "gene"] <- paste0("OC_blood_", names(oc_de)[names(oc_de) != "gene"])

  tissue <- utils::read.csv(file.path(
    project_root, "results", "submission_v32", "supplementary_tables",
    "Table_S2_shared_differentially_expressed_genes.csv"
  ))
  tissue$tissue_concordant <- tissue$logFC_OA * tissue$logFC_OC > 0
  screen <- merge(tissue, oa_de, by = "gene", all.x = TRUE, sort = FALSE)
  screen <- merge(screen, oc_de, by = "gene", all.x = TRUE, sort = FALSE)
  screen$measured_both_blood <- !is.na(screen$OA_blood_logFC) & !is.na(screen$OC_blood_logFC)
  screen$tissue_direction <- sign(screen$logFC_OA)
  screen$all_four_same_direction <- screen$tissue_concordant & screen$measured_both_blood &
    sign(screen$OA_blood_logFC) == screen$tissue_direction &
    sign(screen$OC_blood_logFC) == screen$tissue_direction
  screen$both_blood_nominal <- screen$all_four_same_direction &
    screen$OA_blood_P.Value < 0.05 & screen$OC_blood_P.Value < 0.05
  screen$both_blood_fdr <- screen$all_four_same_direction &
    screen$OA_blood_adj.P.Val < 0.05 & screen$OC_blood_adj.P.Val < 0.05
  screen <- screen[order(!screen$both_blood_fdr, !screen$both_blood_nominal,
    pmax(screen$OA_blood_adj.P.Val, screen$OC_blood_adj.P.Val, na.rm = TRUE)), ]

  positive <- screen[screen$both_blood_fdr %in% TRUE, , drop = FALSE]
  nominal <- screen[screen$both_blood_nominal %in% TRUE, , drop = FALSE]
  counts <- data.frame(
    stage = c(
      "Shared tissue DEGs", "Tissue-concordant DEGs", "Measured in both blood cohorts",
      "All four effects same direction", "Nominal P<0.05 in both blood cohorts",
      "FDR<0.05 in both blood cohorts"
    ),
    genes = c(
      nrow(screen), sum(screen$tissue_concordant),
      sum(screen$tissue_concordant & screen$measured_both_blood),
      sum(screen$all_four_same_direction, na.rm = TRUE),
      nrow(nominal), nrow(positive)
    ),
    stringsAsFactors = FALSE
  )
  dataset_audit <- data.frame(
    dataset_id = c("GSE48556", "GSE31682"),
    disease = c("OA", "OC"),
    biospecimen = c("PBMC", "blood cell fraction"),
    platform = c("GPL6947 Illumina HumanHT-12 V3.0", "GPL2986 ABI Human Genome Survey Microarray V2"),
    disease_n = c(106L, 48L), control_n = c(33L, 20L),
    mapped_genes = c(nrow(oa_gene$expression), nrow(oc_gene$expression)),
    official_url = c(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48556",
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31682"
    ),
    stringsAsFactors = FALSE
  )
  decision <- data.frame(
    rule = "tissue concordant; measured in both blood cohorts; all four effects same sign; FDR<0.05 in each blood cohort",
    positive_gene_count = nrow(positive),
    nominal_gene_count_internal_only = nrow(nominal),
    manuscript_status = if (nrow(positive) > 0L) "ELIGIBLE" else "OMIT_ENTIRE_BLOOD_MODULE",
    rationale = if (nrow(positive) > 0L) {
      "At least one gene met the prespecified four-layer direction and dual-FDR rule."
    } else {
      "No gene met the prespecified dual-cohort FDR rule; blood is not claimed as a systemic component."
    },
    stringsAsFactors = FALSE
  )

  safe_write_csv(dataset_audit, file.path(output, "blood_dataset_audit.csv"))
  safe_write_csv(oa_gene$mapping, file.path(output, "GSE48556_probe_mapping_selected.csv"))
  safe_write_csv(oc_gene$mapping, file.path(output, "GSE31682_probe_mapping_selected.csv"))
  safe_write_csv(oa_de, file.path(output, "GSE48556_blood_DE_all.csv"))
  safe_write_csv(oc_de, file.path(output, "GSE31682_blood_DE_all.csv"))
  safe_write_csv(screen, file.path(output, "tissue_blood_four_layer_screen_full_internal.csv"))
  safe_write_csv(nominal, file.path(output, "nominal_four_layer_genes_internal_only.csv"))
  safe_write_csv(positive, file.path(output, "FDR_supported_systemic_component_genes.csv"))
  safe_write_csv(counts, file.path(output, "blood_screen_attrition.csv"))
  safe_write_csv(decision, file.path(output, "blood_module_decision.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "sessionInfo.txt"))
  message("Blood module decision: ", decision$manuscript_status, "; FDR-supported genes = ", nrow(positive))
  invisible(list(decision = decision, positive = positive, counts = counts))
}
