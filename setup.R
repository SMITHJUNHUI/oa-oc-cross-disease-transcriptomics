#!/usr/bin/env Rscript

options(repos = c(
  MRCIEU = "https://mrcieu.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

args <- commandArgs(trailingOnly = TRUE)
install_missing <- "--install" %in% args
snapshot <- "--snapshot" %in% args
restore <- "--restore" %in% args

cran_packages <- c(
  "data.table", "digest", "ggplot2", "glmnet", "jsonlite", "Matrix", "pheatmap",
  "pROC", "randomForest", "readxl", "renv", "scales", "survival", "survminer",
  "igraph", "remotes",
  "testthat", "timeROC", "uwot", "WGCNA", "yaml"
)

bioc_packages <- c(
  "AnnotationDbi", "Biobase", "BiocParallel", "BiocSingular", "bluster",
  "clusterProfiler", "edgeR", "enrichplot", "GEOquery", "hgu133plus2.db",
  "limma", "org.Hs.eg.db", "rhdf5", "S4Vectors", "scater", "scDblFinder",
  "scran", "scuttle", "SingleCellExperiment", "SingleR", "UCell",
  "SummarizedExperiment"
)

optional_packages <- c("MRPRESSO", "TwoSampleMR")
installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
missing_bioc <- setdiff(bioc_packages, installed)
missing_optional <- setdiff(optional_packages, installed)
missing_github <- setdiff("CellChat", installed)

cat("Missing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
cat("Missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
cat("Missing optional packages:", paste(missing_optional, collapse = ", "), "\n")
cat("Missing GitHub packages:", paste(missing_github, collapse = ", "), "\n")

if (restore) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }
  renv::restore(project = ".", prompt = FALSE)
  writeLines(
    paste("restored", as.character(Sys.time())),
    con = file.path("renv", ".restored"),
    useBytes = TRUE
  )
  cat("renv restore completed. The one-click runner will use it automatically.\n")
  quit(status = 0L)
}

if (!install_missing) {
  cat("\nCheck only. Re-run with --install to install missing packages.\n")
  quit(status = as.integer(length(c(missing_cran, missing_bioc, missing_github)) > 0L))
}

if (length(missing_cran) > 0L) {
  install.packages(missing_cran)
}

if (length(missing_bioc) > 0L) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (length(missing_optional) > 0L) {
  install.packages(
    missing_optional,
    repos = c(
      "https://mrcieu.r-universe.dev",
      "https://cloud.r-project.org"
    )
  )
}

if (length(missing_github) > 0L) {
  cellchat_vendor <- file.path(
    "environment", "vendor", "CellChat_2.2.0.9001.tar.gz"
  )
  if (file.exists(cellchat_vendor)) {
    install.packages(cellchat_vendor, repos = NULL, type = "source")
  } else {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("jinworks/CellChat", upgrade = "never")
  }
}

if (snapshot) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("renv installation failed; cannot write lockfile.", call. = FALSE)
  }
  renv::snapshot(prompt = FALSE)
}

cat("Dependency setup completed.\n")
