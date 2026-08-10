testthat::test_that("single-cell thresholds are batch-specific and non-destructive", {
  metrics <- data.frame(
    cell_id = paste0("C", seq_len(12)),
    batch = rep(c("A", "B"), each = 6),
    nCount = c(100, 105, 110, 115, 120, 1000, 500, 510, 520, 530, 540, 550),
    nFeature = c(20, 21, 22, 23, 24, 100, 50, 51, 52, 53, 54, 55),
    percent_mt = c(1, 1, 2, 2, 3, 40, 2, 2, 3, 3, 4, 4)
  )
  thresholds <- derive_sc_qc_thresholds(metrics)
  audited <- apply_sc_qc_thresholds(metrics, thresholds)

  testthat::expect_equal(nrow(thresholds), 2L)
  testthat::expect_equal(nrow(audited), nrow(metrics))
  testthat::expect_identical(audited$cell_id, metrics$cell_id)
  testthat::expect_true(all(c(
    "passes_core_QC", "passes_QC", "failure_reason"
  ) %in% names(audited)))
  testthat::expect_true(any(!audited$passes_QC))
})

testthat::test_that("missing mitochondrial metric is explicit and does not invent values", {
  metrics <- data.frame(
    batch = rep("A", 20),
    nCount = seq(100, 290, by = 10),
    nFeature = seq(50, 145, by = 5),
    percent_mt = NA_real_
  )
  thresholds <- derive_sc_qc_thresholds(metrics)
  audited <- apply_sc_qc_thresholds(metrics, thresholds)

  testthat::expect_true(is.na(thresholds$percent_mt_upper))
  testthat::expect_true(all(is.na(audited$percent_mt)))
  testthat::expect_false(any(grepl(
    "mitochondrial",
    audited$failure_reason,
    fixed = TRUE
  )))
})

testthat::test_that("enabled single-cell config registers all five adapters", {
  project_root <- getOption("ocoa.test_project_root")
  config <- read_project_config(project_root, "config/config.example.yml")
  adapters <- vapply(
    config$single_cell$datasets,
    `[[`,
    character(1),
    "adapter"
  )

  testthat::expect_true(config$modules$single_cell)
  testthat::expect_equal(length(adapters), 5L)
  testthat::expect_equal(anyDuplicated(adapters), 0L)
  testthat::expect_setequal(
    adapters,
    c(
      "gse255460_wide_counts",
      "gse104782_umi_counts",
      "gse169454_tenx_raw_tar",
      "gse180661_h5_csr",
      "tenx_tar_gse154600"
    )
  )
  testthat::expect_identical(config$single_cell$scope, "qc_gate_only")
  testthat::expect_true(config$single_cell$allow_partial)
})

testthat::test_that("GSE104782 cell metadata are parsed without inventing batches", {
  parsed <- .sc_parse_gse104782_cells(c("OA1_S0.1", "OA10_S4.32"))

  testthat::expect_identical(parsed$donor, c("OA1", "OA10"))
  testthat::expect_identical(parsed$state_code, c("S0", "S4"))
  testthat::expect_identical(parsed$within_group_index, c(1L, 32L))
  testthat::expect_error(
    .sc_parse_gse104782_cells("not_a_valid_cell"),
    "do not all match"
  )
})

testthat::test_that("GSE169454 manifest requires seven paired 10x libraries", {
  samples <- c(
    "GSM1_normal1", "GSM2_normal2", "GSM3_normal3",
    "GSM4_oa1", "GSM5_oa2", "GSM6_oa3", "GSM7_oa4"
  )
  members <- unlist(lapply(samples, function(sample_id) {
    c(
      paste0(sample_id, "_filtered_barcodes.tsv.gz"),
      paste0(sample_id, "_filtered_features.tsv.gz"),
      paste0(sample_id, "_filtered_matrix.mtx.gz"),
      paste0(sample_id, "_raw_barcodes.tsv.gz"),
      paste0(sample_id, "_raw_features.tsv.gz"),
      paste0(sample_id, "_raw_matrix.mtx.gz")
    )
  }))
  manifest <- .sc_gse169454_member_manifest(members)

  testthat::expect_equal(nrow(manifest), 42L)
  testthat::expect_equal(length(unique(manifest$sample_id)), 7L)
  testthat::expect_error(
    .sc_gse169454_member_manifest(members[-1L]),
    "not exactly seven"
  )
})

testthat::test_that("CSR H5 row extraction preserves requested cell order", {
  testthat::skip_if_not_installed("rhdf5")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  rhdf5::h5createFile(path)
  rhdf5::h5createGroup(path, "X")
  rhdf5::h5write(c(1, 2, 3, 4, 5, 6), path, "X/data")
  rhdf5::h5write(as.integer(c(0, 2, 1, 3, 0, 3)), path, "X/indices")
  indptr <- c(0, 2, 2, 4, 6)
  rhdf5::h5write(indptr, path, "X/indptr")

  observed <- .sc_read_h5_csr_rows(
    path,
    h5_rows = c(4L, 1L, 3L),
    indptr = indptr,
    n_features = 4L
  )
  expected <- matrix(
    c(
      5, 0, 0, 6,
      1, 0, 2, 0,
      0, 3, 0, 4
    ),
    nrow = 4L,
    ncol = 3L
  )
  testthat::expect_equal(as.matrix(observed), expected)
})

testthat::test_that("partitioned disk CSR restores GSE255460 counts exactly", {
  testthat::skip_if_not_installed("jsonlite")
  bundle_dir <- tempfile("gse255460-csr-")
  dir.create(bundle_dir, recursive = TRUE)
  on.exit(unlink(bundle_dir, recursive = TRUE, force = TRUE), add = TRUE)
  partition_dir <- file.path(bundle_dir, "P1")
  dir.create(partition_dir)
  writeLines(
    c(
      "gene_index_1based\tgene_id",
      "1\tMT-ND1",
      "2\tG1",
      "3\tG2"
    ),
    file.path(bundle_dir, "features.tsv"),
    useBytes = TRUE
  )
  writeLines(
    c(
      "cell_id\tmatrix_cell_id\tglobal_index_1based",
      "A-1\tA.1\t1",
      "B-1\tB.1\t2"
    ),
    file.path(partition_dir, "barcodes.tsv"),
    useBytes = TRUE
  )
  writeBin(
    as.integer(c(1, 3, 5)),
    file.path(partition_dir, "data.i32"),
    size = 4L,
    endian = "little"
  )
  writeBin(
    as.integer(c(0, 1, 0)),
    file.path(partition_dir, "indices.i32"),
    size = 4L,
    endian = "little"
  )
  writeBin(
    as.integer(c(0, 1, 2, 3)),
    file.path(partition_dir, "indptr.i32"),
    size = 4L,
    endian = "little"
  )
  observed_size <- function(path) as.numeric(file.info(path)$size)
  partition <- list(
    partition_id = "P1",
    safe_id = "P1",
    n_cells = 2,
    n_genes = 3,
    nonzero = 3,
    data_file = "P1/data.i32",
    data_bytes = observed_size(file.path(partition_dir, "data.i32")),
    indices_file = "P1/indices.i32",
    indices_bytes = observed_size(
      file.path(partition_dir, "indices.i32")
    ),
    indptr_file = "P1/indptr.i32",
    indptr_bytes = observed_size(file.path(partition_dir, "indptr.i32")),
    barcodes_file = "P1/barcodes.tsv",
    barcodes_bytes = observed_size(
      file.path(partition_dir, "barcodes.tsv")
    )
  )
  manifest <- list(
    format = "gse255460_partitioned_csr",
    format_version = 1,
    orientation = "genes_by_cells",
    sparse_encoding = "CSR",
    index_base = 0,
    value_dtype = "little_endian_int32",
    index_dtype = "little_endian_int32",
    pointer_dtype = "little_endian_int32",
    all_cells_nCount_exact = TRUE,
    all_cells_nFeature_exact = TRUE,
    header_mapping_exact = TRUE,
    features_unique = TRUE,
    n_genes = 3,
    n_cells = 2,
    n_partitions = 1,
    nonzero = 3,
    features_file = "features.tsv",
    partitions = list(partition)
  )
  manifest_path <- file.path(bundle_dir, "manifest.json")
  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  metadata <- data.frame(
    .cell_id = c("A-1", "B-1"),
    ID = rep("P1", 2),
    nCount_RNA = c(6, 3),
    nFeature_RNA = c(2, 1),
    stringsAsFactors = FALSE
  )
  bundle <- list(directory = bundle_dir, manifest_path = manifest_path)
  validated <- .sc_gse255460_validate_manifest(bundle, metadata)
  imported <- .sc_read_gse255460_partition(
    bundle,
    validated,
    "P1",
    metadata
  )
  expected <- matrix(
    c(1, 0, 5, 0, 3, 0),
    nrow = 3,
    dimnames = list(c("MT-ND1", "G1", "G2"), c("A-1", "B-1"))
  )
  testthat::expect_equal(as.matrix(imported$counts), expected)
})

testthat::test_that("adapter errors become machine-readable blockers", {
  result <- single_cell_error_result(
    "TEST",
    "synthetic_adapter",
    simpleError("synthetic failure")
  )
  result <- normalize_single_cell_result(
    result,
    list(id = "TEST"),
    "synthetic_adapter"
  )

  testthat::expect_identical(result$status, "blocked_adapter_error")
  testthat::expect_false(result$downstream_ready)
  testthat::expect_equal(
    names(result$capability),
    c("dataset_id", "adapter", "capability", "status", "reason")
  )
  testthat::expect_match(result$reason, "synthetic failure")
})
