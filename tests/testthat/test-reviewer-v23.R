testthat::test_that("V2.3 candidate evidence summary exposes every requested tier", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v23",
    "supplementary_tables",
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  required <- c(
    "shared_DEG",
    "direction",
    "WGCNA_support",
    "LASSO_support",
    "random_forest_support",
    "strict_nested_frequency",
    "single_cell_context",
    "ten_gene_set_role"
  )
  testthat::expect_equal(nrow(table), 10L)
  testthat::expect_true(all(required %in% names(table)))
  testthat::expect_true(all(grepl(
    "not an optimized predictive signature",
    table$ten_gene_set_role,
    fixed = TRUE
  )))
})

testthat::test_that("V2.3 pathway direction matrix is complete and bounded", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v23",
    "supplementary_tables",
    "Table_S18_Hallmark_pathway_direction_matrix.csv"
  ))
  testthat::expect_equal(nrow(table), 50L)
  testthat::expect_true(all(is.finite(table$OA_NES)))
  testthat::expect_true(all(is.finite(table$OC_NES)))
  testthat::expect_setequal(
    unique(table$direction_class),
    c("concordant", "discordant")
  )
  testthat::expect_true(all(grepl(
    "does not establish a shared mechanism",
    table$inference_boundary,
    fixed = TRUE
  )))
})

testthat::test_that("V2.3 gene-cell-function matrix covers all candidates", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v23",
    "supplementary_tables",
    "Table_S19_gene_cell_function_context_matrix.csv"
  ))
  testthat::expect_equal(nrow(table), 10L)
  testthat::expect_equal(length(unique(table$gene)), 10L)
  testthat::expect_true(all(nzchar(table$OA_cell_context)))
  testthat::expect_true(all(nzchar(table$OC_cell_context)))
  testthat::expect_true(all(nzchar(table$OA_functional_theme)))
  testthat::expect_true(all(nzchar(table$OC_functional_theme)))
})

testthat::test_that("V2.3 figure pack contains eight supplementary figures", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  figure_dir <- file.path(root, "results", "submission_v23", "figures")
  main <- list.files(figure_dir, pattern = "^Figure[1-6]_.*\\.png$")
  supplementary <- list.files(
    figure_dir,
    pattern = "^SupplementaryFigure[1-8]_.*\\.png$"
  )
  testthat::expect_equal(length(main), 6L)
  testthat::expect_equal(length(supplementary), 8L)
  testthat::expect_true(all(file.info(file.path(
    figure_dir,
    c(main, supplementary)
  ))$size > 10000))
})
