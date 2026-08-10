testthat::test_that("V2.4 candidate-centered pathway context is complete", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v24",
    "supplementary_tables",
    "Table_S20_candidate_centered_Hallmark_context.csv"
  ))
  testthat::expect_equal(nrow(table), 400L)
  testthat::expect_setequal(
    unique(table$candidate),
    c("SOX9", "DDIT3", "BNC1", "AKAP12")
  )
  testthat::expect_true(all(table(table$candidate, table$disease) == 50L))
  testthat::expect_true(all(nzchar(table$calculation_status)))
  testthat::expect_true(all(grepl(
    "not single-gene perturbation",
    table$inference_boundary,
    fixed = TRUE
  )))
})

testthat::test_that("V2.4 external validation distinguishes task context", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v24",
    "supplementary_tables",
    "Table_S21_cross_cohort_molecular_separability_context.csv"
  ))
  testthat::expect_equal(nrow(table), 4L)
  testthat::expect_setequal(
    table$dataset_id,
    c("GSE117999", "GSE82107", "GSE54388", "GSE12470")
  )
  testthat::expect_true(all(nzchar(table$tissue_and_comparator)))
  testthat::expect_true(all(nzchar(table$validation_task_scale)))
  testthat::expect_true(all(grepl(
    "not clinical utility",
    table$inference_boundary,
    fixed = TRUE
  )))
})

testthat::test_that("V2.4 calibration reports instability rather than hiding it", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  table <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v24",
    "supplementary_tables",
    "Table_S22b_cross_fitted_calibration_metrics.csv"
  ))
  testthat::expect_equal(nrow(table), 4L)
  testthat::expect_true(all(is.finite(table$Brier_score)))
  testthat::expect_true(all(table$Brier_score >= 0 & table$Brier_score <= 1))
  testthat::expect_equal(
    sum(table$calibration_status == "estimable descriptive sensitivity"),
    1L
  )
  testthat::expect_true(all(grepl(
    "not external probability calibration",
    table$inference_boundary,
    fixed = TRUE
  )))
})

testthat::test_that("V2.4 regulatory context retains explicit evidence boundaries", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  tf <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v24",
    "supplementary_tables",
    "Table_S23a_KnockTF_candidate_regulatory_context.csv"
  ))
  mirna <- utils::read.csv(file.path(
    root,
    "results",
    "submission_v24",
    "supplementary_tables",
    "Table_S23b_miRTarBase_candidate_regulatory_context.csv"
  ))
  testthat::expect_equal(nrow(tf), 157L)
  testthat::expect_equal(nrow(mirna), 363L)
  testthat::expect_true(all(grepl(
    "does not establish TF activity",
    tf$inference_boundary,
    fixed = TRUE
  )))
  testthat::expect_true(all(
    mirna$disease_specific_evidence == "not available in this study"
  ))
})

testthat::test_that("V2.4 figure pack contains eleven supplementary figures", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  figure_dir <- file.path(root, "results", "submission_v24", "figures")
  main <- list.files(figure_dir, pattern = "^Figure[1-6]_.*\\.png$")
  supplementary <- list.files(
    figure_dir,
    pattern = "^SupplementaryFigure([1-9]|10|11)_.*\\.png$"
  )
  testthat::expect_equal(length(main), 6L)
  testthat::expect_equal(length(supplementary), 11L)
  testthat::expect_true(all(file.info(file.path(
    figure_dir,
    c(main, supplementary)
  ))$size > 10000))
})
