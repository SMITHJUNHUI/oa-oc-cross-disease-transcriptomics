testthat::test_that("training-fold screening ranks a training-only signal", {
  expression <- rbind(
    signal = c(0, 0.1, -0.1, 3, 3.1, 2.9),
    noise_a = c(1, 0, 1, 0, 1, 0),
    noise_b = c(2, 2, 1.9, 2, 2.1, 2)
  )
  status <- c(0L, 0L, 0L, 1L, 1L, 1L)
  screen <- submission_training_fold_screen(
    expression,
    status,
    top_n = 2L
  )
  testthat::expect_equal(screen$gene[[1L]], "signal")
  testthat::expect_equal(screen$screening_rank, 1:2)
})

testthat::test_that("strict nested results declare leakage-safe scope", {
  project_root <- if (file.exists("DESCRIPTION")) "." else file.path("..", "..")
  path <- file.path(
    project_root,
    "results",
    "submission",
    "sensitivity",
    "machine_learning_repeated_cv_summary.csv"
  )
  testthat::expect_true(file.exists(path))
  table <- utils::read.csv(path, stringsAsFactors = FALSE)
  testthat::expect_true(all(table$candidate_space == "all measured genes"))
  testthat::expect_true(all(
    table$feature_selection_scope == "outer training fold only"
  ))
})

testthat::test_that("external robustness and tissue-context boundaries persist", {
  project_root <- if (file.exists("DESCRIPTION")) "." else file.path("..", "..")
  sensitivity <- file.path(
    project_root,
    "results",
    "submission",
    "sensitivity"
  )
  external <- utils::read.csv(
    file.path(sensitivity, "external_validation_signed_composite_score.csv"),
    stringsAsFactors = FALSE
  )
  oc <- external[external$disease == "OC", , drop = FALSE]
  testthat::expect_true(all(oc$permutation_empirical_p <= 1 / 1001))
  testthat::expect_true(all(oc$leave_one_out_auc_minimum >= 0.97))

  hpa <- utils::read.csv(
    file.path(sensitivity, "hpa_normal_tissue_context.csv"),
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(sum(hpa$ovary_listed_as_specific), 1)
  testthat::expect_true(all(!hpa$cartilage_in_reference))
})
