testthat::test_that("V2.2 candidate hierarchy and cell-context annotation are complete", {
  root <- normalizePath(
    file.path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  submission <- file.path(root, "results", "submission_v22")
  candidates <- utils::read.csv(file.path(
    submission,
    "supplementary_tables",
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  cell_go <- utils::read.csv(file.path(
    submission,
    "supplementary_tables",
    "Table_S17_cell_type_marker_GO_annotation.csv"
  ))
  testthat::expect_equal(nrow(candidates), 10L)
  testthat::expect_true(all(candidates$shared_primary_DEG_in_OA))
  testthat::expect_true(all(candidates$shared_primary_DEG_in_OC))
  testthat::expect_equal(sort(unique(cell_go$disease)), c("OA", "OC"))
  testthat::expect_true(all(
    cell_go$analysis_scope == "exploratory descriptive functional annotation"
  ))
})

testthat::test_that("V2.2 figure pack contains the new functional-context figure", {
  root <- normalizePath(
    file.path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  figure_dir <- file.path(root, "results", "submission_v22", "figures")
  expected <- file.path(
    figure_dir,
    paste0(
      "SupplementaryFigure6_cell_type_functional_annotation.",
      c("pdf", "png")
    )
  )
  testthat::expect_true(all(file.exists(expected)))
})
