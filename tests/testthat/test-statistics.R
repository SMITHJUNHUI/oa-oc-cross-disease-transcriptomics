testthat::test_that("limma recovers a strong synthetic disease effect", {
  set.seed(1)
  expression <- matrix(
    stats::rnorm(100 * 12, sd = 0.2),
    nrow = 100,
    dimnames = list(paste0("G", seq_len(100)), paste0("S", seq_len(12)))
  )
  group <- factor(
    rep(c("Normal", "Disease"), each = 6),
    levels = c("Normal", "Disease")
  )
  names(group) <- colnames(expression)
  expression[1:5, group == "Disease"] <-
    expression[1:5, group == "Disease"] + 2

  result <- run_limma_differential(expression, group)
  top <- head(result$gene, 5)
  testthat::expect_true(all(paste0("G", 1:5) %in% top))
  testthat::expect_true(all(result$logFC[match(paste0("G", 1:5), result$gene)] > 1.5))
})

testthat::test_that("rank signature score preserves sample names", {
  expression <- matrix(
    seq_len(30),
    nrow = 10,
    dimnames = list(paste0("G", 1:10), c("S1", "S2", "S3"))
  )
  scores <- rank_signature_scores(
    expression,
    list(SetA = c("G1", "G2"), SetB = c("G9", "G10"))
  )
  testthat::expect_equal(rownames(scores), c("SetA", "SetB"))
  testthat::expect_equal(colnames(scores), c("S1", "S2", "S3"))
  testthat::expect_true(all(scores["SetB", ] > scores["SetA", ]))
})

testthat::test_that("external ROC direction is fixed from training logFC", {
  expression <- rbind(
    UP = c(1, 2, 3, 4, 8, 9, 10, 11),
    DOWN = c(11, 10, 9, 8, 4, 3, 2, 1)
  )
  group <- factor(
    rep(c("Normal", "Disease"), each = 4),
    levels = c("Normal", "Disease")
  )

  up <- evaluate_gene_roc(expression, group, "UP", training_logfc = 2)
  down <- evaluate_gene_roc(expression, group, "DOWN", training_logfc = -2)

  testthat::expect_equal(up$row$auc, 1)
  testthat::expect_equal(down$row$auc, 1)
  testthat::expect_equal(up$row$direction, "<")
  testthat::expect_equal(up$row$expected_expression, "higher_in_disease")
  testthat::expect_equal(down$row$expected_expression, "lower_in_disease")
  testthat::expect_error(
    evaluate_gene_roc(expression, group, "UP", training_logfc = 0),
    "finite, non-zero"
  )
})
