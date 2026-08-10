testthat::test_that("gene normalization and aggregation are deterministic", {
  expression <- matrix(
    c(1, 3, 2, 4, 10, 12),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("p1", "p2", "p3"), c("s1", "s2"))
  )
  aggregated <- aggregate_expression_by_gene(
    expression,
    c("GeneA", "GENEA", "GeneB")
  )
  testthat::expect_equal(rownames(aggregated), c("GENEA", "GENEB"))
  testthat::expect_equal(unname(aggregated["GENEA", ]), c(1.5, 3.5))
  testthat::expect_equal(unname(aggregated["GENEB", ]), c(10, 12))
})

testthat::test_that("group patterns classify every sample exactly once", {
  metadata <- data.frame(
    title = c("Normal sample 1", "Disease sample 1"),
    row.names = c("S1", "S2")
  )
  dataset <- list(
    id = "TEST",
    normal_patterns = "^Normal",
    disease_patterns = "^Disease"
  )
  group <- derive_groups(metadata, dataset)
  testthat::expect_equal(as.character(group), c("Normal", "Disease"))
  testthat::expect_equal(names(group), c("S1", "S2"))
})

testthat::test_that("automatic transformation leaves log-scale data unchanged", {
  expression <- matrix(seq(1, 16, length.out = 40), nrow = 10)
  testthat::expect_equal(
    maybe_transform_expression(expression, "auto"),
    expression
  )
})

