testthat::test_that("V2.1 directional quadrants reproduce the primary overlap", {
  table <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "direction_quadrant_counts.csv"
  ))
  testthat::expect_equal(sum(table$genes), 286)
  testthat::expect_equal(sum(table$genes[table$concordant]), 146)
  testthat::expect_equal(sum(table$genes[!table$concordant]), 140)
})

testthat::test_that("GSE54388 PCA is unsupervised and uses all samples", {
  scores <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "GSE54388_unsupervised_PCA_scores.csv"
  ))
  manifest <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "GSE54388_unsupervised_PCA_manifest.csv"
  ))
  testthat::expect_equal(nrow(scores), 22)
  testthat::expect_equal(unname(table(scores$group)[["Normal"]]), 6)
  testthat::expect_equal(unname(table(scores$group)[["Disease"]]), 16)
  testthat::expect_false(manifest$labels_used_in_PCA[[1L]])
  testthat::expect_equal(manifest$genes_entered[[1L]], 2000)
})

testthat::test_that("MR provenance retains the verified accessions and parameters", {
  table <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "supplementary_tables",
    "Table_S12a_MR_estimates_and_provenance.csv"
  ), check.names = FALSE)
  ids <- unique(c(table$id.exposure, table$id.outcome))
  testthat::expect_setequal(ids, c("ebi-a-GCST007092", "ieu-a-1120"))
  testthat::expect_false(any(grepl("GCST90018888|GCST90038686", ids)))
  testthat::expect_setequal(unique(table$nsnp), c(21, 11))
  testthat::expect_true(all(table$instrument_p_threshold == 5e-8))
  testthat::expect_true(all(table$clump_r2 == 0.001))
  testthat::expect_true(all(table$clump_window_kb == 10000))
})

testthat::test_that("cell-context matrix covers both diseases and all candidates", {
  table <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "single_cell_gene_disease_context_matrix.csv"
  ))
  testthat::expect_equal(nrow(table), 20)
  testthat::expect_setequal(unique(table$disease), c("OA", "OC"))
  testthat::expect_equal(length(unique(table$gene)), 10)
  testthat::expect_true(all(table$fraction_detected >= 0))
  testthat::expect_true(all(table$fraction_detected <= 1))
})

testthat::test_that("TCGA context audit is relative and purity is not reported", {
  scores <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "TCGA_OV_relative_context_scores.csv"
  ))
  correlations <- utils::read.csv(file.path(
    project_root,
    "results",
    "submission_v21",
    "analysis",
    "TCGA_OV_hub_gene_context_correlations.csv"
  ))
  testthat::expect_equal(nrow(scores), 307)
  testthat::expect_false("purity" %in% names(scores))
  testthat::expect_equal(nrow(correlations), 27)
  testthat::expect_setequal(
    unique(correlations$gene),
    c(
      "SOX9", "ELF3", "JUNB", "AKAP12", "BNC1",
      "CFI", "DDIT3", "DIRAS3", "HK2"
    )
  )
  testthat::expect_true(all(is.finite(correlations$spearman_rho)))
})

testthat::test_that("V2.1 figure pack contains six main and five supplementary figures", {
  figure_dir <- file.path(
    project_root,
    "results",
    "submission_v21",
    "figures"
  )
  main <- list.files(figure_dir, pattern = "^Figure[1-6]_.*\\.png$")
  supplementary <- list.files(
    figure_dir,
    pattern = "^SupplementaryFigure[1-5]_.*\\.png$"
  )
  testthat::expect_equal(length(main), 6)
  testthat::expect_equal(length(supplementary), 5)
  testthat::expect_true(all(file.info(file.path(figure_dir, c(
    main,
    supplementary
  )))$size > 10000))
})
