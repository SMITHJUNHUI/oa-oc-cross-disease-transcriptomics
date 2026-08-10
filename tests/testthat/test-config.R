testthat::test_that("project config is valid and resolves paths", {
  project_root <- getOption("ocoa.test_project_root")
  config <- read_project_config(project_root, "config/config.example.yml")
  testthat::expect_equal(config$project$name, "OC_OA_cross_disease_study")
  testthat::expect_true(is_absolute_path(config$drug$dgidb_path))
  raw_config <- readLines(
    file.path(project_root, "config", "config.example.yml"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  user_profile_pattern <- paste0("C:", "/Users/")
  absolute_drive_pattern <- "^[[:space:]]+[^#]+: \"[A-Z]:/"
  testthat::expect_false(any(grepl(
    paste(user_profile_pattern, absolute_drive_pattern, sep = "|"),
    raw_config
  )))
  testthat::expect_true(all(c("OA", "OC") %in% vapply(
    config$datasets,
    `[[`,
    character(1),
    "disease"
  )))
})

testthat::test_that("configured input manifest has unique keys", {
  project_root <- getOption("ocoa.test_project_root")
  config <- read_project_config(project_root, "config/config.example.yml")
  manifest <- configured_input_manifest(config)
  testthat::expect_gt(nrow(manifest), 10L)
  testthat::expect_equal(
    anyDuplicated(paste(manifest$category, manifest$key)),
    0L
  )
  testthat::expect_true(all(manifest$required))
})

testthat::test_that("core and early-stage manifests exclude downstream inputs", {
  project_root <- getOption("ocoa.test_project_root")
  config <- read_project_config(project_root, "config/config.yml")
  registry <- pipeline_stage_registry(config)
  core_stages <- names(registry)[
    seq_len(match("validation", names(registry)))
  ]
  core_manifest <- configured_input_manifest(config, core_stages)
  testthat::expect_false(any(
    core_manifest$module %in% c("tcga", "regulatory", "drug")
  ))
  testthat::expect_true(any(grepl("oa_validation", core_manifest$key)))

  differential_stages <- names(registry)[
    seq_len(match("differential", names(registry)))
  ]
  early_manifest <- configured_input_manifest(config, differential_stages)
  testthat::expect_false(any(grepl("validation", early_manifest$key)))
})

testthat::test_that("report summary ignores disabled dataset placeholders", {
  expression <- matrix(
    1:8,
    nrow = 2,
    dimnames = list(c("G1", "G2"), paste0("S", 1:4))
  )
  dataset <- list(
    id = "SYNTHETIC",
    disease = "OA",
    role = "train",
    expression = expression,
    group = factor(
      c("Normal", "Normal", "Disease", "Disease"),
      levels = c("Normal", "Disease")
    )
  )
  state <- list(
    bulk_training = list(dataset),
    bulk_validation = list(status = "disabled")
  )
  lines <- summarize_pipeline_results(state)
  testthat::expect_true(any(grepl("Loaded datasets: 1", lines, fixed = TRUE)))
})

testthat::test_that("renv lockfile is strict JSON with single-cell dependencies", {
  project_root <- getOption("ocoa.test_project_root")
  lock_path <- file.path(project_root, "renv.lock")
  lock_text <- paste(
    readLines(lock_path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  testthat::expect_true(jsonlite::validate(lock_text))
  lock <- jsonlite::fromJSON(lock_text, simplifyVector = FALSE)
  required <- c(
    "BiocParallel", "bluster", "scater", "scDblFinder", "scran", "scuttle",
    "SingleCellExperiment", "SummarizedExperiment", "xgboost"
  )
  testthat::expect_true(all(required %in% names(lock$Packages)))
})
