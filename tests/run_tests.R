if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required. Run: Rscript setup.R --install", call. = FALSE)
}

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

r_files <- sort(list.files(
  file.path(project_root, "R"),
  pattern = "\\.R$",
  full.names = TRUE
))
if (length(r_files) == 0L) {
  stop("No project R files were found under R/.", call. = FALSE)
}
invisible(lapply(r_files, source, local = .GlobalEnv))

options(ocoa.test_project_root = project_root)

testthat::test_dir(
  file.path(project_root, "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = TRUE,
  stop_on_warning = FALSE
)
