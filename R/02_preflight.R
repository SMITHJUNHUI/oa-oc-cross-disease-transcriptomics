required_packages_for_config <- function(config, active_stages = NULL) {
  stage_is_active <- function(stage) {
    is.null(active_stages) || stage %in% active_stages
  }
  packages <- c(
    "AnnotationDbi", "Biobase", "data.table", "digest", "GEOquery", "ggplot2",
    "hgu133plus2.db", "limma", "org.Hs.eg.db", "pheatmap", "pROC", "yaml"
  )
  packages <- c(packages, "jsonlite")

  module_packages <- list(
    enrichment = c("clusterProfiler", "enrichplot"),
    wgcna = "WGCNA",
    machine_learning = c("glmnet", "randomForest"),
    validation = "pROC",
    immune = character(),
    tcga = c("glmnet", "survival", "survminer", "timeROC"),
    regulatory = "clusterProfiler",
    drug = character(),
    mr = c("MRPRESSO", "TwoSampleMR"),
    single_cell = c(
      "BiocParallel", "Matrix", "S4Vectors", "scater", "scDblFinder",
      "scuttle", "SingleCellExperiment", "SummarizedExperiment", "rhdf5"
    ),
    single_cell_downstream = c(
      "BiocParallel", "BiocSingular", "bluster", "edgeR", "limma",
      "Matrix", "readxl", "rhdf5", "S4Vectors", "scater", "scran",
      "scuttle", "SingleCellExperiment", "SingleR",
      "SummarizedExperiment", "uwot"
    )
  )

  for (module in names(module_packages)) {
    if (
      isTRUE(config$modules[[module]]) &&
        stage_is_active(module)
    ) {
      packages <- c(packages, module_packages[[module]])
    }
  }
  sort(unique(packages))
}

find_single_cell_python <- function(config) {
  configured <- config$single_cell$python_executable %||% ""
  configured_candidate <- if (nzchar(configured)) {
    if (file.exists(configured)) configured else Sys.which(configured)
  } else {
    ""
  }
  candidates <- unique(c(
    configured_candidate,
    Sys.getenv("OC_OA_PYTHON", ""),
    Sys.which("python3"),
    Sys.which("python")
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  for (candidate in candidates) {
    version_text <- tryCatch(
      system2(
        candidate,
        "--version",
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(error) character()
    )
    match <- regexec(
      "Python[[:space:]]+([0-9]+)\\.([0-9]+)(?:\\.([0-9]+))?",
      paste(version_text, collapse = " "),
      perl = TRUE
    )
    parts <- regmatches(paste(version_text, collapse = " "), match)[[1L]]
    if (length(parts) >= 3L) {
      major <- as.integer(parts[[2L]])
      minor <- as.integer(parts[[3L]])
      if (major > 3L || (major == 3L && minor >= 8L)) {
        numpy_text <- tryCatch(
          system2(
            candidate,
            args = c(
              "-c",
              shQuote("import numpy; print(numpy.__version__)")
            ),
            stdout = TRUE,
            stderr = TRUE
          ),
          error = function(error) character()
        )
        numpy_version <- if (
          length(numpy_text) == 1L &&
            grepl("^[0-9]+[.][0-9]+", numpy_text[[1L]])
        ) {
          numpy_text[[1L]]
        } else {
          NA_character_
        }
        if (is.na(numpy_version)) next
        return(data.frame(
          runtime = "python",
          executable = normalizePath(
            candidate,
            winslash = "/",
            mustWork = TRUE
          ),
          version = sub("^Python[[:space:]]+", "", parts[[1L]]),
          numpy_version = numpy_version,
          requirement = ">=3.8 with numpy",
          status = "available",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  data.frame(
    runtime = "python",
    executable = NA_character_,
    version = NA_character_,
    numpy_version = NA_character_,
    requirement = ">=3.8 with numpy",
    status = "missing_too_old_or_numpy_unavailable",
    stringsAsFactors = FALSE
  )
}

inspect_input_files <- function(config, active_stages = NULL) {
  manifest <- configured_input_manifest(config, active_stages)
  if (nrow(manifest) == 0L) {
    return(manifest)
  }

  info <- file.info(manifest$path)
  manifest$exists <- !is.na(info$size)
  manifest$size_bytes <- as.numeric(info$size)
  manifest$size_human <- vapply(
    manifest$size_bytes,
    human_bytes,
    character(1)
  )
  manifest$modified_at <- as.character(info$mtime)
  manifest$checksum_md5 <- NA_character_

  threshold <- config$project$checksum_max_mb * 1024^2
  checksum_indices <- which(
    manifest$exists &
      !is.na(manifest$size_bytes) &
      manifest$size_bytes <= threshold
  )
  if (length(checksum_indices) > 0L) {
    manifest$checksum_md5[checksum_indices] <- unname(
      tools::md5sum(manifest$path[checksum_indices])
    )
  }
  manifest
}

inspect_packages <- function(config, active_stages = NULL) {
  packages <- required_packages_for_config(config, active_stages)
  installed <- rownames(installed.packages())
  data.frame(
    package = packages,
    installed = packages %in% installed,
    version = vapply(
      packages,
      function(package) {
        if (!package %in% installed) return(NA_character_)
        as.character(utils::packageVersion(package))
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

run_preflight <- function(config, active_stages = NULL) {
  output_dir <- ensure_dir(config$project$output_dir)
  ensure_dir(config$project$cache_dir)
  ensure_dir(file.path(output_dir, "figures"))
  ensure_dir(file.path(output_dir, "tables"))
  ensure_dir(file.path(output_dir, "manifests"))
  ensure_dir(file.path(output_dir, "reports"))

  if (getRversion() < "4.5.0") {
    stop("R >= 4.5.0 is required. Detected: ", R.version.string, call. = FALSE)
  }

  input_manifest <- inspect_input_files(config, active_stages)
  package_manifest <- inspect_packages(config, active_stages)
  python_manifest <- if (
    isTRUE(config$modules$single_cell) &&
      (is.null(active_stages) || "single_cell" %in% active_stages)
  ) {
    find_single_cell_python(config)
  } else {
    data.frame()
  }

  safe_write_csv(
    input_manifest,
    file.path(output_dir, "manifests", "input_manifest.csv")
  )
  safe_write_csv(
    package_manifest,
    file.path(output_dir, "manifests", "package_manifest.csv")
  )
  if (nrow(python_manifest) > 0L) {
    safe_write_csv(
      python_manifest,
      file.path(output_dir, "manifests", "single_cell_runtime_manifest.csv")
    )
  }
  yaml::write_yaml(
    config[setdiff(names(config), c(".project_root", ".config_path"))],
    file.path(output_dir, "manifests", "resolved_config.yml")
  )

  missing_inputs <- input_manifest[
    input_manifest$required & !input_manifest$exists,
    ,
    drop = FALSE
  ]
  missing_packages <- package_manifest[
    !package_manifest$installed,
    ,
    drop = FALSE
  ]

  if (nrow(input_manifest) > 0L) {
    log_info(
      "Input files: ", sum(input_manifest$exists), "/",
      nrow(input_manifest), " present."
    )
  }
  log_info(
    "R packages: ", sum(package_manifest$installed), "/",
    nrow(package_manifest), " present."
  )

  if (isTRUE(config$modules$mr)) {
    token_name <- config$mr$opengwas_token_env %||% "OPENGWAS_JWT"
    if (!nzchar(Sys.getenv(token_name))) {
      stop(
        "MR is enabled, but environment variable ", token_name,
        " is empty.",
        call. = FALSE
      )
    }
    if (
      length(config$mr$exposure_ids %||% character()) == 0L ||
        length(config$mr$outcome_ids %||% character()) == 0L
    ) {
      stop(
        "MR is enabled, but exposure_ids or outcome_ids are empty.",
        call. = FALSE
      )
    }
  }

  if (nrow(missing_inputs) > 0L) {
    message <- paste(
      paste0(missing_inputs$key, " -> ", missing_inputs$path),
      collapse = "\n"
    )
    stop("Required input files are missing:\n", message, call. = FALSE)
  }

  if (nrow(missing_packages) > 0L) {
    stop(
      "Required packages are missing: ",
      paste(missing_packages$package, collapse = ", "),
      ". Run: Rscript setup.R --install",
      call. = FALSE
    )
  }
  if (
    nrow(python_manifest) > 0L &&
      !identical(python_manifest$status[[1L]], "available")
  ) {
    stop(
      "The GSE255460 bounded-memory adapter requires Python >= 3.8 with numpy. ",
      "Install Python/numpy or set single_cell.python_executable / OC_OA_PYTHON.",
      call. = FALSE
    )
  }

  list(
    ok = TRUE,
    checked_at = as.character(Sys.time()),
    r_version = R.version.string,
    inputs = input_manifest,
    packages = package_manifest,
    runtimes = python_manifest
  )
}
