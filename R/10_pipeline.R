pipeline_stage_registry <- function(config) {
  list(
    preflight = list(enabled = TRUE),
    bulk_training = list(enabled = TRUE),
    differential = list(enabled = TRUE),
    shared = list(enabled = TRUE),
    enrichment = list(enabled = isTRUE(config$modules$enrichment)),
    wgcna = list(enabled = isTRUE(config$modules$wgcna)),
    machine_learning = list(enabled = isTRUE(config$modules$machine_learning)),
    bulk_validation = list(enabled = isTRUE(config$modules$validation)),
    validation = list(enabled = isTRUE(config$modules$validation)),
    immune = list(enabled = isTRUE(config$modules$immune)),
    tcga = list(enabled = isTRUE(config$modules$tcga)),
    regulatory = list(enabled = isTRUE(config$modules$regulatory)),
    drug = list(enabled = isTRUE(config$modules$drug)),
    mr = list(enabled = isTRUE(config$modules$mr)),
    single_cell = list(enabled = isTRUE(config$modules$single_cell)),
    single_cell_downstream = list(
      enabled = isTRUE(config$modules$single_cell_downstream)
    ),
    report = list(enabled = TRUE)
  )
}

run_named_stage <- function(name, state, config, stage_status) {
  switch(
    name,
    preflight = run_preflight(config, config$.active_stages),
    bulk_training = load_bulk_datasets(config, role = "train"),
    differential = run_differential_stage(state$bulk_training, config),
    shared = derive_shared_candidates(state$differential, config),
    enrichment = run_enrichment_stage(
      state$shared,
      state$differential,
      config
    ),
    wgcna = run_wgcna_stage(state$bulk_training, config),
    machine_learning = run_machine_learning_stage(
      state$bulk_training,
      state$shared,
      state$wgcna %||% list(status = "disabled"),
      config
    ),
    bulk_validation = load_bulk_datasets(config, role = "validation"),
    validation = run_validation_stage(
      state$bulk_validation,
      state$machine_learning %||% list(final_genes = state$shared$genes),
      state$shared,
      config
    ),
    immune = run_immune_stage(state$bulk_training, config),
    tcga = run_tcga_stage(
      state$machine_learning %||% list(final_genes = state$shared$genes),
      state$shared,
      config
    ),
    regulatory = run_regulatory_stage(
      state$machine_learning %||% list(final_genes = state$shared$genes),
      state$shared,
      config
    ),
    drug = run_drug_stage(
      state$machine_learning %||% list(final_genes = state$shared$genes),
      state$shared,
      config
    ),
    mr = run_mr_stage(config),
    single_cell = run_single_cell_stage(config),
    single_cell_downstream = run_single_cell_downstream_stage(
      state$single_cell,
      state$machine_learning %||% list(),
      config
    ),
    report = generate_final_report(state, stage_status, config),
    stop("No implementation for pipeline stage: ", name, call. = FALSE)
  )
}

project_fingerprint <- function(config) {
  require_namespace("digest", "pipeline checkpoint fingerprints")
  code_files <- sort(c(
    list.files(
      file.path(config$.project_root, "R"),
      pattern = "\\.[Rr]$",
      full.names = TRUE
    ),
    list.files(
      file.path(config$.project_root, "tools"),
      pattern = "\\.py$",
      full.names = TRUE
    )
  ))
  code_hashes <- unname(tools::md5sum(code_files))
  code_names <- sub(
    paste0(
      "^",
      gsub(
        "([][{}()+*^$|\\\\?.])",
        "\\\\\\1",
        normalizePath(
          config$.project_root,
          winslash = "/",
          mustWork = TRUE
        )
      ),
      "/"
    ),
    "",
    normalizePath(code_files, winslash = "/", mustWork = TRUE)
  )
  input_manifest <- inspect_input_files(config, config$.active_stages)
  input_identity <- if (nrow(input_manifest) == 0L) {
    data.frame()
  } else {
    input_manifest[, c("path", "exists", "size_bytes", "modified_at"), drop = FALSE]
  }
  digest::digest(
    list(
      config = readLines(config$.config_path, warn = FALSE, encoding = "UTF-8"),
      code = stats::setNames(code_hashes, code_names),
      inputs = input_identity,
      r = paste(R.version$major, R.version$minor, sep = ".")
    ),
    algo = "xxhash64"
  )
}

checkpoint_path <- function(config, stage) {
  file.path(
    ensure_dir(config$project$cache_dir),
    paste0(sprintf("%02d", match(
      stage,
      names(pipeline_stage_registry(config))
    )), "_", stage, ".rds")
  )
}

load_valid_checkpoint <- function(path, fingerprint) {
  if (!file.exists(path)) return(NULL)
  checkpoint <- tryCatch(readRDS(path), error = function(error) NULL)
  if (
    is.null(checkpoint) ||
      !identical(checkpoint$fingerprint, fingerprint) ||
      is.null(checkpoint$value)
  ) {
    return(NULL)
  }
  checkpoint
}

write_pipeline_status <- function(stage_status, config) {
  if (length(stage_status) == 0L) return(invisible(NULL))
  safe_write_csv(
    do.call(rbind, stage_status),
    file.path(config$project$output_dir, "manifests", "pipeline_status.csv")
  )
}

run_pipeline <- function(
    project_root,
    config_path = "config/config.yml",
    mode = "full",
    force = FALSE,
    from = "",
    to = ""
) {
  config <- read_project_config(project_root, config_path)
  set.seed(config$project$seed)
  Sys.setenv(TZ = config$project$timezone %||% "UTC")
  ensure_dir(config$project$output_dir)
  log_file <- initialize_logging(config$project$output_dir)
  log_info("Project: ", config$project$name, ".")
  log_info("Mode: ", mode, "; log: ", log_file, ".")

  registry <- pipeline_stage_registry(config)
  stage_names <- names(registry)
  default_to <- switch(
    mode,
    preflight = "preflight",
    core = "validation",
    full = "report",
    stop("Unsupported pipeline mode: ", mode, call. = FALSE)
  )
  final_stage <- if (nzchar(to)) to else default_to
  if (!final_stage %in% stage_names) {
    stop("Unknown --to stage: ", final_stage, call. = FALSE)
  }
  first_force_stage <- if (nzchar(from)) from else stage_names[[1L]]
  if (!first_force_stage %in% stage_names) {
    stop("Unknown --from stage: ", first_force_stage, call. = FALSE)
  }
  final_index <- match(final_stage, stage_names)
  first_force_index <- match(first_force_stage, stage_names)
  stages <- stage_names[seq_len(final_index)]
  config$.active_stages <- if (identical(mode, "preflight")) {
    stage_names[vapply(registry, function(stage) isTRUE(stage$enabled), logical(1))]
  } else {
    stages
  }

  global_fingerprint <- project_fingerprint(config)
  state <- list()
  stage_status <- list()
  upstream_fingerprints <- character()

  for (index in seq_along(stages)) {
    stage <- stages[[index]]
    enabled <- isTRUE(registry[[stage]]$enabled)
    if (!enabled) {
      log_info("Skipping disabled stage: ", stage, ".")
      state[[stage]] <- list(status = "disabled")
      stage_status[[length(stage_status) + 1L]] <- data.frame(
        stage = stage,
        status = "skipped",
        source = "config",
        started_at = as.character(Sys.time()),
        completed_at = as.character(Sys.time()),
        seconds = 0,
        details = "disabled",
        stringsAsFactors = FALSE
      )
      next
    }

    fingerprint <- digest::digest(
      list(
        global = global_fingerprint,
        stage = stage,
        upstream = upstream_fingerprints
      ),
      algo = "xxhash64"
    )
    path <- checkpoint_path(config, stage)
    should_force <- isTRUE(force) && index >= first_force_index
    checkpoint <- if (!should_force && !identical(stage, "preflight") && !identical(stage, "report")) {
      load_valid_checkpoint(path, fingerprint)
    } else {
      NULL
    }
    started <- Sys.time()

    if (!is.null(checkpoint)) {
      log_info("Loading cached stage: ", stage, ".")
      state[[stage]] <- checkpoint$value
      completed <- Sys.time()
      source <- "checkpoint"
      status <- "completed"
      details <- "cache hit"
    } else {
      log_info("Running stage: ", stage, ".")
      value <- tryCatch(
        run_named_stage(stage, state, config, stage_status),
        error = function(error) {
          completed <- Sys.time()
          stage_status[[length(stage_status) + 1L]] <<- data.frame(
            stage = stage,
            status = "failed",
            source = "computed",
            started_at = as.character(started),
            completed_at = as.character(completed),
            seconds = as.numeric(difftime(completed, started, units = "secs")),
            details = compact_error(error),
            stringsAsFactors = FALSE
          )
          write_pipeline_status(stage_status, config)
          log_error("Stage ", stage, " failed: ", conditionMessage(error))
          stop(error)
        }
      )
      state[[stage]] <- value
      completed <- Sys.time()
      source <- "computed"
      status <- "completed"
      details <- ""
      if (!identical(stage, "preflight") && !identical(stage, "report")) {
        atomic_save_rds(
          list(
            fingerprint = fingerprint,
            stage = stage,
            completed_at = as.character(completed),
            value = value
          ),
          path
        )
      }
    }

    stage_status[[length(stage_status) + 1L]] <- data.frame(
      stage = stage,
      status = status,
      source = source,
      started_at = as.character(started),
      completed_at = as.character(completed),
      seconds = as.numeric(difftime(completed, started, units = "secs")),
      details = details,
      stringsAsFactors = FALSE
    )
    upstream_fingerprints[[stage]] <- fingerprint
    write_pipeline_status(stage_status, config)
  }

  write_session_information(config)
  log_info("Pipeline completed through stage: ", final_stage, ".")
  invisible(list(state = state, status = stage_status, config = config))
}
