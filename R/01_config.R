read_project_config <- function(project_root, config_path = "config/config.yml") {
  require_namespace("yaml", "configuration parsing")

  resolved_config <- resolve_path(config_path, project_root, must_work = TRUE)
  config <- yaml::read_yaml(resolved_config)
  if (!is.list(config)) {
    stop("The configuration file did not decode to a YAML mapping.", call. = FALSE)
  }

  config$.project_root <- normalizePath(
    project_root,
    winslash = "/",
    mustWork = TRUE
  )
  config$.config_path <- resolved_config

  config$project <- config$project %||% list()
  config$project$name <- config$project$name %||% "OC_OA_cross_disease_study"
  config$project$seed <- as.integer(config$project$seed %||% 20260726L)
  config$project$strict <- isTRUE(config$project$strict %||% TRUE)
  config$project$checksum_max_mb <- as.numeric(
    config$project$checksum_max_mb %||% 50
  )
  config$project$output_dir <- resolve_path(
    config$project$output_dir %||% "results",
    config$.project_root
  )
  config$project$cache_dir <- resolve_path(
    config$project$cache_dir %||% "results/cache",
    config$.project_root
  )

  config$modules <- config$modules %||% list()
  module_names <- c(
    "enrichment", "wgcna", "machine_learning", "validation", "immune",
    "tcga", "regulatory", "drug", "mr", "single_cell",
    "single_cell_downstream"
  )
  for (name in module_names) {
    config$modules[[name]] <- isTRUE(config$modules[[name]])
  }

  config$datasets <- config$datasets %||% list()
  if (length(config$datasets) == 0L) {
    stop("No datasets are defined under 'datasets:' in the config.", call. = FALSE)
  }
  if (is.null(names(config$datasets)) || any(!nzchar(names(config$datasets)))) {
    stop("Every dataset must have a unique YAML key.", call. = FALSE)
  }

  allowed_loaders <- c("split_tables", "geo_series")
  allowed_roles <- c("train", "validation")
  allowed_diseases <- c("OA", "OC")

  for (name in names(config$datasets)) {
    dataset <- config$datasets[[name]]
    dataset$id <- dataset$id %||% name
    dataset$required <- isTRUE(dataset$required %||% TRUE)
    if (!dataset$loader %in% allowed_loaders) {
      stop("Dataset ", name, " has an unsupported loader.", call. = FALSE)
    }
    if (!dataset$role %in% allowed_roles) {
      stop("Dataset ", name, " must have role train or validation.", call. = FALSE)
    }
    if (!dataset$disease %in% allowed_diseases) {
      stop("Dataset ", name, " must have disease OA or OC.", call. = FALSE)
    }
    required_dataset_fields <- if (identical(dataset$loader, "split_tables")) {
      c("case_path", "control_path", "gene_column")
    } else {
      c("series_path", "normal_patterns", "disease_patterns")
    }
    missing_dataset_fields <- required_dataset_fields[
      !vapply(
        required_dataset_fields,
        function(field) {
          value <- dataset[[field]]
          !is.null(value) && length(value) > 0L &&
            all(!is.na(value)) && all(nzchar(as.character(value)))
        },
        logical(1)
      )
    ]
    if (length(missing_dataset_fields) > 0L) {
      stop(
        "Dataset ", name, " is missing required fields: ",
        paste(missing_dataset_fields, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (
      identical(dataset$loader, "geo_series") &&
        identical(dataset$annotation_method, "gpl_soft") &&
        (is.null(dataset$annotation_path) || !nzchar(dataset$annotation_path))
    ) {
      stop(
        "Dataset ", name,
        " uses annotation_method=gpl_soft but has no annotation_path.",
        call. = FALSE
      )
    }

    path_fields <- c(
      "case_path", "control_path", "series_path", "annotation_path"
    )
    for (field in intersect(path_fields, names(dataset))) {
      dataset[[field]] <- resolve_path(
        dataset[[field]],
        config$.project_root
      )
    }
    config$datasets[[name]] <- dataset
  }

  train_diseases <- vapply(
    config$datasets,
    function(x) if (identical(x$role, "train")) x$disease else NA_character_,
    character(1)
  )
  train_diseases <- stats::na.omit(train_diseases)
  if (!all(c("OA", "OC") %in% train_diseases)) {
    stop("Exactly one or more OA and OC training datasets are required.", call. = FALSE)
  }
  if (isTRUE(config$modules$validation)) {
    validation_count <- sum(vapply(
      config$datasets,
      function(dataset) identical(dataset$role, "validation"),
      logical(1)
    ))
    if (validation_count == 0L) {
      stop(
        "validation is enabled, but no validation datasets are defined.",
        call. = FALSE
      )
    }
  }

  for (
    section in c(
      "gene_sets", "tcga", "regulatory", "drug",
      "single_cell_downstream"
    )
  ) {
    if (is.list(config[[section]])) {
      for (field in names(config[[section]])) {
        if (endsWith(field, "_path") && is.character(config[[section]][[field]])) {
          config[[section]][[field]] <- resolve_path(
            config[[section]][[field]],
            config$.project_root
          )
        }
      }
    }
  }

  config$single_cell <- config$single_cell %||% list()
  config$single_cell$scope <- config$single_cell$scope %||% "qc_gate_only"
  config$single_cell$allow_partial <- isTRUE(
    config$single_cell$allow_partial %||% FALSE
  )
  config$single_cell$minimum_validated_datasets <- as.integer(
    config$single_cell$minimum_validated_datasets %||% 1L
  )
  config$single_cell$lower_nmads <- as.numeric(
    config$single_cell$lower_nmads %||% 3
  )
  config$single_cell$upper_nmads <- as.numeric(
    config$single_cell$upper_nmads %||% 5
  )
  config$single_cell$mt_nmads <- as.numeric(
    config$single_cell$mt_nmads %||% 3
  )
  config$single_cell$python_executable <- as.character(
    config$single_cell$python_executable %||% ""
  )
  if (is.list(config$single_cell$datasets)) {
    config$single_cell$datasets <- lapply(
      config$single_cell$datasets,
      function(dataset) {
        dataset$required <- isTRUE(dataset$required %||% TRUE)
        for (field in names(dataset)) {
          if (
            grepl("(_path|_archive)$", field) &&
              is.character(dataset[[field]])
          ) {
            dataset[[field]] <- resolve_path(
              dataset[[field]],
              config$.project_root
            )
          }
        }
        dataset
      }
    )
  }

  if (isTRUE(config$modules$single_cell)) {
    single_cell_datasets <- config$single_cell$datasets %||% list()
    if (length(single_cell_datasets) == 0L) {
      stop(
        "single_cell is enabled, but no single_cell.datasets are defined.",
        call. = FALSE
      )
    }
    dataset_ids <- vapply(
      single_cell_datasets,
      function(dataset) as.character(dataset$id %||% ""),
      character(1)
    )
    if (any(!nzchar(dataset_ids)) || anyDuplicated(dataset_ids)) {
      stop(
        "Every single-cell dataset must have a non-empty unique id.",
        call. = FALSE
      )
    }
    adapter_fields <- list(
      gse255460_wide_counts = c("metadata_path", "counts_archive"),
      gse104782_umi_counts = "counts_path",
      gse169454_tenx_raw_tar = "archive_path",
      gse180661_h5_csr = c("cells_path", "matrix_h5_path"),
      tenx_tar_gse154600 = "archive_path"
    )
    for (index in seq_along(single_cell_datasets)) {
      dataset <- single_cell_datasets[[index]]
      adapter <- as.character(dataset$adapter %||% "")
      if (!adapter %in% names(adapter_fields)) {
        stop(
          "Single-cell dataset ", dataset$id,
          " has unsupported adapter '", adapter, "'.",
          call. = FALSE
        )
      }
      required_fields <- adapter_fields[[adapter]]
      missing_fields <- required_fields[!vapply(
        required_fields,
        function(field) {
          value <- dataset[[field]]
          !is.null(value) && length(value) == 1L &&
            !is.na(value) && nzchar(as.character(value))
        },
        logical(1)
      )]
      if (length(missing_fields) > 0L) {
        stop(
          "Single-cell dataset ", dataset$id,
          " is missing fields required by ", adapter, ": ",
          paste(missing_fields, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }
    if (
      !identical(config$single_cell$scope, "qc_gate_only") ||
        config$single_cell$minimum_validated_datasets < 1L
    ) {
      stop(
        "The current single-cell release supports scope=qc_gate_only and ",
        "minimum_validated_datasets >= 1.",
        call. = FALSE
      )
    }
  }

  config$single_cell_downstream <- config$single_cell_downstream %||% list()
  config$single_cell_downstream$hvg_n <- as.integer(
    config$single_cell_downstream$hvg_n %||% 3000L
  )
  config$single_cell_downstream$pca_n <- as.integer(
    config$single_cell_downstream$pca_n %||% 30L
  )
  config$single_cell_downstream$neighbors_k <- as.integer(
    config$single_cell_downstream$neighbors_k %||% 20L
  )
  config$single_cell_downstream$visualization_max_cells <- as.integer(
    config$single_cell_downstream$visualization_max_cells %||% 60000L
  )
  config$single_cell_downstream$minimum_pseudobulk_cells <- as.integer(
    config$single_cell_downstream$minimum_pseudobulk_cells %||% 20L
  )
  config$single_cell_downstream$minimum_group_replicates <- as.integer(
    config$single_cell_downstream$minimum_group_replicates %||% 3L
  )
  if (isTRUE(config$modules$single_cell_downstream)) {
    if (!isTRUE(config$modules$single_cell)) {
      stop(
        "single_cell_downstream requires modules.single_cell=true.",
        call. = FALSE
      )
    }
    required_downstream_path <- as.character(
      config$single_cell_downstream$gse104782_cluster_path %||% ""
    )
    if (!nzchar(required_downstream_path)) {
      stop(
        "single_cell_downstream requires gse104782_cluster_path.",
        call. = FALSE
      )
    }
    numeric_parameters <- c(
      config$single_cell_downstream$hvg_n,
      config$single_cell_downstream$pca_n,
      config$single_cell_downstream$neighbors_k,
      config$single_cell_downstream$visualization_max_cells,
      config$single_cell_downstream$minimum_pseudobulk_cells,
      config$single_cell_downstream$minimum_group_replicates
    )
    if (anyNA(numeric_parameters) || any(numeric_parameters < 1L)) {
      stop(
        "single_cell_downstream numeric parameters must be positive.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(config$project$strict)) {
    module_fields <- list(
      tcga = c("expression_path", "clinical_path"),
      regulatory = c(
        "mirtarbase_path", "knocktf_gmt_path", "knocktf_table_path"
      ),
      drug = c("dgidb_path", "ctd_path")
    )
    for (module in names(module_fields)) {
      if (!isTRUE(config$modules[[module]])) next
      section <- config[[module]] %||% list()
      missing <- module_fields[[module]][
        !vapply(
          module_fields[[module]],
          function(field) {
            value <- section[[field]]
            !is.null(value) && length(value) == 1L &&
              !is.na(value) && nzchar(as.character(value))
          },
          logical(1)
        )
      ]
      if (length(missing) > 0L) {
        stop(
          module, " is enabled, but required fields are missing: ",
          paste(missing, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    }
    if (
      isTRUE(config$modules$enrichment) &&
        length(config$gene_sets %||% list()) == 0L
    ) {
      stop(
        "enrichment is enabled, but no gene_sets are defined.",
        call. = FALSE
      )
    }
  }

  config
}

configured_input_manifest <- function(config, active_stages = NULL) {
  stage_is_active <- function(stage) {
    is.null(active_stages) || stage %in% active_stages
  }
  rows <- list()
  add <- function(category, key, path, required, module) {
    if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
      return(invisible(NULL))
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      category = category,
      key = key,
      path = path,
      required = isTRUE(required),
      module = module,
      stringsAsFactors = FALSE
    )
  }

  for (name in names(config$datasets)) {
    dataset <- config$datasets[[name]]
    dataset_stage <- if (identical(dataset$role, "train")) {
      "bulk_training"
    } else {
      "bulk_validation"
    }
    if (!stage_is_active(dataset_stage)) next
    fields <- intersect(
      c("case_path", "control_path", "series_path", "annotation_path"),
      names(dataset)
    )
    for (field in fields) {
      add(
        "bulk_dataset",
        paste(name, field, sep = "."),
        dataset[[field]],
        dataset$required,
        "bulk"
      )
    }
  }

  if (isTRUE(config$modules$enrichment) && stage_is_active("enrichment")) {
    for (name in names(config$gene_sets %||% list())) {
      add(
        "gene_set",
        name,
        config$gene_sets[[name]],
        TRUE,
        "enrichment"
      )
    }
  }

  if (isTRUE(config$modules$immune) && stage_is_active("immune")) {
    add(
      "configuration",
      "immune_signatures",
      file.path(config$.project_root, "config", "immune_signatures.yml"),
      TRUE,
      "immune"
    )
  }

  if (isTRUE(config$modules$tcga) && stage_is_active("tcga")) {
    add("tcga", "expression", config$tcga$expression_path, TRUE, "tcga")
    add("tcga", "clinical", config$tcga$clinical_path, TRUE, "tcga")
  }

  if (isTRUE(config$modules$regulatory) && stage_is_active("regulatory")) {
    for (field in c(
      "mirtarbase_path", "knocktf_gmt_path", "knocktf_table_path"
    )) {
      add(
        "regulatory",
        field,
        config$regulatory[[field]],
        TRUE,
        "regulatory"
      )
    }
  }

  if (isTRUE(config$modules$drug) && stage_is_active("drug")) {
    add("drug", "dgidb", config$drug$dgidb_path, TRUE, "drug")
    add("drug", "ctd", config$drug$ctd_path, TRUE, "drug")
  }

  if (isTRUE(config$modules$single_cell) && stage_is_active("single_cell")) {
    for (dataset in config$single_cell$datasets %||% list()) {
      for (field in names(dataset)) {
        if (grepl("(_path|_archive)$", field)) {
          add(
            "single_cell",
            paste(dataset$id %||% "dataset", field, sep = "."),
            dataset[[field]],
            dataset$required %||% TRUE,
            "single_cell"
          )
        }
      }
    }
  }

  if (
    isTRUE(config$modules$single_cell_downstream) &&
      stage_is_active("single_cell_downstream")
  ) {
    add(
      "single_cell_reference",
      "GSE104782_published_clusters",
      config$single_cell_downstream$gse104782_cluster_path,
      TRUE,
      "single_cell_downstream"
    )
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}
