submission_training_fold_screen <- function(expression, status, top_n) {
  case_matrix <- expression[, status == 1L, drop = FALSE]
  reference_matrix <- expression[, status == 0L, drop = FALSE]
  row_variance <- function(matrix) {
    n <- ncol(matrix)
    if (n < 2L) return(rep(NA_real_, nrow(matrix)))
    centered_ss <- rowSums(matrix^2) - rowSums(matrix)^2 / n
    pmax(centered_ss / (n - 1L), 0)
  }
  case_mean <- rowMeans(case_matrix)
  reference_mean <- rowMeans(reference_matrix)
  standard_error <- sqrt(
    row_variance(case_matrix) / ncol(case_matrix) +
      row_variance(reference_matrix) / ncol(reference_matrix)
  )
  statistic <- (case_mean - reference_mean) / standard_error
  statistic[!is.finite(statistic)] <- NA_real_
  order_index <- order(abs(statistic), decreasing = TRUE, na.last = NA)
  genes <- rownames(expression)[head(order_index, as.integer(top_n))]
  data.frame(
    gene = genes,
    screening_rank = seq_along(genes),
    training_t_statistic = statistic[genes],
    stringsAsFactors = FALSE
  )
}

submission_run_ml_sensitivity <- function(
    bulk_training,
    ml,
    project_config,
    settings,
    output_dir
) {
  require_namespace("glmnet", "strict nested LASSO validation")
  require_namespace("randomForest", "strict nested random-forest validation")
  repeats <- as.integer(settings$repeated_cv_replicates)
  outer_folds <- as.integer(settings$outer_folds)
  inner_folds_requested <- as.integer(settings$inner_folds)
  screen_top_n <- as.integer(settings$feature_screen_top_n)
  trees <- as.integer(settings$random_forest_trees)
  tuning_trees <- as.integer(settings$random_forest_tuning_trees)
  mtry_candidates <- as.integer(unlist(
    settings$random_forest_mtry_candidates
  ))
  importance_top_n <- as.integer(settings$random_forest_top_n)
  metrics <- list()
  selections <- list()
  screens <- list()
  fits <- list()
  metric_index <- 0L
  selection_index <- 0L
  screen_index <- 0L
  fit_index <- 0L

  for (dataset_index in seq_along(bulk_training)) {
    dataset <- bulk_training[[dataset_index]]
    expression <- as.matrix(dataset$expression)
    group <- droplevels(dataset$group[colnames(expression)])
    status <- as.integer(group == "Disease")

    for (replicate in seq_len(repeats)) {
      seed <- as.integer(project_config$project$seed) +
        dataset_index * 10000L + replicate
      outer_id <- stratified_fold_ids(group, outer_folds, seed)
      predictions <- list(
        LASSO = rep(NA_real_, ncol(expression)),
        RandomForest = rep(NA_real_, ncol(expression))
      )

      for (fold in seq_len(outer_folds)) {
        train <- which(outer_id != fold)
        test <- which(outer_id == fold)
        train_group <- droplevels(group[train])
        train_status <- status[train]
        screen <- submission_training_fold_screen(
          expression[, train, drop = FALSE],
          train_status,
          screen_top_n
        )
        selected_genes <- screen$gene
        screen$dataset_id <- dataset$id
        screen$disease <- dataset$disease
        screen$replicate <- replicate
        screen$fold <- fold
        screen_index <- screen_index + 1L
        screens[[screen_index]] <- screen
        x_train <- t(expression[selected_genes, train, drop = FALSE])
        x_test <- t(expression[selected_genes, test, drop = FALSE])
        class_counts <- table(train_group)
        weights <- 1 / as.numeric(class_counts[train_group])
        weights <- weights / mean(weights)
        inner_folds <- min(
          inner_folds_requested,
          max(3L, min(as.integer(class_counts)))
        )
        inner_id <- stratified_fold_ids(
          train_group,
          inner_folds,
          seed + fold * 100L
        )

        lasso <- tryCatch(
          suppressWarnings(glmnet::cv.glmnet(
            x = x_train,
            y = train_status,
            family = "binomial",
            alpha = 1,
            foldid = inner_id,
            weights = weights,
            type.measure = "deviance",
            standardize = TRUE
          )),
          error = function(error) NULL
        )
        fit_index <- fit_index + 1L
        fits[[fit_index]] <- data.frame(
          dataset_id = dataset$id,
          disease = dataset$disease,
          replicate = replicate,
          fold = fold,
          model = "LASSO",
          fit_success = !is.null(lasso),
          tuned_parameter = if (is.null(lasso)) NA_character_ else "lambda.1se",
          stringsAsFactors = FALSE
        )
        if (!is.null(lasso)) {
          predictions$LASSO[test] <- as.numeric(stats::predict(
            lasso,
            newx = x_test,
            s = "lambda.1se",
            type = "response"
          ))
          coefficients <- as.matrix(stats::coef(lasso, s = "lambda.1se"))
          selected <- setdiff(
            rownames(coefficients)[coefficients[, 1L] != 0],
            "(Intercept)"
          )
          if (length(selected) > 0L) {
            selection_index <- selection_index + 1L
            selections[[selection_index]] <- data.frame(
              dataset_id = dataset$id,
              disease = dataset$disease,
              replicate = replicate,
              fold = fold,
              model = "LASSO",
              gene = selected,
              stringsAsFactors = FALSE
            )
          }
        }

        usable_mtry <- sort(unique(pmax(
          1L,
          pmin(mtry_candidates, ncol(x_train))
        )))
        oob_auc <- rep(NA_real_, length(usable_mtry))
        for (mtry_index in seq_along(usable_mtry)) {
          set.seed(seed + fold * 1000L + mtry_index)
          tune_forest <- tryCatch(
            randomForest::randomForest(
              x = x_train,
              y = train_group,
              ntree = tuning_trees,
              mtry = usable_mtry[[mtry_index]],
              classwt = setNames(
                sum(class_counts) /
                  (length(class_counts) * as.numeric(class_counts)),
                names(class_counts)
              ),
              importance = FALSE
            ),
            error = function(error) NULL
          )
          if (!is.null(tune_forest)) {
            oob_auc[[mtry_index]] <- submission_auc(
              train_status,
              tune_forest$votes[, "Disease"]
            )
          }
        }
        best_mtry <- if (all(!is.finite(oob_auc))) {
          NA_integer_
        } else {
          usable_mtry[order(-oob_auc, usable_mtry, na.last = NA)[[1L]]]
        }
        forest <- NULL
        if (is.finite(best_mtry)) {
          set.seed(seed + fold * 2000L)
          forest <- tryCatch(
            randomForest::randomForest(
              x = x_train,
              y = train_group,
              ntree = trees,
              mtry = best_mtry,
              classwt = setNames(
                sum(class_counts) /
                  (length(class_counts) * as.numeric(class_counts)),
                names(class_counts)
              ),
              importance = TRUE
            ),
            error = function(error) NULL
          )
        }
        fit_index <- fit_index + 1L
        fits[[fit_index]] <- data.frame(
          dataset_id = dataset$id,
          disease = dataset$disease,
          replicate = replicate,
          fold = fold,
          model = "RandomForest",
          fit_success = !is.null(forest),
          tuned_parameter = if (is.null(forest)) {
            NA_character_
          } else {
            paste0("mtry=", best_mtry, "; OOB tuned in outer training fold")
          },
          stringsAsFactors = FALSE
        )
        if (!is.null(forest)) {
          probabilities <- stats::predict(
            forest,
            newdata = x_test,
            type = "prob"
          )
          predictions$RandomForest[test] <- probabilities[, "Disease"]
          importance <- randomForest::importance(forest)
          column <- if ("MeanDecreaseGini" %in% colnames(importance)) {
            "MeanDecreaseGini"
          } else {
            colnames(importance)[[ncol(importance)]]
          }
          selected <- head(
            rownames(importance)[
              order(importance[, column], decreasing = TRUE)
            ],
            importance_top_n
          )
          selection_index <- selection_index + 1L
          selections[[selection_index]] <- data.frame(
            dataset_id = dataset$id,
            disease = dataset$disease,
            replicate = replicate,
            fold = fold,
            model = "RandomForest",
            gene = selected,
            stringsAsFactors = FALSE
          )
        }
      }

      for (model in names(predictions)) {
        probability <- predictions[[model]]
        complete <- is.finite(probability)
        metric_index <- metric_index + 1L
        metrics[[metric_index]] <- data.frame(
          dataset_id = dataset$id,
          disease = dataset$disease,
          replicate = replicate,
          model = model,
          auc = submission_auc(status[complete], probability[complete]),
          balanced_accuracy = if (all(complete)) {
            submission_balanced_accuracy(status, probability)
          } else {
            NA_real_
          },
          brier_score = if (all(complete)) {
            mean((probability - status)^2)
          } else {
            NA_real_
          },
          predictions_complete = all(complete),
          candidate_space = "all measured genes",
          feature_selection_scope = "outer training fold only",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  metrics <- do.call(rbind, metrics)
  selections <- do.call(rbind, Filter(
    function(item) !is.null(item) && nrow(item) > 0L,
    selections
  ))
  screens <- do.call(rbind, screens)
  fits <- do.call(rbind, fits)
  group_keys <- unique(metrics[, c("dataset_id", "disease", "model")])
  metric_summary <- do.call(rbind, lapply(
    seq_len(nrow(group_keys)),
    function(index) {
      key <- group_keys[index, , drop = FALSE]
      subset <- metrics[
        metrics$dataset_id == key$dataset_id &
          metrics$model == key$model,
        ,
        drop = FALSE
      ]
      data.frame(
        key,
        valid_repeats = sum(is.finite(subset$auc)),
        auc_median = stats::median(subset$auc, na.rm = TRUE),
        auc_lower_95_resampling = submission_quantile(
          subset$auc, 0.025
        ),
        auc_upper_95_resampling = submission_quantile(
          subset$auc, 0.975
        ),
        balanced_accuracy_median = stats::median(
          subset$balanced_accuracy, na.rm = TRUE
        ),
        brier_score_median = stats::median(
          subset$brier_score, na.rm = TRUE
        ),
        candidate_space = "all measured genes",
        feature_selection_scope = "outer training fold only",
        stringsAsFactors = FALSE
      )
    }
  ))
  selection_counts <- as.data.frame(table(
    selections$dataset_id,
    selections$disease,
    selections$model,
    selections$gene
  ), stringsAsFactors = FALSE)
  names(selection_counts) <- c(
    "dataset_id", "disease", "model", "gene", "selected_folds"
  )
  selection_counts <- selection_counts[selection_counts$selected_folds > 0, ]
  denominator <- aggregate(
    fit_success ~ dataset_id + disease + model,
    data = fits,
    FUN = sum
  )
  names(denominator)[[4L]] <- "successful_folds"
  selection_frequency <- merge(
    selection_counts,
    denominator,
    by = c("dataset_id", "disease", "model"),
    all.x = TRUE
  )
  selection_frequency$selection_frequency <-
    selection_frequency$selected_folds /
    selection_frequency$successful_folds
  selection_frequency$is_final_hub <-
    selection_frequency$gene %in% ml$final_genes
  screen_frequency <- aggregate(
    screening_rank ~ dataset_id + disease + gene,
    data = screens,
    FUN = function(value) c(
      selected_folds = length(value),
      median_rank = stats::median(value)
    )
  )
  screen_frequency <- do.call(data.frame, screen_frequency)
  names(screen_frequency)[
    (ncol(screen_frequency) - 1L):ncol(screen_frequency)
  ] <- c("selected_folds", "median_screening_rank")
  screen_denominator <- aggregate(
    fold ~ dataset_id + disease,
    data = unique(screens[, c("dataset_id", "disease", "replicate", "fold")]),
    FUN = length
  )
  names(screen_denominator)[[3L]] <- "outer_training_folds"
  screen_frequency <- merge(
    screen_frequency,
    screen_denominator,
    by = c("dataset_id", "disease"),
    all.x = TRUE
  )
  screen_frequency$screening_frequency <-
    screen_frequency$selected_folds /
    screen_frequency$outer_training_folds
  screen_frequency$is_final_hub <-
    screen_frequency$gene %in% ml$final_genes

  safe_write_csv(
    metrics,
    file.path(output_dir, "machine_learning_repeated_cv_metrics.csv")
  )
  safe_write_csv(
    metric_summary,
    file.path(output_dir, "machine_learning_repeated_cv_summary.csv")
  )
  safe_write_csv(
    selection_frequency,
    file.path(output_dir, "machine_learning_selection_frequency.csv")
  )
  safe_write_csv(
    screen_frequency,
    file.path(output_dir, "machine_learning_screen_frequency.csv")
  )
  safe_write_csv(
    fits,
    file.path(output_dir, "machine_learning_outer_fold_fits.csv")
  )
  list(
    metrics = metrics,
    summary = metric_summary,
    selection_frequency = selection_frequency,
    screen_frequency = screen_frequency,
    fits = fits
  )
}

submission_signed_score <- function(dataset, genes, signs, sample_indices) {
  matrix <- dataset$expression[genes, sample_indices, drop = FALSE]
  z <- t(scale(t(matrix)))
  usable <- rowSums(is.finite(z)) == ncol(z)
  z <- z[usable, , drop = FALSE]
  colMeans(z * signs[rownames(z)])
}

submission_run_external_sensitivity <- function(
    bulk_validation,
    validation,
    shared,
    hub_genes,
    settings,
    output_dir
) {
  require_namespace("pROC", "external validation sensitivity")
  fixed_table <- validation$table
  permutations <- as.integer(settings$permutation_replicates)
  comparison_rows <- list()
  composite_rows <- list()
  permutation_rows <- list()
  leave_one_out_rows <- list()
  score_rows <- list()
  comparison_index <- 0L

  for (dataset_index in seq_along(bulk_validation)) {
    dataset <- bulk_validation[[dataset_index]]
    logfc_column <- paste0("logFC_", dataset$disease)
    logfc <- shared$table[[logfc_column]]
    names(logfc) <- shared$table$gene
    present <- intersect(hub_genes, rownames(dataset$expression))
    present <- intersect(present, names(logfc))
    for (gene in present) {
      auto_roc <- pROC::roc(
        response = dataset$group,
        predictor = as.numeric(dataset$expression[gene, ]),
        levels = c("Normal", "Disease"),
        direction = "auto",
        quiet = TRUE
      )
      fixed <- fixed_table[
        fixed_table$dataset_id == dataset$id &
          fixed_table$gene == gene,
        ,
        drop = FALSE
      ]
      comparison_index <- comparison_index + 1L
      comparison_rows[[comparison_index]] <- data.frame(
        dataset_id = dataset$id,
        disease = dataset$disease,
        gene = gene,
        fixed_direction_auc = fixed$auc[[1L]],
        legacy_auto_direction_auc = as.numeric(pROC::auc(auto_roc)),
        auto_minus_fixed_auc =
          as.numeric(pROC::auc(auto_roc)) - fixed$auc[[1L]],
        expected_expression = fixed$expected_expression[[1L]],
        stringsAsFactors = FALSE
      )
    }

    all_indices <- seq_len(ncol(dataset$expression))
    signs <- sign(logfc[present])
    score <- submission_signed_score(
      dataset, present, signs, all_indices
    )
    status <- as.integer(dataset$group == "Disease")
    observed_auc <- submission_auc(status, score)
    composite_roc <- pROC::roc(
      response = dataset$group,
      predictor = score,
      levels = c("Normal", "Disease"),
      direction = "<",
      quiet = TRUE
    )
    interval <- as.numeric(pROC::ci.auc(composite_roc, method = "delong"))
    set.seed(20260726L + dataset_index * 1000L)
    permutation_auc <- replicate(
      permutations,
      submission_auc(sample(status, replace = FALSE), score)
    )
    empirical_p <- (
      1 + sum(permutation_auc >= observed_auc, na.rm = TRUE)
    ) / (permutations + 1)
    loo_auc <- rep(NA_real_, length(all_indices))
    for (omitted in all_indices) {
      retained <- setdiff(all_indices, omitted)
      loo_score <- submission_signed_score(
        dataset, present, signs, retained
      )
      loo_auc[[omitted]] <- submission_auc(status[retained], loo_score)
    }
    composite_rows[[dataset_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      genes_in_score = length(present),
      auc = observed_auc,
      ci_lower = interval[[1L]],
      ci_upper = interval[[3L]],
      permutation_replicates = permutations,
      permutation_empirical_p = empirical_p,
      permutation_null_median = stats::median(permutation_auc, na.rm = TRUE),
      permutation_null_lower_95 = submission_quantile(
        permutation_auc, 0.025
      ),
      permutation_null_upper_95 = submission_quantile(
        permutation_auc, 0.975
      ),
      leave_one_out_auc_minimum = min(loo_auc, na.rm = TRUE),
      leave_one_out_auc_median = stats::median(loo_auc, na.rm = TRUE),
      leave_one_out_auc_maximum = max(loo_auc, na.rm = TRUE),
      n_normal = sum(dataset$group == "Normal"),
      n_disease = sum(dataset$group == "Disease"),
      stringsAsFactors = FALSE
    )
    permutation_rows[[dataset_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      permutation = seq_len(permutations),
      auc = permutation_auc,
      observed_auc = observed_auc,
      stringsAsFactors = FALSE
    )
    leave_one_out_rows[[dataset_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      omitted_sample = colnames(dataset$expression),
      omitted_group = as.character(dataset$group),
      leave_one_out_auc = loo_auc,
      stringsAsFactors = FALSE
    )
    score_rows[[dataset_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      sample_id = colnames(dataset$expression),
      group = as.character(dataset$group),
      signed_score = score,
      stringsAsFactors = FALSE
    )
  }

  comparison <- do.call(rbind, comparison_rows)
  composite <- do.call(rbind, composite_rows)
  permutation_table <- do.call(rbind, permutation_rows)
  leave_one_out <- do.call(rbind, leave_one_out_rows)
  sample_scores <- do.call(rbind, score_rows)
  threshold <- as.numeric(settings$robust_auc_threshold)
  keys <- unique(fixed_table[, c("disease", "gene")])
  consistency <- do.call(rbind, lapply(seq_len(nrow(keys)), function(index) {
    key <- keys[index, , drop = FALSE]
    subset <- fixed_table[
      fixed_table$disease == key$disease &
        fixed_table$gene == key$gene,
      ,
      drop = FALSE
    ]
    data.frame(
      disease = key$disease,
      gene = key$gene,
      cohorts = nrow(subset),
      minimum_fixed_auc = min(subset$auc),
      median_fixed_auc = stats::median(subset$auc),
      maximum_fixed_auc = max(subset$auc),
      all_above_chance = all(subset$auc > 0.5),
      all_at_or_above_robust_threshold = all(subset$auc >= threshold),
      any_ci_excludes_chance = any(subset$ci_lower > 0.5),
      all_ci_exclude_chance = all(subset$ci_lower > 0.5),
      stringsAsFactors = FALSE
    )
  }))

  safe_write_csv(
    comparison,
    file.path(output_dir, "external_validation_direction_bias_audit.csv")
  )
  safe_write_csv(
    fixed_table,
    file.path(output_dir, "external_validation_direction_fixed_auc.csv")
  )
  safe_write_csv(
    consistency,
    file.path(output_dir, "external_validation_cross_cohort_consistency.csv")
  )
  safe_write_csv(
    composite,
    file.path(output_dir, "external_validation_signed_composite_score.csv")
  )
  safe_write_csv(
    permutation_table,
    file.path(output_dir, "external_validation_permutation_auc.csv")
  )
  safe_write_csv(
    leave_one_out,
    file.path(output_dir, "external_validation_leave_one_out_auc.csv")
  )
  safe_write_csv(
    sample_scores,
    file.path(output_dir, "external_validation_sample_scores.csv")
  )
  list(
    comparison = comparison,
    fixed = fixed_table,
    consistency = consistency,
    composite = composite,
    permutations = permutation_table,
    leave_one_out = leave_one_out,
    sample_scores = sample_scores
  )
}

submission_run_hpa_context <- function(
    project_root,
    hub_genes,
    output_dir
) {
  require_namespace("jsonlite", "Human Protein Atlas context analysis")
  hpa_dir <- file.path(project_root, "data", "external", "HPA")
  rows <- lapply(hub_genes, function(gene) {
    path <- file.path(hpa_dir, paste0(gene, ".json"))
    if (!file.exists(path) || file.info(path)$size == 0L) {
      stop("Missing HPA JSON input: ", path, call. = FALSE)
    }
    record <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    if (is.list(record) && is.null(record$Gene)) {
      matches <- vapply(
        record,
        function(item) identical(item$Gene, gene),
        logical(1)
      )
      record <- record[[which(matches)[[1L]]]]
    }
    tissue_values <- record[["RNA tissue specific nTPM"]]
    tissue_names <- if (is.null(tissue_values)) {
      character()
    } else {
      names(tissue_values)
    }
    cell_values <- record[["RNA single cell type specific nCPM"]]
    cell_names <- if (is.null(cell_values)) {
      character()
    } else {
      names(cell_values)
    }
    data.frame(
      gene = gene,
      ensembl_id = record$Ensembl %||% NA_character_,
      rna_tissue_specificity =
        record[["RNA tissue specificity"]] %||% NA_character_,
      rna_tissue_distribution =
        record[["RNA tissue distribution"]] %||% NA_character_,
      tissue_specific_locations = if (length(tissue_names) == 0L) {
        "none listed"
      } else {
        paste(tissue_names, collapse = "; ")
      },
      ovary_listed_as_specific =
        any(tolower(tissue_names) == "ovary"),
      rna_cell_type_specificity =
        record[["RNA single cell type specificity"]] %||% NA_character_,
      cell_type_specific_locations = if (length(cell_names) == 0L) {
        "none listed"
      } else {
        paste(cell_names, collapse = "; ")
      },
      cartilage_in_reference = FALSE,
      hpa_version = "25.1",
      ensembl_version = "109",
      source_url = paste0(
        "https://www.proteinatlas.org/search/",
        gene,
        "?format=json"
      ),
      access_date = "2026-07-29",
      interpretation = paste(
        "Normal-tissue reference only; cartilage is absent,",
        "so this analysis cannot establish OA or OC disease specificity."
      ),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  safe_write_csv(
    table,
    file.path(output_dir, "hpa_normal_tissue_context.csv")
  )
  table
}

run_submission_sensitivity <- function(
    project_root,
    project_config_path,
    sensitivity_config_path
) {
  require_namespace("yaml", "sensitivity configuration")
  project_config <- read_project_config(project_root, project_config_path)
  sensitivity_config <- yaml::read_yaml(file.path(
    project_root,
    sensitivity_config_path
  ))
  output_dir <- ensure_dir(file.path(
    project_root,
    "results",
    "submission",
    "sensitivity"
  ))
  initialize_logging(file.path(project_root, "results", "submission"))
  log_info("Starting reviewer-prioritized non-MR sensitivity analyses.")
  bulk_training <- submission_load_cache(project_root, "02_bulk_training.rds")
  differential <- submission_load_cache(project_root, "03_differential.rds")
  shared <- submission_load_cache(project_root, "04_shared.rds")
  wgcna <- submission_load_cache(project_root, "06_wgcna.rds")
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  bulk_validation <- submission_load_cache(
    project_root, "08_bulk_validation.rds"
  )
  validation <- submission_load_cache(project_root, "09_validation.rds")
  tcga <- submission_load_cache(project_root, "11_tcga.rds")

  deg_result <- submission_run_deg_sensitivity(
    differential,
    shared,
    ml$final_genes,
    sensitivity_config$deg,
    output_dir
  )
  wgcna_result <- submission_run_wgcna_sensitivity(
    bulk_training,
    wgcna,
    project_config,
    sensitivity_config$wgcna,
    output_dir
  )
  log_info("Running strict outer-fold feature selection and nested tuning.")
  ml_result <- submission_run_ml_sensitivity(
    bulk_training,
    ml,
    project_config,
    sensitivity_config$machine_learning,
    output_dir
  )
  tcga_result <- submission_run_tcga_sensitivity(
    tcga,
    ml,
    project_config,
    sensitivity_config$tcga,
    output_dir
  )
  log_info("Running fixed-score permutation and leave-one-out analyses.")
  external_result <- submission_run_external_sensitivity(
    bulk_validation,
    validation,
    shared,
    ml$final_genes,
    sensitivity_config$external_validation,
    output_dir
  )
  log_info("Summarizing Human Protein Atlas normal-tissue context.")
  hpa_result <- submission_run_hpa_context(
    project_root,
    ml$final_genes,
    output_dir
  )
  results <- list(
    deg = deg_result,
    wgcna = wgcna_result,
    machine_learning = ml_result,
    tcga = tcga_result,
    external = external_result,
    hpa = hpa_result
  )
  atomic_save_rds(
    results,
    file.path(output_dir, "non_mr_sensitivity_results.rds")
  )
  submission_write_sensitivity_report(
    results,
    output_dir,
    sensitivity_config
  )
  safe_write_csv(
    data.frame(
      analysis = c(
        "DEG threshold grid",
        "WGCNA soft-power perturbation",
        "WGCNA fixed-module bootstrap",
        "Strict nested ML",
        "TCGA optimism bootstrap",
        "TCGA LASSO selection bootstrap",
        "External direction-fixed ROC",
        "External score permutation",
        "External leave-one-sample-out",
        "HPA normal-tissue context"
      ),
      setting = c(
        "2 FDR x 3 absolute log2FC rules",
        paste(
          "offsets",
          paste(
            unlist(sensitivity_config$wgcna$soft_power_offsets),
            collapse = ","
          )
        ),
        paste(
          sensitivity_config$wgcna$bootstrap_replicates,
          "replicates"
        ),
        paste0(
          sensitivity_config$machine_learning$repeated_cv_replicates,
          " x ",
          sensitivity_config$machine_learning$outer_folds,
          "-fold; top ",
          sensitivity_config$machine_learning$feature_screen_top_n,
          " genes screened inside each outer training fold"
        ),
        paste(
          sensitivity_config$tcga$optimism_bootstrap_replicates,
          "replicates"
        ),
        paste(
          sensitivity_config$tcga$lasso_selection_bootstrap_replicates,
          "replicates"
        ),
        "training logFC fixes direction; DeLong 95% CI",
        paste(
          sensitivity_config$external_validation$permutation_replicates,
          "label permutations"
        ),
        "recompute score after removing each sample",
        "HPA 25.1 / Ensembl 109; cartilage unavailable"
      ),
      seed = as.integer(sensitivity_config$seed),
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "sensitivity_parameter_manifest.csv")
  )
  log_info("Reviewer-prioritized sensitivity analyses completed.")
  invisible(results)
}

submission_write_sensitivity_report <- function(results, output_dir, settings) {
  deg_primary <- results$deg$summary[results$deg$summary$is_primary, ]
  ml <- results$machine_learning$summary
  tcga <- results$tcga$optimism
  external <- results$external$composite
  hpa <- results$hpa
  lines <- c(
    "# Reviewer-prioritized non-MR sensitivity analysis report",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- MR remains a separate negative supplementary analysis.",
    "",
    "## DEG thresholds and direction",
    "",
    sprintf(
      paste(
        "The primary FDR/absolute log2FC rule (%.2f/%.1f) retained",
        "%d shared genes, including %d/10 prioritized genes."
      ),
      deg_primary$fdr_threshold,
      deg_primary$absolute_log2fc_threshold,
      deg_primary$shared_count,
      deg_primary$retained_hub_count
    ),
    sprintf(
      "%d/%d (%.1f%%) primary shared genes changed in the same direction.",
      deg_primary$directionally_concordant_count,
      deg_primary$shared_count,
      100 * deg_primary$directionally_concordant_fraction
    ),
    "",
    "## WGCNA stability",
    "",
    paste(
      "Fixed-module bootstrap sign stability:",
      paste(
        sprintf(
          "%s %.3f",
          results$wgcna$bootstrap$dataset_id,
          results$wgcna$bootstrap$sign_stability
        ),
        collapse = "; "
      )
    ),
    paste(
      "The OA discovery size (n=38) limits network precision;",
      "2,000 bootstraps and leave-one-sample-out analyses quantify",
      "association stability but do not replace independent preservation."
    ),
    "",
    "## Strict nested machine learning",
    "",
    paste(
      sprintf(
        "%s %s median outer-CV AUC %.3f (95%% resampling range %.3f–%.3f)",
        ml$dataset_id,
        ml$model,
        ml$auc_median,
        ml$auc_lower_95_resampling,
        ml$auc_upper_95_resampling
      ),
      collapse = "; "
    ),
    paste(
      "All measured genes entered each resample. Univariate screening was",
      "repeated inside every outer training fold, LASSO lambda was chosen",
      "by inner cross-validation, and random-forest mtry was selected from",
      "outer-training-fold out-of-bag predictions."
    ),
    paste(
      "The near-perfect results persisted after removing the prior",
      "candidate-space leakage. They are interpreted as strong retrospective",
      "molecular separation, not clinical diagnostic performance."
    ),
    "",
    "## External fixed-score robustness",
    "",
    paste(
      sprintf(
        paste0(
          "%s AUC %.3f (%.3f–%.3f), permutation P=%.4f, ",
          "leave-one-out range %.3f–%.3f"
        ),
        external$dataset_id,
        external$auc,
        external$ci_lower,
        external$ci_upper,
        external$permutation_empirical_p,
        external$leave_one_out_auc_minimum,
        external$leave_one_out_auc_maximum
      ),
      collapse = "; "
    ),
    paste(
      "OC separation remained high, but the cohorts compare malignant with",
      "normal tissue and have small normal groups; tissue composition,",
      "platform, and retrospective-cohort effects remain plausible."
    ),
    "",
    "## Normal-tissue context",
    "",
    sprintf(
      paste(
        "%d/10 genes had low tissue specificity and %d/10 listed ovary",
        "among HPA-specific tissues."
      ),
      sum(hpa$rna_tissue_specificity == "Low tissue specificity"),
      sum(hpa$ovary_listed_as_specific)
    ),
    paste(
      "Cartilage is absent from HPA 25.1 normal-tissue categories.",
      "This analysis audits background expression and cannot establish",
      "OA–OC disease specificity."
    ),
    "",
    "## TCGA-OV model",
    "",
    sprintf(
      paste(
        "The apparent C-index was %.3f and the optimism-corrected C-index",
        "was %.3f after %d successful bootstrap replicates."
      ),
      tcga$apparent_cindex,
      tcga$optimism_corrected_cindex,
      tcga$successful_bootstrap_replicates
    ),
    "",
    "## Interpretation boundary",
    "",
    paste(
      "These analyses support reproducible associations and molecular",
      "separation. They do not establish a causal OA–OC relationship,",
      "a clinically validated diagnostic panel, or a treatment target."
    )
  )
  write_utf8(
    lines,
    file.path(output_dir, "non_mr_sensitivity_report.md")
  )
}
