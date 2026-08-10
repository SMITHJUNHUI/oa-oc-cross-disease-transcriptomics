submission_load_cache <- function(project_root, file) {
  path <- file.path(project_root, "results", "cache", file)
  if (!file.exists(path)) {
    stop("Required cache is missing: ", path, call. = FALSE)
  }
  object <- readRDS(path)
  object$value %||% object
}

submission_auc <- function(status, score) {
  keep <- is.finite(score) & !is.na(status)
  status <- as.integer(status[keep])
  score <- as.numeric(score[keep])
  positive <- which(status == 1L)
  negative <- which(status == 0L)
  if (length(positive) == 0L || length(negative) == 0L) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  statistic <- sum(ranks[positive]) -
    length(positive) * (length(positive) + 1) / 2
  statistic / (length(positive) * length(negative))
}

submission_balanced_accuracy <- function(status, probability) {
  predicted <- as.integer(probability >= 0.5)
  sensitivity <- mean(predicted[status == 1L] == 1L)
  specificity <- mean(predicted[status == 0L] == 0L)
  mean(c(sensitivity, specificity))
}

submission_stratified_sample <- function(status, replace = TRUE) {
  unlist(lapply(
    sort(unique(status)),
    function(value) {
      indices <- which(status == value)
      sample(indices, length(indices), replace = replace)
    }
  ), use.names = FALSE)
}

submission_quantile <- function(values, probability) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  as.numeric(stats::quantile(values, probability, names = FALSE, type = 8))
}

submission_run_deg_sensitivity <- function(
    differential,
    primary_shared,
    hub_genes,
    settings,
    output_dir
) {
  fdr_values <- as.numeric(unlist(settings$fdr_thresholds))
  fc_values <- as.numeric(unlist(settings$absolute_log2fc_thresholds))
  oa <- differential$oa_train$table
  oc <- differential$oc_train$table
  primary <- primary_shared$genes
  rows <- list()
  membership <- list()
  index <- 0L

  for (fdr in fdr_values) {
    for (fc in fc_values) {
      oa_genes <- oa$gene[oa$adj.P.Val < fdr & abs(oa$logFC) >= fc]
      oc_genes <- oc$gene[oc$adj.P.Val < fdr & abs(oc$logFC) >= fc]
      shared <- intersect(oa_genes, oc_genes)
      oa_logfc <- setNames(oa$logFC, oa$gene)
      oc_logfc <- setNames(oc$logFC, oc$gene)
      concordant <- shared[
        sign(oa_logfc[shared]) == sign(oc_logfc[shared])
      ]
      retained_hubs <- intersect(hub_genes, shared)
      union_primary <- union(primary, shared)
      jaccard <- if (length(union_primary) > 0L) {
        length(intersect(primary, shared)) / length(union_primary)
      } else {
        NA_real_
      }
      index <- index + 1L
      rows[[index]] <- data.frame(
        fdr_threshold = fdr,
        absolute_log2fc_threshold = fc,
        oa_deg_count = length(oa_genes),
        oc_deg_count = length(oc_genes),
        shared_count = length(shared),
        directionally_concordant_count = length(concordant),
        directionally_concordant_fraction = ifelse(
          length(shared) > 0L,
          length(concordant) / length(shared),
          NA_real_
        ),
        primary_set_jaccard = jaccard,
        retained_hub_count = length(retained_hubs),
        retained_hubs = paste(retained_hubs, collapse = ";"),
        is_primary = fdr == as.numeric(settings$primary_fdr) &&
          fc == as.numeric(settings$primary_absolute_log2fc),
        stringsAsFactors = FALSE
      )
      if (length(shared) > 0L) {
        membership[[index]] <- data.frame(
          gene = shared,
          fdr_threshold = fdr,
          absolute_log2fc_threshold = fc,
          logFC_OA = unname(oa_logfc[shared]),
          logFC_OC = unname(oc_logfc[shared]),
          directionally_concordant =
            sign(oa_logfc[shared]) == sign(oc_logfc[shared]),
          is_hub = shared %in% hub_genes,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary <- do.call(rbind, rows)
  membership <- do.call(rbind, Filter(Negate(is.null), membership))
  safe_write_csv(
    summary,
    file.path(output_dir, "deg_threshold_sensitivity_summary.csv")
  )
  safe_write_csv(
    membership,
    file.path(output_dir, "deg_threshold_sensitivity_membership.csv")
  )
  list(summary = summary, membership = membership)
}

submission_module_scores <- function(dataset, genes, expected_sign) {
  genes <- intersect(genes, rownames(dataset$expression))
  matrix <- t(dataset$expression[genes, , drop = FALSE])
  variable <- apply(matrix, 2L, stats::sd, na.rm = TRUE)
  matrix <- matrix[, is.finite(variable) & variable > 0, drop = FALSE]
  if (ncol(matrix) < 2L) {
    stop(dataset$id, " has too few variable module genes.", call. = FALSE)
  }
  scores <- stats::prcomp(matrix, center = TRUE, scale. = TRUE)$x[, 1L]
  trait <- as.integer(dataset$group[rownames(matrix)] == "Disease")
  observed <- stats::cor(scores, trait)
  if (sign(observed) != sign(expected_sign)) {
    scores <- -scores
  }
  list(scores = scores, trait = trait, genes = colnames(matrix))
}

submission_run_wgcna_sensitivity <- function(
    bulk_training,
    wgcna,
    project_config,
    settings,
    output_dir
) {
  perturbation_rows <- list()
  membership_rows <- list()
  bootstrap_rows <- list()
  loso_rows <- list()
  offsets <- as.numeric(unlist(settings$soft_power_offsets))
  bootstrap_replicates <- as.integer(settings$bootstrap_replicates)
  row_index <- 0L
  bootstrap_index <- 0L
  loso_index <- 0L

  for (dataset_key in names(wgcna)) {
    primary <- wgcna[[dataset_key]]
    dataset <- bulk_training[[dataset_key]]
    primary_correlation <- primary$module_table$correlation[
      match(primary$best_module, primary$module_table$module)
    ]
    powers <- sort(unique(pmax(1, primary$power + offsets)))

    for (power in powers) {
      if (power == primary$power) {
        result <- primary
      } else {
        local_settings <- project_config$wgcna
        local_settings$soft_power <- power
        result <- run_wgcna_dataset(
          dataset,
          settings = local_settings,
          output_dir = file.path(output_dir, "wgcna_recompute")
        )
      }
      overlap <- intersect(primary$selected_genes, result$selected_genes)
      union_genes <- union(primary$selected_genes, result$selected_genes)
      result_correlation <- result$module_table$correlation[
        match(result$best_module, result$module_table$module)
      ]
      row_index <- row_index + 1L
      perturbation_rows[[row_index]] <- data.frame(
        dataset_id = dataset$id,
        disease = dataset$disease,
        soft_power = power,
        primary_soft_power = primary$power,
        best_module = result$best_module,
        best_module_gene_count = length(result$selected_genes),
        module_trait_correlation = result_correlation,
        absolute_module_trait_correlation = abs(result_correlation),
        primary_gene_retention = length(overlap) /
          length(primary$selected_genes),
        jaccard_with_primary = length(overlap) / length(union_genes),
        is_primary = power == primary$power,
        stringsAsFactors = FALSE
      )
      membership_rows[[row_index]] <- data.frame(
        dataset_id = dataset$id,
        disease = dataset$disease,
        soft_power = power,
        gene = result$selected_genes,
        in_primary_module = result$selected_genes %in% primary$selected_genes,
        stringsAsFactors = FALSE
      )
    }

    scores <- submission_module_scores(
      dataset,
      primary$selected_genes,
      primary_correlation
    )
    set.seed(as.integer(project_config$project$seed) +
      match(dataset_key, names(wgcna)) * 1000L)
    correlations <- rep(NA_real_, bootstrap_replicates)
    for (replicate in seq_len(bootstrap_replicates)) {
      indices <- submission_stratified_sample(scores$trait, replace = TRUE)
      correlations[[replicate]] <- stats::cor(
        scores$scores[indices],
        scores$trait[indices]
      )
    }
    bootstrap_index <- bootstrap_index + 1L
    bootstrap_rows[[bootstrap_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      module = primary$best_module,
      module_gene_count = length(scores$genes),
      observed_correlation = stats::cor(scores$scores, scores$trait),
      bootstrap_median = stats::median(correlations, na.rm = TRUE),
      bootstrap_lower_95 = submission_quantile(correlations, 0.025),
      bootstrap_upper_95 = submission_quantile(correlations, 0.975),
      sign_stability = mean(
        sign(correlations) == sign(primary_correlation),
        na.rm = TRUE
      ),
      bootstrap_replicates = sum(is.finite(correlations)),
      stringsAsFactors = FALSE
    )

    for (left_out in seq_along(scores$scores)) {
      keep <- setdiff(seq_along(scores$scores), left_out)
      loso_index <- loso_index + 1L
      loso_rows[[loso_index]] <- data.frame(
        dataset_id = dataset$id,
        disease = dataset$disease,
        left_out_sample = names(scores$scores)[[left_out]],
        correlation = stats::cor(
          scores$scores[keep],
          scores$trait[keep]
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  perturbations <- do.call(rbind, perturbation_rows)
  memberships <- do.call(rbind, membership_rows)
  bootstrap <- do.call(rbind, bootstrap_rows)
  loso <- do.call(rbind, loso_rows)
  safe_write_csv(
    perturbations,
    file.path(output_dir, "wgcna_soft_power_perturbation.csv")
  )
  safe_write_csv(
    memberships,
    file.path(output_dir, "wgcna_soft_power_module_membership.csv")
  )
  safe_write_csv(
    bootstrap,
    file.path(output_dir, "wgcna_module_trait_bootstrap.csv")
  )
  safe_write_csv(
    loso,
    file.path(output_dir, "wgcna_module_trait_leave_one_out.csv")
  )
  list(
    perturbations = perturbations,
    memberships = memberships,
    bootstrap = bootstrap,
    loso = loso
  )
}

submission_run_ml_sensitivity <- function(
    bulk_training,
    ml,
    project_config,
    settings,
    output_dir
) {
  require_namespace("glmnet", "repeated LASSO validation")
  require_namespace("randomForest", "repeated random-forest validation")
  repeats <- as.integer(settings$repeated_cv_replicates)
  outer_folds <- as.integer(settings$outer_folds)
  inner_folds_requested <- as.integer(settings$inner_folds)
  trees <- as.integer(settings$random_forest_trees)
  top_n <- as.integer(settings$random_forest_top_n)
  metrics <- list()
  selections <- list()
  metric_index <- 0L
  selection_index <- 0L

  for (dataset_index in seq_along(bulk_training)) {
    dataset <- bulk_training[[dataset_index]]
    genes <- intersect(ml$candidates, rownames(dataset$expression))
    expression <- dataset$expression[genes, , drop = FALSE]
    variance <- apply(expression, 1L, stats::var, na.rm = TRUE)
    expression <- expression[
      is.finite(variance) & variance > 0,
      ,
      drop = FALSE
    ]
    x <- t(expression)
    group <- droplevels(dataset$group[rownames(x)])
    status <- as.integer(group == "Disease")

    for (replicate in seq_len(repeats)) {
      seed <- as.integer(project_config$project$seed) +
        dataset_index * 10000L + replicate
      outer_id <- stratified_fold_ids(group, outer_folds, seed)
      predictions <- list(
        LASSO = rep(NA_real_, nrow(x)),
        RandomForest = rep(NA_real_, nrow(x))
      )

      for (fold in seq_len(outer_folds)) {
        train <- which(outer_id != fold)
        test <- which(outer_id == fold)
        train_group <- droplevels(group[train])
        train_status <- status[train]
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
            x = x[train, , drop = FALSE],
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
        if (!is.null(lasso)) {
          predictions$LASSO[test] <- as.numeric(stats::predict(
            lasso,
            newx = x[test, , drop = FALSE],
            s = "lambda.1se",
            type = "response"
          ))
          coefficients <- as.matrix(stats::coef(lasso, s = "lambda.1se"))
          selected <- setdiff(
            rownames(coefficients)[coefficients[, 1L] != 0],
            "(Intercept)"
          )
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

        class_weights <- sum(class_counts) /
          (length(class_counts) * as.numeric(class_counts))
        names(class_weights) <- names(class_counts)
        set.seed(seed + fold * 1000L)
        forest <- tryCatch(
          randomForest::randomForest(
            x = x[train, , drop = FALSE],
            y = train_group,
            ntree = trees,
            classwt = class_weights,
            importance = TRUE
          ),
          error = function(error) NULL
        )
        if (!is.null(forest)) {
          probabilities <- stats::predict(
            forest,
            newdata = x[test, , drop = FALSE],
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
            top_n
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
  group_keys <- unique(metrics[, c("dataset_id", "disease", "model")])
  summary_rows <- lapply(seq_len(nrow(group_keys)), function(index) {
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
      auc_lower_95_resampling = submission_quantile(subset$auc, 0.025),
      auc_upper_95_resampling = submission_quantile(subset$auc, 0.975),
      balanced_accuracy_median =
        stats::median(subset$balanced_accuracy, na.rm = TRUE),
      brier_score_median = stats::median(subset$brier_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  metric_summary <- do.call(rbind, summary_rows)

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
    fold ~ dataset_id + disease + model,
    data = unique(selections[, c(
      "dataset_id", "disease", "model", "replicate", "fold"
    )]),
    FUN = length
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
  list(
    metrics = metrics,
    summary = metric_summary,
    selection_frequency = selection_frequency
  )
}

submission_cindex <- function(time, status, score) {
  data <- data.frame(time = time, status = status, score = score)
  result <- survival::concordance(
    survival::Surv(time, status) ~ score,
    data = data,
    reverse = TRUE
  )
  as.numeric(result$concordance)
}

submission_tidy_cox <- function(model, term_label = NULL) {
  summary <- summary(model)
  rows <- data.frame(
    term = rownames(summary$coefficients),
    hazard_ratio = summary$conf.int[, "exp(coef)"],
    ci_lower = summary$conf.int[, "lower .95"],
    ci_upper = summary$conf.int[, "upper .95"],
    p_value = summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  if (!is.null(term_label)) {
    rows$model <- term_label
  }
  rows
}

submission_run_tcga_sensitivity <- function(
    tcga,
    ml,
    project_config,
    settings,
    output_dir
) {
  require_namespace("survival", "TCGA model sensitivity")
  require_namespace("glmnet", "TCGA LASSO stability")
  expression <- read_tcga_expression(project_config$tcga$expression_path)
  clinical <- read_tcga_clinical(project_config$tcga$clinical_path)
  prepared <- prepare_tcga_survival(expression, clinical, ml$final_genes)
  selected <- intersect(tcga$lasso$selected$gene, colnames(prepared$expression))
  if (length(selected) == 0L) {
    stop("TCGA sensitivity requires at least one selected gene.", call. = FALSE)
  }

  model_data <- prepared$survival
  standardized <- scale(prepared$expression[, selected, drop = FALSE])
  colnames(standardized) <- paste0("z_", make.names(selected))
  model_data <- cbind(model_data, as.data.frame(standardized))
  selected_terms <- colnames(standardized)
  selected_formula <- stats::as.formula(paste(
    "survival::Surv(time, status) ~",
    paste(selected_terms, collapse = " + ")
  ))
  selected_model <- survival::coxph(
    selected_formula,
    data = model_data,
    ties = "efron",
    x = TRUE
  )
  selected_model_table <- submission_tidy_cox(
    selected_model,
    "selected_gene_multivariable"
  )

  risk <- tcga$lasso$risk
  risk_model_data <- merge(
    prepared$survival,
    risk[, c("patient_id", "risk_score")],
    by = "patient_id",
    all.x = TRUE,
    sort = FALSE
  )
  risk_model_data$risk_z <- as.numeric(scale(risk_model_data$risk_score))
  risk_model <- survival::coxph(
    survival::Surv(time, status) ~ risk_z,
    data = risk_model_data,
    ties = "efron"
  )
  risk_table <- submission_tidy_cox(
    risk_model,
    "continuous_lasso_risk_per_sd"
  )
  time_varying_model <- survival::coxph(
    survival::Surv(time, status) ~ risk_z + tt(risk_z),
    data = risk_model_data,
    ties = "efron",
    tt = function(value, time, ...) {
      value * log(pmax(time, 1))
    }
  )
  time_varying_table <- submission_tidy_cox(
    time_varying_model,
    "continuous_risk_time_varying"
  )
  ph <- survival::cox.zph(risk_model)
  ph_table <- data.frame(
    term = rownames(ph$table),
    chisq = ph$table[, "chisq"],
    df = ph$table[, "df"],
    p_value = ph$table[, "p"],
    stringsAsFactors = FALSE
  )

  clinical_index <- match(risk_model_data$patient_id, clinical$PATIENT_ID)
  risk_model_data$age <- suppressWarnings(
    as.numeric(clinical$AGE[clinical_index])
  )
  stage <- toupper(as.character(clinical$CLINICAL_STAGE[clinical_index]))
  risk_model_data$stage_group <- ifelse(
    grepl("III|IV", stage),
    "advanced_III_IV",
    ifelse(grepl("I|II", stage), "early_I_II", NA_character_)
  )
  risk_model_data$stage_group <- factor(risk_model_data$stage_group)
  risk_model_data$age_z <- as.numeric(scale(risk_model_data$age))
  adjusted_data <- risk_model_data[
    stats::complete.cases(
      risk_model_data[, c("time", "status", "risk_z", "age_z", "stage_group")]
    ),
    ,
    drop = FALSE
  ]
  adjusted_table <- data.frame()
  adjusted_ph_table <- data.frame()
  if (
    nrow(adjusted_data) >= 100L &&
      length(unique(adjusted_data$stage_group)) == 2L
  ) {
    adjusted_model <- survival::coxph(
      survival::Surv(time, status) ~ risk_z + age_z + stage_group,
      data = adjusted_data,
      ties = "efron"
    )
    adjusted_table <- submission_tidy_cox(
      adjusted_model,
      "risk_adjusted_for_age_and_stage"
    )
    adjusted_ph <- survival::cox.zph(adjusted_model)
    adjusted_ph_table <- data.frame(
      term = rownames(adjusted_ph$table),
      chisq = adjusted_ph$table[, "chisq"],
      df = adjusted_ph$table[, "df"],
      p_value = adjusted_ph$table[, "p"],
      stringsAsFactors = FALSE
    )
  }

  bootstrap_replicates <- as.integer(
    settings$optimism_bootstrap_replicates
  )
  set.seed(as.integer(project_config$project$seed) + 50000L)
  original_prediction <- as.numeric(stats::predict(
    selected_model,
    newdata = model_data,
    type = "lp"
  ))
  apparent_cindex <- submission_cindex(
    model_data$time,
    model_data$status,
    original_prediction
  )
  bootstrap_rows <- vector("list", bootstrap_replicates)
  for (replicate in seq_len(bootstrap_replicates)) {
    indices <- submission_stratified_sample(model_data$status, replace = TRUE)
    bootstrap_data <- model_data[indices, , drop = FALSE]
    fitted <- tryCatch(
      survival::coxph(
        selected_formula,
        data = bootstrap_data,
        ties = "efron"
      ),
      error = function(error) NULL
    )
    if (is.null(fitted)) next
    bootstrap_prediction <- as.numeric(stats::predict(
      fitted,
      newdata = bootstrap_data,
      type = "lp"
    ))
    test_prediction <- as.numeric(stats::predict(
      fitted,
      newdata = model_data,
      type = "lp"
    ))
    bootstrap_rows[[replicate]] <- data.frame(
      replicate = replicate,
      bootstrap_cindex = submission_cindex(
        bootstrap_data$time,
        bootstrap_data$status,
        bootstrap_prediction
      ),
      test_cindex = submission_cindex(
        model_data$time,
        model_data$status,
        test_prediction
      ),
      stringsAsFactors = FALSE
    )
  }
  bootstrap <- do.call(rbind, Filter(Negate(is.null), bootstrap_rows))
  bootstrap$optimism <- bootstrap$bootstrap_cindex - bootstrap$test_cindex
  optimism_summary <- data.frame(
    apparent_cindex = apparent_cindex,
    mean_optimism = mean(bootstrap$optimism, na.rm = TRUE),
    optimism_corrected_cindex =
      apparent_cindex - mean(bootstrap$optimism, na.rm = TRUE),
    test_cindex_median = stats::median(bootstrap$test_cindex, na.rm = TRUE),
    test_cindex_lower_95 = submission_quantile(
      bootstrap$test_cindex,
      0.025
    ),
    test_cindex_upper_95 = submission_quantile(
      bootstrap$test_cindex,
      0.975
    ),
    successful_bootstrap_replicates = nrow(bootstrap),
    stringsAsFactors = FALSE
  )

  lasso_replicates <- as.integer(
    settings$lasso_selection_bootstrap_replicates
  )
  lasso_inner_folds <- as.integer(settings$lasso_inner_folds)
  x <- prepared$expression
  y_time <- prepared$survival$time
  y_status <- prepared$survival$status
  lasso_rows <- list()
  lasso_index <- 0L
  set.seed(as.integer(project_config$project$seed) + 60000L)
  for (replicate in seq_len(lasso_replicates)) {
    indices <- submission_stratified_sample(y_status, replace = TRUE)
    folds <- stratified_fold_ids(
      factor(y_status[indices], levels = c(0, 1)),
      lasso_inner_folds,
      as.integer(project_config$project$seed) + 60000L + replicate
    )
    fitted <- tryCatch(
      suppressWarnings(glmnet::cv.glmnet(
        x = x[indices, , drop = FALSE],
        y = survival::Surv(y_time[indices], y_status[indices]),
        family = "cox",
        alpha = 1,
        foldid = folds,
        standardize = TRUE
      )),
      error = function(error) NULL
    )
    if (is.null(fitted)) next
    for (rule in c("lambda.1se", "lambda.min")) {
      coefficients <- as.matrix(stats::coef(fitted, s = rule))
      selected_genes <- rownames(coefficients)[coefficients[, 1L] != 0]
      if (length(selected_genes) == 0L) next
      lasso_index <- lasso_index + 1L
      lasso_rows[[lasso_index]] <- data.frame(
        replicate = replicate,
        lambda_rule = rule,
        gene = selected_genes,
        stringsAsFactors = FALSE
      )
    }
  }
  lasso_selections <- do.call(rbind, lasso_rows)
  lasso_frequency <- as.data.frame(table(
    lasso_selections$lambda_rule,
    lasso_selections$gene
  ), stringsAsFactors = FALSE)
  names(lasso_frequency) <- c("lambda_rule", "gene", "selected_replicates")
  lasso_frequency <- lasso_frequency[
    lasso_frequency$selected_replicates > 0,
    ,
    drop = FALSE
  ]
  lasso_frequency$selection_frequency <-
    lasso_frequency$selected_replicates / lasso_replicates

  safe_write_csv(
    rbind(
      selected_model_table,
      risk_table,
      time_varying_table,
      adjusted_table
    ),
    file.path(output_dir, "tcga_cox_model_sensitivity.csv")
  )
  safe_write_csv(
    ph_table,
    file.path(output_dir, "tcga_continuous_risk_ph_test.csv")
  )
  safe_write_csv(
    adjusted_ph_table,
    file.path(output_dir, "tcga_adjusted_model_ph_test.csv")
  )
  safe_write_csv(
    bootstrap,
    file.path(output_dir, "tcga_optimism_bootstrap_replicates.csv")
  )
  safe_write_csv(
    optimism_summary,
    file.path(output_dir, "tcga_optimism_bootstrap_summary.csv")
  )
  safe_write_csv(
    lasso_frequency,
    file.path(output_dir, "tcga_lasso_selection_stability.csv")
  )
  list(
    cox = rbind(
      selected_model_table,
      risk_table,
      time_varying_table,
      adjusted_table
    ),
    ph = ph_table,
    adjusted_ph = adjusted_ph_table,
    bootstrap = bootstrap,
    optimism = optimism_summary,
    lasso_frequency = lasso_frequency,
    adjusted_n = nrow(adjusted_data)
  )
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
  comparison_rows <- list()
  composite_rows <- list()
  comparison_index <- 0L
  composite_index <- 0L

  for (dataset in bulk_validation) {
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

    matrix <- dataset$expression[present, , drop = FALSE]
    z <- t(scale(t(matrix)))
    usable <- rowSums(is.finite(z)) == ncol(z)
    z <- z[usable, , drop = FALSE]
    genes <- rownames(z)
    signed_z <- z * sign(logfc[genes])
    score <- colMeans(signed_z)
    composite_roc <- pROC::roc(
      response = dataset$group,
      predictor = score,
      levels = c("Normal", "Disease"),
      direction = "<",
      quiet = TRUE
    )
    interval <- as.numeric(pROC::ci.auc(composite_roc, method = "delong"))
    composite_index <- composite_index + 1L
    composite_rows[[composite_index]] <- data.frame(
      dataset_id = dataset$id,
      disease = dataset$disease,
      genes_in_score = length(genes),
      auc = as.numeric(pROC::auc(composite_roc)),
      ci_lower = interval[[1L]],
      ci_upper = interval[[3L]],
      n_normal = sum(dataset$group == "Normal"),
      n_disease = sum(dataset$group == "Disease"),
      stringsAsFactors = FALSE
    )
  }

  comparison <- do.call(rbind, comparison_rows)
  composite <- do.call(rbind, composite_rows)
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
  list(
    comparison = comparison,
    fixed = fixed_table,
    consistency = consistency,
    composite = composite
  )
}

submission_write_sensitivity_report <- function(results, output_dir, settings) {
  deg_primary <- results$deg$summary[results$deg$summary$is_primary, ]
  wgcna_power <- results$wgcna$perturbations
  wgcna_non_primary <- wgcna_power[!wgcna_power$is_primary, ]
  ml <- results$machine_learning$summary
  tcga <- results$tcga$optimism
  external <- results$external$consistency
  composite <- results$external$composite

  lines <- c(
    "# Non-MR sensitivity analysis report",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- All sensitivity settings were fixed before examining these outputs.",
    "- MR is intentionally excluded here and retained separately as a negative supplementary analysis.",
    "",
    "## DEG thresholds",
    "",
    sprintf(
      "The primary FDR/absolute log2FC rule (%.2f/%.1f) retained %d shared genes, including %d of 10 final hub genes.",
      deg_primary$fdr_threshold,
      deg_primary$absolute_log2fc_threshold,
      deg_primary$shared_count,
      deg_primary$retained_hub_count
    ),
    "The full six-rule grid is reported without selecting a post-hoc preferred threshold.",
    "",
    "## WGCNA stability",
    "",
    paste(
      "Soft-power perturbations:",
      paste(
        sprintf(
          "%s power %d: |r|=%.3f, primary-module retention=%.3f",
          wgcna_non_primary$dataset_id,
          wgcna_non_primary$soft_power,
          wgcna_non_primary$absolute_module_trait_correlation,
          wgcna_non_primary$primary_gene_retention
        ),
        collapse = "; "
      )
    ),
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
    "",
    "## Machine-learning resampling",
    "",
    paste(
      sprintf(
        "%s %s median repeated-CV AUC %.3f (95%% resampling interval %.3f–%.3f)",
        ml$dataset_id,
        ml$model,
        ml$auc_median,
        ml$auc_lower_95_resampling,
        ml$auc_upper_95_resampling
      ),
      collapse = "; "
    ),
    "These are internal resampling estimates within a fixed candidate space and are not independent external performance claims.",
    "",
    "## TCGA-OV model",
    "",
    sprintf(
      "The selected-gene Cox model had apparent C-index %.3f and optimism-corrected C-index %.3f after %d successful bootstrap replicates.",
      tcga$apparent_cindex,
      tcga$optimism_corrected_cindex,
      tcga$successful_bootstrap_replicates
    ),
    paste0(
      "Age/stage-adjusted complete-case sample size: ",
      results$tcga$adjusted_n,
      ". The proportional-hazards diagnostic and a time-varying coefficient model are reported; median-split Kaplan–Meier curves remain illustrative."
    ),
    "",
    "## External validation",
    "",
    sprintf(
      "%d of %d disease–gene pairs had direction-fixed AUC >0.5 in both external cohorts; %d met AUC ≥ %.2f in both.",
      sum(external$all_above_chance),
      nrow(external),
      sum(external$all_at_or_above_robust_threshold),
      as.numeric(settings$external_validation$robust_auc_threshold)
    ),
    paste(
      "Signed 10-gene composite AUCs:",
      paste(
        sprintf(
          "%s %.3f (%.3f–%.3f)",
          composite$dataset_id,
          composite$auc,
          composite$ci_lower,
          composite$ci_upper
        ),
        collapse = "; "
      )
    ),
    "Legacy automatically oriented ROC values are retained only in a bias-audit table; direction-fixed values are used for conclusions.",
    "",
    "## Interpretation boundary",
    "",
    "Sensitivity analyses evaluate robustness of associations and prediction summaries. They do not establish a causal OA–OC relationship.",
    "Null, unstable, and direction-inconsistent findings are retained and will be represented in the manuscript and supplements."
  )
  write_utf8(
    lines,
    file.path(output_dir, "non_mr_sensitivity_report.md")
  )
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
  log_info("Starting predeclared non-MR sensitivity analyses.")

  bulk_training <- submission_load_cache(project_root, "02_bulk_training.rds")
  differential <- submission_load_cache(project_root, "03_differential.rds")
  shared <- submission_load_cache(project_root, "04_shared.rds")
  wgcna <- submission_load_cache(project_root, "06_wgcna.rds")
  ml <- submission_load_cache(project_root, "07_machine_learning.rds")
  bulk_validation <- submission_load_cache(
    project_root,
    "08_bulk_validation.rds"
  )
  validation <- submission_load_cache(project_root, "09_validation.rds")
  tcga <- submission_load_cache(project_root, "11_tcga.rds")

  log_info("Running DEG threshold grid.")
  deg_result <- submission_run_deg_sensitivity(
    differential,
    shared,
    ml$final_genes,
    sensitivity_config$deg,
    output_dir
  )
  log_info("Running WGCNA power perturbation and association stability.")
  wgcna_result <- submission_run_wgcna_sensitivity(
    bulk_training,
    wgcna,
    project_config,
    sensitivity_config$wgcna,
    output_dir
  )
  log_info("Running repeated machine-learning resampling.")
  ml_result <- submission_run_ml_sensitivity(
    bulk_training,
    ml,
    project_config,
    sensitivity_config$machine_learning,
    output_dir
  )
  log_info("Running TCGA model sensitivity.")
  tcga_result <- submission_run_tcga_sensitivity(
    tcga,
    ml,
    project_config,
    sensitivity_config$tcga,
    output_dir
  )
  log_info("Running direction-fixed external validation sensitivity.")
  external_result <- submission_run_external_sensitivity(
    bulk_validation,
    validation,
    shared,
    ml$final_genes,
    sensitivity_config$external_validation,
    output_dir
  )
  results <- list(
    deg = deg_result,
    wgcna = wgcna_result,
    machine_learning = ml_result,
    tcga = tcga_result,
    external = external_result
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
        "ML repeated stratified CV",
        "TCGA optimism bootstrap",
        "TCGA LASSO selection bootstrap",
        "External direction-fixed ROC"
      ),
      setting = c(
        paste(
          length(unlist(sensitivity_config$deg$fdr_thresholds)),
          "FDR x",
          length(unlist(
            sensitivity_config$deg$absolute_log2fc_thresholds
          )),
          "absolute log2FC rules"
        ),
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
        paste(
          sensitivity_config$machine_learning$repeated_cv_replicates,
          "x",
          sensitivity_config$machine_learning$outer_folds,
          "fold"
        ),
        paste(
          sensitivity_config$tcga$optimism_bootstrap_replicates,
          "replicates"
        ),
        paste(
          sensitivity_config$tcga$lasso_selection_bootstrap_replicates,
          "replicates"
        ),
        "training logFC fixes direction; DeLong 95% CI"
      ),
      seed = as.integer(sensitivity_config$seed),
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "sensitivity_parameter_manifest.csv")
  )
  log_info("Non-MR sensitivity analyses completed.")
  invisible(results)
}
