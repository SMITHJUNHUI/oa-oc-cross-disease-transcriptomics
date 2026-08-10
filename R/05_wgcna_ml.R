select_wgcna_power <- function(dat_expr, settings) {
  configured <- as.numeric(settings$soft_power %||% 0)
  if (configured > 0) {
    return(list(power = configured, fit = NULL, rule = "configured"))
  }

  powers <- c(1:10, seq(12, 20, by = 2))
  fit <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = "signed",
    verbose = 0
  )$fitIndices
  target <- as.numeric(settings$scale_free_r2 %||% 0.80)
  acceptable <- which(fit$SFT.R.sq >= target & fit$slope < 0)
  if (length(acceptable) > 0L) {
    index <- acceptable[[1L]]
    rule <- "first_scale_free_target"
  } else {
    index <- which.max(fit$SFT.R.sq)
    rule <- "maximum_observed_scale_free_fit"
  }
  list(power = fit$Power[[index]], fit = fit, rule = rule)
}

run_wgcna_dataset <- function(dataset, settings, output_dir) {
  require_namespace("WGCNA", "co-expression network analysis")
  expression <- dataset$expression
  variances <- apply(expression, 1L, stats::mad, na.rm = TRUE)
  variances <- variances[is.finite(variances) & variances > 0]
  if (length(variances) < 100L) {
    stop(dataset$id, " has too few variable genes for WGCNA.", call. = FALSE)
  }

  keep_n <- min(
    as.integer(settings$top_variable_genes %||% 5000L),
    length(variances)
  )
  selected <- names(sort(variances, decreasing = TRUE))[seq_len(keep_n)]
  dat_expr <- t(expression[selected, , drop = FALSE])
  quality <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  dat_expr <- dat_expr[
    quality$goodSamples,
    quality$goodGenes,
    drop = FALSE
  ]
  trait <- as.numeric(dataset$group[rownames(dat_expr)] == "Disease")
  if (length(unique(trait)) != 2L) {
    stop(dataset$id, " WGCNA trait is not binary after QC.", call. = FALSE)
  }

  threads <- as.integer(settings$threads %||% 1L)
  try(WGCNA::allowWGCNAThreads(nThreads = max(1L, threads)), silent = TRUE)
  power_selection <- select_wgcna_power(dat_expr, settings)
  if (!is.null(power_selection$fit)) {
    safe_write_csv(
      power_selection$fit,
      file.path(
        output_dir,
        "tables",
        paste0("WGCNA_", dataset$disease, "_", dataset$id, "_power_fit.csv")
      )
    )
  }

  # WGCNA 1.74 resolves a weighted `cor` helper by name inside nested calls.
  # When WGCNA is used only through `::`, match.fun() can otherwise select
  # stats::cor(), which does not accept weights.x/weights.y.
  had_global_cor <- exists("cor", envir = .GlobalEnv, inherits = FALSE)
  previous_global_cor <- if (had_global_cor) {
    get("cor", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  assign("cor", WGCNA::cor, envir = .GlobalEnv)
  on.exit(
    {
      if (had_global_cor) {
        assign("cor", previous_global_cor, envir = .GlobalEnv)
      } else if (exists("cor", envir = .GlobalEnv, inherits = FALSE)) {
        rm("cor", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )

  model <- WGCNA::blockwiseModules(
    dat_expr,
    power = power_selection$power,
    TOMType = "signed",
    networkType = "signed",
    minModuleSize = as.integer(settings$min_module_size %||% 30L),
    mergeCutHeight = as.numeric(settings$merge_cut_height %||% 0.25),
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    maxBlockSize = as.integer(settings$max_block_size %||% 6000L),
    randomSeed = 20260726L,
    nThreads = max(1L, threads),
    verbose = 0
  )

  eigengenes <- WGCNA::orderMEs(model$MEs)
  correlations <- stats::cor(
    eigengenes,
    trait,
    use = "pairwise.complete.obs"
  )
  p_values <- WGCNA::corPvalueStudent(correlations, nrow(dat_expr))
  module_names <- sub("^ME", "", rownames(correlations))
  module_table <- data.frame(
    module = module_names,
    correlation = as.numeric(correlations[, 1L]),
    p_value = as.numeric(p_values[, 1L]),
    stringsAsFactors = FALSE
  )
  module_table$fdr <- stats::p.adjust(module_table$p_value, method = "BH")
  module_table <- module_table[
    module_table$module != "grey",
    ,
    drop = FALSE
  ]
  if (nrow(module_table) == 0L) {
    stop(dataset$id, " produced no non-grey WGCNA modules.", call. = FALSE)
  }
  module_table <- module_table[
    order(abs(module_table$correlation), decreasing = TRUE),
    ,
    drop = FALSE
  ]
  best_module <- module_table$module[[1L]]
  gene_modules <- data.frame(
    gene = colnames(dat_expr),
    module = model$colors,
    stringsAsFactors = FALSE
  )
  selected_genes <- gene_modules$gene[gene_modules$module == best_module]

  prefix <- paste(dataset$disease, dataset$id, sep = "_")
  safe_write_csv(
    module_table,
    file.path(output_dir, "tables", paste0("WGCNA_", prefix, "_module_trait.csv"))
  )
  safe_write_csv(
    gene_modules,
    file.path(output_dir, "tables", paste0("WGCNA_", prefix, "_gene_modules.csv"))
  )
  write_utf8(
    selected_genes,
    file.path(output_dir, "tables", paste0("WGCNA_", prefix, "_best_module_genes.txt"))
  )

  log_info(
    "WGCNA ", dataset$id, ": power ", power_selection$power,
    ", best module ", best_module, " (", length(selected_genes), " genes, r=",
    sprintf("%.3f", module_table$correlation[[1L]]), ")."
  )

  list(
    dataset_id = dataset$id,
    disease = dataset$disease,
    power = power_selection$power,
    power_rule = power_selection$rule,
    best_module = best_module,
    selected_genes = selected_genes,
    module_table = module_table,
    gene_modules = gene_modules
  )
}

run_wgcna_stage <- function(bulk_datasets, config) {
  if (!isTRUE(config$modules$wgcna)) {
    return(list(status = "disabled"))
  }
  training <- bulk_datasets[
    vapply(bulk_datasets, function(x) identical(x$role, "train"), logical(1))
  ]
  results <- lapply(
    training,
    run_wgcna_dataset,
    settings = config$wgcna,
    output_dir = config$project$output_dir
  )
  names(results) <- names(training)
  results
}

prepare_ml_candidates <- function(shared, wgcna_results, settings) {
  table <- shared$table
  table$evidence_score <-
    -log10(pmax(table$adj.P.Val_OA, .Machine$double.xmin)) +
    -log10(pmax(table$adj.P.Val_OC, .Machine$double.xmin)) +
    abs(table$logFC_OA) + abs(table$logFC_OC)
  table <- table[order(table$evidence_score, decreasing = TRUE), , drop = FALSE]

  candidates <- table$gene
  if (!identical(wgcna_results$status %||% "", "disabled")) {
    module_genes <- unique(unlist(lapply(
      wgcna_results,
      function(x) x$selected_genes %||% character()
    )))
    filtered <- intersect(candidates, module_genes)
    minimum <- as.integer(settings$minimum_candidate_genes %||% 5L)
    if (length(filtered) >= minimum) {
      candidates <- filtered
    } else {
      log_warn(
        "WGCNA filtering retained only ", length(filtered),
        " shared genes; using ranked shared DEGs instead."
      )
    }
  }

  maximum <- as.integer(settings$maximum_candidate_genes %||% 100L)
  candidates <- head(candidates, maximum)
  unique(candidates)
}

stratified_fold_ids <- function(group, folds, seed) {
  set.seed(seed)
  fold_id <- integer(length(group))
  for (level in levels(group)) {
    indices <- which(group == level)
    fold_id[indices] <- sample(rep(seq_len(folds), length.out = length(indices)))
  }
  fold_id
}

run_ml_for_dataset <- function(dataset, candidates, settings, seed) {
  require_namespace("glmnet", "LASSO feature selection")
  require_namespace("randomForest", "random-forest feature selection")

  genes <- intersect(candidates, rownames(dataset$expression))
  variances <- apply(
    dataset$expression[genes, , drop = FALSE],
    1L,
    stats::var,
    na.rm = TRUE
  )
  genes <- genes[is.finite(variances) & variances > 0]
  if (length(genes) < 2L) {
    stop(dataset$id, " has fewer than two ML candidate genes.", call. = FALSE)
  }

  x <- t(dataset$expression[genes, , drop = FALSE])
  y_factor <- droplevels(dataset$group[colnames(dataset$expression)])
  y <- as.integer(y_factor == "Disease")
  class_counts <- table(y_factor)
  balance_rule <- settings$class_balance %||% "inverse_frequency"
  observation_weights <- rep(1, length(y_factor))
  forest_class_weights <- NULL
  if (identical(balance_rule, "inverse_frequency")) {
    observation_weights <- 1 / as.numeric(class_counts[y_factor])
    observation_weights <- observation_weights / mean(observation_weights)
    forest_class_weights <- sum(class_counts) /
      (length(class_counts) * as.numeric(class_counts))
    names(forest_class_weights) <- names(class_counts)
  } else if (!identical(balance_rule, "none")) {
    stop("class_balance must be inverse_frequency or none.", call. = FALSE)
  }
  minimum_class <- min(table(y_factor))
  folds <- min(10L, max(3L, as.integer(minimum_class)))
  fold_id <- stratified_fold_ids(y_factor, folds, seed)

  set.seed(seed)
  lasso <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    foldid = fold_id,
    weights = observation_weights,
    type.measure = "deviance",
    standardize = TRUE
  )
  lambda_rule <- settings$lasso_rule %||% "lambda.1se"
  if (!lambda_rule %in% c("lambda.1se", "lambda.min")) {
    stop("lasso_rule must be lambda.1se or lambda.min.", call. = FALSE)
  }
  coefficients <- as.matrix(stats::coef(lasso, s = lambda_rule))
  coefficient_table <- data.frame(
    gene = rownames(coefficients),
    coefficient = as.numeric(coefficients[, 1L]),
    stringsAsFactors = FALSE
  )
  coefficient_table <- coefficient_table[
    coefficient_table$gene != "(Intercept)",
    ,
    drop = FALSE
  ]
  lasso_genes <- coefficient_table$gene[coefficient_table$coefficient != 0]

  set.seed(seed)
  forest <- randomForest::randomForest(
    x = x,
    y = y_factor,
    ntree = as.integer(settings$random_forest_trees %||% 500L),
    classwt = forest_class_weights,
    importance = TRUE
  )
  importance <- randomForest::importance(forest)
  importance_column <- if ("MeanDecreaseGini" %in% colnames(importance)) {
    "MeanDecreaseGini"
  } else {
    colnames(importance)[[ncol(importance)]]
  }
  importance_table <- data.frame(
    gene = rownames(importance),
    importance = as.numeric(importance[, importance_column]),
    stringsAsFactors = FALSE
  )
  importance_table <- importance_table[
    order(importance_table$importance, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rf_genes <- head(
    importance_table$gene,
    as.integer(settings$random_forest_top_n %||% 20L)
  )
  consensus <- intersect(lasso_genes, rf_genes)
  if (length(consensus) == 0L) {
    consensus <- unique(c(lasso_genes, head(rf_genes, 5L)))
  }

  list(
    disease = dataset$disease,
    dataset_id = dataset$id,
    candidates = genes,
    lasso_genes = lasso_genes,
    rf_genes = rf_genes,
    consensus = consensus,
    lasso_coefficients = coefficient_table,
    rf_importance = importance_table
  )
}

run_machine_learning_stage <- function(
    bulk_datasets,
    shared,
    wgcna_results,
    config
) {
  if (!isTRUE(config$modules$machine_learning)) {
    return(list(status = "disabled", final_genes = shared$genes))
  }
  settings <- config$machine_learning
  candidates <- prepare_ml_candidates(shared, wgcna_results, settings)
  training <- bulk_datasets[
    vapply(bulk_datasets, function(x) identical(x$role, "train"), logical(1))
  ]
  results <- lapply(
    seq_along(training),
    function(index) {
      run_ml_for_dataset(
        training[[index]],
        candidates,
        settings,
        seed = config$project$seed + index
      )
    }
  )
  names(results) <- names(training)

  for (result in results) {
    prefix <- paste(result$disease, result$dataset_id, sep = "_")
    safe_write_csv(
      result$lasso_coefficients,
      file.path(
        config$project$output_dir,
        "tables",
        paste0("ML_", prefix, "_lasso_coefficients.csv")
      )
    )
    safe_write_csv(
      result$rf_importance,
      file.path(
        config$project$output_dir,
        "tables",
        paste0("ML_", prefix, "_random_forest_importance.csv")
      )
    )
  }

  disease_consensus <- split(
    results,
    vapply(results, `[[`, character(1), "disease")
  )
  final <- Reduce(
    intersect,
    lapply(disease_consensus, function(x) unique(unlist(lapply(x, `[[`, "consensus"))))
  )
  minimum_final <- as.integer(settings$minimum_final_genes %||% 2L)
  if (length(final) < minimum_final) {
    votes <- table(unlist(lapply(results, function(x) {
      c(x$lasso_genes, x$rf_genes)
    })))
    ranked <- names(sort(votes, decreasing = TRUE))
    final <- unique(c(final, ranked))
    final <- final[final %in% candidates]
    final <- head(final, max(minimum_final, min(10L, length(final))))
  }
  if (length(final) == 0L) {
    stop("Machine-learning selection returned no genes.", call. = FALSE)
  }

  write_utf8(
    final,
    file.path(config$project$output_dir, "tables", "final_hub_genes.txt")
  )
  summary <- data.frame(
    gene = final,
    selected_by_OA = vapply(
      final,
      function(gene) any(vapply(
        results,
        function(x) x$disease == "OA" && gene %in% x$consensus,
        logical(1)
      )),
      logical(1)
    ),
    selected_by_OC = vapply(
      final,
      function(gene) any(vapply(
        results,
        function(x) x$disease == "OC" && gene %in% x$consensus,
        logical(1)
      )),
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    summary,
    file.path(config$project$output_dir, "tables", "final_hub_gene_evidence.csv")
  )
  log_info("Machine-learning final genes: ", paste(final, collapse = ", "), ".")

  list(
    candidates = candidates,
    datasets = results,
    final_genes = final,
    evidence = summary
  )
}
