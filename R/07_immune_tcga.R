rank_signature_scores <- function(expression, signatures) {
  ranked <- apply(
    expression,
    2L,
    rank,
    ties.method = "average",
    na.last = "keep"
  )
  ranked <- (ranked - 1) / max(1, nrow(expression) - 1)
  scores <- lapply(
    signatures,
    function(genes) {
      present <- intersect(normalize_gene_symbols(genes), rownames(expression))
      if (length(present) == 0L) {
        return(rep(NA_real_, ncol(expression)))
      }
      colMeans(ranked[present, , drop = FALSE], na.rm = TRUE)
    }
  )
  scores <- do.call(rbind, scores)
  rownames(scores) <- names(signatures)
  colnames(scores) <- colnames(expression)
  scores
}

analyse_immune_dataset <- function(dataset, signatures, output_dir) {
  scores <- rank_signature_scores(dataset$expression, signatures)
  keep <- rowSums(is.finite(scores)) == ncol(scores)
  scores <- scores[keep, , drop = FALSE]
  if (nrow(scores) == 0L) {
    stop(dataset$id, " has no mappable immune signatures.", call. = FALSE)
  }

  group <- dataset$group[colnames(scores)]
  statistics <- do.call(
    rbind,
    lapply(
      rownames(scores),
      function(signature) {
        values <- scores[signature, ]
        normal <- values[group == "Normal"]
        disease <- values[group == "Disease"]
        test <- stats::wilcox.test(disease, normal, exact = FALSE)
        data.frame(
          signature = signature,
          median_normal = stats::median(normal, na.rm = TRUE),
          median_disease = stats::median(disease, na.rm = TRUE),
          median_difference = stats::median(disease, na.rm = TRUE) -
            stats::median(normal, na.rm = TRUE),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  statistics$fdr <- stats::p.adjust(statistics$p_value, method = "BH")
  statistics$dataset_id <- dataset$id
  statistics$disease <- dataset$disease
  statistics <- statistics[
    order(statistics$fdr, -abs(statistics$median_difference)),
    ,
    drop = FALSE
  ]

  prefix <- paste(dataset$disease, dataset$id, sep = "_")
  safe_write_csv(
    statistics,
    file.path(output_dir, "tables", paste0("immune_", prefix, "_statistics.csv"))
  )
  safe_write_csv(
    as.data.frame(scores),
    file.path(output_dir, "tables", paste0("immune_", prefix, "_scores.csv")),
    row.names = TRUE
  )

  annotation <- data.frame(
    Group = group,
    row.names = names(group)
  )
  grDevices::pdf(
    file.path(output_dir, "figures", paste0("immune_", prefix, "_heatmap.pdf")),
    width = 9,
    height = max(5, 0.35 * nrow(scores) + 2)
  )
  pheatmap::pheatmap(
    scores,
    scale = "row",
    annotation_col = annotation,
    show_colnames = FALSE,
    main = paste("Rank-based immune signatures:", dataset$disease, dataset$id)
  )
  grDevices::dev.off()

  top <- head(statistics$signature, min(6L, nrow(statistics)))
  long <- data.frame(
    signature = rep(top, each = ncol(scores)),
    sample = rep(colnames(scores), times = length(top)),
    score = as.numeric(t(scores[top, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
  long$group <- as.character(group[long$sample])
  figure <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = group, y = score, fill = group)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA) +
    ggplot2::facet_wrap(~ signature, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Normal = "#4DAF4A", Disease = "#E41A1C")) +
    ggplot2::labs(x = NULL, y = "Mean within-sample percentile rank", fill = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  ggplot2::ggsave(
    file.path(output_dir, "figures", paste0("immune_", prefix, "_top_signatures.pdf")),
    figure,
    width = 10,
    height = 7
  )

  list(scores = scores, statistics = statistics)
}

run_immune_stage <- function(bulk_datasets, config) {
  if (!isTRUE(config$modules$immune)) {
    return(list(status = "disabled"))
  }
  signature_path <- file.path(
    config$.project_root,
    "config",
    "immune_signatures.yml"
  )
  signatures <- yaml::read_yaml(signature_path)
  training <- bulk_datasets[
    vapply(bulk_datasets, function(x) identical(x$role, "train"), logical(1))
  ]
  results <- lapply(
    training,
    analyse_immune_dataset,
    signatures = signatures,
    output_dir = config$project$output_dir
  )
  names(results) <- names(training)
  results
}

read_tcga_expression <- function(path) {
  require_namespace("data.table", "TCGA expression matrix")
  raw <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  gene_column <- if ("Hugo_Symbol" %in% names(raw)) {
    "Hugo_Symbol"
  } else {
    names(raw)[[1L]]
  }
  non_expression <- intersect(
    c(gene_column, "Entrez_Gene_Id", "Entrez_Gene_ID"),
    names(raw)
  )
  genes <- raw[[gene_column]]
  expression <- safe_numeric_matrix(
    raw[setdiff(names(raw), non_expression)],
    "TCGA expression"
  )
  expression <- aggregate_expression_by_gene(expression, genes)
  maybe_transform_expression(expression, "auto")
}

read_tcga_clinical <- function(path) {
  require_namespace("data.table", "TCGA clinical data")
  clinical <- data.table::fread(
    path,
    skip = "OTHER_PATIENT_ID",
    data.table = FALSE,
    check.names = FALSE,
    fill = TRUE
  )
  required <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
  missing <- setdiff(required, names(clinical))
  if (length(missing) > 0L) {
    stop(
      "TCGA clinical file is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  clinical
}

prepare_tcga_survival <- function(expression, clinical, genes) {
  present <- intersect(genes, rownames(expression))
  if (length(present) == 0L) {
    stop("No selected genes were present in TCGA expression.", call. = FALSE)
  }
  samples <- colnames(expression)
  patient_ids <- substr(samples, 1L, 12L)
  clinical_index <- match(patient_ids, clinical$PATIENT_ID)
  time <- suppressWarnings(as.numeric(clinical$OS_MONTHS[clinical_index]))
  status_text <- toupper(as.character(clinical$OS_STATUS[clinical_index]))
  status <- as.integer(grepl("(^1)|DECEASED", status_text))
  keep <- !is.na(time) & time > 0 & !is.na(status)
  if (sum(keep) < 20L) {
    stop("Fewer than 20 TCGA samples had usable survival data.", call. = FALSE)
  }

  survival <- data.frame(
    sample = samples[keep],
    patient_id = patient_ids[keep],
    time = time[keep],
    status = status[keep],
    stringsAsFactors = FALSE
  )
  survival <- survival[!duplicated(survival$patient_id), , drop = FALSE]
  matrix <- t(expression[present, survival$sample, drop = FALSE])
  list(survival = survival, expression = matrix, genes = present)
}

run_univariate_cox <- function(prepared) {
  rows <- lapply(
    prepared$genes,
    function(gene) {
      expression <- as.numeric(scale(prepared$expression[, gene]))
      if (any(!is.finite(expression))) {
        return(NULL)
      }
      data <- transform(prepared$survival, expression = expression)
      model <- tryCatch(
        survival::coxph(
          survival::Surv(time, status) ~ expression,
          data = data,
          ties = "efron"
        ),
        error = function(error) NULL
      )
      if (is.null(model)) return(NULL)
      summary <- summary(model)
      interval <- as.numeric(exp(stats::confint(model))[1L, ])
      data.frame(
        gene = gene,
        hazard_ratio_per_sd = exp(stats::coef(model))[[1L]],
        ci_lower = interval[1L],
        ci_upper = interval[2L],
        p_value = summary$coefficients[1L, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  )
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame())
  }
  table <- do.call(rbind, rows)
  table$fdr <- stats::p.adjust(table$p_value, method = "BH")
  table[order(table$fdr), , drop = FALSE]
}

run_lasso_cox <- function(prepared, config) {
  require_namespace("glmnet", "TCGA LASSO Cox model")
  x <- prepared$expression
  y <- survival::Surv(prepared$survival$time, prepared$survival$status)
  if (sum(prepared$survival$status) < 10L || ncol(x) < 2L) {
    return(list(status = "insufficient_events_or_genes"))
  }

  set.seed(config$project$seed)
  folds <- min(10L, max(3L, floor(nrow(x) / 10L)))
  model <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = folds,
    standardize = TRUE
  )
  coefficients <- as.matrix(stats::coef(model, s = "lambda.1se"))
  selected <- rownames(coefficients)[coefficients[, 1L] != 0]
  lambda <- "lambda.1se"
  if (length(selected) == 0L) {
    coefficients <- as.matrix(stats::coef(model, s = "lambda.min"))
    selected <- rownames(coefficients)[coefficients[, 1L] != 0]
    lambda <- "lambda.min"
  }
  if (length(selected) == 0L) {
    return(list(status = "no_selected_genes", model = model))
  }

  risk <- as.numeric(stats::predict(model, newx = x, s = lambda, type = "link"))
  risk_table <- prepared$survival
  risk_table$risk_score <- risk
  risk_table$risk_group <- ifelse(
    risk > stats::median(risk, na.rm = TRUE),
    "High",
    "Low"
  )
  selected_table <- data.frame(
    gene = selected,
    coefficient = coefficients[selected, 1L],
    lambda_rule = lambda,
    stringsAsFactors = FALSE
  )

  list(
    status = "ok",
    model = model,
    selected = selected_table,
    risk = risk_table,
    lambda_rule = lambda
  )
}

plot_tcga_survival <- function(lasso, config) {
  if (!identical(lasso$status, "ok")) return(invisible(NULL))
  survival_was_attached <- "package:survival" %in% search()
  if (!survival_was_attached) {
    suppressPackageStartupMessages(
      library("survival", character.only = TRUE)
    )
    on.exit(
      detach("package:survival", character.only = TRUE, unload = FALSE),
      add = TRUE
    )
  }
  output_dir <- config$project$output_dir
  fit <- survival::survfit(
    survival::Surv(time, status) ~ risk_group,
    data = lasso$risk
  )
  plot <- survminer::ggsurvplot(
    fit,
    data = lasso$risk,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Overall survival (months)",
    title = "TCGA-OV multigene risk model"
  )
  grDevices::pdf(
    file.path(output_dir, "figures", "TCGA_OV_risk_KM.pdf"),
    width = 8,
    height = 8
  )
  print(plot)
  grDevices::dev.off()

  times <- as.numeric(config$tcga$time_points_months %||% c(12, 36, 60))
  times <- times[times < max(lasso$risk$time, na.rm = TRUE)]
  if (length(times) > 0L) {
    roc <- timeROC::timeROC(
      T = lasso$risk$time,
      delta = lasso$risk$status,
      marker = lasso$risk$risk_score,
      cause = 1,
      times = times,
      iid = FALSE
    )
    colours <- grDevices::hcl.colors(length(times), palette = "Dark 3")
    grDevices::pdf(
      file.path(output_dir, "figures", "TCGA_OV_risk_timeROC.pdf"),
      width = 7,
      height = 7
    )
    plot(roc, time = times[[1L]], col = colours[[1L]], title = FALSE)
    if (length(times) > 1L) {
      for (index in 2:length(times)) {
        plot(roc, time = times[[index]], col = colours[[index]], add = TRUE)
      }
    }
    graphics::abline(0, 1, lty = 2, col = "grey60")
    graphics::legend(
      "bottomright",
      legend = paste0(times, " months (AUC=", sprintf("%.2f", roc$AUC), ")"),
      col = colours,
      lty = 1,
      bty = "n"
    )
    grDevices::dev.off()
  }
  invisible(NULL)
}

run_tcga_stage <- function(ml_results, shared, config) {
  if (!isTRUE(config$modules$tcga)) {
    return(list(status = "disabled"))
  }
  require_namespace("survival", "TCGA survival analysis")
  require_namespace("survminer", "Kaplan-Meier visualization")
  require_namespace("timeROC", "time-dependent ROC")

  genes <- ml_results$final_genes %||% shared$genes
  expression <- read_tcga_expression(config$tcga$expression_path)
  clinical <- read_tcga_clinical(config$tcga$clinical_path)
  prepared <- prepare_tcga_survival(expression, clinical, genes)
  safe_write_csv(
    prepared$survival,
    file.path(config$project$output_dir, "tables", "TCGA_OV_survival_cohort.csv")
  )

  univariate <- run_univariate_cox(prepared)
  safe_write_csv(
    univariate,
    file.path(config$project$output_dir, "tables", "TCGA_OV_univariate_Cox.csv")
  )
  lasso <- run_lasso_cox(prepared, config)
  if (identical(lasso$status, "ok")) {
    safe_write_csv(
      lasso$selected,
      file.path(config$project$output_dir, "tables", "TCGA_OV_LASSO_Cox_genes.csv")
    )
    safe_write_csv(
      lasso$risk,
      file.path(config$project$output_dir, "tables", "TCGA_OV_risk_scores.csv")
    )
    plot_tcga_survival(lasso, config)
  } else {
    log_warn("TCGA LASSO Cox model was not fitted: ", lasso$status, ".")
  }

  list(
    genes = prepared$genes,
    samples = nrow(prepared$survival),
    events = sum(prepared$survival$status),
    univariate = univariate,
    lasso = lasso
  )
}
