summarize_pipeline_results <- function(state) {
  lines <- c()

  bulk_datasets <- c(
    state$bulk_training %||% list(),
    state$bulk_validation %||% list()
  )
  bulk_datasets <- Filter(
    function(dataset) {
      is.list(dataset) &&
        all(c("id", "disease", "role", "expression", "group") %in% names(dataset))
    },
    bulk_datasets
  )
  if (length(bulk_datasets) > 0L) {
    datasets <- bulk_datasets
    lines <- c(
      lines,
      "## Bulk transcriptomic datasets",
      "",
      sprintf("- Loaded datasets: %d", length(datasets))
    )
    for (dataset in datasets) {
      counts <- table(dataset$group)
      lines <- c(
        lines,
        sprintf(
          "- %s (%s, %s): %d genes; %d normal and %d disease samples",
          dataset$id,
          dataset$disease,
          dataset$role,
          nrow(dataset$expression),
          unname(counts["Normal"]),
          unname(counts["Disease"])
        )
      )
    }
    lines <- c(lines, "")
  }

  if (!is.null(state$shared)) {
    lines <- c(
      lines,
      "## Cross-disease candidates",
      "",
      sprintf("- Shared OA-OC DEGs: %d", length(state$shared$genes)),
      sprintf("- Selection rule: %s", state$shared$rule),
      sprintf("- |log2FC| threshold: %s", state$shared$threshold),
      ""
    )
  }

  if (!is.null(state$machine_learning)) {
    genes <- state$machine_learning$final_genes %||% character()
    lines <- c(
      lines,
      "## Feature selection",
      "",
      sprintf("- Final hub genes (%d): %s", length(genes), paste(genes, collapse = ", ")),
      ""
    )
  }

  if (!is.null(state$validation$table) && nrow(state$validation$table) > 0L) {
    validation <- state$validation$table
    lines <- c(
      lines,
      "## External validation",
      "",
      sprintf(
        "- Evaluated %d gene-dataset combinations across %d validation cohorts.",
        nrow(validation),
        length(unique(validation$dataset_id))
      ),
      sprintf(
        "- AUC range: %.3f–%.3f",
        min(validation$auc, na.rm = TRUE),
        max(validation$auc, na.rm = TRUE)
      ),
      ""
    )
  }

  if (!is.null(state$tcga) && identical(state$tcga$status %||% "ok", "ok")) {
    lines <- c(
      lines,
      "## TCGA-OV",
      "",
      sprintf("- Survival cohort: %d patients, %d events.", state$tcga$samples, state$tcga$events),
      ""
    )
  }

  if (!is.null(state$mr) && identical(state$mr$status %||% "disabled", "ok")) {
    mr_direction_lines <- if (
      !is.null(state$mr$directions) &&
        nrow(state$mr$directions) > 0L
    ) {
      paste0(
        "- Directions: ",
        paste0(
          state$mr$directions$exposure_id,
          " -> ",
          state$mr$directions$outcome_id,
          collapse = "; "
        ),
        "."
      )
    } else {
      c(
        sprintf(
          "- Exposure IDs: %s.",
          paste(state$mr$exposure_ids, collapse = ", ")
        ),
        sprintf(
          "- Outcome IDs: %s.",
          paste(state$mr$outcome_ids, collapse = ", ")
        )
      )
    }
    lines <- c(
      lines,
      "## Mendelian randomization",
      "",
      mr_direction_lines,
      sprintf(
        "- MR estimate rows: %d.",
        nrow(state$mr$estimates %||% data.frame())
      ),
      sprintf(
        "- Sensitivity audit rows: %d.",
        nrow(state$mr$sensitivity_audit %||% data.frame())
      ),
      ""
    )
  }

  if (
    !is.null(state$single_cell) &&
      !identical(state$single_cell$status %||% "disabled", "disabled")
  ) {
    single_cell <- state$single_cell
    lines <- c(
      lines,
      "## Single-cell QC gate",
      "",
      sprintf("- Gate status: %s.", single_cell$status),
      sprintf(
        "- Fully validated local datasets: %s.",
        paste(single_cell$validated_datasets %||% character(), collapse = ", ")
      ),
      sprintf(
        "- Blocked/limited local datasets: %s.",
        if (length(single_cell$blocked_datasets %||% character()) == 0L) {
          "none"
        } else {
          paste(single_cell$blocked_datasets, collapse = ", ")
        }
      ),
      sprintf(
        "- Downstream single-cell inference enabled: %s.",
        ifelse(isTRUE(single_cell$downstream_enabled), "yes", "no")
      ),
      ""
    )
  }

  if (
    !is.null(state$single_cell_downstream) &&
      identical(
        state$single_cell_downstream$status %||% "disabled",
        "completed"
      )
  ) {
    downstream <- state$single_cell_downstream
    lines <- c(
      lines,
      "## Single-cell downstream analysis",
      "",
      sprintf(
        "- Completed datasets: %s.",
        paste(downstream$datasets_completed, collapse = ", ")
      ),
      sprintf(
        "- QC-pass cells analysed: %s.",
        format(downstream$cells_analysed, big.mark = ",", scientific = FALSE)
      ),
      "- OA and OC were analysed separately; no cross-disease latent integration was performed.",
      "- Replicate-aware pseudobulk tests were emitted only for cell types meeting the pre-specified sample-count gate.",
      ""
    )
  }

  if (!is.null(state$validation$table) && nrow(state$validation$table) > 0L) {
    validation <- state$validation$table
    auc_line <- grep("^- AUC range:", lines)
    if (length(auc_line) == 1L) {
      lines[auc_line] <- sprintf(
        "- Direction-fixed AUC range: %.3f–%.3f.",
        min(validation$auc, na.rm = TRUE),
        max(validation$auc, na.rm = TRUE)
      )
      lines <- append(
        lines,
        "- Expected expression directions were fixed from the corresponding training dataset before evaluation.",
        after = auc_line
      )
    }
  }

  lines
}

write_session_information <- function(config) {
  path <- file.path(config$project$output_dir, "manifests", "session_info.txt")
  output <- capture.output(utils::sessionInfo())
  write_utf8(output, path)
  path
}

generate_final_report <- function(state, stage_status, config) {
  status_table <- do.call(rbind, stage_status)
  safe_write_csv(
    status_table,
    file.path(config$project$output_dir, "manifests", "pipeline_status.csv")
  )
  write_session_information(config)

  lines <- c(
    paste0("# ", config$project$name, " — pipeline report"),
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- R: ", R.version.string),
    paste0("- Seed: ", config$project$seed),
    "",
    "## Stage status",
    "",
    "| Stage | Status | Source | Seconds |",
    "|---|---|---:|---:|"
  )
  lines[1] <- paste0("# ", config$project$name, " — pipeline report")
  for (index in seq_len(nrow(status_table))) {
    row <- status_table[index, ]
    lines <- c(
      lines,
      sprintf(
        "| %s | %s | %s | %.2f |",
        row$stage,
        row$status,
        row$source,
        as.numeric(row$seconds)
      )
    )
  }
  lines <- c(
    lines,
    "",
    summarize_pipeline_results(state),
    "## Interpretation boundary",
    "",
    "- Outputs are computational evidence and require biological interpretation.",
    "- External validation cohorts must remain independent of feature selection.",
    "- Completed bidirectional MR is retained as a negative supplementary analysis: it did not detect causal evidence and does not prove absence of an effect.",
    "- Shared molecular features, prognostic associations, and single-cell localization must not be described as evidence that OA causes OC or that OC causes OA.",
    "- Drug and regulatory-network tables are hypothesis-generating supplementary evidence, not treatment recommendations.",
    "- Single-cell QC remains a mandatory gate; downstream annotation and pseudobulk inference run only for datasets that passed it, and OA/OC are not integrated across diseases.",
    ""
  )

  report_path <- file.path(
    config$project$output_dir,
    "reports",
    "pipeline_report.md"
  )
  write_utf8(lines, report_path)
  report_path
}
