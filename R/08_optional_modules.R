run_regulatory_stage <- function(ml_results, shared, config) {
  if (!isTRUE(config$modules$regulatory)) {
    return(list(status = "disabled"))
  }
  genes <- normalize_gene_symbols(ml_results$final_genes %||% shared$genes)
  output_dir <- config$project$output_dir

  log_info("Reading miRTarBase interactions.")
  mirtar <- data.table::fread(
    config$regulatory$mirtarbase_path,
    select = c("miRNA", "Target Gene"),
    data.table = FALSE
  )
  mirtar$target_normalized <- normalize_gene_symbols(mirtar[["Target Gene"]])
  mirtar <- unique(mirtar[
    mirtar$target_normalized %in% genes,
    c("miRNA", "Target Gene", "target_normalized")
  ])
  safe_write_csv(
    mirtar,
    file.path(output_dir, "tables", "supplementary_miRTarBase_hub_interactions.csv")
  )

  log_info("Reading KnockTF GMT interactions.")
  knocktf_gmt <- clusterProfiler::read.gmt(
    config$regulatory$knocktf_gmt_path
  )
  names(knocktf_gmt) <- c("tf", "target")
  knocktf_gmt$target_normalized <- normalize_gene_symbols(knocktf_gmt$target)
  knocktf_gmt <- unique(knocktf_gmt[
    knocktf_gmt$target_normalized %in% genes,
    ,
    drop = FALSE
  ])
  safe_write_csv(
    knocktf_gmt,
    file.path(output_dir, "tables", "supplementary_KnockTF_GMT_hub_interactions.csv")
  )

  log_info("Reading selected columns from detailed KnockTF table.")
  knocktf <- data.table::fread(
    config$regulatory$knocktf_table_path,
    select = c("TF", "Gene", "Log2FC", "P_value", "Corrected_P", "up_down"),
    data.table = FALSE,
    showProgress = FALSE
  )
  knocktf$target_normalized <- normalize_gene_symbols(knocktf$Gene)
  knocktf <- knocktf[
    knocktf$target_normalized %in% genes,
    ,
    drop = FALSE
  ]
  safe_write_csv(
    knocktf,
    file.path(output_dir, "tables", "supplementary_KnockTF_hub_interactions.csv")
  )

  list(
    status = "ok",
    mirtarbase = mirtar,
    knocktf_gmt = knocktf_gmt,
    knocktf = knocktf
  )
}

read_ctd_interactions <- function(path, genes) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header_lines <- readLines(connection, n = 100L, warn = FALSE)
  fields_index <- grep("^# Fields:", header_lines)[1L]
  if (is.na(fields_index)) {
    stop("Could not locate CTD field declaration.", call. = FALSE)
  }
  field_line_index <- which(
    seq_along(header_lines) > fields_index &
      grepl("^#\\s*ChemicalName,", header_lines, perl = TRUE)
  )[1L]
  if (is.na(field_line_index)) {
    stop("Could not locate CTD field names.", call. = FALSE)
  }
  field_names <- strsplit(
    sub("^#\\s*", "", header_lines[[field_line_index]], perl = TRUE),
    ",",
    fixed = TRUE
  )[[1L]]
  candidate_data_lines <- which(
    seq_along(header_lines) > field_line_index &
      nzchar(trimws(header_lines)) &
      !startsWith(header_lines, "#")
  )
  if (length(candidate_data_lines) == 0L) {
    stop("Could not locate the first CTD data row.", call. = FALSE)
  }
  first_data_line <- candidate_data_lines[[1L]]

  table <- data.table::fread(
    path,
    skip = first_data_line - 1L,
    header = FALSE,
    col.names = field_names,
    data.table = FALSE,
    showProgress = FALSE
  )
  keep_columns <- intersect(
    c(
      "ChemicalName", "ChemicalID", "GeneSymbol", "GeneID",
      "Organism", "OrganismID", "Interaction", "InteractionActions",
      "PubMedIDs"
    ),
    names(table)
  )
  table <- table[, keep_columns, drop = FALSE]
  table$gene_normalized <- normalize_gene_symbols(table$GeneSymbol)
  table <- table[table$gene_normalized %in% genes, , drop = FALSE]
  if ("OrganismID" %in% names(table)) {
    table <- table[is.na(table$OrganismID) | table$OrganismID == 9606, , drop = FALSE]
  }
  table
}

run_drug_stage <- function(ml_results, shared, config) {
  if (!isTRUE(config$modules$drug)) {
    return(list(status = "disabled"))
  }
  genes <- normalize_gene_symbols(ml_results$final_genes %||% shared$genes)
  output_dir <- config$project$output_dir

  log_info("Reading DGIdb interactions from ", config$drug$dgidb_path, ".")
  dgidb <- data.table::fread(
    config$drug$dgidb_path,
    skip = "gene_claim_name",
    data.table = FALSE,
    showProgress = FALSE
  )
  required <- c("gene_name", "drug_name")
  if (!all(required %in% names(dgidb))) {
    stop("DGIdb table is missing gene_name or drug_name.", call. = FALSE)
  }
  dgidb$gene_normalized <- normalize_gene_symbols(dgidb$gene_name)
  dgidb <- dgidb[dgidb$gene_normalized %in% genes, , drop = FALSE]
  if (
    isTRUE(config$drug$approved_only) &&
      "drug_is_approved" %in% names(dgidb)
  ) {
    approved <- tolower(as.character(dgidb$drug_is_approved)) == "true"
    dgidb <- dgidb[approved, , drop = FALSE]
  }
  dgidb <- dgidb[!is.na(dgidb$drug_name) & nzchar(dgidb$drug_name), , drop = FALSE]
  safe_write_csv(
    dgidb,
    file.path(output_dir, "tables", "supplementary_DGIdb_hub_interactions.csv")
  )

  log_info("Reading CTD chemical-gene interactions.")
  ctd <- read_ctd_interactions(config$drug$ctd_path, genes)
  safe_write_csv(
    ctd,
    file.path(output_dir, "tables", "supplementary_CTD_hub_interactions.csv")
  )

  common <- character()
  if (nrow(dgidb) > 0L && nrow(ctd) > 0L) {
    common_upper <- intersect(
      toupper(unique(dgidb$drug_name)),
      toupper(unique(ctd$ChemicalName))
    )
    common <- unique(c(
      dgidb$drug_name[toupper(dgidb$drug_name) %in% common_upper],
      ctd$ChemicalName[toupper(ctd$ChemicalName) %in% common_upper]
    ))
  }
  write_utf8(
    sort(common),
    file.path(output_dir, "tables", "supplementary_common_DGIdb_CTD_compounds.txt")
  )

  list(status = "ok", dgidb = dgidb, ctd = ctd, common_compounds = common)
}

.mr_capture_table <- function(label, function_call) {
  tryCatch(
    {
      value <- function_call()
      list(
        label = label,
        status = "completed",
        data = as.data.frame(value),
        error = ""
      )
    },
    error = function(error) {
      list(
        label = label,
        status = "unavailable",
        data = data.frame(),
        error = conditionMessage(error)
      )
    }
  )
}

.mr_skipped_table <- function(label, reason) {
  list(
    label = label,
    status = "not_applicable",
    data = data.frame(),
    error = reason
  )
}

.mr_write_capture <- function(capture, path) {
  if (nrow(capture$data) > 0L) {
    safe_write_csv(capture$data, path)
  }
  data.frame(
    analysis = capture$label,
    status = capture$status,
    rows = nrow(capture$data),
    details = capture$error,
    stringsAsFactors = FALSE
  )
}

run_mr_stage <- function(config) {
  if (!isTRUE(config$modules$mr)) {
    return(list(
      status = "disabled",
      reason = "Enable only after exposure/outcome IDs are pre-specified."
    ))
  }
  require_namespace("TwoSampleMR", "Mendelian randomization")
  require_namespace("ieugwasr", "OpenGWAS metadata and authentication")
  output_dir <- config$project$output_dir
  exposures <- unlist(config$mr$exposure_ids, use.names = FALSE)
  outcomes <- unlist(config$mr$outcome_ids, use.names = FALSE)
  token_name <- config$mr$opengwas_token_env %||% "OPENGWAS_JWT"
  opengwas_jwt <- trimws(Sys.getenv(token_name, unset = ""))
  if (!nzchar(opengwas_jwt)) {
    stop("OpenGWAS JWT is empty in ", token_name, ".", call. = FALSE)
  }
  if (length(exposures) == 0L || length(outcomes) == 0L) {
    stop("MR exposure_ids and outcome_ids must be non-empty.", call. = FALSE)
  }
  mr_dir <- ensure_dir(file.path(output_dir, "mr"))
  metadata <- ieugwasr::gwasinfo(
    unique(c(exposures, outcomes)),
    opengwas_jwt = opengwas_jwt
  )
  safe_write_csv(
    as.data.frame(metadata),
    file.path(mr_dir, "OpenGWAS_dataset_metadata.csv")
  )
  results <- list()
  combined_estimates <- list()
  combined_audit <- list()

  for (exposure_id in exposures) {
    instruments <- TwoSampleMR::extract_instruments(
      outcomes = exposure_id,
      p1 = as.numeric(config$mr$pvalue_threshold %||% 5e-8),
      clump = TRUE,
      r2 = as.numeric(config$mr$clump_r2 %||% 0.001),
      kb = as.numeric(config$mr$clump_kb %||% 10000),
      opengwas_jwt = opengwas_jwt
    )
    if (is.null(instruments) || nrow(instruments) == 0L) {
      stop("No instruments were extracted for ", exposure_id, ".", call. = FALSE)
    }
    instruments$F_statistic <- with(
      instruments,
      (beta.exposure / se.exposure)^2
    )
    for (outcome_id in outcomes) {
      key <- paste(exposure_id, outcome_id, sep = "__")
      prefix <- file.path(mr_dir, paste0("MR_", clean_filename(key)))
      outcome_data <- TwoSampleMR::extract_outcome_data(
        snps = instruments$SNP,
        outcomes = outcome_id,
        opengwas_jwt = opengwas_jwt
      )
      if (is.null(outcome_data) || nrow(outcome_data) == 0L) {
        stop("No outcome associations were found for ", key, ".", call. = FALSE)
      }
      harmonised <- TwoSampleMR::harmonise_data(
        exposure_dat = instruments,
        outcome_dat = outcome_data,
        action = 2
      )
      analysis_data <- harmonised[
        !is.na(harmonised$mr_keep) & harmonised$mr_keep,
        ,
        drop = FALSE
      ]
      if (nrow(analysis_data) == 0L) {
        stop("No harmonised instruments passed mr_keep for ", key, ".", call. = FALSE)
      }
      estimate <- TwoSampleMR::mr(analysis_data)
      odds_ratio <- tryCatch(
        TwoSampleMR::generate_odds_ratios(estimate),
        error = function(error) data.frame()
      )
      safe_write_csv(
        instruments,
        paste0(prefix, "_instruments.csv")
      )
      safe_write_csv(
        outcome_data,
        paste0(prefix, "_outcome_associations.csv")
      )
      safe_write_csv(
        harmonised,
        paste0(prefix, "_harmonised.csv")
      )
      safe_write_csv(
        estimate,
        paste0(prefix, "_estimates.csv")
      )
      if (nrow(odds_ratio) > 0L) {
        safe_write_csv(
          odds_ratio,
          paste0(prefix, "_odds_ratios.csv")
        )
      }

      nsnp <- nrow(analysis_data)
      sensitivity <- list(
        heterogeneity = if (nsnp >= 2L) {
          .mr_capture_table(
            "Cochran_Q_heterogeneity",
            function() TwoSampleMR::mr_heterogeneity(analysis_data)
          )
        } else {
          .mr_skipped_table("Cochran_Q_heterogeneity", "Fewer than 2 SNPs.")
        },
        pleiotropy = if (nsnp >= 3L) {
          .mr_capture_table(
            "MR_Egger_intercept",
            function() TwoSampleMR::mr_pleiotropy_test(analysis_data)
          )
        } else {
          .mr_skipped_table("MR_Egger_intercept", "Fewer than 3 SNPs.")
        },
        single_snp = .mr_capture_table(
          "single_SNP_and_all_methods",
          function() TwoSampleMR::mr_singlesnp(analysis_data)
        ),
        leave_one_out = if (nsnp >= 3L) {
          .mr_capture_table(
            "leave_one_out",
            function() TwoSampleMR::mr_leaveoneout(analysis_data)
          )
        } else {
          .mr_skipped_table("leave_one_out", "Fewer than 3 SNPs.")
        },
        directionality = .mr_capture_table(
          "Steiger_directionality_test",
          function() TwoSampleMR::directionality_test(analysis_data)
        ),
        steiger_filter = .mr_capture_table(
          "Steiger_SNP_filtering",
          function() TwoSampleMR::steiger_filtering(analysis_data)
        )
      )
      audit <- do.call(
        rbind,
        Map(
          function(capture, name) {
            .mr_write_capture(
              capture,
              paste0(prefix, "_", name, ".csv")
            )
          },
          sensitivity,
          names(sensitivity)
        )
      )

      presso <- if (nsnp >= 4L && requireNamespace("MRPRESSO", quietly = TRUE)) {
        tryCatch(
          {
            value <- MRPRESSO::mr_presso(
              BetaOutcome = "beta.outcome",
              BetaExposure = "beta.exposure",
              SdOutcome = "se.outcome",
              SdExposure = "se.exposure",
              data = analysis_data,
              OUTLIERtest = TRUE,
              DISTORTIONtest = TRUE,
              SignifThreshold = 0.05,
              NbDistribution = 1000,
              seed = as.integer(config$project$seed)
            )
            atomic_save_rds(value, paste0(prefix, "_MRPRESSO.rds"))
            write_utf8(
              capture.output(print(value)),
              paste0(prefix, "_MRPRESSO.txt")
            )
            list(status = "completed", value = value, error = "")
          },
          error = function(error) {
            list(
              status = "unavailable",
              value = NULL,
              error = conditionMessage(error)
            )
          }
        )
      } else {
        list(
          status = "not_applicable",
          value = NULL,
          error = "MR-PRESSO requires at least 4 SNPs and the MRPRESSO package."
        )
      }
      audit <- rbind(
        audit,
        data.frame(
          analysis = "MR_PRESSO",
          status = presso$status,
          rows = if (is.null(presso$value)) 0L else 1L,
          details = presso$error,
          stringsAsFactors = FALSE
        )
      )
      audit$exposure_id <- exposure_id
      audit$outcome_id <- outcome_id
      safe_write_csv(audit, paste0(prefix, "_sensitivity_audit.csv"))

      estimate$exposure_id <- exposure_id
      estimate$outcome_id <- outcome_id
      combined_estimates[[key]] <- estimate
      combined_audit[[key]] <- audit
      results[[key]] <- list(
        instruments = instruments,
        outcome_data = outcome_data,
        harmonised = harmonised,
        analysis_data = analysis_data,
        estimates = estimate,
        odds_ratios = odds_ratio,
        sensitivity = sensitivity,
        mr_presso = presso,
        sensitivity_audit = audit
      )
    }
  }
  estimate_table <- data.table::rbindlist(
    combined_estimates,
    use.names = TRUE,
    fill = TRUE
  )
  audit_table <- data.table::rbindlist(
    combined_audit,
    use.names = TRUE,
    fill = TRUE
  )
  estimate_path <- file.path(mr_dir, "MR_combined_estimates.csv")
  audit_path <- file.path(mr_dir, "MR_sensitivity_audit.csv")
  safe_write_csv(estimate_table, estimate_path)
  safe_write_csv(audit_table, audit_path)
  report_path <- file.path(output_dir, "reports", "mr_report.md")
  report_lines <- c(
    "# Mendelian randomization report",
    "",
    paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- Exposure IDs: ", paste(exposures, collapse = ", ")),
    paste0("- Outcome IDs: ", paste(outcomes, collapse = ", ")),
    paste0(
      "- Instrument rule: p<",
      config$mr$pvalue_threshold,
      "; clump r2<",
      config$mr$clump_r2,
      "; window ",
      config$mr$clump_kb,
      " kb."
    ),
    "",
    "## Output",
    "",
    paste0("- Combined estimates: `", estimate_path, "`"),
    paste0("- Sensitivity audit: `", audit_path, "`"),
    "- Per-analysis instruments, harmonised data, heterogeneity, Egger intercept, single-SNP, leave-one-out, Steiger and MR-PRESSO outputs are stored in `results/mr/`.",
    "",
    "## Interpretation boundary",
    "",
    "- Instrument strength is reported using per-SNP F statistics.",
    "- Ancestry compatibility and sample overlap must be confirmed from source studies before causal interpretation.",
    "- MR-PRESSO and other sensitivity tests are diagnostic and do not by themselves prove causality.",
    ""
  )
  write_utf8(report_lines, report_path)
  list(
    status = "ok",
    exposure_ids = exposures,
    outcome_ids = outcomes,
    analyses = results,
    estimates = estimate_table,
    sensitivity_audit = audit_table,
    dataset_metadata = metadata,
    report = report_path
  )
}
