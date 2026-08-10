submission_build_data_sources_original <- submission_build_data_sources
submission_build_supplementary_tables_original <-
  submission_build_supplementary_tables
submission_build_claim_registry_original <- submission_build_claim_registry
submission_build_reproducibility_checklist_original <-
  submission_build_reproducibility_checklist

submission_metadata_field_summary <- function(metadata, field) {
  names_lower <- tolower(names(metadata))
  matches <- switch(
    field,
    age = grepl("(^|[^a-z])age([^a-z]|$)", names_lower),
    sex = grepl("sex|gender", names_lower),
    stage = grepl("stage", names_lower),
    grade = grepl("grade", names_lower),
    treatment = grepl("treatment|therapy|chemotherapy", names_lower) &
      !grepl("protocol", names_lower),
    rep(FALSE, length(names_lower))
  )
  columns <- names(metadata)[matches]
  if (length(columns) == 0L) {
    return("not reported in the local processed metadata")
  }
  observed <- Reduce(
    `|`,
    lapply(columns, function(column) {
      value <- trimws(as.character(metadata[[column]]))
      !is.na(value) & nzchar(value) &
        !tolower(value) %in% c("na", "n/a", "unknown", "not available")
    })
  )
  sprintf(
    "%d/%d records; field(s): %s",
    sum(observed),
    nrow(metadata),
    paste(columns, collapse = "; ")
  )
}

submission_build_data_sources <- function(project_root, output_dir) {
  sources <- submission_build_data_sources_original(
    project_root,
    output_dir
  )
  for (field in c("age", "sex", "stage", "grade", "treatment")) {
    sources[[paste0(field, "_metadata")]] <- "not applicable"
  }
  bulk_training <- submission_load_cache(project_root, "02_bulk_training.rds")
  bulk_validation <- submission_load_cache(
    project_root, "08_bulk_validation.rds"
  )
  bulk <- c(bulk_training, bulk_validation)
  for (dataset in bulk) {
    row <- which(sources$source_id == dataset$id)
    if (length(row) != 1L) next
    for (field in c("age", "sex", "stage", "grade", "treatment")) {
      sources[row, paste0(field, "_metadata")] <-
        submission_metadata_field_summary(dataset$metadata, field)
    }
  }
  tcga_row <- which(sources$source_id == "TCGA-OV")
  if (length(tcga_row) == 1L) {
    sources$age_metadata[tcga_row] <-
      "available; included in age/stage-adjusted complete-case model (n=303)"
    sources$stage_metadata[tcga_row] <-
      "available; included in age/stage-adjusted complete-case model (n=303)"
    sources$sex_metadata[tcga_row] <- "not modeled"
    sources$grade_metadata[tcga_row] <- "not modeled"
    sources$treatment_metadata[tcga_row] <- "not modeled"
  }
  sc_rows <- sources$modality == "single-cell RNA-seq"
  for (field in c("age", "sex", "stage", "grade", "treatment")) {
    sources[sc_rows, paste0(field, "_metadata")] <-
      "not consistently available or used for cross-disease inference"
  }
  hpa <- data.frame(
    source_id = "Human Protein Atlas",
    disease = "not applicable",
    modality = "normal-tissue and normal-cell expression annotation",
    analysis_role = "tissue-background confounding audit",
    source_design = paste(
      "HPA/GTEx consensus categories and HPA single-cell-type categories"
    ),
    analysis_units = "10 prioritized genes",
    public_url = "https://www.proteinatlas.org/humanproteome/tissue/data",
    accession_or_version = "HPA 25.1; Ensembl 109",
    access_date = "2026-07-29",
    notes = paste(
      "Cartilage is absent from the normal-tissue reference;",
      "annotations cannot establish OA or OC disease specificity."
    ),
    age_metadata = "not applicable",
    sex_metadata = "not applicable",
    stage_metadata = "not applicable",
    grade_metadata = "not applicable",
    treatment_metadata = "not applicable",
    stringsAsFactors = FALSE
  )
  sources <- rbind(sources, hpa)
  safe_write_csv(sources, file.path(output_dir, "data_source_manifest.csv"))
  sources
}

submission_build_claim_registry <- function(project_root, output_dir) {
  claims <- submission_build_claim_registry_original(
    project_root,
    output_dir
  )
  claims$manuscript_claim[claims$claim_id == "C05"] <- paste(
    "Strict nested resampling repeated feature screening across all measured",
    "genes inside each outer training fold and still showed near-perfect",
    "internal molecular separation."
  )
  claims$primary_data[claims$claim_id == "C05"] <- paste0(
    "results/submission/sensitivity/",
    "machine_learning_repeated_cv_summary.csv"
  )
  claims$allowed_wording[claims$claim_id == "C05"] <-
    "internal molecular separability under strict nested resampling"
  claims$prohibited_wording[claims$claim_id == "C05"] <-
    "clinical diagnostic accuracy or universally validated biomarkers"
  claims$manuscript_claim[claims$claim_id == "C06"] <- paste(
    "The fixed ten-gene score separated OC tumor and normal samples in two",
    "external cohorts and was stable to 1,000 label permutations and",
    "leave-one-sample-out analysis, but is not a clinical diagnostic model."
  )
  claims$primary_data[claims$claim_id == "C06"] <- paste0(
    "results/submission/sensitivity/external_validation_",
    "signed_composite_score.csv; external_validation_permutation_auc.csv; ",
    "external_validation_leave_one_out_auc.csv"
  )
  claims$allowed_wording[claims$claim_id == "C06"] <-
    "cross-cohort molecular separation with explicit robustness limits"
  claims$manuscript_claim[claims$claim_id == "C08"] <- paste(
    "Single-cell results place candidate expression in distinct cartilage",
    "and ovarian-tumor cellular contexts; eligible OA pseudobulk contrasts",
    "support cell-state-specific associations."
  )
  claims$allowed_wording[claims$claim_id == "C08"] <-
    "distinct cellular contexts and replicated pseudobulk association"
  hpa_claim <- data.frame(
    claim_id = "C13",
    manuscript_claim = paste(
      "HPA normal-reference annotations do not support describing the",
      "ten candidates as uniformly cartilage- or ovary-specific; only",
      "DIRAS3 listed ovary among its group-enriched tissues, and cartilage",
      "was unavailable in the reference."
    ),
    primary_data = paste0(
      "results/submission/sensitivity/hpa_normal_tissue_context.csv"
    ),
    figure_or_table = "Figure S4; Table S14",
    allowed_wording = "normal-tissue background-expression context",
    prohibited_wording = "proof of OA–OC disease specificity",
    status = "verified against current outputs",
    stringsAsFactors = FALSE
  )
  claims <- rbind(claims, hpa_claim)
  safe_write_csv(claims, file.path(output_dir, "claim_evidence_registry.csv"))
  claims
}

submission_build_reproducibility_checklist <- function(
    project_root,
    output_dir
) {
  checklist <- submission_build_reproducibility_checklist_original(
    project_root,
    output_dir
  )
  repository_row <- data.frame(
    item_id = sprintf("R%02d", nrow(checklist) + 1L),
    domain = "repository deposition",
    item = "Public repository URL and archival DOI are finalized.",
    status = "pending author action",
    evidence = paste(
      "The local one-command project is complete;",
      "GitHub/Zenodo deposition has not been performed."
    ),
    stringsAsFactors = FALSE
  )
  checklist <- rbind(checklist, repository_row)
  safe_write_csv(
    checklist,
    file.path(output_dir, "reproducibility_checklist.csv")
  )
  checklist
}

submission_build_supplementary_tables <- function(project_root, output_dir) {
  index <- submission_build_supplementary_tables_original(
    project_root,
    output_dir
  )
  table_dir <- file.path(output_dir, "supplementary_tables")
  sensitivity_dir <- file.path(output_dir, "sensitivity")
  ml <- submission_bind_fill(
    lapply(
      c(
        "machine_learning_repeated_cv_summary.csv",
        "machine_learning_selection_frequency.csv",
        "machine_learning_screen_frequency.csv",
        "machine_learning_outer_fold_fits.csv"
      ),
      function(name) submission_read_csv(file.path(sensitivity_dir, name))
    ),
    c(
      "nested_cv_summary",
      "nested_model_selection_frequency",
      "outer_fold_screen_frequency",
      "outer_fold_fit_audit"
    )
  )
  safe_write_csv(
    ml,
    file.path(table_dir, "Table_S5_machine_learning_resampling.csv")
  )
  external <- submission_bind_fill(
    lapply(
      c(
        "external_validation_direction_fixed_auc.csv",
        "external_validation_cross_cohort_consistency.csv",
        "external_validation_signed_composite_score.csv",
        "external_validation_direction_bias_audit.csv",
        "external_validation_permutation_auc.csv",
        "external_validation_leave_one_out_auc.csv",
        "external_validation_sample_scores.csv"
      ),
      function(name) submission_read_csv(file.path(sensitivity_dir, name))
    ),
    c(
      "gene_auc",
      "cross_cohort_consistency",
      "signed_score",
      "direction_bias_audit",
      "permutation_auc",
      "leave_one_out_auc",
      "sample_scores"
    )
  )
  safe_write_csv(
    external,
    file.path(table_dir, "Table_S6_external_validation.csv")
  )
  submission_copy_table(
    project_root,
    "results/submission/sensitivity/hpa_normal_tissue_context.csv",
    table_dir,
    "Table_S14_HPA_normal_tissue_context.csv"
  )
  index$contents[index$table_id == "Table S1"] <- paste(
    "Public accessions, study roles, analysis units, and availability of",
    "age, sex, stage, grade, and treatment metadata."
  )
  index$contents[index$table_id == "Table S5"] <- paste(
    "Strict nested resampling, outer-fold feature screening, model",
    "selection frequencies, and fit audit."
  )
  index$contents[index$table_id == "Table S6"] <- paste(
    "Direction-fixed gene AUCs and score AUCs, 1,000 label permutations,",
    "leave-one-sample-out influence, and sample-level scores."
  )
  index <- rbind(
    index,
    data.frame(
      table_id = "Table S14",
      filename = "Table_S14_HPA_normal_tissue_context.csv",
      title = "Human Protein Atlas normal-tissue context",
      contents = paste(
        "HPA 25.1 tissue- and cell-type-specificity categories for the",
        "ten candidates, with explicit absence of cartilage."
      ),
      source = paste0(
        "results/submission/sensitivity/hpa_normal_tissue_context.csv"
      ),
      stringsAsFactors = FALSE
    )
  )
  safe_write_csv(
    index,
    file.path(table_dir, "supplementary_table_index.csv")
  )
  index_lines <- c(
    "# Supplementary table index",
    "",
    paste(
      "All tables are UTF-8 CSV files. Interaction tables remain",
      "hypothesis-generating and are not treatment recommendations."
    ),
    "",
    "| Table | File | Title | Contents |",
    "|---|---|---|---|",
    apply(index, 1L, function(row) {
      sprintf(
        "| %s | `%s` | %s | %s |",
        row[["table_id"]],
        row[["filename"]],
        row[["title"]],
        row[["contents"]]
      )
    })
  )
  write_utf8(index_lines, file.path(table_dir, "README.md"))
  index
}
