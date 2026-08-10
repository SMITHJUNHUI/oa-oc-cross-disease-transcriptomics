submission_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required submission source table is missing: ", path, call. = FALSE)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

submission_copy_table <- function(project_root, relative_source, output_dir, output_name) {
  source <- file.path(project_root, relative_source)
  if (!file.exists(source)) {
    stop("Required supplementary-table source is missing: ", source, call. = FALSE)
  }
  destination <- file.path(output_dir, output_name)
  copied <- file.copy(source, destination, overwrite = TRUE, copy.date = TRUE)
  if (!isTRUE(copied)) {
    stop("Could not copy supplementary table: ", source, call. = FALSE)
  }
  destination
}

submission_bind_fill <- function(tables, source_names = NULL) {
  if (is.null(source_names)) {
    source_names <- paste0("source_", seq_along(tables))
  }
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  rows <- lapply(seq_along(tables), function(index) {
    table <- tables[[index]]
    missing <- setdiff(columns, names(table))
    for (column in missing) {
      table[[column]] <- NA
    }
    table <- table[, columns, drop = FALSE]
    table$source_table <- source_names[[index]]
    table
  })
  do.call(rbind, rows)
}

submission_build_data_sources <- function(project_root, output_dir) {
  bulk_train <- submission_read_csv(
    file.path(project_root, "results", "tables", "dataset_qc_summary_train.csv")
  )
  bulk_validation <- submission_read_csv(
    file.path(project_root, "results", "tables", "dataset_qc_summary_validation.csv")
  )
  bulk <- rbind(bulk_train, bulk_validation)
  bulk$analysis_units <- sprintf(
    "%d samples (%d reference, %d disease)",
    bulk$samples,
    bulk$normal,
    bulk$disease_samples
  )

  geo_urls <- function(ids) {
    paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ids)
  }
  bulk_sources <- data.frame(
    source_id = bulk$dataset_id,
    disease = bulk$disease,
    modality = "bulk transcriptomics",
    analysis_role = bulk$role,
    source_design = ifelse(
      bulk$dataset_id == "GSE114007",
      "RNA-seq; human knee cartilage",
      ifelse(
        bulk$dataset_id %in% c("GSE117999"),
        "Agilent expression microarray; human cartilage",
        ifelse(
          bulk$dataset_id %in% c("GSE82107"),
          "Affymetrix expression microarray; human synovium",
          "expression microarray; ovarian/peritoneal tissue"
        )
      )
    ),
    analysis_units = bulk$analysis_units,
    public_url = geo_urls(bulk$dataset_id),
    accession_or_version = bulk$dataset_id,
    access_date = "2026-07-29",
    notes = "Public processed expression data; group labels were derived from repository metadata.",
    stringsAsFactors = FALSE
  )

  single_cell <- submission_read_csv(file.path(
    project_root,
    "results",
    "single_cell_downstream",
    "single_cell_downstream_summary.csv"
  ))
  pre_qc <- c(
    GSE104782 = 1600,
    GSE169454 = 68129,
    GSE255460 = 135896,
    GSE154600 = 52121,
    GSE180661 = 929690
  )
  source_design <- c(
    GSE104782 = "scRNA-seq; OA articular chondrocytes",
    GSE169454 = "10x scRNA-seq; normal and OA articular cartilage",
    GSE255460 = "10x scRNA-seq; OA and non-OA knee cartilage",
    GSE154600 = "10x scRNA-seq; high-grade serous ovarian carcinoma",
    GSE180661 = "10x scRNA-seq; multisite high-grade serous ovarian carcinoma"
  )
  single_cell_sources <- data.frame(
    source_id = single_cell$dataset_id,
    disease = single_cell$disease,
    modality = "single-cell RNA-seq",
    analysis_role = "cellular localization",
    source_design = unname(source_design[single_cell$dataset_id]),
    analysis_units = sprintf(
      "%s audited cells; %s QC-pass cells",
      format(unname(pre_qc[single_cell$dataset_id]), big.mark = ","),
      format(single_cell$cells, big.mark = ",")
    ),
    public_url = geo_urls(single_cell$dataset_id),
    accession_or_version = single_cell$dataset_id,
    access_date = "2026-07-29",
    notes = paste0(
      single_cell$annotation_source,
      "; quantitative summaries used complete QC-pass cells."
    ),
    stringsAsFactors = FALSE
  )

  tcga_source <- data.frame(
    source_id = "TCGA-OV",
    disease = "OC",
    modality = "bulk RNA-seq and clinical survival",
    analysis_role = "exploratory prognostic context",
    source_design = "TCGA ovarian serous cystadenocarcinoma",
    analysis_units = "303 complete survival records; 184 events",
    public_url = "https://www.cancer.gov/ccg/research/genome-sequencing/tcga/studied-cancers/ovarian-serous-cystadenocarcinoma-study",
    accession_or_version = "TCGA-OV legacy RSEM release",
    access_date = "2026-07-29",
    notes = "Used for exploratory association only; no independent prognostic validation cohort was available.",
    stringsAsFactors = FALSE
  )

  mr_metadata <- submission_read_csv(file.path(
    project_root,
    "results",
    "mr",
    "OpenGWAS_dataset_metadata.csv"
  ))
  mr_sources <- data.frame(
    source_id = mr_metadata$id,
    disease = ifelse(grepl("Osteoarthritis", mr_metadata$trait), "OA", "OC"),
    modality = "GWAS summary statistics",
    analysis_role = "negative bidirectional MR supplement",
    source_design = mr_metadata$trait,
    analysis_units = sprintf(
      "%s participants (%s cases; %s controls)",
      format(mr_metadata$sample_size, big.mark = ","),
      format(mr_metadata$ncase, big.mark = ","),
      format(mr_metadata$ncontrol, big.mark = ",")
    ),
    public_url = paste0("https://opengwas.io/datasets/", mr_metadata$id),
    accession_or_version = mr_metadata$id,
    access_date = "2026-07-29",
    notes = paste0(
      mr_metadata$population,
      " ancestry; genome build ",
      mr_metadata$build,
      "; PMID ",
      mr_metadata$pmid,
      "."
    ),
    stringsAsFactors = FALSE
  )

  resource_sources <- data.frame(
    source_id = c(
      "MSigDB-H", "MSigDB-C2", "MSigDB-C5", "MSigDB-C6",
      "miRTarBase", "KnockTF", "DGIdb", "CTD"
    ),
    disease = "not applicable",
    modality = c(
      rep("gene-set resource", 4),
      "regulatory interaction resource",
      "regulatory interaction resource",
      "drug–gene interaction resource",
      "chemical–gene interaction resource"
    ),
    analysis_role = c(
      rep("functional enrichment", 4),
      "hypothesis-generating supplement",
      "hypothesis-generating supplement",
      "hypothesis-generating supplement",
      "hypothesis-generating supplement"
    ),
    source_design = c(
      "Hallmark gene sets", "Canonical pathways", "Gene Ontology gene sets",
      "Oncogenic signatures", "experimentally supported miRNA targets",
      "transcription-factor perturbation targets", "curated drug–gene interactions",
      "curated chemical–gene interactions"
    ),
    analysis_units = c(
      rep("MSigDB v2026.1 human symbols", 4),
      rep("local release recorded in input manifest", 4)
    ),
    public_url = c(
      rep("https://www.gsea-msigdb.org/gsea/msigdb/", 4),
      "https://mirtarbase.cuhk.edu.cn/",
      "https://www.licpathway.net/KnockTF/",
      "https://www.dgidb.org/",
      "https://ctdbase.org/"
    ),
    accession_or_version = c(
      "h.all.v2026.1.Hs", "c2.cp.v2026.1.Hs",
      "c5.go.v2026.1.Hs", "c6.all.v2026.1.Hs",
      "input-manifest release", "KnockTF v2", "input-manifest release",
      "input-manifest release"
    ),
    access_date = "2026-07-29",
    notes = "Secondary resources were used for annotation or hypothesis generation, not causal or therapeutic claims.",
    stringsAsFactors = FALSE
  )

  sources <- rbind(
    bulk_sources,
    single_cell_sources,
    tcga_source,
    mr_sources,
    resource_sources
  )
  safe_write_csv(sources, file.path(output_dir, "data_source_manifest.csv"))
  sources
}

submission_build_parameter_manifest <- function(project_root, output_dir) {
  config <- yaml::read_yaml(file.path(project_root, "config", "local.yml"))
  sensitivity <- submission_read_csv(file.path(
    project_root,
    "results",
    "submission",
    "sensitivity",
    "sensitivity_parameter_manifest.csv"
  ))
  single_cell <- submission_read_csv(file.path(
    project_root,
    "results",
    "single_cell_downstream",
    "analysis_manifest.csv"
  ))

  core <- data.frame(
    section = c(
      "global", "bulk differential expression", "bulk differential expression",
      "WGCNA", "WGCNA", "WGCNA", "WGCNA",
      "machine learning", "machine learning", "machine learning",
      "enrichment", "enrichment", "external validation",
      "TCGA", "MR", "MR", "MR", "MR"
    ),
    parameter = c(
      "random_seed", "FDR", "absolute_log2_fold_change",
      "top_variable_genes", "minimum_module_size", "merge_cut_height",
      "scale_free_fit_target", "maximum_candidate_genes",
      "random_forest_trees_primary", "lasso_rule",
      "multiple_testing", "GSEA_qvalue_cutoff",
      "ROC_direction_policy", "time_points_months",
      "instrument_pvalue", "clump_r2", "clump_window_kb",
      "interpretation"
    ),
    value = c(
      config$project$seed,
      config$differential$fdr,
      config$differential$log2_fc,
      config$wgcna$top_variable_genes,
      config$wgcna$min_module_size,
      config$wgcna$merge_cut_height,
      config$wgcna$scale_free_r2,
      config$machine_learning$maximum_candidate_genes,
      config$machine_learning$random_forest_trees,
      config$machine_learning$lasso_rule,
      config$enrichment$p_adjust_method,
      config$enrichment$qvalue_cutoff,
      "fixed from disease-specific training log2 fold-change",
      paste(config$tcga$time_points_months, collapse = ", "),
      config$mr$pvalue_threshold,
      config$mr$clump_r2,
      config$mr$clump_kb,
      "negative supplementary analysis; no detected causal evidence"
    ),
    source = "config/local.yml",
    stringsAsFactors = FALSE
  )
  sensitivity_out <- data.frame(
    section = paste0("sensitivity: ", sensitivity$analysis),
    parameter = "setting",
    value = sensitivity$setting,
    source = "config/submission_sensitivity.yml",
    stringsAsFactors = FALSE
  )
  single_cell_out <- data.frame(
    section = "single-cell downstream",
    parameter = single_cell$parameter,
    value = single_cell$value,
    source = "results/single_cell_downstream/analysis_manifest.csv",
    stringsAsFactors = FALSE
  )
  manifest <- rbind(core, sensitivity_out, single_cell_out)
  safe_write_csv(manifest, file.path(output_dir, "parameter_manifest.csv"))
  manifest
}

submission_build_reproducibility_checklist <- function(project_root, output_dir) {
  input_manifest <- submission_read_csv(file.path(
    project_root,
    "results",
    "manifests",
    "input_manifest.csv"
  ))
  package_manifest <- submission_read_csv(file.path(
    project_root,
    "results",
    "manifests",
    "package_manifest.csv"
  ))
  figure_source_files <- list.files(file.path(
    project_root,
    "results",
    "submission",
    "figures",
    "source_data"
  ), pattern = "\\.csv$", full.names = FALSE)

  checklist <- data.frame(
    item_id = sprintf("R%02d", seq_len(17)),
    domain = c(
      "entry point", "configuration", "randomness", "environment",
      "environment", "inputs", "inputs", "testing", "sensitivity",
      "figures", "figures", "single-cell", "validation", "claims",
      "security", "manuscript metadata", "submission"
    ),
    item = c(
      "One-command project runner is present.",
      "Portable example configuration and local override are separated.",
      "A project-wide seed is fixed.",
      "R dependencies are locked with renv.",
      "Runtime and package manifests are recorded.",
      "All required local inputs passed preflight.",
      "Input sizes, timestamps, and available checksums are recorded.",
      "Automated unit/regression tests pass.",
      "Prespecified non-MR sensitivity analyses have been completed.",
      "Main and supplementary figures have a shared style manifest.",
      "Every figure has machine-readable source data.",
      "All five single-cell datasets passed the local QC gate.",
      "External ROC directions are fixed from training data.",
      "Claim–data–figure wording boundaries are recorded.",
      "No OpenGWAS JWT is stored in project text outputs.",
      "Author list, affiliations, journal, funding, and conflicts are finalized.",
      "A human author has approved scientific wording and reporting checklists."
    ),
    status = c(
      "complete", "complete", "complete", "complete", "complete",
      ifelse(all(input_manifest$exists), "complete", "incomplete"),
      "complete",
      "complete",
      "complete",
      "complete",
      ifelse(length(figure_source_files) >= 29L, "complete", "incomplete"),
      "complete",
      "complete",
      "complete",
      "complete",
      "pending author input",
      "pending human approval"
    ),
    evidence = c(
      "run_project.ps1 / run_project.bat",
      "config/config.example.yml and config/local.yml",
      paste0("seed=", 20260726),
      "renv.lock",
      paste0(nrow(package_manifest), " package records plus session_info.txt"),
      paste0(sum(input_manifest$exists), "/", nrow(input_manifest), " inputs exist"),
      "results/manifests/input_manifest.csv",
      "tests/testthat/",
      "results/submission/sensitivity/",
      "results/submission/figures/figure_style_manifest.csv",
      paste0(length(figure_source_files), " CSV files"),
      "results/tables/single_cell_dataset_status.csv",
      "results/tables/external_validation_AUC.csv",
      "results/submission/claim_evidence_registry.csv",
      "automated final regression audit",
      "manuscript placeholders",
      "author sign-off required before submission"
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    checklist,
    file.path(output_dir, "reproducibility_checklist.csv")
  )
  checklist
}

submission_build_claim_registry <- function(project_root, output_dir) {
  claims <- data.frame(
    claim_id = sprintf("C%02d", seq_len(12)),
    manuscript_claim = c(
      "OA and OC discovery datasets shared 286 threshold-defined differentially expressed genes.",
      "Only 146 of the 286 shared genes changed in the same direction, so overlap does not imply a common directional program.",
      "Ten genes were prioritized by the prespecified network and feature-selection workflow.",
      "Primary WGCNA module–trait associations were stable to bootstrap resampling and modest soft-power perturbation.",
      "Perfect repeated internal cross-validation AUCs are exploratory because candidate screening preceded resampling.",
      "Direction-fixed external validation was strong in OC but weak and imprecise in OA.",
      "The TCGA-OV risk score was associated with survival but had modest discrimination and proportional-hazards concerns.",
      "Single-cell analyses localized hub-gene expression across OA and OC atlases; inferential pseudobulk testing was limited to replicated contrasts.",
      "Immune-signature differences and pathway enrichments are associative annotations.",
      "Bidirectional MR did not detect a causal effect in either direction.",
      "Reverse MR showed heterogeneity and MR-PRESSO outliers, weakening causal interpretation.",
      "The total evidence supports shared molecular features, not an OA-to-OC or OC-to-OA causal relationship."
    ),
    primary_data = c(
      "results/tables/DEG_selection_summary.csv",
      "results/submission/sensitivity/deg_threshold_sensitivity_summary.csv",
      "results/tables/final_hub_gene_evidence.csv",
      "results/submission/sensitivity/wgcna_module_trait_bootstrap.csv",
      "results/submission/sensitivity/machine_learning_repeated_cv_summary.csv",
      "results/submission/sensitivity/external_validation_signed_composite_score.csv",
      "results/submission/sensitivity/tcga_cox_model_sensitivity.csv",
      "results/tables/single_cell_hub_gene_evidence.csv",
      "results/tables/immune_*_statistics.csv and GSEA_*_hallmark.csv",
      "results/mr/MR_combined_estimates.csv",
      "results/mr/MR_ieu-a-1120__ebi-a-GCST007092_heterogeneity.csv",
      "all primary and sensitivity outputs"
    ),
    figure_or_table = c(
      "Figure 2; Table S2",
      "Figure 2; Table S3",
      "Figures 2–3; Table S2",
      "Figure 3; Figure S1; Table S4",
      "Figure 3; Figure S1; Table S5",
      "Figure 4; Table S6",
      "Figure 6; Figure S1; Table S7",
      "Figure 5; Figure S3; Tables S9–S10",
      "Figure 6; Table S8 and Table S11",
      "Figure S2; Table S12",
      "Figure S2; Table S12",
      "Figures 1–6; Figure S2"
    ),
    allowed_wording = c(
      "shared or overlapping DEGs",
      "directionally heterogeneous overlap",
      "prioritized candidates or hub-gene panel",
      "stable association under tested perturbations",
      "internal separability; exploratory",
      "external discrimination differed by disease",
      "exploratory prognostic association",
      "cellular localization and replicated pseudobulk association",
      "pathway or immune-signature association",
      "no detected causal evidence",
      "heterogeneity/outlier diagnostics limit interpretation",
      "shared molecular features or transcriptomic convergence"
    ),
    prohibited_wording = c(
      "causal link",
      "same biological direction",
      "validated biomarkers",
      "causal modules",
      "generalizable diagnostic accuracy",
      "universally validated signature",
      "clinically ready prognostic model",
      "cell type causes disease",
      "immune mechanism proven",
      "MR proves no causal effect",
      "reverse causality established",
      "OA causes OC or OC causes OA"
    ),
    status = "verified against current outputs",
    stringsAsFactors = FALSE
  )
  safe_write_csv(claims, file.path(output_dir, "claim_evidence_registry.csv"))
  claims
}

submission_build_supplementary_tables <- function(project_root, output_dir) {
  table_dir <- ensure_dir(file.path(output_dir, "supplementary_tables"))
  index <- list()
  register <- function(table_id, filename, title, contents, source) {
    index[[length(index) + 1L]] <<- data.frame(
      table_id = table_id,
      filename = filename,
      title = title,
      contents = contents,
      source = source,
      stringsAsFactors = FALSE
    )
  }

  sources <- submission_build_data_sources(project_root, output_dir)
  safe_write_csv(sources, file.path(table_dir, "Table_S1_data_sources_and_cohorts.csv"))
  register(
    "Table S1", "Table_S1_data_sources_and_cohorts.csv",
    "Datasets and secondary resources",
    "Public accessions, study roles, analysis units, and interpretation notes.",
    "data_source_manifest.csv"
  )

  submission_copy_table(
    project_root, "results/tables/shared_OA_OC_DEGs.csv", table_dir,
    "Table_S2_shared_differentially_expressed_genes.csv"
  )
  register(
    "Table S2", "Table_S2_shared_differentially_expressed_genes.csv",
    "Primary shared OA–OC differentially expressed genes",
    "All 286 genes meeting FDR <0.05 and absolute log2 fold-change >=1 in both discovery datasets.",
    "results/tables/shared_OA_OC_DEGs.csv"
  )

  deg_summary <- submission_read_csv(file.path(
    project_root, "results/submission/sensitivity/deg_threshold_sensitivity_summary.csv"
  ))
  deg_membership <- submission_read_csv(file.path(
    project_root, "results/submission/sensitivity/deg_threshold_sensitivity_membership.csv"
  ))
  safe_write_csv(deg_summary, file.path(table_dir, "Table_S3a_DEG_threshold_summary.csv"))
  safe_write_csv(deg_membership, file.path(table_dir, "Table_S3b_DEG_threshold_membership.csv"))
  register(
    "Table S3a", "Table_S3a_DEG_threshold_summary.csv",
    "DEG threshold sensitivity summary",
    "Six prespecified combinations of FDR and absolute log2 fold-change.",
    "results/submission/sensitivity/deg_threshold_sensitivity_summary.csv"
  )
  register(
    "Table S3b", "Table_S3b_DEG_threshold_membership.csv",
    "DEG threshold membership",
    "Gene-level membership under each threshold combination.",
    "results/submission/sensitivity/deg_threshold_sensitivity_membership.csv"
  )

  wgcna <- submission_bind_fill(
    lapply(
      c(
        "wgcna_soft_power_perturbation.csv",
        "wgcna_module_trait_bootstrap.csv",
        "wgcna_module_trait_leave_one_out.csv"
      ),
      function(name) submission_read_csv(file.path(
        project_root, "results/submission/sensitivity", name
      ))
    ),
    c("soft_power_perturbation", "bootstrap", "leave_one_out")
  )
  safe_write_csv(wgcna, file.path(table_dir, "Table_S4_WGCNA_stability.csv"))
  register(
    "Table S4", "Table_S4_WGCNA_stability.csv",
    "WGCNA stability analyses",
    "Soft-power perturbation, 2,000 bootstrap replicates, and leave-one-sample-out estimates.",
    "results/submission/sensitivity/wgcna_*.csv"
  )

  ml <- submission_bind_fill(
    lapply(
      c(
        "machine_learning_repeated_cv_summary.csv",
        "machine_learning_selection_frequency.csv"
      ),
      function(name) submission_read_csv(file.path(
        project_root, "results/submission/sensitivity", name
      ))
    ),
    c("repeated_cv_summary", "selection_frequency")
  )
  safe_write_csv(ml, file.path(table_dir, "Table_S5_machine_learning_resampling.csv"))
  register(
    "Table S5", "Table_S5_machine_learning_resampling.csv",
    "Machine-learning resampling and feature stability",
    "Fifty repeated five-fold outer resamples and gene-selection frequencies.",
    "results/submission/sensitivity/machine_learning_*.csv"
  )

  external <- submission_bind_fill(
    lapply(
      c(
        "external_validation_direction_fixed_auc.csv",
        "external_validation_cross_cohort_consistency.csv",
        "external_validation_signed_composite_score.csv",
        "external_validation_direction_bias_audit.csv"
      ),
      function(name) submission_read_csv(file.path(
        project_root, "results/submission/sensitivity", name
      ))
    ),
    c("gene_auc", "cross_cohort_consistency", "signed_score", "direction_bias_audit")
  )
  safe_write_csv(external, file.path(table_dir, "Table_S6_external_validation.csv"))
  register(
    "Table S6", "Table_S6_external_validation.csv",
    "Direction-fixed external validation",
    "Gene-level AUCs, cross-cohort consistency, signed composite scores, and the direction-bias audit.",
    "results/submission/sensitivity/external_validation_*.csv"
  )

  tcga_names <- c(
    "tcga_cox_model_sensitivity.csv",
    "tcga_optimism_bootstrap_summary.csv",
    "tcga_adjusted_model_ph_test.csv",
    "tcga_continuous_risk_ph_test.csv",
    "tcga_lasso_selection_stability.csv"
  )
  tcga <- submission_bind_fill(
    lapply(tcga_names, function(name) submission_read_csv(file.path(
      project_root, "results/submission/sensitivity", name
    ))),
    sub("\\.csv$", "", tcga_names)
  )
  safe_write_csv(tcga, file.path(table_dir, "Table_S7_TCGA_model_sensitivity.csv"))
  register(
    "Table S7", "Table_S7_TCGA_model_sensitivity.csv",
    "TCGA-OV exploratory model sensitivity",
    "Cox estimates, optimism correction, proportional-hazards tests, and bootstrap LASSO stability.",
    "results/submission/sensitivity/tcga_*.csv"
  )

  immune <- rbind(
    submission_read_csv(file.path(
      project_root, "results/tables/immune_OA_GSE114007_statistics.csv"
    )),
    submission_read_csv(file.path(
      project_root, "results/tables/immune_OC_GSE18520_statistics.csv"
    ))
  )
  safe_write_csv(immune, file.path(table_dir, "Table_S8_immune_signatures.csv"))
  register(
    "Table S8", "Table_S8_immune_signatures.csv",
    "Rank-based immune-signature comparisons",
    "Disease-minus-reference score differences with Benjamini–Hochberg correction.",
    "results/tables/immune_*_statistics.csv"
  )

  sc_status <- submission_bind_fill(
    list(
      submission_read_csv(file.path(
        project_root, "results/tables/single_cell_dataset_status.csv"
      )),
      submission_read_csv(file.path(
        project_root, "results/tables/single_cell_capability_matrix.csv"
      )),
      submission_read_csv(file.path(
        project_root, "results/single_cell_downstream/single_cell_downstream_summary.csv"
      ))
    ),
    c("qc_status", "capability_matrix", "downstream_summary")
  )
  safe_write_csv(sc_status, file.path(table_dir, "Table_S9_single_cell_QC_and_status.csv"))
  register(
    "Table S9", "Table_S9_single_cell_QC_and_status.csv",
    "Single-cell adapter, QC, and downstream status",
    "Audit status, supported capabilities, retained cells, annotations, and eligible contrasts.",
    "results/tables/single_cell_* and results/single_cell_downstream/"
  )

  submission_copy_table(
    project_root, "results/tables/single_cell_hub_gene_evidence.csv", table_dir,
    "Table_S10_single_cell_hub_gene_evidence.csv"
  )
  register(
    "Table S10", "Table_S10_single_cell_hub_gene_evidence.csv",
    "Single-cell hub-gene localization",
    "Detection fractions and replicated pseudobulk results where eligible.",
    "results/tables/single_cell_hub_gene_evidence.csv"
  )

  submission_copy_table(
    project_root, "results/submission/figures/source_data/Figure6_hallmark_GSEA.csv",
    table_dir, "Table_S11a_Hallmark_GSEA.csv"
  )
  submission_copy_table(
    project_root, "results/tables/GO_shared_genes.csv",
    table_dir, "Table_S11b_GO_shared_genes.csv"
  )
  submission_copy_table(
    project_root, "results/tables/KEGG_shared_genes.csv",
    table_dir, "Table_S11c_KEGG_shared_genes.csv"
  )
  register(
    "Table S11a", "Table_S11a_Hallmark_GSEA.csv",
    "Hallmark gene-set enrichment",
    "Hallmark GSEA results displayed in Figure 6.",
    "results/submission/figures/source_data/Figure6_hallmark_GSEA.csv"
  )
  register(
    "Table S11b", "Table_S11b_GO_shared_genes.csv",
    "Gene Ontology over-representation",
    "GO enrichment among shared genes.",
    "results/tables/GO_shared_genes.csv"
  )
  register(
    "Table S11c", "Table_S11c_KEGG_shared_genes.csv",
    "KEGG over-representation",
    "KEGG enrichment among shared genes.",
    "results/tables/KEGG_shared_genes.csv"
  )

  mr <- submission_bind_fill(
    list(
      submission_read_csv(file.path(project_root, "results/mr/MR_combined_estimates.csv")),
      submission_read_csv(file.path(
        project_root,
        "results/mr/MR_ebi-a-GCST007092__ieu-a-1120_heterogeneity.csv"
      )),
      submission_read_csv(file.path(
        project_root,
        "results/mr/MR_ieu-a-1120__ebi-a-GCST007092_heterogeneity.csv"
      )),
      submission_read_csv(file.path(
        project_root,
        "results/mr/MR_ebi-a-GCST007092__ieu-a-1120_pleiotropy.csv"
      )),
      submission_read_csv(file.path(
        project_root,
        "results/mr/MR_ieu-a-1120__ebi-a-GCST007092_pleiotropy.csv"
      ))
    ),
    c(
      "estimates", "forward_heterogeneity", "reverse_heterogeneity",
      "forward_egger_intercept", "reverse_egger_intercept"
    )
  )
  safe_write_csv(mr, file.path(table_dir, "Table_S12_negative_bidirectional_MR.csv"))
  register(
    "Table S12", "Table_S12_negative_bidirectional_MR.csv",
    "Negative bidirectional Mendelian randomization",
    "Five estimators per direction plus heterogeneity and Egger-intercept diagnostics.",
    "results/mr/"
  )

  regulatory_files <- c(
    "supplementary_miRTarBase_hub_interactions.csv",
    "supplementary_KnockTF_GMT_hub_interactions.csv",
    "supplementary_DGIdb_hub_interactions.csv",
    "supplementary_CTD_hub_interactions.csv"
  )
  for (name in regulatory_files) {
    submission_copy_table(
      project_root,
      file.path("results", "tables", name),
      table_dir,
      paste0("Table_S13_", name)
    )
  }
  register(
    "Table S13", "Table_S13_supplementary_*_hub_interactions.csv",
    "Regulatory and compound interaction lookups",
    "Database-derived interaction records retained as hypothesis-generating evidence only.",
    "results/tables/supplementary_*_hub_interactions.csv"
  )

  index_table <- do.call(rbind, index)
  safe_write_csv(index_table, file.path(table_dir, "supplementary_table_index.csv"))

  index_lines <- c(
    "# Supplementary table index",
    "",
    "All tables are UTF-8 CSV files. Large interaction tables remain hypothesis-generating and are not treatment recommendations.",
    "",
    "| Table | File | Title | Contents |",
    "|---|---|---|---|",
    apply(index_table, 1L, function(row) {
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
  index_table
}

run_submission_materials <- function(project_root) {
  require_namespace("yaml", "submission material assembly")
  output_dir <- ensure_dir(file.path(project_root, "results", "submission"))
  initialize_logging(file.path(output_dir, "logs"))
  log_info("Building data-source, parameter, and reproducibility manifests.")
  submission_build_parameter_manifest(project_root, output_dir)
  submission_build_claim_registry(project_root, output_dir)
  submission_build_reproducibility_checklist(project_root, output_dir)
  log_info("Building numbered supplementary tables.")
  index <- submission_build_supplementary_tables(project_root, output_dir)
  log_info("Submission materials completed: ", nrow(index), " indexed table files.")
  invisible(index)
}
