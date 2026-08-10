v30_communication_base_runner <- run_reviewer_v30

v30_cc_seed <- function(value) {
  values <- utf8ToInt(as.character(value))
  as.integer((sum(values * seq_along(values)) + 20260726) %% .Machine$integer.max)
}

v30_cc_select_cells <- function(annotation, cell_column, type_column, maximum = 150L) {
  annotation <- as.data.frame(annotation)
  annotation[[cell_column]] <- as.character(annotation[[cell_column]])
  annotation[[type_column]] <- as.character(annotation[[type_column]])
  annotation <- annotation[
    nzchar(annotation[[cell_column]]) &
      nzchar(annotation[[type_column]]) &
      !annotation[[type_column]] %in% c("Unassigned", "Ambiguous", "Other"),
    ,
    drop = FALSE
  ]
  selected <- lapply(split(seq_len(nrow(annotation)), annotation[[type_column]]), function(index) {
    set.seed(v30_cc_seed(paste(annotation[[type_column]][index[[1L]]], nrow(annotation), sep = "|")))
    if (length(index) > maximum) sample(index, maximum) else index
  })
  annotation[sort(unlist(selected, use.names = FALSE)), , drop = FALSE]
}

v30_cc_signature <- function(files, parameters) {
  files <- sort(unique(files[file.exists(files)]))
  info <- file.info(files)
  digest::digest(list(
    files = data.frame(
      path = normalizePath(files, winslash = "/", mustWork = TRUE),
      size = info$size,
      modified = as.numeric(info$mtime),
      stringsAsFactors = FALSE
    ),
    parameters = parameters,
    CellChat = as.character(utils::packageVersion("CellChat"))
  ), algo = "sha256")
}

v30_cc_cache_path <- function(cache_dir, dataset_id, sample_id) {
  file.path(
    cache_dir,
    paste0("cellchat_", dataset_id, "_", gsub("[^A-Za-z0-9]+", "_", sample_id), ".rds")
  )
}

v30_cc_read_cache <- function(cache_dir, dataset_id, sample_id, signature) {
  path <- v30_cc_cache_path(cache_dir, dataset_id, sample_id)
  if (!file.exists(path)) return(NULL)
  cached <- readRDS(path)
  if (identical(cached$signature, signature)) cached$value else NULL
}

v30_cc_run_sample <- function(counts, metadata, dataset_id, sample_id, condition) {
  require_namespace("CellChat", "sample-resolved cell-cell communication")
  keep <- table(metadata$cell_type)
  keep <- names(keep)[keep >= 20L]
  selected <- metadata$cell_type %in% keep
  counts <- counts[, selected, drop = FALSE]
  metadata <- metadata[selected, , drop = FALSE]
  if (ncol(counts) < 40L || length(unique(metadata$cell_type)) < 2L) {
    return(list(
      interactions = data.frame(),
      audit = data.frame(
        dataset_id = dataset_id,
        sample = sample_id,
        condition = condition,
        cells = ncol(counts),
        cell_types = length(unique(metadata$cell_type)),
        interactions = 0L,
        total_probability = 0,
        status = "insufficient cells or cell types",
        stringsAsFactors = FALSE
      )
    ))
  }
  rownames(metadata) <- metadata$cell_id
  metadata$samples <- factor(rep(sample_id, nrow(metadata)))
  normalized <- CellChat::normalizeData(counts, scale.factor = 10000, do.log = TRUE)
  object <- CellChat::createCellChat(
    object = normalized,
    meta = metadata,
    group.by = "cell_type",
    datatype = "RNA",
    do.sparse = TRUE
  )
  database_environment <- new.env(parent = emptyenv())
  utils::data("CellChatDB.human", package = "CellChat", envir = database_environment)
  object@DB <- get("CellChatDB.human", envir = database_environment)
  object <- CellChat::subsetData(object)
  object <- CellChat::identifyOverExpressedGenes(
    object,
    thresh.pc = 0,
    thresh.fc = 0,
    thresh.p = 0.05,
    min.cells = 10,
    do.fast = FALSE
  )
  object <- CellChat::identifyOverExpressedInteractions(object)
  object <- CellChat::computeCommunProb(
    object,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE,
    distance.use = FALSE,
    nboot = 50,
    seed.use = v30_cc_seed(paste(dataset_id, sample_id, sep = "|"))
  )
  object <- CellChat::filterCommunication(object, min.cells = 20)
  object <- CellChat::computeCommunProbPathway(object, thresh = 0.05)
  object <- CellChat::aggregateNet(object, thresh = 0.05)
  interactions <- as.data.frame(CellChat::subsetCommunication(object, thresh = 0.05))
  if (nrow(interactions) > 0L) {
    interactions$dataset_id <- dataset_id
    interactions$sample <- sample_id
    interactions$condition <- condition
  }
  list(
    interactions = interactions,
    audit = data.frame(
      dataset_id = dataset_id,
      sample = sample_id,
      condition = condition,
      cells = ncol(counts),
      cell_types = length(unique(metadata$cell_type)),
      interactions = nrow(interactions),
      total_probability = if (nrow(interactions) > 0L) sum(interactions$prob) else 0,
      status = "completed",
      stringsAsFactors = FALSE
    )
  )
}

v30_cc_cached_sample <- function(
    counts,
    metadata,
    dataset_id,
    sample_id,
    condition,
    signature,
    cache_dir
) {
  cache_path <- v30_cc_cache_path(cache_dir, dataset_id, sample_id)
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$signature, signature)) return(cached$value)
  }
  value <- tryCatch(
    v30_cc_run_sample(counts, metadata, dataset_id, sample_id, condition),
    error = function(error) list(
      interactions = data.frame(),
      audit = data.frame(
        dataset_id = dataset_id,
        sample = sample_id,
        condition = condition,
        cells = ncol(counts),
        cell_types = length(unique(metadata$cell_type)),
        interactions = 0L,
        total_probability = 0,
        status = paste0("failed: ", conditionMessage(error)),
        stringsAsFactors = FALSE
      )
    )
  )
  saveRDS(list(signature = signature, value = value), cache_path, compress = "xz")
  value
}

v30_cc_oc_samples <- function(project_root, cache_dir) {
  annotation_path <- file.path(
    project_root, "results", "single_cell_downstream", "GSE154600", "cell_annotations.tsv.gz"
  )
  annotation <- data.table::fread(annotation_path, showProgress = FALSE)
  sce_paths <- sort(list.files(
    file.path(project_root, "results", "cache", "single_cell", "GSE154600"),
    pattern = "_qc_sce\\.rds$",
    full.names = TRUE
  ))
  results <- vector("list", length(sce_paths))
  for (index in seq_along(sce_paths)) {
    path <- sce_paths[[index]]
    sample_id <- sub("_qc_sce\\.rds$", "", basename(path))
    log_info("V3.0 CellChat: GSE154600 sample ", index, "/", length(sce_paths), ".")
    signature <- v30_cc_signature(
      c(annotation_path, path),
      list(dataset = "GSE154600", sample = sample_id, maximum = 150L, nboot = 50L)
    )
    cached <- v30_cc_read_cache(cache_dir, "GSE154600", sample_id, signature)
    if (!is.null(cached)) {
      results[[index]] <- cached
      next
    }
    sample_annotation <- annotation[batch == sample_id]
    sample_annotation <- v30_cc_select_cells(sample_annotation, "cell_id", "cell_type", 150L)
    sce <- readRDS(path)
    raw_counts <- SummarizedExperiment::assay(sce, "counts")
    cell_index <- match(sample_annotation$cell_id, colnames(raw_counts))
    if (anyNA(cell_index)) stop("GSE154600 CellChat annotation mismatch: ", sample_id, call. = FALSE)
    raw_counts <- raw_counts[, cell_index, drop = FALSE]
    counts <- .scd_collapse_gene_symbols(
      raw_counts,
      SummarizedExperiment::rowData(sce)$gene_symbol,
      SummarizedExperiment::rowData(sce)$gene_id
    )
    metadata <- data.frame(
      cell_id = colnames(counts),
      cell_type = sample_annotation$cell_type,
      stringsAsFactors = FALSE
    )
    results[[index]] <- v30_cc_cached_sample(
      counts, metadata, "GSE154600", sample_id, "OC", signature, cache_dir
    )
    rm(sce, counts)
    invisible(gc())
  }
  results
}

v30_cc_oa_samples <- function(project_root, cache_dir) {
  local_config <- file.path(project_root, "config", "local.yml")
  config <- read_project_config(
    project_root,
    if (file.exists(local_config)) "config/local.yml" else "config/config.yml"
  )
  dataset <- .scd_dataset_config(config, "GSE255460")
  source_metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
  bundle <- .sc_gse255460_ensure_sparse_bundle(dataset, config)
  validated <- .sc_gse255460_validate_manifest(bundle, source_metadata)
  annotation_path <- file.path(
    project_root, "results", "single_cell_downstream", "GSE255460",
    "cell_annotations_all_QC_pass.tsv.gz"
  )
  annotation <- data.table::fread(annotation_path, showProgress = FALSE)
  donors <- sort(unique(as.character(annotation$donor)))
  results <- vector("list", length(donors))
  for (index in seq_along(donors)) {
    donor_id <- donors[[index]]
    log_info("V3.0 CellChat: GSE255460 donor ", index, "/", length(donors), ".")
    donor_annotation <- annotation[annotation$donor == donor_id, , drop = FALSE]
    condition <- unique(stats::na.omit(as.character(donor_annotation$trait)))
    condition <- if (length(condition) == 1L) condition else "mixed"
    signature <- v30_cc_signature(
      c(annotation_path, dataset$metadata_path),
      list(
        dataset = "GSE255460", sample = donor_id, maximum = 150L,
        nboot = 50L, donor_filter = "strict_v2"
      )
    )
    cached <- v30_cc_read_cache(cache_dir, "GSE255460", donor_id, signature)
    if (!is.null(cached)) {
      results[[index]] <- cached
      next
    }
    donor_annotation <- v30_cc_select_cells(donor_annotation, "cell_id", "celltype", 150L)
    partitions <- sort(unique(as.character(donor_annotation$ID)))
    count_parts <- list()
    meta_parts <- list()
    for (partition in partitions) {
      imported <- .sc_read_gse255460_partition(bundle, validated, partition, source_metadata)
      part_annotation <- donor_annotation[donor_annotation$ID == partition, , drop = FALSE]
      cell_index <- match(part_annotation$cell_id, colnames(imported$counts))
      if (anyNA(cell_index)) stop("GSE255460 CellChat annotation mismatch: ", donor_id, call. = FALSE)
      count_parts[[partition]] <- imported$counts[, cell_index, drop = FALSE]
      meta_parts[[partition]] <- data.frame(
        cell_id = colnames(imported$counts)[cell_index],
        cell_type = part_annotation$celltype,
        stringsAsFactors = FALSE
      )
    }
    counts <- do.call(cbind, count_parts)
    metadata <- do.call(rbind, meta_parts)
    metadata <- metadata[match(colnames(counts), metadata$cell_id), , drop = FALSE]
    results[[index]] <- v30_cc_cached_sample(
      counts, metadata, "GSE255460", donor_id, condition, signature, cache_dir
    )
    rm(counts, count_parts)
    invisible(gc())
  }
  results
}

v30_cc_gene_tokens <- function(value) {
  value <- toupper(as.character(value))
  unique(unlist(strsplit(value, "[^A-Z0-9.-]+"), use.names = FALSE))
}

v30_cc_consensus <- function(interactions, audit, shared_genes, candidates) {
  if (nrow(interactions) == 0L) stop("No CellChat interactions passed the sample-level gate.", call. = FALSE)
  interactions <- data.table::as.data.table(interactions)
  audit <- data.table::as.data.table(audit)
  interactions[, interaction_key := paste(source, target, interaction_name_2, sep = "|")]
  interactions[, context := paste(dataset_id, condition, sep = ":")]
  audit[, context := paste(dataset_id, condition, sep = ":")]
  totals <- audit[status == "completed", .(samples_total = data.table::uniqueN(sample)), by = context]
  thresholds <- totals[, .(
    context,
    samples_total,
    samples_required = pmax(2L, ceiling(samples_total * 0.5))
  )]
  consensus <- interactions[, .(
    samples_present = data.table::uniqueN(sample),
    median_probability = stats::median(prob),
    total_probability = sum(prob),
    pathway_name = paste(sort(unique(pathway_name)), collapse = ";"),
    ligand = paste(sort(unique(ligand)), collapse = ";"),
    receptor = paste(sort(unique(receptor)), collapse = ";")
  ), by = .(context, dataset_id, condition, source, target, interaction_name_2, interaction_key)]
  consensus <- merge(consensus, thresholds, by = "context", all.x = TRUE)
  consensus[, prevalence := samples_present / samples_total]
  consensus[, consensus_status := samples_present >= samples_required]
  consensus[, inference_boundary := paste0(
    "CellChat v", as.character(utils::packageVersion("CellChat")),
    " sample-level consensus from balanced subsamples; no OA-OC probability comparison, causal signaling, or spatial contact is inferred."
  )]
  pathway <- interactions[, .(
    edge_count = .N,
    pathway_probability = sum(prob)
  ), by = .(context, dataset_id, condition, sample, pathway_name)]
  pathway <- pathway[, .(
    samples_present = data.table::uniqueN(sample),
    median_probability = stats::median(pathway_probability),
    median_edges = as.numeric(stats::median(edge_count))
  ), by = .(context, dataset_id, condition, pathway_name)]
  pathway <- merge(pathway, thresholds, by = "context", all.x = TRUE)
  pathway[, prevalence := samples_present / samples_total]
  pathway[, consensus_status := samples_present >= samples_required]
  anchored <- consensus[consensus_status == TRUE]
  anchored[, ligand_tokens := lapply(ligand, v30_cc_gene_tokens)]
  anchored[, receptor_tokens := lapply(receptor, v30_cc_gene_tokens)]
  anchored[, shared_gene_overlap := vapply(seq_len(.N), function(index) {
    paste(sort(intersect(c(ligand_tokens[[index]], receptor_tokens[[index]]), shared_genes)), collapse = ";")
  }, character(1))]
  anchored[, candidate_overlap := vapply(seq_len(.N), function(index) {
    paste(sort(intersect(c(ligand_tokens[[index]], receptor_tokens[[index]]), candidates)), collapse = ";")
  }, character(1))]
  anchored <- anchored[nzchar(shared_gene_overlap) | nzchar(candidate_overlap)]
  anchored[, c("ligand_tokens", "receptor_tokens") := NULL]
  list(consensus = consensus, pathway = pathway, anchored = anchored, thresholds = thresholds)
}

v30_nichenet_overlay <- function(project_root, consensus, candidates) {
  prior_path <- file.path(
    project_root, "results", "api_cache", "nichenet_v2",
    "ligand_target_matrix_nsga2r_final.rds"
  )
  if (!file.exists(prior_path)) {
    return(list(status = "not run: official NicheNet v2 human prior matrix unavailable", edges = data.frame()))
  }
  prior <- readRDS(prior_path)
  target_index <- match(candidates, rownames(prior))
  target_index <- target_index[!is.na(target_index)]
  if (length(target_index) == 0L) {
    return(list(status = "not run: no candidate target mapped to NicheNet prior", edges = data.frame()))
  }
  consensus <- consensus[consensus$consensus_status, , drop = FALSE]
  ligand_tokens <- unique(unlist(lapply(consensus$ligand, v30_cc_gene_tokens), use.names = FALSE))
  ligand_index <- match(ligand_tokens, colnames(prior))
  ligand_index <- ligand_index[!is.na(ligand_index)]
  if (length(ligand_index) == 0L) {
    return(list(status = "not run: no consensus CellChat ligand mapped to NicheNet prior", edges = data.frame()))
  }
  subset <- as.matrix(prior[target_index, ligand_index, drop = FALSE])
  values <- which(is.finite(subset) & subset > 0, arr.ind = TRUE)
  if (nrow(values) == 0L) {
    return(list(status = "completed: no positive prior weights", edges = data.frame()))
  }
  edges <- data.frame(
    target = rownames(subset)[values[, 1L]],
    ligand = colnames(subset)[values[, 2L]],
    regulatory_potential = subset[values],
    stringsAsFactors = FALSE
  )
  contexts <- unique(consensus[, c("context", "dataset_id", "condition", "source", "ligand")])
  contexts$ligand_token <- lapply(contexts$ligand, v30_cc_gene_tokens)
  contexts <- do.call(rbind, lapply(seq_len(nrow(contexts)), function(index) {
    data.frame(
      context = contexts$context[[index]],
      dataset_id = contexts$dataset_id[[index]],
      condition = contexts$condition[[index]],
      source = contexts$source[[index]],
      ligand = contexts$ligand_token[[index]],
      stringsAsFactors = FALSE
    )
  }))
  edges <- merge(edges, contexts, by = "ligand")
  edges$prior_source <- "NicheNet v2 human ligand-target matrix; Zenodo record 7074291"
  edges$inference_boundary <- paste0(
    "Prior-consistency overlay only: CellChat consensus ligands were linked to the fixed ten-gene set through ",
    "NicheNet regulatory-potential weights; ligand activity, differential regulation, mediation, and causality were not inferred."
  )
  list(status = "completed prior-consistency overlay", edges = edges)
}

v30_cc_build_figures <- function(paths, audit, consensus, pathway, anchored, nichenet) {
  audit <- audit[audit$status == "completed", , drop = FALSE]
  audit$context <- paste(audit$dataset_id, audit$condition, sep = ":")
  display_context <- c(
    "GSE154600:OC" = "OC atlas",
    "GSE255460:Control" = "OA control",
    "GSE255460:OA" = "OA"
  )
  audit$context_display <- unname(display_context[audit$context])
  audit$context_display[is.na(audit$context_display)] <- audit$context[is.na(audit$context_display)]
  p1 <- ggplot2::ggplot(audit, ggplot2::aes(x = interactions, y = total_probability, colour = context_display)) +
    ggplot2::geom_point(size = 2.2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 2, max.overlaps = Inf) +
    ggplot2::labs(
      title = "Sample-resolved CellChat",
      subtitle = "Balanced cell-type subsamples; one point per biological sample",
      x = "Significant CellChat interactions",
      y = "Summed communication probability",
      colour = NULL
    ) + submission_theme(7) + ggplot2::theme(legend.position = "bottom")
  pathway_plot <- pathway[pathway$consensus_status, , drop = FALSE]
  pathway_plot <- pathway_plot[order(-pathway_plot$prevalence, -pathway_plot$median_probability), ]
  keep <- unique(head(pathway_plot$pathway_name, 15L))
  pathway_plot <- pathway_plot[pathway_plot$pathway_name %in% keep, , drop = FALSE]
  pathway_plot$context_display <- unname(display_context[pathway_plot$context])
  pathway_plot$context_display[is.na(pathway_plot$context_display)] <- pathway_plot$context[is.na(pathway_plot$context_display)]
  p2 <- ggplot2::ggplot(
    pathway_plot,
    ggplot2::aes(x = context_display, y = pathway_name, size = median_probability, fill = prevalence)
  ) +
    ggplot2::geom_point(shape = 21, colour = "#344054") +
    ggplot2::scale_fill_gradient(low = "#E8F1F8", high = "#1261A0", limits = c(0, 1)) +
    ggplot2::labs(
      title = "Consensus pathways",
      subtitle = "Context-specific sample-consensus rule",
      x = NULL, y = NULL, fill = "Prevalence", size = "Median probability"
    ) + submission_theme(7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))
  pair <- consensus[consensus$consensus_status, , drop = FALSE]
  pair_probability <- aggregate(
    total_probability ~ context + source + target,
    data = pair,
    FUN = sum
  )
  pair_support <- aggregate(
    samples_present ~ context + source + target,
    data = pair,
    FUN = max
  )
  pair <- merge(
    pair_probability,
    pair_support,
    by = c("context", "source", "target"),
    all = TRUE
  )
  pair <- pair[order(-pair$total_probability), , drop = FALSE]
  pair <- do.call(rbind, lapply(split(pair, pair$context), head, 8L))
  pair$context_display <- unname(display_context[pair$context])
  pair$context_display[is.na(pair$context_display)] <- pair$context[is.na(pair$context_display)]
  pair$pair_label <- paste(pair$source, "\u2192", pair$target)
  pair$pair_context <- paste(pair$pair_label, pair$context_display, sep = " ||| ")
  pair$pair_context <- factor(pair$pair_context, levels = rev(unique(pair$pair_context)))
  p3 <- ggplot2::ggplot(
    pair,
    ggplot2::aes(x = total_probability, y = pair_context, size = samples_present, fill = context_display)
  ) +
    ggplot2::geom_point(shape = 21, colour = "#344054") +
    ggplot2::facet_wrap(~context_display, scales = "free_y", nrow = 1) +
    ggplot2::scale_y_discrete(labels = function(x) sub(" \\|\\|\\| .*$", "", x)) +
    ggplot2::scale_fill_manual(values = c("OC atlas" = "#CC79A7", "OA control" = "#56B4E9", "OA" = "#D55E00")) +
    ggplot2::labs(
      title = "Consensus sender-receiver pairs",
      subtitle = "Top pairs within each context; probabilities are descriptive within atlas",
      x = "Summed communication probability", y = NULL, size = "Sample support", fill = NULL
    ) + submission_theme(6.6) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5.4),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", size = 6.5)
    )
  figure <- submission_panel_tag((p1 | p2) / p3 + patchwork::plot_layout(heights = c(0.9, 1.35)))
  submission_save_plot(
    figure,
    "SupplementaryFigure13_sample_consensus_CellChat",
    paths$figures,
    width_mm = 185,
    height_mm = 185
  )
  if (nrow(nichenet$edges) > 0L) {
    edges <- nichenet$edges
    edges <- edges[order(-edges$regulatory_potential), , drop = FALSE]
    edges <- do.call(rbind, lapply(split(edges, edges$context), head, 30L))
    edges$ligand_context <- paste(edges$context, edges$ligand, sep = " | ")
    p5 <- ggplot2::ggplot(
      edges,
      ggplot2::aes(x = target, y = ligand_context, fill = regulatory_potential)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient(low = "#FFF7EC", high = "#8C2D04") +
      ggplot2::labs(
        title = "CellChat-NicheNet prior-consistency overlay",
        subtitle = "Consensus ligands mapped to the fixed ten-gene set; top prior weights per context",
        x = "Fixed candidate target", y = "Atlas/condition and consensus ligand",
        fill = "Regulatory potential"
      ) + submission_theme(7) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    submission_save_plot(
      submission_panel_tag(p5),
      "SupplementaryFigure14_NicheNet_prior_overlay",
      paths$figures,
      width_mm = 185,
      height_mm = 145
    )
  }
}

v30_cc_update_documentation <- function(paths, nichenet) {
  legends_path <- file.path(paths$figures, "figure_legends.md")
  legends <- readLines(legends_path, warn = FALSE, encoding = "UTF-8")
  existing_communication <- grep(
    "^## Supplementary Figure 13\\.",
    legends
  )
  if (length(existing_communication) > 0L) {
    legends <- legends[seq_len(existing_communication[[1L]] - 1L)]
  }
  additions <- c(
    "",
    "## Supplementary Figure 13. Sample-consensus CellChat communication context",
    "",
    paste0(
      "CellChat was run separately for each biological sample after deterministic cell-type-balanced subsampling. ",
      "Panels show sample-level network size/strength, cross-sample pathway prevalence, and the strongest consensus sender-receiver contexts. ",
      "OA and OC communication probabilities were not compared numerically because tissues, cell labels, and comparator availability differ."
    ),
    "",
    "## Supplementary Figure 14. NicheNet prior-consistency overlay",
    "",
    paste0(
      "CellChat consensus ligands were projected onto the official NicheNet v2 human ligand-target prior for the fixed ten-gene set. ",
      "This is a prior-consistency overlay, not ligand-activity inference, differential regulation, mediation, or causal signaling. Status: ",
      nichenet$status, "."
    )
  )
  writeLines(c(legends, additions), legends_path, useBytes = TRUE)
  index_path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE)
  new_rows <- data.frame(
    table_id = c("Table S26a", "Table S26b", "Table S26c", "Table S26d", "Table S27", "Table S28"),
    filename = c(
      "Table_S26a_CellChat_sample_audit.csv",
      "Table_S26b_CellChat_sample_interactions.csv",
      "Table_S26c_CellChat_consensus_pathways.csv",
      "Table_S26d_shared_DEG_anchored_CellChat.csv",
      "Table_S27_NicheNet_prior_overlay.csv",
      "Table_S28_communication_feasibility_boundaries.csv"
    ),
    title = c(
      "CellChat sample audit", "CellChat sample-level interactions",
      "CellChat consensus pathways", "Shared-DEG-anchored CellChat interactions",
      "NicheNet prior-consistency overlay", "Communication-analysis feasibility and boundaries"
    ),
    contents = c(
      "Cells, labels, interaction counts, probability sums, and status per biological sample.",
      "Significant ligand-receptor interactions retained from each sample-specific CellChat model.",
      "Pathway prevalence and median strength under context-specific sample-consensus thresholds.",
      "Consensus interactions whose ligand or receptor overlaps a shared DEG or fixed candidate.",
      "NicheNet v2 prior weights linking CellChat consensus ligands to the fixed ten-gene set.",
      "Dataset eligibility, design constraints, parameter choices, and inference boundaries."
    ),
    source = c(
      "sample-specific CellChat models", "sample-specific CellChat models",
      "sample-consensus aggregation", "sample-consensus aggregation plus shared DEG set",
      "CellChat consensus plus NicheNet v2 human prior", "prespecified communication-analysis gate"
    ),
    stringsAsFactors = FALSE
  )
  normalized_existing <- sub("^Table ", "", index$table_id)
  normalized_new <- sub("^Table ", "", new_rows$table_id)
  index <- index[!normalized_existing %in% normalized_new, , drop = FALSE]
  safe_write_csv(rbind(index, new_rows), index_path)
}

v30_run_communication_context <- function(project_root, paths) {
  require_namespace("CellChat", "sample-resolved exploratory communication analysis")
  require_namespace("data.table", "communication result aggregation")
  cache_dir <- ensure_dir(file.path(project_root, "results", "cache", "submission_v30", "cellchat"))
  results <- c(v30_cc_oc_samples(project_root, cache_dir), v30_cc_oa_samples(project_root, cache_dir))
  audit <- data.table::rbindlist(lapply(results, `[[`, "audit"), fill = TRUE)
  interactions <- data.table::rbindlist(lapply(results, `[[`, "interactions"), fill = TRUE)
  candidate_table <- v30_candidate_definition(paths)
  candidates <- candidate_table$gene
  shared <- utils::read.csv(file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv"))
  shared_genes <- toupper(as.character(shared$gene[as.logical(shared$primary_shared)]))
  consensus <- v30_cc_consensus(interactions, audit, shared_genes, candidates)
  nichenet <- v30_nichenet_overlay(project_root, consensus$consensus, candidates)
  safe_write_csv(as.data.frame(audit), file.path(paths$tables, "Table_S26a_CellChat_sample_audit.csv"))
  safe_write_csv(as.data.frame(interactions), file.path(paths$tables, "Table_S26b_CellChat_sample_interactions.csv"))
  safe_write_csv(as.data.frame(consensus$pathway), file.path(paths$tables, "Table_S26c_CellChat_consensus_pathways.csv"))
  safe_write_csv(as.data.frame(consensus$anchored), file.path(paths$tables, "Table_S26d_shared_DEG_anchored_CellChat.csv"))
  safe_write_csv(as.data.frame(nichenet$edges), file.path(paths$tables, "Table_S27_NicheNet_prior_overlay.csv"))
  feasibility <- data.frame(
    method = c("CellChat", "NicheNet ligand activity", "NicheNet prior overlay"),
    status = c(
      "included as sample-consensus exploratory analysis",
      "not performed",
      nichenet$status
    ),
    rationale = c(
      "GSE255460 has author labels and 11 biological donors; GSE154600 has five tumor samples with transferred labels. Models are fitted per sample.",
      "OC lacks a symmetric disease/reference receiver contrast, and the fixed ten-gene target set is too small to support an unbiased activity test.",
      "Consensus CellChat ligands can be checked against an external ligand-target prior without claiming activity."
    ),
    inference_boundary = c(
      "No direct OA-OC probability comparison, spatial contact, signaling flux, or causality.",
      "Exclusion prevents a database-driven activity claim unsupported by the available design.",
      "External prior consistency only; no regulation, mediation, or mechanism."
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(feasibility, file.path(paths$tables, "Table_S28_communication_feasibility_boundaries.csv"))
  v30_cc_build_figures(
    paths, as.data.frame(audit), as.data.frame(consensus$consensus),
    as.data.frame(consensus$pathway), as.data.frame(consensus$anchored), nichenet
  )
  v30_cc_update_documentation(paths, nichenet)
  list(audit = audit, interactions = interactions, consensus = consensus, nichenet = nichenet)
}

run_reviewer_v30 <- function(project_root) {
  paths <- v30_communication_base_runner(project_root)
  log_info("Adding sample-consensus CellChat and bounded NicheNet prior context.")
  communication <- v30_run_communication_context(project_root, paths)
  log_info(
    "V3.0 communication context completed: ", nrow(communication$interactions),
    " sample-level CellChat interactions; NicheNet status = ", communication$nichenet$status, "."
  )
  invisible(paths)
}
