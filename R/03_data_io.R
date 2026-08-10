normalize_gene_symbols <- function(symbols) {
  symbols <- trimws(as.character(symbols))
  symbols <- sub("\\s+///.*$", "", symbols)
  symbols <- sub(";.*$", "", symbols)
  symbols <- sub("\\|.*$", "", symbols)
  symbols <- toupper(trimws(symbols))
  invalid <- is.na(symbols) |
    !nzchar(symbols) |
    symbols %in% c("---", "NA", "N/A", "NULL") |
    grepl("^[0-9]+$", symbols)
  symbols[invalid] <- NA_character_
  symbols
}

aggregate_expression_by_gene <- function(expression, symbols) {
  expression <- safe_numeric_matrix(expression, "expression matrix")
  symbols <- normalize_gene_symbols(symbols)
  keep <- !is.na(symbols)
  expression <- expression[keep, , drop = FALSE]
  symbols <- symbols[keep]

  if (nrow(expression) == 0L) {
    stop("No valid gene symbols remained after annotation.", call. = FALSE)
  }

  missing <- is.na(expression)
  expression_zeroed <- expression
  expression_zeroed[missing] <- 0
  sums <- rowsum(
    expression_zeroed,
    group = symbols,
    reorder = FALSE,
    na.rm = TRUE
  )
  observed <- rowsum(
    1L * !missing,
    group = symbols,
    reorder = FALSE,
    na.rm = TRUE
  )
  averages <- sums / observed
  averages[observed == 0] <- NA_real_
  storage.mode(averages) <- "double"
  averages
}

maybe_transform_expression <- function(expression, method = "auto") {
  expression <- safe_numeric_matrix(expression, "expression matrix")
  method <- tolower(method %||% "auto")
  if (identical(method, "none")) {
    return(expression)
  }
  if (identical(method, "log2")) {
    if (any(expression < 0, na.rm = TRUE)) {
      stop("Cannot log2-transform a matrix containing negative values.", call. = FALSE)
    }
    return(log2(expression + 1))
  }
  if (identical(method, "log2_shift")) {
    minimum <- min(expression, na.rm = TRUE)
    offset <- if (minimum <= 0) 1 - minimum else 0
    log_info(
      "Applying log2(x + offset) transformation with offset ",
      sprintf("%.6g", offset), "."
    )
    return(log2(expression + offset))
  }
  if (!identical(method, "auto")) {
    stop("Unsupported expression transform: ", method, call. = FALSE)
  }

  quantiles <- stats::quantile(
    expression,
    probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE,
    names = FALSE
  )
  should_log <- quantiles[[5L]] > 100 ||
    (quantiles[[6L]] - quantiles[[1L]] > 50 && quantiles[[2L]] > 0)
  if (should_log) {
    if (quantiles[[1L]] < 0) {
      log_warn("Auto-transform detected a large range but negative values exist; no log2 applied.")
      return(expression)
    }
    log_info("Applying automatic log2(x + 1) transformation.")
    return(log2(expression + 1))
  }
  expression
}

match_any_pattern <- function(values, patterns) {
  if (length(patterns %||% character()) == 0L) {
    return(rep(FALSE, length(values)))
  }
  Reduce(
    `|`,
    lapply(
      patterns,
      function(pattern) grepl(pattern, values, ignore.case = TRUE, perl = TRUE)
    )
  )
}

derive_groups <- function(metadata, dataset) {
  candidate_columns <- c("title", "description", "source_name_ch1")
  group_column <- candidate_columns[candidate_columns %in% names(metadata)][1L]
  if (is.na(group_column)) {
    stop(dataset$id, " has no supported metadata column for grouping.", call. = FALSE)
  }

  labels <- trimws(as.character(metadata[[group_column]]))
  normal <- match_any_pattern(labels, dataset$normal_patterns)
  disease <- match_any_pattern(labels, dataset$disease_patterns)

  if (any(normal & disease)) {
    stop(
      dataset$id, " has samples matching both normal and disease patterns: ",
      paste(unique(labels[normal & disease]), collapse = "; "),
      call. = FALSE
    )
  }
  if (any(!normal & !disease)) {
    stop(
      dataset$id, " has unclassified samples: ",
      paste(unique(labels[!normal & !disease]), collapse = "; "),
      call. = FALSE
    )
  }

  group <- ifelse(disease, "Disease", "Normal")
  group <- factor(group, levels = c("Normal", "Disease"))
  names(group) <- rownames(metadata)
  attr(group, "metadata_column") <- group_column
  group
}

extract_symbol_from_annotation <- function(values) {
  values <- as.character(values)
  vapply(
    values,
    function(value) {
      if (is.na(value) || !nzchar(trimws(value))) return(NA_character_)
      first_gene <- strsplit(value, "\\s+///\\s+", perl = TRUE)[[1L]][[1L]]
      tokens <- trimws(strsplit(first_gene, "\\s+//\\s+|;", perl = TRUE)[[1L]])
      candidates <- tokens[
        grepl("^[A-Za-z][A-Za-z0-9._-]{1,30}$", tokens) &
          !grepl("^(NM_|NR_|XM_|XR_|ENSG|ENST|ILMN_|A_[0-9])", tokens) &
          !tolower(tokens) %in% c("na", "null", "unknown", "control")
      ]
      if (length(candidates) == 0L) {
        return(normalize_gene_symbols(first_gene)[[1L]])
      }
      normalize_gene_symbols(candidates[[1L]])[[1L]]
    },
    character(1)
  )
}

read_gpl_soft_mapping <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  in_table <- FALSE
  field_names <- NULL
  probe_chunks <- list()
  symbol_chunks <- list()
  control_chunks <- list()

  repeat {
    lines <- readLines(connection, n = 10000L, warn = FALSE)
    if (length(lines) == 0L) break

    if (!in_table) {
      begin <- grep("^!platform_table_begin", lines, ignore.case = TRUE)
      if (length(begin) == 0L) next
      in_table <- TRUE
      lines <- if (begin[[1L]] < length(lines)) {
        lines[(begin[[1L]] + 1L):length(lines)]
      } else {
        character()
      }
    }

    if (is.null(field_names)) {
      if (length(lines) == 0L) next
      field_names <- strsplit(lines[[1L]], "\t", fixed = TRUE)[[1L]]
      lines <- lines[-1L]
    }
    if (length(lines) == 0L) next

    end <- grep("^!platform_table_end", lines, ignore.case = TRUE)
    finished <- length(end) > 0L
    if (finished) {
      lines <- if (end[[1L]] > 1L) lines[seq_len(end[[1L]] - 1L)] else character()
    }

    if (length(lines) > 0L) {
      probe_index <- match("ID", toupper(field_names))
      gene_index <- match("GENE_SYMBOL", toupper(field_names))
      control_index <- match("CONTROL_TYPE", toupper(field_names))
      if (is.na(probe_index) || is.na(gene_index)) {
        stop(
          "GPL table must contain ID and GENE_SYMBOL columns: ", path,
          call. = FALSE
        )
      }
      parts <- strsplit(lines, "\t", fixed = TRUE)
      extract_column <- function(index, default = "") {
        vapply(
          parts,
          function(values) {
            if (is.na(index) || length(values) < index) default else values[[index]]
          },
          character(1)
        )
      }
      probe_chunks[[length(probe_chunks) + 1L]] <- extract_column(probe_index)
      symbol_chunks[[length(symbol_chunks) + 1L]] <- extract_column(gene_index)
      control_chunks[[length(control_chunks) + 1L]] <- extract_column(
        control_index,
        default = "FALSE"
      )
    }
    if (finished) break
  }

  if (is.null(field_names) || length(probe_chunks) == 0L) {
    stop("No platform table was found in: ", path, call. = FALSE)
  }
  data.frame(
    probe = unlist(probe_chunks, use.names = FALSE),
    symbol = extract_symbol_from_annotation(
      unlist(symbol_chunks, use.names = FALSE)
    ),
    control = unlist(control_chunks, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

build_gpl_mapping <- function(dataset, cache_dir) {
  cache_path <- file.path(
    ensure_dir(file.path(cache_dir, "annotations")),
    paste0(clean_filename(dataset$id), "_probe_to_symbol.rds")
  )
  if (file.exists(cache_path)) {
    log_info("Using cached annotation for ", dataset$id, ".")
    return(readRDS(cache_path))
  }

  log_info(
    "Streaming platform annotation for ", dataset$id,
    "; only the platform table is read."
  )
  mapping <- read_gpl_soft_mapping(dataset$annotation_path)
  control_values <- tolower(trimws(mapping$control))
  keep <- !control_values %in% c("true", "1", "yes", "pos", "positive")
  mapping <- mapping[, c("probe", "symbol"), drop = FALSE]
  mapping <- mapping[keep & !is.na(mapping$symbol) & nzchar(mapping$probe), ]
  mapping <- mapping[!duplicated(mapping$probe), ]
  atomic_save_rds(mapping, cache_path)
  mapping
}

annotate_expression <- function(expression, dataset, cache_dir) {
  method <- tolower(dataset$annotation_method %||% "none")
  if (identical(method, "none")) {
    return(aggregate_expression_by_gene(expression, rownames(expression)))
  }

  if (identical(method, "hgu133plus2")) {
    require_namespace("AnnotationDbi", "Affymetrix probe annotation")
    require_namespace("hgu133plus2.db", "Affymetrix GPL570 annotation")
    probes <- rownames(expression)
    mapping <- AnnotationDbi::select(
      hgu133plus2.db::hgu133plus2.db,
      keys = unique(probes),
      columns = "SYMBOL",
      keytype = "PROBEID"
    )
    mapping <- mapping[
      !is.na(mapping$SYMBOL) & !duplicated(mapping$PROBEID),
      ,
      drop = FALSE
    ]
    symbols <- mapping$SYMBOL[match(probes, mapping$PROBEID)]
    return(aggregate_expression_by_gene(expression, symbols))
  }

  if (identical(method, "gpl_soft")) {
    mapping <- build_gpl_mapping(dataset, cache_dir)
    symbols <- mapping$symbol[match(rownames(expression), mapping$probe)]
    return(aggregate_expression_by_gene(expression, symbols))
  }

  stop(
    "Unsupported annotation_method '", method, "' for ", dataset$id,
    call. = FALSE
  )
}

validate_expression_dataset <- function(object) {
  expression <- object$expression
  group <- object$group
  if (!is.matrix(expression) || nrow(expression) < 100L || ncol(expression) < 4L) {
    stop(object$id, " expression matrix has implausible dimensions.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression))) {
    stop(object$id, " still has duplicated gene symbols.", call. = FALSE)
  }
  if (!identical(colnames(expression), names(group))) {
    stop(object$id, " sample names are not aligned with group labels.", call. = FALSE)
  }
  if (any(table(group) < 2L)) {
    stop(object$id, " has fewer than two samples in a group.", call. = FALSE)
  }
  if (anyNA(expression)) {
    fraction <- mean(is.na(expression))
    if (fraction > 0.01) {
      stop(object$id, " contains more than 1% missing values.", call. = FALSE)
    }
    row_medians <- apply(expression, 1L, stats::median, na.rm = TRUE)
    indices <- which(is.na(expression), arr.ind = TRUE)
    expression[indices] <- row_medians[indices[, "row"]]
    object$expression <- expression
  }
  object
}

load_split_tables <- function(dataset) {
  require_namespace("data.table", "delimited expression matrices")
  cases <- data.table::fread(dataset$case_path, data.table = FALSE)
  controls <- data.table::fread(dataset$control_path, data.table = FALSE)
  gene_column <- dataset$gene_column %||% names(cases)[[1L]]
  if (!gene_column %in% names(cases) || !gene_column %in% names(controls)) {
    stop(dataset$id, " gene column was not found in both split tables.", call. = FALSE)
  }

  exclude_patterns <- dataset$exclude_column_patterns %||%
    c("^(Average|Mean|Median|Min|Max)(\\s|$)")
  is_summary_column <- function(columns) {
    match_any_pattern(columns, exclude_patterns)
  }
  case_columns <- setdiff(names(cases), gene_column)
  control_columns <- setdiff(names(controls), gene_column)
  case_columns <- case_columns[!is_summary_column(case_columns)]
  control_columns <- control_columns[!is_summary_column(control_columns)]
  if (length(case_columns) < 2L || length(control_columns) < 2L) {
    stop(dataset$id, " has too few sample columns after summary-column filtering.", call. = FALSE)
  }
  log_info(
    dataset$id, ": excluding non-sample columns: ",
    paste(
      c(
        setdiff(setdiff(names(cases), gene_column), case_columns),
        setdiff(setdiff(names(controls), gene_column), control_columns)
      ),
      collapse = ", "
    ),
    "."
  )

  case_genes <- normalize_gene_symbols(cases[[gene_column]])
  control_genes <- normalize_gene_symbols(controls[[gene_column]])
  case_matrix <- safe_numeric_matrix(
    cases[case_columns],
    paste(dataset$id, "case matrix")
  )
  control_matrix <- safe_numeric_matrix(
    controls[control_columns],
    paste(dataset$id, "control matrix")
  )
  rownames(case_matrix) <- case_genes
  rownames(control_matrix) <- control_genes
  case_matrix <- aggregate_expression_by_gene(case_matrix, rownames(case_matrix))
  control_matrix <- aggregate_expression_by_gene(control_matrix, rownames(control_matrix))
  genes <- intersect(rownames(case_matrix), rownames(control_matrix))
  expression <- cbind(
    case_matrix[genes, , drop = FALSE],
    control_matrix[genes, , drop = FALSE]
  )
  expression <- maybe_transform_expression(expression, dataset$transform)
  group <- factor(
    c(rep("Disease", ncol(case_matrix)), rep("Normal", ncol(control_matrix))),
    levels = c("Normal", "Disease")
  )
  names(group) <- colnames(expression)
  metadata <- data.frame(
    sample = colnames(expression),
    title = colnames(expression),
    group = as.character(group),
    row.names = colnames(expression),
    stringsAsFactors = FALSE
  )

  list(
    id = dataset$id,
    disease = dataset$disease,
    role = dataset$role,
    expression = expression,
    group = group,
    metadata = metadata
  )
}

load_geo_series <- function(dataset, cache_dir) {
  require_namespace("GEOquery", "GEO Series Matrix files")
  geo <- GEOquery::getGEO(
    filename = dataset$series_path,
    getGPL = FALSE
  )
  if (is.list(geo)) geo <- geo[[1L]]
  expression <- Biobase::exprs(geo)
  metadata <- Biobase::pData(geo)
  group <- derive_groups(metadata, dataset)
  group <- group[colnames(expression)]
  missing_fraction <- colMeans(!is.finite(expression))
  maximum_missing <- as.numeric(dataset$max_sample_missing_fraction %||% 0.05)
  usable <- missing_fraction <= maximum_missing
  dropped_samples <- data.frame(
    sample = colnames(expression)[!usable],
    missing_fraction = missing_fraction[!usable],
    group = as.character(group[!usable]),
    stringsAsFactors = FALSE
  )
  if (nrow(dropped_samples) > 0L) {
    log_warn(
      dataset$id, ": dropping ", nrow(dropped_samples),
      " samples above the missing-value threshold: ",
      paste(dropped_samples$sample, collapse = ", "), "."
    )
    expression <- expression[, usable, drop = FALSE]
    group <- droplevels(group[usable])
    metadata <- metadata[colnames(expression), , drop = FALSE]
  }
  expression <- maybe_transform_expression(expression, dataset$transform)
  expression <- annotate_expression(expression, dataset, cache_dir)
  metadata$group <- as.character(group[rownames(metadata)])

  list(
    id = dataset$id,
    disease = dataset$disease,
    role = dataset$role,
    expression = expression,
    group = group,
    metadata = metadata,
    dropped_samples = dropped_samples
  )
}

load_bulk_dataset <- function(dataset, cache_dir) {
  log_info("Loading ", dataset$id, " (", dataset$disease, ", ", dataset$role, ").")
  object <- switch(
    dataset$loader,
    split_tables = load_split_tables(dataset),
    geo_series = load_geo_series(dataset, cache_dir),
    stop("Unsupported loader: ", dataset$loader, call. = FALSE)
  )
  object <- validate_expression_dataset(object)
  log_info(
    dataset$id, ": ", nrow(object$expression), " genes x ",
    ncol(object$expression), " samples; groups ",
    paste(names(table(object$group)), table(object$group), collapse = ", "), "."
  )
  object
}

load_bulk_datasets <- function(config, role = NULL) {
  dataset_configs <- config$datasets
  if (!is.null(role)) {
    dataset_configs <- dataset_configs[
      vapply(
        dataset_configs,
        function(dataset) identical(dataset$role, role),
        logical(1)
      )
    ]
  }
  datasets <- lapply(
    dataset_configs,
    load_bulk_dataset,
    cache_dir = config$project$cache_dir
  )
  names(datasets) <- names(dataset_configs)

  qc <- do.call(
    rbind,
    lapply(
      datasets,
      function(object) {
        counts <- table(object$group)
        data.frame(
          dataset_id = object$id,
          disease = object$disease,
          role = object$role,
          genes = nrow(object$expression),
          samples = ncol(object$expression),
          normal = unname(counts["Normal"]),
          disease_samples = unname(counts["Disease"]),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  safe_write_csv(
    qc,
    file.path(
      config$project$output_dir,
      "tables",
      paste0(
        "dataset_qc_summary",
        if (is.null(role)) "" else paste0("_", role),
        ".csv"
      )
    )
  )
  dropped <- do.call(
    rbind,
    lapply(
      datasets,
      function(object) {
        table <- object$dropped_samples %||% data.frame()
        if (nrow(table) == 0L) return(NULL)
        table$dataset_id <- object$id
        table$disease <- object$disease
        table
      }
    )
  )
  if (!is.null(dropped) && nrow(dropped) > 0L) {
    safe_write_csv(
      dropped,
      file.path(
        config$project$output_dir,
        "tables",
        paste0(
          "dropped_samples",
          if (is.null(role)) "" else paste0("_", role),
          ".csv"
        )
      )
    )
  }
  datasets
}

load_all_bulk_datasets <- function(config) {
  load_bulk_datasets(config, role = NULL)
}
