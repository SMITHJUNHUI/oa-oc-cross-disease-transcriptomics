v30_advanced_base_runner <- run_reviewer_v30

v30_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probability, names = FALSE, type = 8))
}

v30_string_inputs <- function(project_root) {
  base <- file.path(project_root, "results", "api_cache", "string_v12")
  list(
    mapping = file.path(base, "shared_deg_string_ids.json"),
    physical = file.path(base, "shared_deg_physical_score700.json"),
    functional = file.path(base, "shared_deg_functional_score700.json")
  )
}

v30_read_string_json <- function(path) {
  require_namespace("jsonlite", "STRING API cache import")
  if (!file.exists(path)) {
    stop("Missing cached STRING response: ", path, call. = FALSE)
  }
  value <- jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
  if (is.data.frame(value)) {
    return(value)
  }
  if (is.list(value) && length(value) > 0L && is.data.frame(value[[1L]])) {
    return(value[[1L]])
  }
  stop("Unexpected STRING response structure: ", path, call. = FALSE)
}

v30_ppi_topology_row <- function(graph, network_type, subset_name) {
  require_namespace("igraph", "PPI topology")
  node_count <- igraph::vcount(graph)
  edge_count <- igraph::ecount(graph)
  degrees <- if (node_count > 0L) igraph::degree(graph) else numeric()
  components <- if (node_count > 0L) igraph::components(graph) else NULL
  largest <- if (is.null(components)) 0L else max(components$csize)
  data.frame(
    network_type = network_type,
    subset = subset_name,
    mapped_nodes = node_count,
    connected_nodes = sum(degrees > 0),
    edges = edge_count,
    density = if (node_count > 1L) igraph::edge_density(graph, loops = FALSE) else NA_real_,
    mean_degree = if (node_count > 0L) mean(degrees) else NA_real_,
    median_degree = if (node_count > 0L) stats::median(degrees) else NA_real_,
    components = if (is.null(components)) 0L else components$no,
    largest_component_nodes = largest,
    global_clustering = if (node_count > 2L) {
      igraph::transitivity(graph, type = "globalundirected", isolates = "zero")
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

v30_ppi_permutation <- function(graph, node_table, network_type, iterations = 10000L) {
  require_namespace("igraph", "PPI label permutation")
  edge_frame <- igraph::as_data_frame(graph, what = "edges")
  nodes <- as.character(node_table$gene)
  observed_labels <- stats::setNames(as.character(node_table$direction_binary), nodes)
  edge_a <- match(edge_frame$from, nodes)
  edge_b <- match(edge_frame$to, nodes)
  if (anyNA(c(edge_a, edge_b))) {
    stop("PPI permutation edges do not map to the audited node table.", call. = FALSE)
  }
  statistic <- function(labels) {
    concordant <- labels == "concordant"
    n_concordant <- sum(concordant)
    n_discordant <- length(labels) - n_concordant
    within_concordant <- sum(concordant[edge_a] & concordant[edge_b])
    within_discordant <- sum(!concordant[edge_a] & !concordant[edge_b])
    density_concordant <- if (n_concordant > 1L) {
      2 * within_concordant / (n_concordant * (n_concordant - 1))
    } else {
      NA_real_
    }
    density_discordant <- if (n_discordant > 1L) {
      2 * within_discordant / (n_discordant * (n_discordant - 1))
    } else {
      NA_real_
    }
    c(
      concordant_edges = within_concordant,
      discordant_edges = within_discordant,
      density_difference = density_concordant - density_discordant
    )
  }
  observed <- statistic(observed_labels)
  set.seed(20260726)
  null <- replicate(iterations, statistic(sample(observed_labels, replace = FALSE)))
  null <- as.data.frame(t(null))
  null$iteration <- seq_len(nrow(null))
  null$network_type <- network_type
  summary <- do.call(rbind, lapply(names(observed), function(metric) {
    values <- null[[metric]]
    center <- mean(values)
    spread <- stats::sd(values)
    data.frame(
      network_type = network_type,
      metric = metric,
      observed = unname(observed[[metric]]),
      null_mean = center,
      null_sd = spread,
      z_score = if (is.finite(spread) && spread > 0) {
        (unname(observed[[metric]]) - center) / spread
      } else {
        NA_real_
      },
      empirical_p_two_sided = (
        1 + sum(abs(values - center) >= abs(unname(observed[[metric]]) - center))
      ) / (length(values) + 1),
      permutations = iterations,
      stringsAsFactors = FALSE
    )
  }))
  list(summary = summary, null = null)
}

v30_prepare_ppi <- function(project_root, paths) {
  require_namespace("igraph", "direction-aware PPI analysis")
  inputs <- v30_string_inputs(project_root)
  if (any(!vapply(inputs, file.exists, logical(1)))) {
    stop(
      "STRING v12.0 cache is incomplete. Refresh the audited API responses first.",
      call. = FALSE
    )
  }
  shared_path <- file.path(
    paths$tables,
    "Table_S2_shared_differentially_expressed_genes.csv"
  )
  shared <- utils::read.csv(shared_path, check.names = FALSE)
  shared$gene <- toupper(as.character(shared$gene))
  shared$direction_quadrant <- ifelse(
    shared$logFC_OA > 0 & shared$logFC_OC > 0,
    "OA higher / OC higher",
    ifelse(
      shared$logFC_OA < 0 & shared$logFC_OC < 0,
      "OA lower / OC lower",
      ifelse(
        shared$logFC_OA > 0 & shared$logFC_OC < 0,
        "OA higher / OC lower",
        "OA lower / OC higher"
      )
    )
  )
  shared$direction_binary <- ifelse(
    shared$directionally_concordant,
    "concordant",
    "discordant"
  )
  candidates <- utils::read.csv(file.path(
    paths$tables,
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  candidate_genes <- toupper(as.character(
    candidates$gene[as.logical(candidates$fixed_score_member)]
  ))

  mapping <- v30_read_string_json(inputs$mapping)
  mapping$query_gene <- toupper(as.character(mapping$queryItem))
  mapping$preferredName <- toupper(as.character(mapping$preferredName))
  mapping$mapping_status <- ifelse(
    mapping$query_gene == mapping$preferredName,
    "exact symbol",
    "alias resolved"
  )
  missing <- setdiff(shared$gene, mapping$query_gene)
  if (length(missing) > 0L) {
    mapping <- rbind(
      mapping,
      data.frame(
        queryIndex = NA_integer_,
        queryItem = missing,
        stringId = NA_character_,
        ncbiTaxonId = 9606,
        taxonName = "Homo sapiens",
        preferredName = NA_character_,
        annotation = NA_character_,
        query_gene = missing,
        mapping_status = "unmapped",
        stringsAsFactors = FALSE
      )
    )
  }
  mapping <- merge(
    mapping,
    shared[, c(
      "gene", "logFC_OA", "logFC_OC", "direction_quadrant", "direction_binary"
    )],
    by.x = "query_gene",
    by.y = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  mapping$candidate <- mapping$query_gene %in% candidate_genes

  build_network <- function(path, network_type) {
    edge_raw <- v30_read_string_json(path)
    id_to_gene <- stats::setNames(mapping$query_gene, mapping$stringId)
    edge_raw$gene_A <- unname(id_to_gene[as.character(edge_raw$stringId_A)])
    edge_raw$gene_B <- unname(id_to_gene[as.character(edge_raw$stringId_B)])
    edge_raw <- edge_raw[
      !is.na(edge_raw$gene_A) &
        !is.na(edge_raw$gene_B) &
        edge_raw$gene_A != edge_raw$gene_B,
      ,
      drop = FALSE
    ]
    edge_raw$pair_a <- pmin(edge_raw$gene_A, edge_raw$gene_B)
    edge_raw$pair_b <- pmax(edge_raw$gene_A, edge_raw$gene_B)
    edge_raw <- edge_raw[
      order(edge_raw$pair_a, edge_raw$pair_b, -edge_raw$score),
      ,
      drop = FALSE
    ]
    edge_raw <- edge_raw[
      !duplicated(paste(edge_raw$pair_a, edge_raw$pair_b, sep = "||")),
      ,
      drop = FALSE
    ]
    edge_raw$network_type <- network_type
    edge_raw$required_score <- 700L
    edge_raw$string_version <- "12.0"
    vertices <- mapping[
      mapping$mapping_status != "unmapped",
      c(
        "query_gene", "stringId", "preferredName", "mapping_status",
        "logFC_OA", "logFC_OC", "direction_quadrant", "direction_binary",
        "candidate"
      ),
      drop = FALSE
    ]
    names(vertices)[names(vertices) == "query_gene"] <- "gene"
    graph <- igraph::graph_from_data_frame(
      edge_raw[, c("pair_a", "pair_b", "score"), drop = FALSE],
      directed = FALSE,
      vertices = vertices
    )
    igraph::E(graph)$score <- as.numeric(igraph::E(graph)$score)
    degree <- igraph::degree(graph)
    betweenness <- igraph::betweenness(graph, normalized = TRUE, weights = NULL)
    clustering <- igraph::transitivity(
      graph,
      type = "localundirected",
      isolates = "zero"
    )
    node_metrics <- data.frame(
      gene = igraph::V(graph)$name,
      network_type = network_type,
      degree = as.numeric(degree),
      betweenness = as.numeric(betweenness),
      clustering_coefficient = as.numeric(clustering),
      stringsAsFactors = FALSE
    )
    node_metrics <- merge(node_metrics, vertices, by = "gene", sort = FALSE)
    topology <- rbind(
      v30_ppi_topology_row(graph, network_type, "all_mapped_shared_DEGs"),
      v30_ppi_topology_row(
        igraph::induced_subgraph(
          graph,
          vids = vertices$gene[vertices$direction_binary == "concordant"]
        ),
        network_type,
        "concordant"
      ),
      v30_ppi_topology_row(
        igraph::induced_subgraph(
          graph,
          vids = vertices$gene[vertices$direction_binary == "discordant"]
        ),
        network_type,
        "discordant"
      )
    )
    permutation <- v30_ppi_permutation(
      graph,
      vertices[, c("gene", "direction_binary")],
      network_type
    )
    list(
      graph = graph,
      edges = edge_raw,
      nodes = node_metrics,
      topology = topology,
      permutation = permutation
    )
  }

  physical <- build_network(inputs$physical, "high-confidence physical")
  functional <- build_network(inputs$functional, "high-confidence functional")
  topology <- rbind(physical$topology, functional$topology)
  permutation <- rbind(
    physical$permutation$summary,
    functional$permutation$summary
  )
  mapping$accessed_on <- "2026-08-01"
  mapping$string_version <- "12.0"
  mapping$species <- "Homo sapiens (NCBI taxon 9606)"
  mapping$inference_boundary <- paste0(
    "STRING mapping/association evidence; absence of an edge is not evidence of ",
    "no biological interaction."
  )
  physical$edges$inference_boundary <- paste0(
    "High-confidence STRING physical association among submitted shared DEG products; ",
    "not direct experimental validation in OA or OC."
  )
  functional$edges$inference_boundary <- paste0(
    "High-confidence STRING functional association sensitivity network; evidence ",
    "channels may include indirect and predicted associations."
  )
  physical$nodes$hub_definition <- "Descriptive within-network topology; not used for candidate selection"
  functional$nodes$hub_definition <- "Sensitivity topology; not used for candidate selection"
  topology$comparison_boundary <- paste0(
    "Descriptive induced-subgraph metrics include mapped isolates. Group-size effects ",
    "are addressed by fixed-size label permutation."
  )
  permutation$interpretation_boundary <- paste0(
    "Permutation tests direction-label organization within the fixed STRING graph; ",
    "it does not test disease mechanism or causal protein interaction."
  )

  safe_write_csv(mapping, file.path(paths$tables, "Table_S25a_STRING_mapping_audit.csv"))
  safe_write_csv(
    rbind(physical$edges, functional$edges),
    file.path(paths$tables, "Table_S25b_direction_aware_STRING_edges.csv")
  )
  safe_write_csv(
    rbind(physical$nodes, functional$nodes),
    file.path(paths$tables, "Table_S25c_STRING_node_topology.csv")
  )
  safe_write_csv(
    topology,
    file.path(paths$tables, "Table_S25d_STRING_network_topology.csv")
  )
  safe_write_csv(
    permutation,
    file.path(paths$tables, "Table_S25e_STRING_direction_label_permutation.csv")
  )
  safe_write_csv(
    physical$permutation$null,
    file.path(paths$source, "SupplementaryFigure12_physical_label_permutation_null.csv")
  )
  safe_write_csv(mapping, file.path(paths$source, "Figure2_STRING_mapping.csv"))
  safe_write_csv(physical$edges, file.path(paths$source, "Figure2_STRING_physical_edges.csv"))
  safe_write_csv(physical$nodes, file.path(paths$source, "Figure2_STRING_physical_nodes.csv"))
  list(
    mapping = mapping,
    physical = physical,
    functional = functional,
    topology = topology,
    permutation = permutation
  )
}

v30_ppi_plot <- function(ppi, compact = TRUE) {
  require_namespace("igraph", "PPI network layout")
  graph <- ppi$physical$graph
  graph <- igraph::induced_subgraph(graph, vids = names(igraph::degree(graph))[igraph::degree(graph) > 0])
  set.seed(20260726)
  coordinates <- igraph::layout_with_fr(
    graph,
    weights = igraph::E(graph)$score,
    niter = 1500,
    grid = "nogrid"
  )
  nodes <- ppi$physical$nodes[
    match(igraph::V(graph)$name, ppi$physical$nodes$gene),
    ,
    drop = FALSE
  ]
  nodes$x <- coordinates[, 1L]
  nodes$y <- coordinates[, 2L]
  edge_frame <- igraph::as_data_frame(graph, what = "edges")
  edge_frame$x <- nodes$x[match(edge_frame$from, nodes$gene)]
  edge_frame$y <- nodes$y[match(edge_frame$from, nodes$gene)]
  edge_frame$xend <- nodes$x[match(edge_frame$to, nodes$gene)]
  edge_frame$yend <- nodes$y[match(edge_frame$to, nodes$gene)]
  top_hubs <- head(nodes$gene[order(-nodes$degree, nodes$gene)], if (compact) 7L else 12L)
  labels <- nodes[nodes$gene %in% unique(c(top_hubs, nodes$gene[nodes$candidate])), , drop = FALSE]
  colors <- c(concordant = "#009E73", discordant = "#D55E00")
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_frame,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, alpha = score),
      linewidth = if (compact) 0.35 else 0.48,
      colour = "#98A2B3"
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(
        x = x,
        y = y,
        fill = direction_binary,
        size = degree,
        shape = candidate
      ),
      colour = "#344054",
      stroke = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(x = x, y = y, label = gene),
      size = if (compact) 1.9 else 2.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.16,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23)) +
    ggplot2::scale_size_continuous(range = if (compact) c(1.2, 4.2) else c(1.7, 6)) +
    ggplot2::scale_alpha_continuous(range = c(0.25, 0.8), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = if (compact) "Physical STRING context" else "Direction-aware physical STRING network",
      subtitle = paste0(
        igraph::vcount(graph), " connected products; ", igraph::ecount(graph),
        if (compact) " edges (score >=0.700)" else " STRING v12.0 edges (score >=0.700)"
      ),
      fill = NULL,
      shape = "Ten-gene set",
      size = "Degree",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::guides(
      fill = if (compact) "none" else ggplot2::guide_legend(order = 1),
      shape = if (compact) "none" else ggplot2::guide_legend(order = 2),
      size = if (compact) "none" else ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = if (compact) 8 else 10),
      plot.subtitle = ggplot2::element_text(size = if (compact) 6 else 7),
      legend.position = if (compact) "none" else "bottom",
      legend.text = ggplot2::element_text(size = 5.5),
      legend.title = ggplot2::element_text(size = 5.5),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
  plot
}

v30_build_figure2_with_ppi <- function(project_root, paths, ppi) {
  common <- utils::read.csv(file.path(paths$source, "Figure2_common_gene_effects_quadrants.csv"))
  oa <- utils::read.csv(file.path(paths$source, "Figure2_OA_DEG.csv"))
  oc <- utils::read.csv(file.path(paths$source, "Figure2_OC_DEG.csv"))
  hub_genes <- submission_load_cache(project_root, "07_machine_learning.rds")$final_genes
  p1 <- submission_volcano_plot(oa, "OA", hub_genes)
  p2 <- submission_volcano_plot(oc, "OC", hub_genes)
  quadrant_levels <- c(
    "OA higher / OC higher", "OA lower / OC lower",
    "OA higher / OC lower", "OA lower / OC higher"
  )
  common$direction_quadrant <- factor(
    common$direction_quadrant,
    levels = c("Not primary shared", quadrant_levels)
  )
  shared_only <- common[common$primary_shared, , drop = FALSE]
  candidates <- common[common$hub, , drop = FALSE]
  quadrant_colors <- c(
    "OA higher / OC higher" = "#009E73",
    "OA lower / OC lower" = "#6A3D9A",
    "OA higher / OC lower" = "#0072B2",
    "OA lower / OC higher" = "#D55E00"
  )
  p3_base <- ggplot2::ggplot(common, ggplot2::aes(x = logFC_OA, y = logFC_OC)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.3) +
    ggplot2::geom_point(colour = "#D4D7DB", size = 0.31, alpha = 0.30) +
    ggplot2::geom_point(
      data = shared_only,
      ggplot2::aes(colour = direction_quadrant),
      size = 0.82,
      alpha = 0.72
    ) +
    ggplot2::scale_colour_manual(
      values = quadrant_colors,
      guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::labs(
      title = "Shared membership does not imply shared direction",
      subtitle = "146/286 concordant and 140/286 discordant",
      x = "OA log2 fold change",
      y = "OC log2 fold change",
      colour = NULL
    ) +
    submission_theme(7.3) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 5.4),
      legend.key.width = grid::unit(7, "mm")
    )
  inset <- ggplot2::ggplot(
    candidates,
    ggplot2::aes(x = logFC_OA, y = logFC_OC, colour = direction_quadrant)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#9CA3AF", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#9CA3AF", linewidth = 0.25) +
    ggplot2::geom_point(size = 1.8) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = gene),
      size = 2,
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = quadrant_colors) +
    ggplot2::labs(title = "Ten-gene evidence summary", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 6.2, base_family = "Arial") +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", colour = "#6B7280"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 6.8)
    )
  p3 <- p3_base + patchwork::inset_element(
    inset,
    left = 0.50,
    bottom = 0.04,
    right = 0.99,
    top = 0.55
  )
  p4 <- v30_ppi_plot(ppi, compact = TRUE)
  figure <- submission_panel_tag((p1 | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "Figure2_gene_direction",
    paths$figures,
    width_mm = 185,
    height_mm = 165
  )
}

v30_build_ppi_supplement <- function(paths, ppi) {
  network <- v30_ppi_plot(ppi, compact = FALSE) +
    ggplot2::theme(legend.position = "none")
  nodes <- ppi$physical$nodes
  nodes$label <- ifelse(nodes$degree >= v30_quantile(nodes$degree, 0.95), nodes$gene, "")
  p2 <- ggplot2::ggplot(
    nodes,
    ggplot2::aes(x = degree, y = betweenness, colour = direction_binary, shape = candidate)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggrepel::geom_text_repel(
      data = nodes[nzchar(nodes$label), , drop = FALSE],
      ggplot2::aes(label = label),
      size = 2.1,
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(concordant = "#009E73", discordant = "#D55E00"),
      labels = c(concordant = "Concordant", discordant = "Discordant")
    ) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17),
      labels = c(`FALSE` = "Other", `TRUE` = "Ten-gene member")
    ) +
    ggplot2::labs(
      title = "Descriptive node topology",
      subtitle = "Topology did not determine the ten-gene set",
      x = "Degree",
      y = "Normalized betweenness",
      colour = NULL,
      shape = "Ten-gene set"
    ) +
    submission_theme(7) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1),
      shape = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )
  topology <- ppi$topology[
    ppi$topology$subset %in% c("concordant", "discordant"),
    ,
    drop = FALSE
  ]
  topology$network_type <- factor(
    topology$network_type,
    levels = c("high-confidence physical", "high-confidence functional")
  )
  p3 <- ggplot2::ggplot(
    topology,
    ggplot2::aes(x = subset, y = density, fill = network_type)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75), width = 0.66, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c(
      `high-confidence physical` = "#4E79A7",
      `high-confidence functional` = "#F28E2B"
    )) +
    ggplot2::labs(
      title = "Induced-subgraph density",
      subtitle = "Mapped isolates retained",
      x = NULL,
      y = "Edge density",
      fill = NULL
    ) +
    submission_theme(7) +
    ggplot2::theme(legend.position = "top")
  null <- ppi$physical$permutation$null
  observed <- ppi$physical$permutation$summary
  observed <- observed$observed[observed$metric == "density_difference"]
  p4 <- ggplot2::ggplot(null, ggplot2::aes(x = density_difference)) +
    ggplot2::geom_histogram(bins = 40, fill = "#DCE6F1", colour = "white") +
    ggplot2::geom_vline(xintercept = observed, colour = "#C00000", linewidth = 0.8) +
    ggplot2::labs(
      title = "Direction-label permutation",
      subtitle = "10,000 fixed-size permutations; red = observed",
      x = "Concordant density minus discordant density",
      y = "Permutations"
    ) +
    submission_theme(7)
  figure <- submission_panel_tag((network | p2) / (p3 | p4))
  submission_save_plot(
    figure,
    "SupplementaryFigure12_direction_aware_STRING_network",
    paths$figures,
    width_mm = 185,
    height_mm = 185
  )
}

v30_candidate_definition <- function(paths) {
  table <- utils::read.csv(file.path(
    paths$tables,
    "Table_S16_candidate_prioritization_matrix.csv"
  ))
  table <- table[as.logical(table$fixed_score_member), , drop = FALSE]
  table$gene <- toupper(as.character(table$gene))
  table
}

v30_ucell_features <- function(candidate_table) {
  signed <- function(effect) {
    ifelse(
      effect >= 0,
      paste0(candidate_table$gene, "+"),
      paste0(candidate_table$gene, "-")
    )
  }
  list(
    candidate_abundance = candidate_table$gene,
    OA_direction = signed(as.numeric(candidate_table$log2FC_OA)),
    OC_direction = signed(as.numeric(candidate_table$log2FC_OC))
  )
}

v30_candidate_count_summary <- function(counts, metadata, genes, dataset_id, disease, scope) {
  present <- intersect(genes, rownames(counts))
  groups <- unique(metadata[, c("context", "cell_type"), drop = FALSE])
  rows <- vector("list", nrow(groups))
  for (index in seq_len(nrow(groups))) {
    selected <- metadata$context == groups$context[[index]] &
      metadata$cell_type == groups$cell_type[[index]]
    selected[is.na(selected)] <- FALSE
    values <- matrix(0, nrow = length(genes), ncol = sum(selected))
    rownames(values) <- genes
    if (length(present) > 0L && sum(selected) > 0L) {
      values[match(present, genes), ] <- as.matrix(counts[present, selected, drop = FALSE])
    }
    rows[[index]] <- data.frame(
      dataset_id = dataset_id,
      disease = disease,
      context = groups$context[[index]],
      cell_type = groups$cell_type[[index]],
      gene = genes,
      cells = sum(selected),
      samples = length(unique(metadata$sample[selected])),
      expression_metric = "mean_umi_per_cell",
      expression_value = if (sum(selected) > 0L) rowSums(values) / sum(selected) else NA_real_,
      fraction_detected = if (sum(selected) > 0L) rowSums(values > 0) / sum(selected) else NA_real_,
      sampling_scope = scope,
      genes_present = length(present),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

v30_score_ucell_chunk <- function(
    counts,
    metadata,
    candidate_table,
    dataset_id,
    disease,
    scope
) {
  require_namespace("UCell", "rank-based candidate-set scoring")
  require_namespace("BiocParallel", "serial UCell scoring")
  genes <- candidate_table$gene
  if (!identical(colnames(counts), metadata$cell_id)) {
    index <- match(colnames(counts), metadata$cell_id)
    if (anyNA(index)) {
      stop("UCell metadata do not cover all count columns in ", dataset_id, call. = FALSE)
    }
    metadata <- metadata[index, , drop = FALSE]
  }
  signatures <- v30_ucell_features(candidate_table)
  score <- UCell::ScoreSignatures_UCell(
    matrix = counts,
    features = signatures,
    maxRank = 1500,
    chunk.size = 500,
    missing_genes = "impute",
    BPPARAM = BiocParallel::SerialParam(),
    force.gc = TRUE
  )
  score <- as.data.frame(score)
  score$cell_id <- rownames(score)
  names(score) <- sub("_UCell$", "", names(score))
  metadata <- metadata[match(score$cell_id, metadata$cell_id), , drop = FALSE]
  result <- data.frame(
    cell_id = score$cell_id,
    dataset_id = dataset_id,
    disease = disease,
    sample = metadata$sample,
    context = metadata$context,
    cell_type = metadata$cell_type,
    sampling_scope = scope,
    candidate_abundance_UCell = score$candidate_abundance,
    OA_direction_UCell = score$OA_direction,
    OC_direction_UCell = score$OC_direction,
    stringsAsFactors = FALSE
  )
  result$matched_direction_UCell <- ifelse(
    disease == "OA",
    result$OA_direction_UCell,
    result$OC_direction_UCell
  )
  evidence <- v30_candidate_count_summary(
    counts,
    metadata,
    genes,
    dataset_id,
    disease,
    scope
  )
  list(scores = result, evidence = evidence)
}

v30_ucell_cache_signature <- function(project_root, paths) {
  files <- c(
    file.path(paths$tables, "Table_S16_candidate_prioritization_matrix.csv"),
    file.path(project_root, "results", "single_cell_downstream", "GSE104782", "cell_annotations.tsv.gz"),
    file.path(project_root, "results", "single_cell_downstream", "GSE154600", "cell_annotations.tsv.gz"),
    file.path(project_root, "results", "single_cell_downstream", "GSE169454", "cell_annotations.tsv.gz"),
    file.path(project_root, "results", "single_cell_downstream", "GSE180661", "annotation_reference_sampling.csv"),
    file.path(project_root, "results", "single_cell_downstream", "GSE255460", "cell_annotations_all_QC_pass.tsv.gz"),
    list.files(
      file.path(project_root, "results", "cache", "single_cell", "GSE154600"),
      pattern = "_qc_sce\\.rds$",
      full.names = TRUE
    ),
    list.files(
      file.path(project_root, "results", "cache", "single_cell", "GSE169454"),
      pattern = "_qc_sce\\.rds$",
      full.names = TRUE
    ),
    file.path(
      project_root,
      "results", "cache", "single_cell", "GSE104782", "GSE104782_umi_qc_sce.rds"
    ),
    file.path(
      project_root,
      "results", "cache", "single_cell_downstream", "GSE180661",
      "balanced_reference_sce.rds"
    ),
    file.path(
      project_root,
      "results", "cache", "single_cell", "GSE255460", "partitioned_csr", "manifest.json"
    )
  )
  files <- sort(unique(files[file.exists(files)]))
  info <- file.info(files)
  digest::digest(list(
    files = data.frame(
      path = normalizePath(files, winslash = "/", mustWork = TRUE),
      size = info$size,
      mtime = as.numeric(info$mtime),
      stringsAsFactors = FALSE
    ),
    UCell = as.character(utils::packageVersion("UCell")),
    maxRank = 1500L,
    score_design = "candidate_abundance_plus_disease_direction_v1"
  ))
}

v30_prepare_ucell <- function(project_root, paths) {
  require_namespace("data.table", "single-cell score tables")
  require_namespace("UCell", "rank-based candidate-set scoring")
  local_config <- file.path(project_root, "config", "local.yml")
  config <- read_project_config(
    project_root,
    if (file.exists(local_config)) "config/local.yml" else "config/config.yml"
  )
  candidate_table <- v30_candidate_definition(paths)
  cache_dir <- ensure_dir(file.path(project_root, "results", "cache", "submission_v30"))
  cache_path <- file.path(cache_dir, "candidate_ucell_scores.rds")
  signature <- v30_ucell_cache_signature(project_root, paths)
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$signature, signature)) {
      return(cached$value)
    }
  }
  score_parts <- list()
  evidence_parts <- list()
  append_result <- function(value) {
    score_parts[[length(score_parts) + 1L]] <<- value$scores
    evidence_parts[[length(evidence_parts) + 1L]] <<- value$evidence
  }
  annotations <- list(
    GSE104782 = data.table::fread(file.path(
      project_root, "results", "single_cell_downstream", "GSE104782",
      "cell_annotations.tsv.gz"
    ), showProgress = FALSE),
    GSE154600 = data.table::fread(file.path(
      project_root, "results", "single_cell_downstream", "GSE154600",
      "cell_annotations.tsv.gz"
    ), showProgress = FALSE),
    GSE169454 = data.table::fread(file.path(
      project_root, "results", "single_cell_downstream", "GSE169454",
      "cell_annotations.tsv.gz"
    ), showProgress = FALSE),
    GSE255460 = data.table::fread(file.path(
      project_root, "results", "single_cell_downstream", "GSE255460",
      "cell_annotations_all_QC_pass.tsv.gz"
    ), showProgress = FALSE)
  )

  log_info("V3.0 UCell: GSE104782.")
  sce <- readRDS(file.path(
    project_root,
    "results", "cache", "single_cell", "GSE104782", "GSE104782_umi_qc_sce.rds"
  ))
  pass <- as.logical(SummarizedExperiment::colData(sce)$passes_QC)
  pass[is.na(pass)] <- FALSE
  counts <- SummarizedExperiment::assay(sce, "counts")[, pass, drop = FALSE]
  cell_ids <- colnames(counts)
  annotation <- annotations$GSE104782[match(cell_ids, cell_id)]
  metadata <- data.frame(
    cell_id = cell_ids,
    sample = annotation$donor,
    context = "OA_all",
    cell_type = annotation$published_cell_type,
    stringsAsFactors = FALSE
  )
  append_result(v30_score_ucell_chunk(
    counts, metadata, candidate_table, "GSE104782", "OA", "all QC-pass OA cells"
  ))
  rm(sce, counts)
  invisible(gc())

  gse154_paths <- sort(list.files(
    file.path(project_root, "results", "cache", "single_cell", "GSE154600"),
    pattern = "_qc_sce\\.rds$",
    full.names = TRUE
  ))
  for (index in seq_along(gse154_paths)) {
    log_info("V3.0 UCell: GSE154600 sample ", index, "/", length(gse154_paths), ".")
    sce <- readRDS(gse154_paths[[index]])
    pass <- as.logical(SummarizedExperiment::colData(sce)$passes_QC)
    pass[is.na(pass)] <- FALSE
    raw_counts <- SummarizedExperiment::assay(sce, "counts")[, pass, drop = FALSE]
    counts <- .scd_collapse_gene_symbols(
      raw_counts,
      SummarizedExperiment::rowData(sce)$gene_symbol,
      SummarizedExperiment::rowData(sce)$gene_id
    )
    cell_ids <- colnames(counts)
    annotation <- annotations$GSE154600[match(cell_ids, cell_id)]
    metadata <- data.frame(
      cell_id = cell_ids,
      sample = annotation$batch,
      context = "HGSOC_all",
      cell_type = annotation$cell_type,
      stringsAsFactors = FALSE
    )
    append_result(v30_score_ucell_chunk(
      counts, metadata, candidate_table, "GSE154600", "OC", "all QC-pass tumor cells"
    ))
    rm(sce, raw_counts, counts)
    invisible(gc())
  }

  gse169_paths <- sort(list.files(
    file.path(project_root, "results", "cache", "single_cell", "GSE169454"),
    pattern = "_oa[1-4]_qc_sce\\.rds$",
    full.names = TRUE,
    ignore.case = TRUE
  ))
  for (index in seq_along(gse169_paths)) {
    log_info("V3.0 UCell: GSE169454 OA sample ", index, "/", length(gse169_paths), ".")
    sce <- readRDS(gse169_paths[[index]])
    pass <- as.logical(SummarizedExperiment::colData(sce)$passes_QC)
    pass[is.na(pass)] <- FALSE
    raw_counts <- SummarizedExperiment::assay(sce, "counts")[, pass, drop = FALSE]
    counts <- .scd_collapse_gene_symbols(
      raw_counts,
      SummarizedExperiment::rowData(sce)$gene_symbol,
      SummarizedExperiment::rowData(sce)$gene_id
    )
    cell_ids <- colnames(counts)
    annotation <- annotations$GSE169454[match(cell_ids, cell_id)]
    metadata <- data.frame(
      cell_id = cell_ids,
      sample = annotation$batch,
      context = "OA_only",
      cell_type = annotation$cell_type,
      stringsAsFactors = FALSE
    )
    append_result(v30_score_ucell_chunk(
      counts, metadata, candidate_table, "GSE169454", "OA", "four OA samples only"
    ))
    rm(sce, raw_counts, counts)
    invisible(gc())
  }

  log_info("V3.0 UCell: GSE180661 balanced reference sample.")
  sce <- readRDS(file.path(
    project_root,
    "results", "cache", "single_cell_downstream", "GSE180661",
    "balanced_reference_sce.rds"
  ))
  counts <- SummarizedExperiment::assay(sce, "counts")
  metadata <- data.frame(
    cell_id = colnames(counts),
    sample = as.character(SummarizedExperiment::colData(sce)$patient_id),
    context = "balanced_reference",
    cell_type = as.character(SummarizedExperiment::colData(sce)$reference_label),
    stringsAsFactors = FALSE
  )
  append_result(v30_score_ucell_chunk(
    counts, metadata, candidate_table, "GSE180661", "OC",
    "deterministic patient-balanced reference sample"
  ))
  rm(sce, counts)
  invisible(gc())

  log_info("V3.0 UCell: GSE255460 OA partitions.")
  dataset <- .scd_dataset_config(config, "GSE255460")
  source_metadata <- .oa_sc_read_gse255460_metadata(dataset$metadata_path)
  bundle <- .sc_gse255460_ensure_sparse_bundle(dataset, config)
  validated <- .sc_gse255460_validate_manifest(bundle, source_metadata)
  oa_annotation <- annotations$GSE255460[trait == "OA"]
  partitions <- sort(unique(as.character(oa_annotation$ID)))
  for (index in seq_along(partitions)) {
    partition_id <- partitions[[index]]
    log_info("V3.0 UCell: GSE255460 partition ", index, "/", length(partitions), ".")
    imported <- .sc_read_gse255460_partition(
      bundle, validated, partition_id, source_metadata
    )
    part_annotation <- oa_annotation[ID == partition_id]
    cell_index <- match(part_annotation$cell_id, colnames(imported$counts))
    if (anyNA(cell_index)) {
      stop("GSE255460 UCell cells do not map to partition ", partition_id, call. = FALSE)
    }
    counts <- imported$counts[, cell_index, drop = FALSE]
    metadata <- data.frame(
      cell_id = colnames(counts),
      sample = part_annotation$donor,
      context = paste0("OA:", part_annotation$group),
      cell_type = part_annotation$celltype,
      stringsAsFactors = FALSE
    )
    append_result(v30_score_ucell_chunk(
      counts, metadata, candidate_table, "GSE255460", "OA", "all QC-pass OA cells"
    ))
    rm(imported, counts)
    invisible(gc())
  }
  value <- list(
    scores = data.table::rbindlist(score_parts, use.names = TRUE, fill = TRUE),
    evidence = data.table::rbindlist(evidence_parts, use.names = TRUE, fill = TRUE),
    parameters = data.frame(
      method = "UCell",
      version = as.character(utils::packageVersion("UCell")),
      maxRank = 1500L,
      candidate_genes = paste(candidate_table$gene, collapse = ";"),
      score_interpretation = paste0(
        "Rank-based candidate-set expression score and disease-direction compatibility; ",
        "not pathway activity or mechanism."
      ),
      stringsAsFactors = FALSE
    )
  )
  atomic_save_rds(list(signature = signature, value = value), cache_path)
  value
}

v30_compute_ccss <- function(evidence, thresholds = c(50L, 100L, 200L)) {
  require_namespace("data.table", "candidate cell-context specificity")
  evidence <- data.table::as.data.table(evidence)
  evidence[, eligible_label := !cell_type %in% c("Unassigned", "Ambiguous", "Other")]
  rows <- list()
  for (threshold in thresholds) {
    part <- data.table::copy(evidence)
    part[, minimum_cells := threshold]
    part[, eligible_primary := eligible_label & cells >= threshold]
    groups <- unique(part[, .(dataset_id, disease, context, gene)])
    for (index in seq_len(nrow(groups))) {
      key <- groups[index]
      selected <- part[
        dataset_id == key$dataset_id & disease == key$disease &
          context == key$context & gene == key$gene
      ]
      eligible <- selected$eligible_primary
      cell_type_count <- sum(eligible)
      selected[, eligible_cell_types := cell_type_count]
      selected[, relative_detection := NA_real_]
      selected[, tau_specificity := NA_real_]
      selected[, candidate_cell_specificity_score := NA_real_]
      selected[, score_status := "insufficient eligible cell types"]
      if (cell_type_count >= 2L) {
        maximum <- max(selected$fraction_detected[eligible], na.rm = TRUE)
        if (is.finite(maximum) && maximum > 0) {
          relative <- selected$fraction_detected[eligible] / maximum
          tau <- sum(1 - relative) / (cell_type_count - 1)
          selected$relative_detection[eligible] <- relative
          selected$tau_specificity[eligible] <- tau
          selected$candidate_cell_specificity_score[eligible] <- tau * relative
          selected$score_status[eligible] <- "estimated"
        } else {
          selected$relative_detection[eligible] <- 0
          selected$tau_specificity[eligible] <- 0
          selected$candidate_cell_specificity_score[eligible] <- 0
          selected$score_status[eligible] <- "not detected"
        }
      }
      rows[[length(rows) + 1L]] <- selected
    }
  }
  context_scores <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  eligible_scores <- context_scores[
    eligible_primary == TRUE & is.finite(candidate_cell_specificity_score)
  ]
  dataset_scores <- eligible_scores[, .(
    median_CCSS = stats::median(candidate_cell_specificity_score),
    median_fraction_detected = stats::median(fraction_detected),
    contexts_contributing = data.table::uniqueN(context),
    context_values = paste(sort(unique(context)), collapse = ";"),
    cells_min = min(cells),
    cells_max = max(cells)
  ), by = .(dataset_id, disease, gene, cell_type, minimum_cells)]
  consensus <- dataset_scores[, .(
    median_CCSS = stats::median(median_CCSS),
    minimum_CCSS = min(median_CCSS),
    maximum_CCSS = max(median_CCSS),
    median_fraction_detected = stats::median(median_fraction_detected),
    n_datasets_contributing = data.table::uniqueN(dataset_id),
    datasets = paste(sort(unique(dataset_id)), collapse = ";"),
    evidence_status = ifelse(
      data.table::uniqueN(dataset_id) >= 2L,
      "replicated exact label",
      "dataset-specific exploratory"
    )
  ), by = .(disease, gene, cell_type, minimum_cells)]
  context_scores$score_definition <- paste0(
    "r = detection fraction / maximum eligible detection fraction; ",
    "tau = sum(1-r)/(K-1); CCSS = tau*r."
  )
  context_scores$inference_boundary <- paste0(
    "Dataset/context-specific descriptive localization score; not differential ",
    "expression, gene-set activity, cell function, or causality."
  )
  consensus$inference_boundary <- paste0(
    "Median across datasets only for identical source labels; missing labels remain ",
    "missing and OA/OC values are not integrated."
  )
  list(
    context = as.data.frame(context_scores),
    dataset = as.data.frame(dataset_scores),
    consensus = as.data.frame(consensus)
  )
}

v30_summarize_ucell <- function(scores) {
  require_namespace("data.table", "UCell sample-aware summaries")
  scores <- data.table::as.data.table(scores)
  scores[, eligible_label := !cell_type %in% c("Unassigned", "Ambiguous", "Other")]
  cell_counts <- scores[, .(cells = .N), by = .(
    dataset_id, disease, sampling_scope, sample, cell_type, eligible_label
  )]
  sample_medians <- scores[, .(
    candidate_abundance_UCell = stats::median(candidate_abundance_UCell),
    matched_direction_UCell = stats::median(matched_direction_UCell)
  ), by = .(
    dataset_id, disease, sampling_scope, sample, cell_type, eligible_label
  )]
  sample_medians <- merge(
    sample_medians,
    cell_counts,
    by = c(
      "dataset_id", "disease", "sampling_scope", "sample", "cell_type",
      "eligible_label"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  summary <- sample_medians[, .(
    n_samples = data.table::uniqueN(sample),
    n_cells = sum(cells),
    candidate_UCell_median = stats::median(candidate_abundance_UCell),
    candidate_UCell_q1 = v30_quantile(candidate_abundance_UCell, 0.25),
    candidate_UCell_q3 = v30_quantile(candidate_abundance_UCell, 0.75),
    matched_direction_UCell_median = stats::median(matched_direction_UCell),
    matched_direction_UCell_q1 = v30_quantile(matched_direction_UCell, 0.25),
    matched_direction_UCell_q3 = v30_quantile(matched_direction_UCell, 0.75)
  ), by = .(dataset_id, disease, sampling_scope, cell_type, eligible_label)]
  summary[, relative_candidate_rank := if (.N > 1L) {
    (rank(candidate_UCell_median, ties.method = "average") - 1) / (.N - 1)
  } else {
    1
  }, by = dataset_id]
  summary[, relative_direction_rank := if (.N > 1L) {
    (rank(matched_direction_UCell_median, ties.method = "average") - 1) / (.N - 1)
  } else {
    1
  }, by = dataset_id]
  summary[, score_definition := paste0(
    "UCell v", as.character(utils::packageVersion("UCell")),
    "; maxRank 1500; cell-level ranks summarized first within biological sample, ",
    "then across samples."
  )]
  summary[, inference_boundary := paste0(
    "The unsigned ten-gene score represents joint rank abundance, not pathway activity. ",
    "The signed score represents compatibility with discovery directions, not a ",
    "single-cell differential-expression effect."
  )]
  list(
    sample = as.data.frame(sample_medians),
    cell_type = as.data.frame(summary)
  )
}

v30_single_cell_advanced <- function(project_root, paths) {
  ucell <- v30_prepare_ucell(project_root, paths)
  evidence <- as.data.frame(ucell$evidence)
  evidence$annotation_boundary <- paste0(
    "Source labels are retained verbatim; Fibroblast is not relabeled CAF and ",
    "Ovarian.cancer.cell is not relabeled tumor epithelial."
  )
  evidence$comparison_boundary <- paste0(
    "Raw UMI means are not compared across datasets; fraction_detected supports ",
    "within-dataset CCSS only."
  )
  safe_write_csv(
    evidence,
    file.path(paths$tables, "Table_S10_single_cell_hub_gene_evidence.csv")
  )
  ccss <- v30_compute_ccss(evidence)
  ucell_summary <- v30_summarize_ucell(ucell$scores)
  safe_write_csv(
    ccss$context,
    file.path(paths$tables, "Table_S24a_dataset_context_CCSS.csv")
  )
  safe_write_csv(
    ccss$consensus,
    file.path(paths$tables, "Table_S24b_disease_consensus_CCSS.csv")
  )
  safe_write_csv(
    ucell_summary$cell_type,
    file.path(paths$tables, "Table_S24c_sample_aware_UCell_summary.csv")
  )
  write_sc_table(
    as.data.frame(ucell$scores),
    file.path(paths$analysis, "single_cell_candidate_UCell_per_cell.tsv.gz")
  )
  safe_write_csv(
    ucell$parameters,
    file.path(paths$analysis, "single_cell_UCell_parameters.csv")
  )
  safe_write_csv(
    ccss$consensus[ccss$consensus$minimum_cells == 100L, , drop = FALSE],
    file.path(paths$source, "Figure4_CCSS_consensus.csv")
  )
  safe_write_csv(
    ucell_summary$cell_type,
    file.path(paths$source, "Figure4_UCell_cell_type_summary.csv")
  )
  list(
    ucell = ucell,
    ccss = ccss,
    ucell_summary = ucell_summary,
    evidence = evidence
  )
}

v30_build_figure4_advanced <- function(paths, single_cell) {
  genes <- c("SOX9", "ELF3", "JUNB", "AKAP12", "BNC1", "CFI", "DDIT3", "DIRAS3", "EFEMP1", "HK2")
  consensus <- single_cell$ccss$consensus
  consensus <- consensus[
    consensus$minimum_cells == 100L & is.finite(consensus$median_CCSS),
    ,
    drop = FALSE
  ]
  consensus$gene <- factor(consensus$gene, levels = rev(genes))
  consensus$cell_type <- factor(
    consensus$cell_type,
    levels = sort(unique(consensus$cell_type))
  )
  p1 <- ggplot2::ggplot(
    consensus,
    ggplot2::aes(
      x = cell_type,
      y = gene,
      size = median_fraction_detected,
      fill = median_CCSS,
      alpha = evidence_status
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "#475467", stroke = 0.25) +
    ggplot2::facet_grid(. ~ disease, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradientn(
      colours = c("#F7FBFF", "#6BAED6", "#08306B"),
      limits = c(0, 1),
      oob = scales::squish
    ) +
    ggplot2::scale_size_continuous(range = c(0.7, 4.8), limits = c(0, 1), labels = scales::percent) +
    ggplot2::scale_alpha_manual(values = c(
      `replicated exact label` = 1,
      `dataset-specific exploratory` = 0.58
    )) +
    ggplot2::labs(
      title = "Candidate cell-context specificity",
      subtitle = "CCSS from detection fractions; >=100 cells per eligible label",
      x = NULL,
      y = NULL,
      fill = "Median CCSS",
      size = "Detection",
      alpha = NULL
    ) +
    submission_theme(6.6) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, size = 5.5),
      axis.text.y = ggplot2::element_text(size = 6),
      legend.position = "bottom"
    )

  ucell <- single_cell$ucell_summary$cell_type
  ucell <- ucell[ucell$eligible_label, , drop = FALSE]
  ucell$dataset_id <- factor(
    ucell$dataset_id,
    levels = rev(c("GSE104782", "GSE169454", "GSE255460", "GSE154600", "GSE180661"))
  )
  ucell$cell_type <- factor(ucell$cell_type, levels = sort(unique(ucell$cell_type)))
  p2 <- ggplot2::ggplot(
    ucell,
    ggplot2::aes(
      x = cell_type,
      y = dataset_id,
      fill = relative_candidate_rank,
      size = pmin(n_samples, 10)
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "#475467", stroke = 0.3) +
    ggplot2::facet_grid(. ~ disease, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFF7EC", "#FC8D59", "#7F0000"),
      limits = c(0, 1)
    ) +
    ggplot2::scale_size_continuous(range = c(1.1, 4.6), breaks = c(1, 5, 10)) +
    ggplot2::labs(
      title = "Sample-aware ten-gene UCell localization",
      subtitle = "Color is the within-atlas rank of sample-median UCell; values are not integrated across diseases",
      x = NULL,
      y = NULL,
      fill = "Within-atlas rank",
      size = "Samples\n(capped at 10)"
    ) +
    submission_theme(6.6) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, size = 5.5),
      axis.text.y = ggplot2::element_text(size = 6.2),
      legend.position = "bottom"
    )

  pseudo <- utils::read.csv(file.path(
    paths$source,
    "Figure5_hub_pseudobulk_evidence.csv"
  ))
  pseudo <- pseudo[is.finite(pseudo$logFC) & pseudo$FDR < 0.05, , drop = FALSE]
  if (nrow(pseudo) == 0L) {
    p3 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No eligible FDR-significant pseudobulk effects") +
      ggplot2::theme_void()
  } else {
    pseudo$label <- factor(
      paste(pseudo$dataset_id, pseudo$cell_type, pseudo$gene, sep = " | "),
      levels = rev(unique(paste(pseudo$dataset_id, pseudo$cell_type, pseudo$gene, sep = " | ")))
    )
    p3 <- ggplot2::ggplot(pseudo, ggplot2::aes(x = logFC, y = label, colour = gene)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#98A2B3") +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        title = "Eligible sample-level pseudobulk contrasts",
        subtitle = "Only contrasts with biological replication and FDR <0.05",
        x = "Pseudobulk log2 fold change",
        y = NULL,
        colour = "Gene"
      ) +
      submission_theme(6.6) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.3), legend.position = "right")
  }
  figure <- submission_panel_tag(p1 / p2 / p3 + patchwork::plot_layout(heights = c(1.15, 0.95, 1)))
  submission_save_plot(
    figure,
    "Figure4_cellular_context",
    paths$figures,
    height_mm = 225
  )
}

v30_replace_legend_section <- function(text, heading, next_heading, body) {
  pattern <- paste0(
    "(?s)",
    gsub("([\\W])", "\\\\\\1", heading, perl = TRUE),
    ".*?(?=",
    gsub("([\\W])", "\\\\\\1", next_heading, perl = TRUE),
    ")"
  )
  replacement <- paste0(heading, "\n\n", body, "\n\n")
  sub(pattern, replacement, text, perl = TRUE)
}

v30_update_advanced_legends <- function(paths) {
  path <- file.path(paths$figures, "figure_legends.md")
  text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  text <- v30_replace_legend_section(
    text,
    "## Figure 2. Gene-level convergence with an approximately balanced directional structure",
    "## Figure 3.",
    paste0(
      "**A-B,** OA and OC differential-expression volcano plots. **C,** The 286 primary shared DEGs ",
      "occupied four directional quadrants (146 concordant and 140 discordant). **D,** The inset isolates the ",
      "ten-gene interpretable evidence summary. **E,** High-confidence physical associations among ",
      "mapped shared-DEG products from STRING v12.0 (score >=0.700), colored by direction class. ",
      "Network topology was descriptive, was not used to select the ten genes, and does not establish ",
      "physical interaction in OA or OC tissue. Threshold sensitivity remains in Figure S1 and full ",
      "network audits are in Figure S12."
    )
  )
  text <- v30_replace_legend_section(
    text,
    "## Figure 4. Within-atlas candidate detection and cell-context localization",
    "## Figure 5.",
    paste0(
      "**A,** Candidate Cell Context Specificity Score (CCSS) based on within-dataset detection fractions ",
      "after excluding Unassigned, Ambiguous, Other, and strata with fewer than 100 cells. Exact source ",
      "labels were retained; disease summaries used medians only across identical labels. **B,** UCell ",
      "rank-based scores for the unsigned ten-gene evidence set, summarized within biological sample ",
      "before cell type and displayed as within-atlas ranks across all five datasets. GSE180661 used the ",
      "predeclared deterministic patient-balanced reference sample. **C,** Eligible FDR-significant ",
      "sample-level pseudobulk effects. UCell denotes joint rank abundance rather than pathway activity, ",
      "and neither CCSS nor UCell establishes cell-type-specific function or mechanism."
    )
  )
  text <- sub(
    "(?s)\\n*## Supplementary Figure 12\\..*$",
    "",
    text,
    perl = TRUE
  )
  text <- paste0(
    text,
    "\n\n## Supplementary Figure 12. Direction-aware STRING interaction context\n\n",
    paste0(
      "**A,** Complete connected high-confidence physical network among shared-DEG products. **B,** ",
      "Descriptive node degree and normalized betweenness; the ten-gene evidence set is marked but network ",
      "topology did not determine candidate selection. **C,** Concordant and discordant induced-subgraph ",
      "density in physical and functional high-confidence STRING networks, with mapped isolates retained. ",
      "**D,** Fixed-size direction-label permutation (10,000 iterations) for the physical-network density ",
      "difference. STRING associations are database evidence, not direct interaction measurements in OA or OC."
    ),
    "\n"
  )
  writeLines(text, path, useBytes = TRUE)
}

v30_update_advanced_index <- function(paths) {
  path <- file.path(paths$tables, "supplementary_table_index.csv")
  index <- utils::read.csv(path, check.names = FALSE)
  index <- index[!index$table_id %in% c("Table S24a-c", "Table S25a-e"), , drop = FALSE]
  additions <- data.frame(
    table_id = c("Table S24a-c", "Table S25a-e"),
    title = c(
      "Candidate cell-context specificity and sample-aware UCell summaries",
      "Direction-aware STRING mapping, edges, topology, and permutation"
    ),
    contents = c(
      paste0(
        "Dataset/context CCSS at 50/100/200-cell thresholds, exact-label disease consensus, ",
        "and sample-aware unsigned/direction-compatible UCell summaries."
      ),
      paste0(
        "STRING v12.0 mapping audit; high-confidence physical and functional edges; node and ",
        "subgraph topology; 10,000 fixed-size direction-label permutations."
      )
    ),
    stringsAsFactors = FALSE
  )
  missing_columns <- setdiff(names(index), names(additions))
  for (column in missing_columns) {
    additions[[column]] <- NA_character_
  }
  additions <- additions[, names(index), drop = FALSE]
  index <- rbind(index, additions)
  safe_write_csv(index, path)
}

v30_update_advanced_registry <- function(paths) {
  path <- file.path(paths$root, "claim_evidence_registry_v30.csv")
  registry <- utils::read.csv(path, check.names = FALSE)
  row <- registry$claim_id == "C08"
  registry$manuscript_claim[row] <- paste0(
    "Across five disease-specific atlases, CCSS and sample-aware UCell summaries placed the ",
    "ten-gene evidence set in distinct annotated cellular contexts; eligible pseudobulk ",
    "contrasts supplied sample-level association evidence."
  )
  registry$figure_or_table[row] <- paste0(
    "Figure 4; Figures S3 and S6; Tables S9-S10, S19, and S24"
  )
  registry$allowed_wording[row] <- paste0(
    "within-atlas CCSS; rank-based UCell localization; sample-level pseudobulk association"
  )
  registry$prohibited_wording[row] <- paste0(
    "cross-disease integrated cell state; pathway activity; cell type causes disease"
  )
  if (!"C25" %in% registry$claim_id) {
    template <- registry[1L, , drop = FALSE]
    template[1, ] <- NA
    template$claim_id <- "C25"
    template$manuscript_claim <- paste0(
      "High-confidence STRING networks contained structured associations among a subset of ",
      "mapped shared-DEG products, with direction-aware topology assessed by fixed-size label permutation."
    )
    template$primary_data <- "STRING v12.0 API cache and direction-aware graph analysis"
    template$figure_or_table <- "Figure 2; Figure S12; Tables S25a-e"
    template$allowed_wording <- "interaction context; database-supported association; descriptive topology"
    template$prohibited_wording <- "shared pathogenic network; validated protein interaction; PPI-defined disease mechanism"
    template$status <- "audited"
    registry <- rbind(registry, template)
  }
  safe_write_csv(registry, path)
}

run_reviewer_v30 <- function(project_root) {
  paths <- v30_advanced_base_runner(project_root)
  log_info("Adding direction-aware STRING interaction context.")
  ppi <- v30_prepare_ppi(project_root, paths)
  v30_build_figure2_with_ppi(project_root, paths, ppi)
  v30_build_ppi_supplement(paths, ppi)
  log_info("Adding CCSS, UCell, and sample-aware single-cell localization.")
  single_cell <- v30_single_cell_advanced(project_root, paths)
  v30_build_figure4_advanced(paths, single_cell)
  v30_update_advanced_legends(paths)
  v30_update_advanced_index(paths)
  v30_update_advanced_registry(paths)
  log_info(
    "V3.0 advanced layers completed: bounded PPI context and five-atlas single-cell scoring."
  )
  invisible(paths)
}
