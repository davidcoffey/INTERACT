# Integrated Immune-Tumor Signaling Viewer
# Global — loads precomputed data only; run precompute.R first
# David Coffey davidcoffey@miami.edu

# 0. Libraries ----
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinycssloaders)
  library(shinyjs)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(ggsci)
  library(visNetwork)
  library(igraph)
  library(survival)
  library(survminer)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(readr)
})

# Load results from shinyOptions
results_opt <- shiny::getShinyOption("results")
if (is.null(results_opt)) {
  # Fallback for testing
  if (file.exists("Data/INTERACT results.Rda")) {
    load("Data/INTERACT results.Rda")
    results_opt <- results

  }
}

precomputed_data <- results_opt$precompute
if (is.null(precomputed_data)) {
  stop("Precomputed data not found in results object. Please run RunInteract() first.")
}

message("Loading precomputed app data...")
list2env(precomputed_data, environment())

# Backward compatibility: map old 'patient' variables to 'sample' if they exist
if (exists("patient.ids.vec") && !exists("sample.ids.vec")) sample.ids.vec <- patient.ids.vec
if (exists("patient.coords.base") && !exists("sample.coords.base")) sample.coords.base <- patient.coords.base
if (exists("patient.stats.summary") && !exists("sample.stats.summary")) sample.stats.summary <- patient.stats.summary
if (exists("patient.interaction.counts") && !exists("sample.interaction.counts")) sample.interaction.counts <- patient.interaction.counts
if (exists("pt.profile.long") && !exists("sample.profile.long")) sample.profile.long <- pt.profile.long
if (exists("pt.pca.x") && !exists("sample.pca.x")) sample.pca.x <- pt.pca.x

# Fallbacks for removed clinical variables to prevent app crashes
if (!exists("survival.df")) {
  survival.df <- data.frame(
    sample = sample.ids.vec, cluster = "Unknown",
    SPECTRUM_SEQ = sample.ids.vec,
    DEMOG_GENDER = "Unknown", DEMOG_SAMPLEAGE = NA,
    chr17p_Loss = "Unknown", chr1q_Gain = "Unknown",
    t4_14_MMSET = "Unknown", t14_16_MAF = "Unknown",
    t14_20_MAFB = "Unknown", t11_14_CCND1 = "Unknown",
    t_MYC = "Unknown", best.response = "Unknown",
    derived.pfs.day = NA, censpfs = NA,
    derived.os.day = NA, censos = NA,
    pfs.months = NA, os.months = NA,
    stringsAsFactors = FALSE
  )
}
if (!exists("clinical.sample")) clinical.sample <- survival.df
if (!exists("clinical.var.labels")) clinical.var.labels <- c(DEMOG_GENDER = "Sex", cluster = "Cluster")
if (!exists("ligand.axes.available")) ligand.axes.available <- character(0)
if (!exists("driver.genes")) driver.genes <- character(0)
if (!exists("circos.data")) circos.data <- data.frame()

message("Done.")

message("Indexing interactions...")
full_interactions = data.table::as.data.table(results_opt$interactions)
if (!data.table::haskey(full_interactions)) {
  data.table::setkey(full_interactions, sample)
}

if (!exists("all.populations")) {
  all.populations <- sort(unique(c(full_interactions$sender, full_interactions$receiver)))
}
if (!exists("pop.colors")) {
  pop.colors <- setNames(rep("#333333", length(all.populations)), all.populations)
}
if (!exists("all.gene.symbols")) {
  all.gene.symbols <- sort(unique(c(
    unique(unlist(strsplit(as.character(full_interactions$ligand), "_"))),
    unique(unlist(strsplit(as.character(full_interactions$receptor), "_")))
  )))
}

if (!exists("edge.gradient.summary") || is.null(edge.gradient.summary)) {
  if (!is.null(results_opt$edge.gradients)) {
    edge.gradient.summary <- results_opt$edge.gradients
  } else {
    edge.gradient.summary <- data.frame()
  }
}

if (!exists("mean.layer.weights.summary") || is.null(mean.layer.weights.summary) || nrow(mean.layer.weights.summary) == 0) {
  if (exists("full_interactions") && nrow(full_interactions) > 0 &&
      all(c("receiver", "weight.pathways", "weight.programs", "weight.genomic") %in% colnames(full_interactions))) {
    mean.layer.weights.summary <- full_interactions %>%
      dplyr::group_by(receiver) %>%
      dplyr::summarize(
        Pathways = mean(weight.pathways, na.rm = TRUE),
        Programs = mean(weight.programs, na.rm = TRUE),
        Genomics = mean(weight.genomic, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    mean.layer.weights.summary <- data.frame()
  }
}

if (!exists("rank.shift.summary") || is.null(rank.shift.summary) || nrow(rank.shift.summary) == 0) {
  if (exists("full_interactions") && nrow(full_interactions) > 0 &&
      all(c("sample", "receiver", "base.score", "final.score") %in% colnames(full_interactions))) {
    rank.shift.summary <- full_interactions %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        base.rank = dplyr::min_rank(dplyr::desc(base.score)),
        final.rank = dplyr::min_rank(dplyr::desc(final.score)),
        rank.shift = base.rank - final.rank
      ) %>%
      dplyr::group_by(receiver) %>%
      dplyr::summarize(
        mean.rank.shift = mean(rank.shift, na.rm = TRUE),
        mean.abs.rank.shift = round(mean(abs(rank.shift), na.rm = TRUE), 2),
        pct.shifted = round(mean(abs(rank.shift) > 0, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      )
  }
}

# 2. Helper Functions ----

# URL builders
genecards.url = function(gene) {
  sprintf("https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s", gene)
}
reactome.url = function(pathway) {
  sprintf("https://reactome.org/content/query?q=%s&types=Pathway",
          utils::URLencode(pathway, reserved = TRUE))
}

gene.link.html = function(gene) {
  sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>',
          genecards.url(gene), gene)
}

# Load program.gene.lists from results
program.gene.lists <- results_opt$program.gene.lists
if (is.null(program.gene.lists)) {
  # Fallback for old results
  program.gene.lists <- data.frame(program.name=character(), url=character(), stringsAsFactors=FALSE)
}
immport.links = setNames(program.gene.lists$url, program.gene.lists$program.name)

pathway.link.html = function(pathway) {
  ifelse(is.na(pathway) | pathway == "",
    "—",
    ifelse(pathway %in% names(immport.links),
      sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', immport.links[pathway], pathway),
      sprintf('<a href="https://saezlab.github.io/progeny/reference/model_human_full.html" target="_blank" rel="noopener">%s</a>', pathway)
    )
  )
}
omnipath.link.html = function(ligand, receptor, label) {
  ifelse(is.na(ligand) | is.na(receptor) | ligand == "" | receptor == "",
    "—",
    sprintf('<a href="https://explore.omnipathdb.org/search?q=%s%%2C+%s%%2C&tab=interactions&species=9606" target="_blank" rel="noopener">%s</a>',
            utils::URLencode(ligand, reserved = TRUE), 
            utils::URLencode(receptor, reserved = TRUE), 
            label)
  )
}

# Load a single sample's interaction data from memory
load.sample.dt = function(pt.id) {
  pt.data <- full_interactions[sample == pt.id]
  if (nrow(pt.data) == 0) return(data.table::data.table())
  return(pt.data)
}

# Build sample visNetwork graph from a per-sample data.table
build.sample.graph = function(pt.dt, top.n = 150, node.fill.by = "role", node.outline.by = "role", edge.color.by = "pathway", highlight.cell = "None", search.ligand = NULL, search.receptor = NULL) {
  if (is.null(pt.dt) || nrow(pt.dt) == 0) return(NULL)

  # Network Input & Size - Prunes the network to include only the top 25% of all cell-cell interactions per patient based on sum.score to filter background noise
  cell.network.edges = pt.dt[sender != receiver, .(sum.score = sum(final.score, na.rm = TRUE)), by = .(sender, receiver)]
  if (nrow(cell.network.edges) > 0) {
    score.threshold = stats::quantile(cell.network.edges$sum.score, probs = 0.75, na.rm = TRUE)
    valid.pairs = cell.network.edges[sum.score >= score.threshold]
    
    # Filter molecular interactions to only those matching valid cross-population cell-cell pairs
    pt.filtered = pt.dt[sender != receiver]
    pt.filtered = pt.filtered[paste(sender, receiver, sep = "_") %in% paste(valid.pairs$sender, valid.pairs$receiver, sep = "_")]
  } else {
    pt.filtered = pt.dt[sender != receiver]
  }

  pt.data = pt.filtered[order(-final.score)][1:min(.N, top.n)]

  edges.agg = pt.data[, .(
    weight      = sum(final.score, na.rm = TRUE),
    n.pairs     = .N,
    top.lr      = paste(head(paste0(source_genesymbol, "→", target_genesymbol), 5), collapse = ", "),
    top.lr.1    = head(paste0(source_genesymbol, "→", target_genesymbol), 1),
    pct.genomic = mean(as.numeric(abs(weight.genomic) > 1e-10), na.rm = TRUE) * 100,
    has.searched.gene = if (!is.null(search.ligand) || !is.null(search.receptor)) {
                          any((length(search.ligand) == 0 | source_genesymbol %in% search.ligand) & 
                              (length(search.receptor) == 0 | target_genesymbol %in% search.receptor))
                        } else FALSE,
    top.axis    = if (any(!is.na(ligand.axis)))
                    names(which.max(table(ligand.axis[!is.na(ligand.axis)])))
                  else NA_character_
  ), by = .(sender, receiver)]
  data.table::setnames(edges.agg, c("sender", "receiver"), c("from", "to"))

  g = igraph::graph_from_data_frame(as.data.frame(edges.agg), directed = TRUE)

  in.str  = igraph::strength(g, mode = "in",  weights = abs(igraph::E(g)$weight))
  out.str = igraph::strength(g, mode = "out", weights = abs(igraph::E(g)$weight))

  # Betweenness Weights - Weighted: Inverts strength weights to compute biological distance (D = 1/W)
  # A strong interaction (high score) means a "short" biological distance
  distance.weights = 1 / igraph::E(g)$weight
  between = igraph::betweenness(g, directed = TRUE, weights = distance.weights)
  # Normalize to [0, 1] for robust role-thresholding comparison
  if (max(between, na.rm = TRUE) > 0) {
    between = between / max(between, na.rm = TRUE)
  }

  # Use hits_scores which replaces deprecated hub_score and authority_score
  hits    = igraph::hits_scores(g, weights = abs(igraph::E(g)$weight))
  hub.sc  = hits$hub
  auth.sc = hits$authority

  role = dplyr::case_when(
    between >= 0.3               ~ "Mediator",
    hub.sc  > 0.5 & auth.sc > 0.5 ~ "Influencer",
    out.str  > in.str  * 1.5    ~ "Sender",
    in.str   > out.str * 1.5    ~ "Receiver",
    TRUE                         ~ "Influencer"
  )

  role.colors = c(
    Sender     = "#E64B35FF",
    Receiver   = "#4DBBD5FF",
    Influencer = "#00A087FF",
    Mediator   = "#3C5488FF"
  )

  node.names = igraph::V(g)$name
  cell.fill  = pop.colors[node.names]
  cell.fill[is.na(cell.fill)] = "#AAAAAA"

  get_node_color = function(color_by, role, hub.sc, auth.sc, cell.fill) {
    if (color_by == "cell_type") {
      return(cell.fill)
    } else if (color_by == "hub") {
      pal = colorRampPalette(c("white", "red"))(100)
      return(pal[as.numeric(cut(hub.sc, breaks = 100))])
    } else if (color_by == "authority") {
      pal = colorRampPalette(c("white", "blue"))(100)
      return(pal[as.numeric(cut(auth.sc, breaks = 100))])
    } else if (color_by == "none") {
      return("transparent")
    } else {
      return(role.colors[role])
    }
  }

  bg.col = get_node_color(node.fill.by, role, hub.sc, auth.sc, cell.fill)
  border.col = get_node_color(node.outline.by, role, hub.sc, auth.sc, cell.fill)

  border.width = ifelse(node.outline.by == "none", 0, 4)

  # Cell highlighting logic
  node.opacity = rep(1, length(node.names))
  edge.opacity = rep(0.5, nrow(edges.agg)) # default edge opacity

  is.gene.search = length(search.ligand) > 0 || length(search.receptor) > 0

  if (highlight.cell != "None" && highlight.cell %in% all.populations) {
    short.highlight = gsub("scRNAseq ", "", highlight.cell)
    # nodes to highlight
    is.hi.node = (node.names == short.highlight) | 
                 (node.names %in% edges.agg$from[edges.agg$to == highlight.cell]) | 
                 (node.names %in% edges.agg$to[edges.agg$from == highlight.cell])
    node.opacity[!is.hi.node] = 0.15
    
    # edges to highlight
    is.hi.edge = (edges.agg$from == highlight.cell) | (edges.agg$to == highlight.cell)
    edge.opacity[!is.hi.edge] = 0.1
    edge.opacity[is.hi.edge] = 0.9
  } else if (is.gene.search) {
    # edges to highlight
    is.hi.edge = edges.agg$has.searched.gene
    edge.opacity[!is.hi.edge] = 0.1
    edge.opacity[is.hi.edge] = 0.9
    
    # nodes to highlight
    hi.nodes = unique(c(edges.agg$from[is.hi.edge], edges.agg$to[is.hi.edge]))
    is.hi.node = (paste0("scRNAseq ", node.names) %in% hi.nodes | node.names %in% hi.nodes)
    node.opacity[!is.hi.node] = 0.15
  }

  apply_opacity = function(col_vec, opac_vec) {
    sapply(seq_along(col_vec), function(i) {
      if (col_vec[i] == "transparent") return("transparent")
      if (grepl("^rgba", col_vec[i])) {
        # Extract RGB components and replace alpha
        rgb_vals = gsub("rgba\\(|\\)", "", col_vec[i])
        rgb_vals = strsplit(rgb_vals, ",")[[1]]
        return(sprintf("rgba(%s,%s,%s,%.2f)", trimws(rgb_vals[1]), trimws(rgb_vals[2]), trimws(rgb_vals[3]), opac_vec[i]))
      }
      col = col2rgb(col_vec[i])
      sprintf("rgba(%d,%d,%d,%.2f)", col[1], col[2], col[3], opac_vec[i])
    })
  }

  nodes = data.frame(
    id                       = node.names,
    label                    = gsub("scRNAseq ", "", node.names),
    color.background         = apply_opacity(bg.col, node.opacity),
    color.border             = apply_opacity(border.col, node.opacity),
    color.highlight.background = bg.col,
    color.highlight.border   = border.col,
    borderWidth              = border.width,
    size                     = pmax(22, 15 + 35 * (hub.sc + auth.sc) / 2),
    shadow                   = TRUE,
    title                    = paste0(
      "<b>", node.names, "</b><br>",
      "Role: <b>", role, "</b><br>",
      "Out-strength: ", round(out.str, 2), "<br>",
      "In-strength: ", round(in.str, 2), "<br>",
      "Betweenness: ", round(between, 4), "<br>",
      "Hub score: ", round(hub.sc, 4), "<br>",
      "Authority score: ", round(auth.sc, 4)
    ),
    stringsAsFactors = FALSE
  )

  leg_info = list(type = "none")
  edge.col.vec = rep("rgba(150,150,150,0.5)", nrow(edges.agg))
  if (edge.color.by == "pathway") {
    unq.axes = unique(edges.agg$top.axis[!is.na(edges.agg$top.axis)])
    if(length(unq.axes) > 0) {
      axis.pal = setNames(ggsci::pal_igv()(length(unq.axes)), unq.axes)
      edge.col.vec = ifelse(!is.na(edges.agg$top.axis), axis.pal[edges.agg$top.axis], "rgba(150,150,150,0.5)")
      leg_info = list(type = "discrete", pal = axis.pal)
    }
  } else if (edge.color.by == "lr") {
    unq.lr = unique(edges.agg$top.lr.1)
    lr.pal = setNames(ggsci::pal_igv()(length(unq.lr)), unq.lr)
    edge.col.vec = lr.pal[edges.agg$top.lr.1]
    leg_info = list(type = "discrete", pal = lr.pal)
  } else if (edge.color.by == "weight") {
    pal = colorRampPalette(c("lightgrey", "black"))(100)
    edge.col.vec = pal[as.numeric(cut(edges.agg$weight, breaks = 100))]
    leg_info = list(type = "continuous", pal = pal)
  } else if (edge.color.by == "genomic") {
    has.genomic = edges.agg$pct.genomic > 0
    edge.col.vec = ifelse(has.genomic, "#E64B35FF", "rgba(150,150,150,0.5)")
    leg_info = list(type = "discrete", pal = c("Genomic Hit" = "#E64B35FF", "No Hit" = "rgba(150,150,150,0.5)"))
  }

  # Apply edge opacity
  edge.col.vec = apply_opacity(edge.col.vec, edge.opacity)

  edges.vis = as.data.frame(edges.agg) %>%
    dplyr::mutate(
      id = paste0(from, "_", to),
      value = weight,
      color = edge.col.vec,
      title = paste0(
        "<b>", from, " → ", to, "</b><br>",
        "Total weight: ", round(weight, 3), "<br>",
        "L-R pairs: ", n.pairs, "<br>",
        "Genomically weighted: ", round(pct.genomic, 1), "%<br>",
        "<i>Top pairs:</i> ", top.lr
      )
    )

  metrics.df = data.frame(
    node            = node.names,
    role            = role,
    out.strength    = round(out.str, 3),
    in.strength     = round(in.str, 3),
    betweenness     = round(between, 4),
    hub.score       = round(hub.sc, 4),
    authority.score = round(auth.sc, 4),
    stringsAsFactors = FALSE
  )

  list(nodes = nodes, edges = edges.vis, metrics = metrics.df, legend = leg_info)
}



# Cluster enrichment from reactive cluster vector
compute.cluster.enrichment = function(cluster.vec) {
  cluster.df = data.frame(sample = names(cluster.vec), cluster = as.character(cluster.vec),
                           stringsAsFactors = FALSE)
  pt.profile.long %>%
    as.data.frame() %>%
    dplyr::inner_join(cluster.df, by = "sample") %>%
    dplyr::group_by(cluster, sender, receiver) %>%
    dplyr::summarize(mean.score = mean(total.score, na.rm = TRUE), .groups = "drop")
}

# Null-coalescing operator
`%||%` = function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# Null-coalescing operator (used in server.R)
`%||%` = function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1]) && a[1] != "") a else b

message("Global setup complete.")
