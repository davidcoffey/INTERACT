# Integrated Immune-Tumor Signaling Viewer
# Server Logic

server = function(input, output, session) {

  # 0. Initialization ----
  # Pre-register plotly shiny event IDs to prevent initialization warnings before tabs render
  session$userData$plotlyShinyEventIDs = unique(c(session$userData$plotlyShinyEventIDs, "plotly_click-cluster.scatter"))

  updateSelectizeInput(session, "sample.id",
                       choices  = sort(sample.ids.vec),
                       selected = sort(sample.ids.vec)[1],
                       server   = TRUE)

  updateSelectInput(session, "home.table.axis",
                    choices  = c("All" = "all", setNames(ligand.axes.available, ligand.axes.available)))
  updateSelectInput(session, "home.table.sender",
                    choices  = c("All" = "all", setNames(all.populations, all.populations)))
  updateSelectInput(session, "home.table.receiver",
                    choices  = c("All" = "all", setNames(all.populations, all.populations)))

  updateSelectizeInput(session, "genomic.highlight.genes",
                       choices  = all.gene.symbols,
                       server   = TRUE)

  # Update sample network search genes
  updateSelectizeInput(session, "sample.ligand.search", choices = all.gene.symbols, server = TRUE)
  updateSelectizeInput(session, "sample.receptor.search", choices = all.gene.symbols, server = TRUE)

  # Update compare tab inputs
  updateSelectizeInput(session, "compare.sample.ids", choices = sort(sample.ids.vec), server = TRUE)
  updateSelectizeInput(session, "compare.ligand.search", choices = all.gene.symbols, server = TRUE)
  updateSelectizeInput(session, "compare.receptor.search", choices = all.gene.symbols, server = TRUE)

  # Sample L-R filter choices updated when sample changes
  observeEvent(input$sample.id, {
    pt.senders   = sort(unique(all.populations))
    pt.receivers = sort(unique(all.populations))
    updateSelectInput(session, "sample.lr.sender",
                      choices = c("All" = "all", setNames(pt.senders, pt.senders)))
    updateSelectInput(session, "sample.lr.receiver",
                      choices = c("All" = "all", setNames(pt.receivers, pt.receivers)))
  })

  # 1. Home Tab ----

  # Value boxes
  output$vbox.samples = renderValueBox({
    valueBox(
      value    = format(length(sample.ids.vec), big.mark = ","),
      subtitle = "Samples",
      icon     = icon("users"),
      color    = "blue"
    )
  })
  output$vbox.interactions = renderValueBox({
    valueBox(
      value    = format(nrow(top.lr.global), big.mark = ","),
      subtitle = "Unique Interactions",
      icon     = icon("project-diagram"),
      color    = "purple"
    )
  })
  output$vbox.celltypes = renderValueBox({
    valueBox(
      value    = length(all.populations),
      subtitle = "Cell Populations",
      icon     = icon("share-nodes"),
      color    = "green"
    )
  })
  output$vbox.genomic.pct = renderValueBox({
    valueBox(
      value    = paste0(round(weight.layer.summary$pct.genomic.weighted, 1), "%"),
      subtitle = "Genomically Weighted Interactions",
      icon     = icon("dna"),
      color    = "orange"
    )
  })

  # Global heatmap
  output$home.heatmap = renderPlotly({
    mat.wide = sender.receiver.agg %>%
      dplyr::select(sender, receiver, mean.score) %>%
      tidyr::pivot_wider(names_from = receiver, values_from = mean.score, values_fill = 0)
    mat.vals = as.matrix(mat.wide[, -1])
    rownames(mat.vals) = mat.wide$sender

    plotly::plot_ly(
      x         = colnames(mat.vals),
      y         = rownames(mat.vals),
      z         = mat.vals,
      type      = "heatmap",
      colorscale = "Blues",
      reversescale = TRUE,
      hovertemplate = "Sender: %{y}<br>Receiver: %{x}<br>Mean Score: %{z:.4f}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis  = list(title = "Receiver", tickangle = -45, tickfont = list(size = 9)),
        yaxis  = list(title = "Sender",   tickfont  = list(size = 9)),
        margin = list(b = 140, l = 200)
      )
  })

  # Sample Interaction Counts
  output$home.interaction.counts = renderPlotly({
    df = sample.interaction.counts %>%
      dplyr::arrange(dplyr::desc(interaction.count))

    top10 = head(df, 10) %>% dplyr::mutate(group = "Top 10")
    bot10 = tail(df, 10) %>% dplyr::mutate(group = "Bottom 10")

    plot_df = rbind(top10, bot10) %>%
      dplyr::filter(!is.na(interaction.count))

    plotly::plot_ly(plot_df,
      x = ~interaction.count,
      y = ~reorder(sample, interaction.count),
      type = "bar", orientation = "h",
      color = ~group,
      colors = c("Top 10" = "#E64B35FF", "Bottom 10" = "#4DBBD5FF"),
      hovertemplate = "Sample: %{y}<br>Interactions: %{x}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis = list(title = "Total Interactions"),
        yaxis = list(title = ""),
        margin = list(l = 100),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15)
      )
  })

  # Top sender-receiver bar
  output$home.bar.axes = renderPlotly({
    top.axes = sender.receiver.agg %>%
      dplyr::arrange(dplyr::desc(mean.score)) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::mutate(pair = paste0(gsub("scRNAseq ", "", sender), " → ",
                                   gsub("scRNAseq ", "", receiver)))

    top.axes = top.axes %>%
      dplyr::filter(!is.na(mean.score))

    plotly::plot_ly(top.axes,
      x = ~mean.score, y = ~reorder(pair, mean.score),
      type = "bar", orientation = "h",
      marker = list(color = ggsci::pal_npg()(1)),
      hovertemplate = "%{y}<br>Mean Score: %{x:.4f}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis  = list(title = "Mean Interaction Score"),
        yaxis  = list(title = "", tickfont = list(size = 9)),
        margin = list(l = 180)
      )
  })

  # Top Ligand-Receptor Sankey
  output$home.sankey = renderPlotly({
    df = top.ecosystem.events %>%
      dplyr::group_by(source_genesymbol, target_genesymbol) %>%
      dplyr::summarise(mean.score = abs(sum(mean.final.score)), .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(mean.score)) %>%
      dplyr::slice_head(n = 15)

    nodes = data.frame(name = unique(c(df$source_genesymbol, df$target_genesymbol)))
    df$source = match(df$source_genesymbol, nodes$name) - 1
    df$target = match(df$target_genesymbol, nodes$name) - 1

    plotly::plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes$name,
        pad = 15,
        thickness = 20,
        line = list(color = "black", width = 0.5)
      ),
      link = list(
        source = df$source,
        target = df$target,
        value  = df$mean.score,
        hovertemplate = "%{source.label} &rarr; %{target.label}<br>Mean Score: %{value:.2f}<extra></extra>"
      )
    ) %>%
      plotly::layout(margin = list(l = 10, r = 10, t = 10, b = 10))
  })
  # Ecosystem table
  output$home.programs.dt = renderDT({
    req(program.gene.lists)
    if(nrow(program.gene.lists) == 0) return(NULL)
    df = program.gene.lists %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      dplyr::mutate(
        Program = sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', url, program.name),
        Genes = count
      ) %>%
      dplyr::select(Program, Genes)

    DT::datatable(df,
      escape = FALSE,
      rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 10,
        dom = "tp",
        scrollY = "320px",
        paging = FALSE
      )
    )
  })

  # Ecosystem table
  output$home.eco.table = renderDT({
    n.show = input$home.table.n
    df = top.ecosystem.events %>% as.data.frame()

    if (input$home.table.axis    != "all") df = dplyr::filter(df, ligand.axis == input$home.table.axis)
    if (input$home.table.sender  != "all") df = dplyr::filter(df, sender == input$home.table.sender)
    if (input$home.table.receiver!= "all") df = dplyr::filter(df, receiver == input$home.table.receiver)

    df = df %>%
      dplyr::slice_head(n = n.show) %>%
      dplyr::mutate(
        Ligand        = gene.link.html(source_genesymbol),
        Receptor      = gene.link.html(target_genesymbol),
        `Ligand Axis` = omnipath.link.html(source_genesymbol, target_genesymbol, ligand.axis),
        `Sender`      = gsub("scRNAseq ", "", sender),
        `Receiver`    = gsub("scRNAseq ", "", receiver),
        `Mean Score`  = round(mean.final.score, 4),
        `Base Score`  = round(mean.base.score, 4),
        `N Samples`  = n.samples,
        `% Genomic`   = round(genomic.rate, 1)
      ) %>%
      dplyr::select(Ligand, Receptor, `Ligand Axis`, Sender, Receiver,
                    `Mean Score`, `Base Score`, `N Samples`, `% Genomic`)

    DT::datatable(df,
      escape    = FALSE,
      rownames  = FALSE,
      selection = "none",
      options   = list(
        pageLength = 15, scrollX = TRUE,
        order      = list(list(5, "desc")),
        dom        = "ftip"
      )
    ) %>%
      DT::formatStyle("% Genomic",
        background = DT::styleColorBar(c(0, 100), "#4DBBD5AA"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      ) %>%
      DT::formatStyle("Mean Score",
        background = DT::styleColorBar(c(0, max(df$`Mean Score`, na.rm=TRUE)), "#E64B35AA"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
      ) %>%
      DT::formatStyle("Base Score",
        background = DT::styleColorBar(c(0, max(df$`Base Score`, na.rm=TRUE)), "#00A087AA"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
      ) %>%
      DT::formatStyle("N Samples",
        background = DT::styleColorBar(c(0, max(df$`N Samples`, na.rm=TRUE)), "#3C5488AA"),
        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat", backgroundPosition = "center"
      )
  })

  # 2. Sample Clusters Tab ----

  # Reactive cluster assignment — supports Louvain resolution, k-means, and precomputed clusters
  reactive.clusters = reactive({
    # 1. Louvain clustering by resolution if available
    if (!is.null(input$cluster.res) && exists("all.louvain") && is.list(all.louvain) && length(all.louvain) > 0) {
      res.val = as.numeric(input$cluster.res)
      res.key = paste0("res_", sprintf("%.1f", res.val))
      if (!is.null(all.louvain[[res.key]])) {
        cl = all.louvain[[res.key]]
        return(if (all(grepl("^Cluster", cl))) cl else paste0("Cluster ", cl))
      }
    }

    # 2. Precomputed k-means if input$cluster.k is available
    if (!is.null(input$cluster.k) && length(input$cluster.k) > 0 && !is.na(input$cluster.k)) {
      k = as.integer(input$cluster.k)
      if (!is.na(k) && length(k) > 0) {
        key = paste0("k", k)
        if (exists("all.kmeans") && is.list(all.kmeans) && !is.null(all.kmeans[[key]])) {
          return(paste0("Cluster ", all.kmeans[[key]]$cluster))
        }
      }
    }

    # 3. Base coords precomputed clusters
    if (exists("sample.coords.base") && "cluster" %in% colnames(sample.coords.base)) {
      cl = as.character(sample.coords.base$cluster)
      return(if (all(grepl("^Cluster", cl))) cl else paste0("Cluster ", cl))
    }

    # Fallback
    rep("Cluster 1", length(sample.ids.vec))
  })

  # Sample coords with reactive clusters
  filtered.coords = reactive({
    clusters = reactive.clusters()
    df = sample.coords.base %>%
      dplyr::mutate(cluster = clusters)
    df
  })

  # Consistent color mapping across cluster scatter and network graph
  cluster.color.mapping = reactive({
    df = filtered.coords()
    col.by = input$cluster.color %||% "cluster"

    vals = sort(unique(as.character(df[[col.by]])))
    vals[is.na(vals) | vals == ""] = "Unknown"
    n.vals = length(vals)

    if (n.vals <= 10) {
      pal = ggsci::pal_npg("nrc")(max(3, n.vals))[1:n.vals]
    } else if (requireNamespace("pals", quietly = TRUE)) {
      pal = as.character(pals::polychrome(n.vals))
    } else {
      pal = grDevices::rainbow(n.vals)
    }
    stats::setNames(pal, vals)
  })

  # Dynamic Projection Control
  output$cluster.proj.ui = renderUI({
    has.gnn = (exists("gnn.embeddings") && !is.null(gnn.embeddings) && is.data.frame(gnn.embeddings) && nrow(gnn.embeddings) > 0) ||
              (exists("embed.label") && any(grepl("GNN", embed.label)))

    if (has.gnn) {
      radioButtons("cluster.proj", "Projection:",
                   choices = c("GNN UMAP" = "umap", "PCA" = "pca"),
                   selected = "umap", inline = TRUE)
    } else {
      radioButtons("cluster.proj", "Projection:",
                   choices = c("PCA" = "pca"),
                   selected = "pca", inline = TRUE)
    }
  })

  output$cluster.map.subtitle = renderUI({
    p(class = "text-muted", style = "font-size: 12px; margin-bottom: 5px;",
      "Click any sample node to inspect its microenvironment.")
  })

  output$cluster.network.subtitle = renderUI({
    p(class = "text-muted", style = "font-size: 12px; margin-bottom: 5px;",
      "K-nearest neighbor patient similarity graph.")
  })

  output$cluster.scatter = renderPlotly({
    df = filtered.coords()
    col.by = input$cluster.color %||% "cluster"
    df[[col.by]][is.na(df[[col.by]])] = "Unknown"

    pal.map = cluster.color.mapping()

    proj = input$cluster.proj %||% (if (exists("embed.label") && any(grepl("GNN", embed.label))) "umap" else "pca")
    if (proj == "umap") {
      x.col = "Dim1"; y.col = "Dim2"
      x.lab = embed.label[1]; y.lab = embed.label[2]
    } else {
      x.col = "PC1"; y.col = "PC2"
      x.lab = "PC1"; y.lab = "PC2"
    }

    leg.title = if (col.by == "cluster") "Cluster" else unname(clinical.var.labels[col.by])
    if (is.null(leg.title) || is.na(leg.title)) leg.title = col.by

    df.plot = df %>%
      dplyr::filter(!is.na(.data[[x.col]]), !is.na(.data[[y.col]]))

    plotly::plot_ly(df.plot,
      x           = ~get(x.col),
      y           = ~get(y.col),
      color       = ~get(col.by),
      colors      = pal.map,
      customdata  = ~sample,
      type        = "scatter", mode = "markers",
      text        = ~paste0("Sample: ", sample,
                            "<br>Cluster: ", cluster),
      hoverinfo   = "text",
      marker      = list(size = 8, opacity = 0.8,
                         line = list(color = "white", width = 0.5)),
      source      = "cluster.scatter"
    ) %>%
      plotly::layout(
        xaxis  = list(title = x.lab, titlefont = list(color = "#000000")),
        yaxis  = list(title = y.lab, titlefont = list(color = "#000000")),
        legend = list(title = list(text = leg.title, font = list(color = "#000000")),
                      font = list(color = "#000000"))
      ) %>%
      plotly::event_register("plotly_click")
  })

  outputOptions(output, "cluster.scatter", suspendWhenHidden = FALSE)

  # Tab switch on cluster scatter click
  observeEvent(suppressWarnings(plotly::event_data("plotly_click", source = "cluster.scatter")), {
    click.data = suppressWarnings(plotly::event_data("plotly_click", source = "cluster.scatter"))
    req(click.data)
    pt = click.data$customdata
    if (!is.null(pt)) {
      updateSelectizeInput(session, "sample.id", selected = pt)
      session$sendCustomMessage("switch-to-tab", "sample")
    }
  }, ignoreInit = TRUE)

  # KNN network computation
  sample.knn.edges = reactive({
    if (exists("sample.knn.edges.df") && is.data.frame(sample.knn.edges.df) && nrow(sample.knn.edges.df) > 0) {
      return(sample.knn.edges.df)
    }
    if (exists("sample.knn.edges") && is.data.frame(sample.knn.edges) && nrow(sample.knn.edges) > 0) {
      return(sample.knn.edges)
    }
    k.nn = 5
    dist.mat = as.matrix(dist(pt.pca.x[, 1:n.pcs.use]))
    edge.list = lapply(1:nrow(dist.mat), function(i) {
      nn = order(dist.mat[i, ])[2:(k.nn + 1)]
      data.frame(
        from = rownames(dist.mat)[i],
        to = rownames(dist.mat)[nn],
        weight = 1 / (1 + dist.mat[i, nn]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, edge.list)
  })

  output$cluster.network = renderVisNetwork({
    df = filtered.coords()
    col.by = input$cluster.color %||% "cluster"
    df[[col.by]][is.na(df[[col.by]])] = "Unknown"

    pal.map = cluster.color.mapping()

    nodes = data.frame(
      id = df$sample,
      label = df$sample,
      color = unname(pal.map[as.character(df[[col.by]])]),
      size = 15,
      title = paste0("Sample: ", df$sample, "<br>", col.by, ": ", df[[col.by]]),
      stringsAsFactors = FALSE
    )

    edges = sample.knn.edges()

    visNetwork::visNetwork(nodes, edges, width = "100%") %>%
      visNetwork::visPhysics(
        solver = "forceAtlas2Based",
        stabilization = list(enabled = TRUE, iterations = 150),
        forceAtlas2Based = list(gravitationalConstant = -50, springLength = 100)
      ) %>%
      visNetwork::visEdges(color = "rgba(150,150,150,0.5)") %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1)) %>%
      visNetwork::visInteraction(hover = TRUE, tooltipDelay = 100) %>%
      visNetwork::visEvents(select = "function(properties) {
        Shiny.setInputValue('cluster.network_selectedNodes', properties.nodes);
      }")
  })

  observeEvent(input$cluster.network_selectedNodes, {
    sel = input$cluster.network_selectedNodes
    if (length(sel) > 0) {
      pt = sel[[1]]
      updateSelectizeInput(session, "sample.id", selected = pt)
      session$sendCustomMessage("switch-to-tab", "sample")
    }
  })

  # Cluster-Specific Edge Attributions (GNN Integrated Gradients)
  get.cluster.edge.attributions = reactive({
    target.cluster = input$cluster.ig.select

    # 1. Use precomputed edge.gradient.summary if available
    if (exists("edge.gradient.summary") && is.data.frame(edge.gradient.summary) && nrow(edge.gradient.summary) > 0 &&
        "cluster" %in% colnames(edge.gradient.summary)) {
      if (is.null(target.cluster) || !target.cluster %in% edge.gradient.summary$cluster) {
        target.cluster = unique(edge.gradient.summary$cluster)[1]
      }
      sub.df = edge.gradient.summary[edge.gradient.summary$cluster == target.cluster, , drop = FALSE]
      if (nrow(sub.df) > 0) {
        return(sub.df)
      }
    }

    # 2. Dynamic multi-channel attribution fallback using full_interactions
    if (!exists("full_interactions") || nrow(full_interactions) == 0) {
      return(data.frame())
    }

    cluster.df = data.frame(sample = sample.ids.vec, cluster = reactive.clusters(), stringsAsFactors = FALSE)
    if (is.null(target.cluster) || !target.cluster %in% cluster.df$cluster) {
      target.cluster = unique(cluster.df$cluster)[1]
    }

    dt.merged = merge(as.data.frame(full_interactions), cluster.df, by = "sample")
    if (nrow(dt.merged) == 0) return(data.frame())

    dt.merged$edge.label = paste0(
      gsub("scRNAseq ", "", dt.merged$sender), " [", dt.merged$source_genesymbol, "] → ",
      gsub("scRNAseq ", "", dt.merged$receiver), " [", dt.merged$target_genesymbol, "]"
    )

    global.means = dt.merged %>%
      dplyr::group_by(edge.label, sender, receiver) %>%
      dplyr::summarize(
        glob.rna  = mean(base.score, na.rm = TRUE),
        glob.path = mean(weight.pathways, na.rm = TRUE),
        glob.prog = mean(weight.programs, na.rm = TRUE),
        glob.gen  = mean(weight.genomic, na.rm = TRUE),
        .groups = "drop"
      )

    clust.dt = dt.merged %>% dplyr::filter(cluster == target.cluster)
    if (nrow(clust.dt) == 0) return(data.frame())

    clust.means = clust.dt %>%
      dplyr::group_by(edge.label, sender, receiver) %>%
      dplyr::summarize(
        c.rna  = mean(base.score, na.rm = TRUE),
        c.path = mean(weight.pathways, na.rm = TRUE),
        c.prog = mean(weight.programs, na.rm = TRUE),
        c.gen  = mean(weight.genomic, na.rm = TRUE),
        .groups = "drop"
      )

    diff.df = merge(clust.means, global.means, by = c("edge.label", "sender", "receiver"))
    diff.df = diff.df %>%
      dplyr::mutate(
        cluster    = target.cluster,
        edge_name  = edge.label,
        ig_rna     = c.rna - glob.rna,
        ig_pathway = c.path - glob.path,
        ig_program = c.prog - glob.prog,
        ig_snv     = (c.gen - glob.gen) * 0.5,
        ig_cnv     = (c.gen - glob.gen) * 0.5,
        ig_score   = abs(ig_rna) + abs(ig_pathway) + abs(ig_program) + abs(c.gen - glob.gen)
      ) %>%
      dplyr::arrange(dplyr::desc(ig_score))

    return(diff.df)
  })

  output$cluster.ig.select.ui = renderUI({
    choices = if (exists("edge.gradient.summary") && is.data.frame(edge.gradient.summary) && nrow(edge.gradient.summary) > 0 && "cluster" %in% colnames(edge.gradient.summary)) {
      sort(unique(edge.gradient.summary$cluster))
    } else {
      sort(unique(reactive.clusters()))
    }
    if (length(choices) == 0) choices = "Cluster 1"
    selectInput("cluster.ig.select", "Patient Cluster:",
                choices = choices, selected = choices[1])
  })

  output$cluster.ig.plot = renderPlotly({
    df = get.cluster.edge.attributions()
    if (is.null(df) || nrow(df) == 0) {
      return(plotly::plot_ly() %>% plotly::layout(title = "No edge attribution data available"))
    }

    top.n = if (!is.null(input$cluster.ig.top.n)) input$cluster.ig.top.n else 15
    barmode = if (!is.null(input$cluster.ig.barmode)) input$cluster.ig.barmode else "relative"

    df.plot = df %>%
      dplyr::arrange(dplyr::desc(ig_score)) %>%
      dplyr::slice_head(n = top.n)

    if (!"edge_label" %in% colnames(df.plot)) {
      df.plot$edge_label = df.plot$edge_name
    }
    df.plot$edge_label = gsub("scRNAseq ", "", df.plot$edge_label)

    # Ascending order for vertical display
    df.plot = df.plot %>%
      dplyr::arrange(ig_score) %>%
      dplyr::mutate(edge_label = factor(edge_label, levels = unique(edge_label)))

    plotly::plot_ly() %>%
      plotly::add_trace(
        data = df.plot, y = ~edge_label, x = ~ig_rna, name = "RNA Expression",
        type = "bar", orientation = "h", marker = list(color = "#2CA02C"),
        hovertemplate = "%{y}<br>RNA: %{x:.4f}<extra></extra>"
      ) %>%
      plotly::add_trace(
        data = df.plot, y = ~edge_label, x = ~ig_pathway, name = "Pathway Activation",
        type = "bar", orientation = "h", marker = list(color = "#FF7F0E"),
        hovertemplate = "%{y}<br>Pathway: %{x:.4f}<extra></extra>"
      ) %>%
      plotly::add_trace(
        data = df.plot, y = ~edge_label, x = ~ig_program, name = "Phenotypic Programs",
        type = "bar", orientation = "h", marker = list(color = "#1F77B4"),
        hovertemplate = "%{y}<br>Program: %{x:.4f}<extra></extra>"
      ) %>%
      plotly::add_trace(
        data = df.plot, y = ~edge_label, x = ~ig_snv, name = "Somatic SNVs",
        type = "bar", orientation = "h", marker = list(color = "#D62728"),
        hovertemplate = "%{y}<br>SNV: %{x:.4f}<extra></extra>"
      ) %>%
      plotly::add_trace(
        data = df.plot, y = ~edge_label, x = ~ig_cnv, name = "Copy Number (CNVs)",
        type = "bar", orientation = "h", marker = list(color = "#9467BD"),
        hovertemplate = "%{y}<br>CNV: %{x:.4f}<extra></extra>"
      ) %>%
      plotly::layout(
        barmode = barmode,
        xaxis   = list(title = "Attribution Score (GNN Integrated Gradients)", titlefont = list(color = "#000000")),
        yaxis   = list(title = "", tickfont = list(size = 9, color = "#000000")),
        margin  = list(l = 280, r = 50, t = 40, b = 50),
        legend  = list(orientation = "h", x = 0, y = 1.12, font = list(color = "#000000"))
      )
  })

  # 3. Sample Overview Tab ----

  # Load per-sample data.table from disk (cached within session)
  sample.dt.cache = reactiveValues(id = NULL, dt = NULL)

  get.sample.dt = reactive({
    req(input$sample.id)
    if (is.null(sample.dt.cache$id) || sample.dt.cache$id != input$sample.id) {
      sample.dt.cache$dt = load.sample.dt(input$sample.id)
      sample.dt.cache$id = input$sample.id
    }
    sample.dt.cache$dt
  })

  observe({
    updateSelectInput(session, "sample.cell.highlight", choices = c("None" = "None", all.populations))
    updateSelectInput(session, "compare.cell.highlight", choices = c("None" = "None", all.populations))
  })

  sample.graph.data = reactive({
    req(input$sample.id)
    withProgress(message = "Loading sample data...", value = 0.3, {
      pt.dt = get.sample.dt()

      lig.search = input$sample.ligand.search
      rec.search = input$sample.receptor.search

      setProgress(0.7, message = "Building network...")
      build.sample.graph(pt.dt, top.n = input$sample.top.n,
                          node.fill.by = input$sample.node.fill,
                          node.outline.by = input$sample.node.outline,
                          edge.color.by = input$sample.edge.color,
                          highlight.cell = input$sample.cell.highlight,
                          search.ligand = lig.search, search.receptor = rec.search)
    })
  })

  output$sample.network = renderVisNetwork({
    gdata = sample.graph.data()
    req(gdata)

    layout.choice = input$sample.network.layout %||% "layout_in_circle"
    layout.choice = input$sample.network.layout %||% "layout_nicely"

    vn = visNetwork::visNetwork(gdata$nodes, gdata$edges,
                                main  = paste("Signaling Network:", input$sample.id),
                                width = "100%")

    if (layout.choice == "physics") {
      vn = vn %>%
        visNetwork::visPhysics(
          solver           = "forceAtlas2Based",
          stabilization    = list(enabled = TRUE, iterations = 150),
          forceAtlas2Based = list(gravitationalConstant = -60, springLength = 120)
        )
    } else {
      vn = vn %>%
        visNetwork::visIgraphLayout(layout = layout.choice, physics = FALSE, smooth = TRUE)
    }

    vn %>%
      visNetwork::visEdges(
        arrows  = list(to = list(enabled = TRUE, scaleFactor = 0.6)),
        smooth  = list(type = "curvedCW", roundness = 0.2),
        scaling = list(min = 1, max = 10),
        color   = list(color = "rgba(150,150,150,0.5)", highlight = "#E64B35FF")
      ) %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
        nodesIdSelection = FALSE
      ) %>%
      visNetwork::visInteraction(hover = TRUE, tooltipDelay = 100, dragNodes = TRUE) %>%
      visNetwork::visLayout(randomSeed = 42) %>%
      visNetwork::visEvents(
        select = "function(properties) {
          Shiny.setInputValue('sample.network_selectedNodes', properties.nodes);
          Shiny.setInputValue('sample.network_selectedEdges', properties.edges);
        }",
        deselectNode = "function(properties) {
          Shiny.setInputValue('sample.network_selectedNodes', properties.nodes);
        }",
        deselectEdge = "function(properties) {
          Shiny.setInputValue('sample.network_selectedEdges', properties.edges);
        }"
      )
  })
  # Node legend helper
  build.node.legend.ui = function(node.fill.input) {
    if (node.fill.input == "hub") {
      grad = "linear-gradient(to right, white, red)"
      tagList(
        div(style = paste0("height: 12px; width: 100%; border-radius: 2px; background: ", grad, "; border: 1px solid #ccc;")),
        div(style = "display: flex; justify-content: space-between; font-size: 11px; margin-top: 2px;",
            span("Low Hub Score"), span("High Hub Score"))
      )
    } else if (node.fill.input == "authority") {
      grad = "linear-gradient(to right, white, blue)"
      tagList(
        div(style = paste0("height: 12px; width: 100%; border-radius: 2px; background: ", grad, "; border: 1px solid #ccc;")),
        div(style = "display: flex; justify-content: space-between; font-size: 11px; margin-top: 2px;",
            span("Low Auth Score"), span("High Auth Score"))
      )
    } else if (node.fill.input == "role") {
      div(style = "font-size: 13px; line-height: 1.6; margin-bottom: 10px;",
        div(tags$span(class="role-bullet", style="background:#E64B35FF;"),
            tags$b("Sender"), " — high out-strength, low in-strength"),
        div(tags$span(class="role-bullet", style="background:#4DBBD5FF;"),
            tags$b("Receiver"), " — high in-strength, low out-strength"),
        div(tags$span(class="role-bullet", style="background:#00A087FF;"),
            tags$b("Influencer"), " — high hub + authority"),
        div(tags$span(class="role-bullet", style="background:#3C5488FF;"),
            tags$b("Mediator"), " — high betweenness")
      )
    } else {
      cell_html = lapply(names(pop.colors), function(p) {
        div(style = "display: inline-block; margin-right: 10px; margin-bottom: 5px; font-size: 11px;",
            span(style = paste0("display: inline-block; width: 10px; height: 10px; background-color: ", pop.colors[p], "; border-radius: 50%; margin-right: 4px; vertical-align: middle;")),
            span(style = "vertical-align: middle;", p))
      })
      div(style = "display: flex; flex-wrap: wrap; max-height: 250px; overflow-y: auto;", cell_html)
    }
  }

  output$sample.node.legend = renderUI({
    build.node.legend.ui(input$sample.node.fill)
  })

  output$sample.node.outline.legend = renderUI({
    req(input$sample.node.outline)
    build.node.legend.ui(input$sample.node.outline)
  })

  output$compare.node.legend = renderUI({
    req(input$compare.node.fill)
    build.node.legend.ui(input$compare.node.fill)
  })

  output$compare.node.outline.legend = renderUI({
    req(input$compare.node.outline)
    build.node.legend.ui(input$compare.node.outline)
  })

  output$sample.edge.legend = renderUI({
    gdata = sample.graph.data()
    req(gdata, gdata$legend)
    leg = gdata$legend

    if (leg$type == "discrete") {
      # Render list of color blocks
      items = lapply(names(leg$pal), function(nm) {
        tags$div(style = "margin-bottom: 3px; font-size: 13px;",
          tags$span(class="role-bullet", style=paste0("background:", leg$pal[[nm]], ";")),
          nm
        )
      })
      do.call(tagList, items)
    } else if (leg$type == "continuous") {
      # Render a gradient bar
      grad = paste0("linear-gradient(to right, ", leg$pal[1], ", ", leg$pal[length(leg$pal)], ")")
      tagList(
        div(style = paste0("height: 12px; width: 100%; border-radius: 2px; background: ", grad, ";")),
        div(style = "display: flex; justify-content: space-between; font-size: 11px; margin-top: 2px;",
            span("Low"), span("High"))
      )
    } else {
      tags$i(style = "font-size:12px; color:#888;", "No legend available.")
    }
  })

  output$sample.metrics.dt = renderDT({
    gdata = sample.graph.data()
    req(gdata)

    df = gdata$metrics %>%
      dplyr::mutate(node = gsub("scRNAseq ", "", node)) %>%
      dplyr::arrange(dplyr::desc(hub.score))

    hlt = input$sample.cell.highlight
    if (!is.null(hlt) && hlt != "None") {
      df = df[df$node == hlt, , drop = FALSE]
    }

    sel_nodes = input$sample.network_selectedNodes
    if (length(sel_nodes) > 0) {
      df = df[df$node %in% gsub("scRNAseq ", "", sel_nodes), , drop = FALSE]
    }

    DT::datatable(df,
      rownames  = FALSE,
      selection = "none",
      colnames  = c("Cell Type", "Role", "Out-Strength", "In-Strength",
                    "Betweenness", "Hub Score", "Authority Score"),
      options   = list(pageLength = 10, scrollX = TRUE, dom = "ftp")
    ) %>%
      DT::formatStyle("role",
        backgroundColor = DT::styleEqual(
          c("Sender", "Receiver", "Influencer", "Mediator"),
          c("#E64B3530", "#4DBBD530", "#00A08730", "#3C548830")
        )
      ) %>%
      DT::formatStyle(
        c("out.strength", "in.strength", "betweenness", "hub.score", "authority.score"),
        background = DT::styleColorBar(range(c(0, 1)), '#3c8dbc50'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })

  output$sample.lr.dt = renderDT({
    req(input$sample.id)
    n.show = input$sample.lr.n

    pt.dt = get.sample.dt()[order(-final.score)]
    if (input$sample.lr.sender   != "all") pt.dt = pt.dt[sender   == input$sample.lr.sender]
    if (input$sample.lr.receiver != "all") pt.dt = pt.dt[receiver == input$sample.lr.receiver]

    hlt = input$sample.cell.highlight
    if (!is.null(hlt) && hlt != "None") {
      pt.dt = pt.dt[sender == hlt | receiver == hlt]
    }

    sel_edges = input$sample.network_selectedEdges
    if (length(sel_edges) > 0) {
      pt.dt = pt.dt[paste0(sender, "_", receiver) %in% sel_edges]
    }

    df = as.data.frame(pt.dt[1:min(.N, n.show)]) %>%
      dplyr::mutate(
        Ligand     = gene.link.html(source_genesymbol),
        Receptor   = gene.link.html(target_genesymbol),
        Sender     = gsub("scRNAseq ", "", sender),
        Receiver   = gsub("scRNAseq ", "", receiver),
        `Final Score` = round(final.score, 4),
        `Base Score`  = round(base.score, 4),
        `Gen. Weight` = round(weight.genomic, 2),
        `Path. Weight`= round(weight.pathways, 2),
        `Prog. Weight`= round(weight.programs, 2),
        `Ligand Axis` = omnipath.link.html(source_genesymbol, target_genesymbol, ligand.axis),
        `Receptor Axis` = pathway.link.html(receptor.axis)
      ) %>%
      dplyr::select(Ligand, Receptor, Sender, Receiver,
                    `Final Score`, `Base Score`, `Gen. Weight`, `Path. Weight`, `Prog. Weight`,
                    `Ligand Axis`, `Receptor Axis`)

    DT::datatable(df,
      escape    = FALSE,
      rownames  = FALSE,
      selection = "none",
      options   = list(
        pageLength = 15, scrollX = TRUE,
        order = list(list(4, "desc")), dom = "ftip"
      )
    ) %>%
      DT::formatStyle("Gen. Weight",
        backgroundColor = DT::styleInterval(1, c("white", "#E64B3520"))
      ) %>%
      DT::formatStyle(
        c("Final Score", "Base Score", "Gen. Weight", "Path. Weight", "Prog. Weight"),
        background = DT::styleColorBar(range(c(0, max(df$`Final Score`, na.rm=TRUE))), '#00a65a50'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })

  output$sample.sankey.ui = renderUI({
    hi_cell = input$sample.cell.highlight
    sel_nodes = input$sample.network_selectedNodes
    has_selection = (hi_cell != "None") || (length(sel_nodes) > 0)

    if (!has_selection) {
      return(p(style = "color:#888; font-size:13px; text-align:center; margin-top:200px;",
               "Select a cell type in the controls to the left or click a node in the network to view its specific ligand-receptor flow."))
    }
    withSpinner(plotlyOutput("sample.sankey", height = "520px"))
  })

  output$sample.sankey = renderPlotly({
    hi_cell = input$sample.cell.highlight
    sel_nodes = input$sample.network_selectedNodes
    has_selection = (hi_cell != "None") || (length(sel_nodes) > 0)
    req(has_selection)

    pt.dt = get.sample.dt()[order(-final.score)]

    targets = c()
    if (hi_cell != "None") targets = c(targets, hi_cell)
    if (length(sel_nodes) > 0) targets = c(targets, sel_nodes)

    # Filter only interactions involving the highlighted cell types
    pt.dt = pt.dt[sender %in% targets | receiver %in% targets]
    if (nrow(pt.dt) == 0) return(plotly::plot_ly() %>% plotly::layout(annotations = list(text="No interactions", showarrow=FALSE)))

    df_sankey = head(pt.dt, 50) %>%
      dplyr::mutate(
        lig_full = paste0(gsub("scRNAseq ", "", sender), " : ", source_genesymbol),
        rec_full = paste0(target_genesymbol, " : ", gsub("scRNAseq ", "", receiver))
      )

    nodes = unique(c(df_sankey$lig_full, df_sankey$rec_full))

    plotly::plot_ly(
      type = "sankey",
      orientation = "h",
      node = list(
        label = nodes,
        pad = 10,
        thickness = 20,
        line = list(color = "black", width = 0.5)
      ),
      link = list(
        source = match(df_sankey$lig_full, nodes) - 1,
        target = match(df_sankey$rec_full, nodes) - 1,
        value  = abs(df_sankey$final.score),
        customdata = paste0("Score: ", round(df_sankey$final.score, 3)),
        hovertemplate = "%{source.label} &rarr; %{target.label}<br>%{customdata}<extra></extra>"
      )
    ) %>%
      plotly::layout(margin = list(t = 30, b = 20, l = 10, r = 10))
  })

  output$sample.clinical.ui = renderUI({
    req(input$sample.id)
    clin = clinical.sample %>%
      dplyr::filter(SPECTRUM_SEQ == input$sample.id)

    if (nrow(clin) == 0) {
      return(p(style = "color:#888; font-size:12px;",
               "No clinical data available for this sample."))
    }
    clin = clin[1, ]

    tags$table(style = "width:100%; font-size:12.5px;",
      tags$tr(tags$td(tags$b("Sample ID")),   tags$td(clin$PUBLIC_ID)),
      tags$tr(tags$td(tags$b("Sex")),           tags$td(clin$DEMOG_GENDER)),
      tags$tr(tags$td(tags$b("Age")),           tags$td(clin$DEMOG_SAMPLEAGE)),
      tags$tr(tags$td(tags$b("Best Response")), tags$td(clin$best.response)),
      tags$tr(tags$td(tags$b("PFS")),
              tags$td(paste0(clin$derived.pfs.day, "d (event=",
                              clin$censpfs, ")"))),
      tags$tr(tags$td(tags$b("del(17p)")),   tags$td(as.character(clin$chr17p_Loss))),
      tags$tr(tags$td(tags$b("gain(1q)")),   tags$td(as.character(clin$chr1q_Gain))),
      tags$tr(tags$td(tags$b("t(4;14)")),    tags$td(as.character(clin$t4_14_MMSET))),
      tags$tr(tags$td(tags$b("t(11;14)")),   tags$td(as.character(clin$t11_14_CCND1))),
      tags$tr(tags$td(tags$b("MYC")),        tags$td(as.character(clin$t_MYC)))
    )
  })

  # 3b. Sample Comparison Tab ----

  output$compare.cohort.cluster.ui = renderUI({
    k = input$k.clusters
    choices = c("All", paste("Cluster", 1:k))
    clusters = reactive.clusters()
    u.clust = if (!is.null(clusters)) sort(unique(clusters)) else NULL
    choices = if (length(u.clust) > 0) c("All", u.clust) else "All"
    selectInput("cohort.cluster", "Cluster:", choices = choices, selected = "All")
  })

  output$cohort.filter.ligand.ui = renderUI({
    choices = sort(unique(sample.stats.summary$`Top Ligand`))
    selectizeInput("cohort.filter.ligand", "Top Ligand:", choices = choices, multiple = TRUE, options = list(placeholder = "Search..."))
  })

  output$cohort.filter.receptor.ui = renderUI({
    choices = sort(unique(sample.stats.summary$`Top Receptor`))
    selectizeInput("cohort.filter.receptor", "Top Receptor:", choices = choices, multiple = TRUE, options = list(placeholder = "Search..."))
  })

  output$cohort.filter.pathway.ui = renderUI({
    choices = sort(unique(sample.stats.summary$`Top Pathway`))
    selectizeInput("cohort.filter.pathway", "Top Pathway:", choices = choices, multiple = TRUE, options = list(placeholder = "Search..."))
  })

  output$cohort.filter.program.ui = renderUI({
    choices = sort(unique(sample.stats.summary$`Top Program`))
    selectizeInput("cohort.filter.program", "Top Program:", choices = choices, multiple = TRUE, options = list(placeholder = "Search..."))
  })

  output$cohort.filter.celltype.ui = renderUI({
    choices = sort(unique(c(full_interactions$sender, full_interactions$receiver)))
    choices = gsub("scRNAseq ", "", choices)
    selectizeInput("cohort.filter.celltype", "Cell Type:", choices = choices, multiple = TRUE, options = list(placeholder = "Search..."))
  })

  compare.cohort.data = reactive({
    if (!exists("sample.stats.summary")) return(NULL)

    df = sample.stats.summary %>%
      dplyr::rename(Sample = sample) %>%
      dplyr::select(
        Sample,
        `Number of interactions`,
        `Top Ligand`,
        `Top Receptor`,
        `Top Cell-Cell Interaction`,
        `Top Program`,
        `Top Pathway`,
        `Mean Edge Weight`
      )

    # Filter by Top Ligand
    if (length(input$cohort.filter.ligand) > 0) {
      df = df %>% dplyr::filter(`Top Ligand` %in% input$cohort.filter.ligand)
    }
    # Filter by Top Receptor
    if (length(input$cohort.filter.receptor) > 0) {
      df = df %>% dplyr::filter(`Top Receptor` %in% input$cohort.filter.receptor)
    }
    # Filter by Top Pathway
    if (length(input$cohort.filter.pathway) > 0) {
      df = df %>% dplyr::filter(`Top Pathway` %in% input$cohort.filter.pathway)
    }
    # Filter by Top Program
    if (length(input$cohort.filter.program) > 0) {
      df = df %>% dplyr::filter(`Top Program` %in% input$cohort.filter.program)
    }
    # Filter by Cell Type (matching sender or receiver in Top Cell-Cell Interaction)
    if (length(input$cohort.filter.celltype) > 0) {
      df = df %>% dplyr::filter(grepl(paste(input$cohort.filter.celltype, collapse="|"), `Top Cell-Cell Interaction`))
    }

    df %>%
      dplyr::mutate(
        `Top Ligand` = gene.link.html(`Top Ligand`),
        `Top Receptor` = gene.link.html(`Top Receptor`),
        `Top Program` = pathway.link.html(`Top Program`)
      )
  })


  output$compare.cohort.dt = renderDT({
    df = compare.cohort.data()
    if (is.null(df)) return(NULL)

    # Add two empty rows to eliminate excess white space underneath
    empty_rows = as.data.frame(matrix(NA, nrow = 2, ncol = ncol(df)))
    colnames(empty_rows) = colnames(df)
    df = rbind(df, empty_rows)
    DT::datatable(df,
      escape = FALSE,
      rownames = FALSE,
      selection = "multiple",
      options = list(pageLength = 12, scrollX = TRUE, dom = "ftp")
    ) %>%
      DT::formatStyle("Mean Edge Weight",
        background = DT::styleColorBar(range(c(0, max(df$`Mean Edge Weight`, na.rm=TRUE))), '#00a65a50'),
        backgroundSize = '98% 88%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })

  output$compare.legend = renderUI({
    pts = input$compare.sample.ids
    if (is.null(pts) || length(pts) == 0) return(tags$i(style = "font-size:12px; color:#888;", "Select samples to view legend."))

    pt = pts[1]
    pt.dt = load.sample.dt(pt)
    if (is.null(pt.dt) || nrow(pt.dt) == 0) return(NULL)

    lig.search = input$compare.ligand.search
    rec.search = input$compare.receptor.search

    gdata = build.sample.graph(pt.dt, top.n = input$compare.top.n,
                                node.fill.by = input$compare.node.fill,
                                node.outline.by = input$compare.node.outline,
                                edge.color.by = input$compare.edge.color,
                                highlight.cell = input$compare.cell.highlight,
                                search.ligand = lig.search, search.receptor = rec.search)
    req(gdata, gdata$legend)
    leg = gdata$legend

    if (leg$type == "discrete") {
      items = lapply(names(leg$pal), function(nm) {
        tags$div(style = "margin-bottom: 3px; font-size: 13px;",
          tags$span(class="role-bullet", style=paste0("background:", leg$pal[[nm]], ";")),
          nm
        )
      })
      do.call(tagList, items)
    } else if (leg$type == "continuous") {
      grad = paste0("linear-gradient(to right, ", leg$pal[1], ", ", leg$pal[length(leg$pal)], ")")
      tagList(
        div(style = paste0("height: 12px; width: 100%; border-radius: 2px; background: ", grad, ";")),
        div(style = "display: flex; justify-content: space-between; font-size: 11px; margin-top: 2px;",
            span("Low"), span("High"))
      )
    } else {
      tags$i(style = "font-size:12px; color:#888;", "No legend available.")
    }
  })

  observeEvent(input$compare.cohort.dt_rows_selected, {
    df = compare.cohort.data()
    sel_rows = input$compare.cohort.dt_rows_selected
    if (length(sel_rows) > 0) {
      sel_pts = df$Sample[sel_rows]
      if (length(sel_pts) > 12) sel_pts = sel_pts[1:12]
      updateSelectizeInput(session, "compare.sample.ids", selected = sel_pts)
    } else {
      updateSelectizeInput(session, "compare.sample.ids", selected = character(0))
    }
  }, ignoreNULL = FALSE)

  output$compare.network.ui = renderUI({
    pts = input$compare.sample.ids
    if (is.null(pts) || length(pts) == 0) {
      return(tags$div(style = "padding: 50px; text-align: center; color: #888;",
                      "Select at least one sample to begin comparison."))
    }

    n_pts = length(pts)
    col_width = if (n_pts == 1) {
      12
    } else if (n_pts == 2 || n_pts == 4) {
      6
    } else {
      4
    }

    # Generate UI elements
    plot_ui_list = lapply(seq_along(pts), function(i) {
      pt = pts[i]
      ns_id = paste0("compare_net_", i)
      column(col_width,
        box(title = paste("Sample:", pt), status = "primary", solidHeader = TRUE, width = 12,
            withSpinner(visNetworkOutput(ns_id, height = "400px"))
        )
      )
    })

    do.call(fluidRow, plot_ui_list)
  })

  observe({
    pts = input$compare.sample.ids
    if (is.null(pts) || length(pts) == 0) return()

    lapply(seq_along(pts), function(i) {
      local({
        my_i = i
        pt = pts[my_i]
        ns_id = paste0("compare_net_", my_i)

        output[[ns_id]] = renderVisNetwork({
          pt.dt = load.sample.dt(pt)
          if (is.null(pt.dt) || nrow(pt.dt) == 0) return(NULL)

          lig.search = input$compare.ligand.search
          rec.search = input$compare.receptor.search

          gdata = build.sample.graph(pt.dt, top.n = input$compare.top.n,
                                      node.fill.by = input$compare.node.fill,
                                      node.outline.by = input$compare.node.outline,
                                      edge.color.by = input$compare.edge.color,
                                      highlight.cell = input$compare.cell.highlight,
                                      search.ligand = lig.search, search.receptor = rec.search)
          if (is.null(gdata)) return(NULL)

          layout.choice = input$compare.network.layout %||% "layout_in_circle"
          layout.choice = input$compare.network.layout %||% "layout_nicely"

          vn = visNetwork::visNetwork(gdata$nodes, gdata$edges, width = "100%")

          if (layout.choice == "physics") {
            vn = vn %>%
              visNetwork::visPhysics(
                solver           = "forceAtlas2Based",
                stabilization    = list(enabled = TRUE, iterations = 150),
                forceAtlas2Based = list(gravitationalConstant = -60, springLength = 120)
              )
          } else {
            vn = vn %>%
              visNetwork::visIgraphLayout(layout = layout.choice, physics = FALSE, smooth = TRUE)
          }

          vn %>%
            visNetwork::visEdges(
              arrows  = list(to = list(enabled = TRUE, scaleFactor = 0.6)),
              smooth  = list(type = "curvedCW", roundness = 0.2),
              scaling = list(min = 1, max = 10)
            ) %>%
            visNetwork::visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1), nodesIdSelection = FALSE) %>%
            visNetwork::visInteraction(hover = TRUE, tooltipDelay = 100, dragNodes = TRUE) %>%
            visNetwork::visLayout(randomSeed = 42)
        })
      })
    })
  })

  # 4. Genomic Drivers Tab ----

  output$genomic.circos.grid = renderUI({
    axes = ligand.axes.available
    req(axes)

    plot.ids = paste0("circos_plot_", gsub("[^a-zA-Z0-9]", "_", axes))
    n.axes   = length(axes)
    n.cols   = if (n.axes <= 2) 2 else if (n.axes <= 4) 2 else 3
    col.width = 12 %/% n.cols

    panels = lapply(seq_along(axes), function(i) {
      column(col.width,
        div(class = "circos-panel",
          plotOutput(plot.ids[i], height = "380px")
        )
      )
    })

    # Organize into rows of n.cols
    n.rows = ceiling(n.axes / n.cols)
    rows   = lapply(seq_len(n.rows), function(r) {
      idx = ((r - 1) * n.cols + 1):min(r * n.cols, n.axes)
      fluidRow(panels[idx])
    })
    tagList(rows)
  })

  # Create renderPlot for each axis dynamically
  observe({
    axes     = ligand.axes.available
    top.n    = input$genomic.top.n
    hi.genes = input$genomic.highlight.genes
    req(axes)

    lapply(axes, function(ax) {
      local({
        local.ax    = ax
        local.top.n = top.n
        local.genes = hi.genes
        plot.id     = paste0("circos_plot_", gsub("[^a-zA-Z0-9]", "_", local.ax))

        output[[plot.id]] = renderPlot({
          draw.circos.chord(local.ax,
                            top.n          = local.top.n,
                            highlight.genes = local.genes)
        }, bg = "white")
      })
    })
  })

  # 5. Pathway Activity Map Tab ----

  # Update cluster choices dynamically
  observe({
    k = as.integer(input$cluster.k)
    clusters = reactive.clusters()
    req(clusters)
    u.clust = sort(unique(clusters))
    updateSelectInput(session, "pathway.cluster.subset",
                      choices = c("All" = "All", setNames(u.clust, u.clust)))
  })

  output$pathway.heatmap = renderPlotly({
    score_type = input$pathway.score.type
    subset_cluster = input$pathway.cluster.subset

    score_list = if (score_type == "pathway") modality.progeny.scores else functional.program.scores
    mean_data = if (score_type == "pathway") progeny.mean.data else program.mean.data

    if (!is.null(subset_cluster) && subset_cluster != "All") {
      clusters = reactive.clusters()
      pt_ids = sample.ids.vec[clusters == subset_cluster]

      plot.df = lapply(names(score_list), function(pop) {
        mat = score_list[[pop]]
        if (is.null(mat) || nrow(mat) == 0) return(NULL)
        mat_sub = mat[rownames(mat) %in% pt_ids, , drop = FALSE]
        if (nrow(mat_sub) == 0) return(NULL)
        means = colMeans(mat_sub, na.rm = TRUE)
        tibble::tibble(population = pop, pathway = names(means), mean.score = as.numeric(means))
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::filter(!population %in% c("CNV", "SNV"))
    } else {
      plot.df = mean_data %>% dplyr::rename(pathway = if(score_type == "pathway") "pathway" else "program")
    }

    mat.wide = plot.df %>%
      tidyr::pivot_wider(names_from = pathway, values_from = mean.score, values_fill = 0)
    mat.vals = as.matrix(mat.wide[, -1])
    rownames(mat.vals) = mat.wide$population
    
    # Safeguard against NA/NaN for hclust
    mat.vals[is.na(mat.vals) | is.nan(mat.vals)] = 0

    if (nrow(mat.vals) > 2) {
      row.order = hclust(dist(mat.vals))$order
      mat.vals  = mat.vals[row.order, , drop = FALSE]
    }
    if (ncol(mat.vals) > 2) {
      col.order = hclust(dist(t(mat.vals)))$order
      mat.vals  = mat.vals[, col.order, drop = FALSE]
    }

    plotly::plot_ly(
      x         = colnames(mat.vals),
      y         = rownames(mat.vals),
      z         = mat.vals,
      type      = "heatmap",
      colorscale = "RdBu",
      reversescale = TRUE,
      colorbar  = list(title = "Mean Activity<br>Score"),
      hovertemplate = "Population: %{y}<br>Pathway: %{x}<br>Mean Score: %{z:.3f}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis  = list(title = "Pathway", tickangle = -40),
        yaxis  = list(title = "", tickfont = list(size = 9)),
        margin = list(l = 230, b = 80)
      )
  })

  output$pathway.violin = renderPlotly({
    score_type = input$pathway.score.type
    subset_cluster = input$pathway.cluster.subset

    score_list = if (score_type == "pathway") modality.progeny.scores else functional.program.scores
    mean_data = if (score_type == "pathway") progeny.mean.data else program.mean.data

    if (!is.null(subset_cluster) && subset_cluster != "All") {
      clusters = reactive.clusters()
      pt_ids = sample.ids.vec[clusters == subset_cluster]
      plot.df = lapply(names(score_list), function(pop) {
        mat = score_list[[pop]]
        if (is.null(mat) || nrow(mat) == 0) return(NULL)
        mat_sub = mat[rownames(mat) %in% pt_ids, , drop = FALSE]
        if (nrow(mat_sub) == 0) return(NULL)
        means = colMeans(mat_sub, na.rm = TRUE)
        tibble::tibble(population = pop, pathway = names(means), mean.score = as.numeric(means))
      }) %>%
        dplyr::bind_rows() %>%
        dplyr::filter(!population %in% c("CNV", "SNV"))
    } else {
      plot.df = mean_data %>% dplyr::rename(pathway = if(score_type == "pathway") "pathway" else "program")
    }

    pathway_medians = plot.df %>%
      dplyr::group_by(pathway) %>%
      dplyr::summarize(med = median(mean.score, na.rm = TRUE), .groups = 'drop') %>%
      dplyr::arrange(med)

    # Map medians to RdBu colors
    med_vals = pathway_medians$med
    pal = colorRampPalette(c("royalblue", "firebrick"))(100)
    med_norm = as.numeric(cut(med_vals, breaks = 100))
    pathway_colors = setNames(pal[med_norm], pathway_medians$pathway)

    plotly::plot_ly(plot.df,
      x     = ~pathway,
      y     = ~mean.score,
      type  = "box",
      color = ~pathway,
      colors = pathway_colors,
      boxpoints = "all", jitter = 0.3, pointpos = 0,
      hovertemplate = "Pathway: %{x}<br>Score: %{y:.3f}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis    = list(title = if(score_type == "pathway") "PROGENy Pathway" else "Biological Program"),
        yaxis    = list(title = "Mean Score"),
        showlegend = FALSE
      )
  })

  # 6. Clinical Results Tab ----

  # Clinical data joined with reactive clusters
  clinical.surv.data = reactive({
    clusters = reactive.clusters()
    clust.df = data.frame(sample = sample.ids.vec, cluster = clusters,
                           stringsAsFactors = FALSE)
    survival.df %>%
      dplyr::select(-dplyr::any_of("cluster")) %>%
      dplyr::left_join(clust.df, by = "sample")
  })

  output$clinical.km = renderPlot({
    df        = clinical.surv.data()
    endpoint  = input$clinical.endpoint
    group.col = input$clinical.group

    if (endpoint == "pfs") {
      df = df %>%
        dplyr::filter(!is.na(pfs.months), !is.na(censpfs)) %>%
        dplyr::mutate(surv.time = pfs.months, surv.event = as.integer(censpfs))
      time.lab = "Progression-Free Survival (months)"
    } else {
      df = df %>%
        dplyr::filter(!is.na(os.months), !is.na(censos)) %>%
        dplyr::mutate(surv.time = os.months, surv.event = as.integer(censos))
      time.lab = "Overall Survival (months)"
    }

    if (!group.col %in% colnames(df)) {
      plot.new(); text(0.5, 0.5, "Grouping variable not found.", cex = 1.2); return()
    }

    df$strata  = factor(df[[group.col]])
    surv.fit   = survival::survfit(Surv(surv.time, surv.event) ~ strata, data = df)

    n.groups   = nlevels(df$strata)
    pal.colors = if (n.groups <= 10) ggsci::pal_npg()(10)[1:n.groups] else grDevices::rainbow(n.groups)

    survminer::ggsurvplot(
      surv.fit,
      data         = df,
      palette      = pal.colors,
      conf.int     = input$clinical.km.ci,
      risk.table   = input$clinical.km.risk,
      pval         = TRUE,
      pval.method  = TRUE,
      xlab         = time.lab,
      ylab         = "Survival Probability",
      legend.title = clinical.var.labels[group.col] %||% group.col,
      ggtheme      = ggplot2::theme_classic(base_size = 13),
      tables.theme = survminer::theme_cleantable(),
      risk.table.y.text.col = TRUE
    )
  })

  cox.model.result = eventReactive(input$clinical.run.cox, {
    df       = clinical.surv.data()
    endpoint = input$clinical.endpoint
    covars   = input$clinical.cox.covars

    if (endpoint == "pfs") {
      df = df %>%
        dplyr::filter(!is.na(pfs.months), !is.na(censpfs)) %>%
        dplyr::mutate(surv.time = pfs.months, surv.event = as.integer(censpfs))
    } else {
      df = df %>%
        dplyr::filter(!is.na(os.months), !is.na(censos)) %>%
        dplyr::mutate(surv.time = os.months, surv.event = as.integer(censos))
    }

    df$Cluster = relevel(factor(df$cluster), ref = "Cluster 1")
    df$Age = df$DEMOG_SAMPLEAGE
    df$Sex = df$DEMOG_GENDER
    df$`del(17p)` = df$chr17p_Loss
    df$`gain(1q)` = df$chr1q_Gain
    df$`t(4;14)` = df$t4_14_MMSET
    df$`t(14;16)` = df$t14_16_MAF
    df$`t(14;20)` = df$t14_20_MAFB
    df$`t(11;14)` = df$t11_14_CCND1
    df$MYC = df$t_MYC

    clean_map = c(
      "DEMOG_SAMPLEAGE" = "Age",
      "DEMOG_GENDER" = "Sex",
      "chr17p_Loss" = "del(17p)",
      "chr1q_Gain" = "gain(1q)",
      "t4_14_MMSET" = "t(4;14)",
      "t14_16_MAF" = "t(14;16)",
      "t14_20_MAFB" = "t(14;20)",
      "t11_14_CCND1" = "t(11;14)",
      "t_MYC" = "MYC"
    )
    mapped_covars = unname(clean_map[covars])
    mapped_covars = mapped_covars[!is.na(mapped_covars)]

    rhs.terms = c("Cluster", mapped_covars)

    # Use plain Surv() — survival:: prefix in formula strings breaks ggforest
    formula = stats::as.formula(
      paste("Surv(surv.time, surv.event) ~", paste(paste0("`", rhs.terms, "`"), collapse = " + "))
    )

    model = tryCatch(survival::coxph(formula, data = df), error = function(e) e$message)
    if (is.character(model)) {
      list(model = NULL, data = df, error = model)
    } else {
      list(model = model, data = df, error = NULL)
    }
  }, ignoreNULL = FALSE)

  output$clinical.forest = renderPlot({
    res = cox.model.result()
    if (!is.null(res$error)) {
      plot.new()
      text(0.5, 0.5, paste("Cox Model failed:\n\n", res$error, "\n\nToo many covariates selected for the number of samples/events."), cex = 1.1)
      return()
    }
    if (is.null(res$model)) {
      plot.new(); text(0.5, 0.5, "Click 'Fit Cox Model' to run analysis.", cex = 1.2); return()
    }
    survminer::ggforest(res$model, data = res$data,
                        main = "Hazard Ratios — Cox PH Model", fontsize = 1.1)
  })



  # 7. Signal Diagnostics Tab ----

  total.n   = if (exists("sender.receiver.agg") && "n.pairs" %in% colnames(sender.receiver.agg)) sum(as.numeric(sender.receiver.agg$n.pairs), na.rm = TRUE) else 1
  if (total.n == 0 || is.na(total.n)) total.n = 1
  n.boosted = if (exists("weight.direction.summary") && "Boosted" %in% colnames(weight.direction.summary)) sum(weight.direction.summary$Boosted, na.rm = TRUE) else 0
  n.supp    = if (exists("weight.direction.summary") && "Suppressed" %in% colnames(weight.direction.summary)) sum(weight.direction.summary$Suppressed, na.rm = TRUE) else 0
  n.neutral = if (exists("weight.direction.summary") && "Neutral" %in% colnames(weight.direction.summary)) sum(weight.direction.summary$Neutral, na.rm = TRUE) else 0

  output$diag.vbox.boosted = renderValueBox({
    valueBox(
      value    = paste0(round(n.boosted / total.n * 100, 1), "%"),
      subtitle = "Boosted Interactions",
      icon     = icon("arrow-up"),
      color    = "red"
    )
  })
  output$diag.vbox.suppressed = renderValueBox({
    valueBox(
      value    = paste0(round(n.supp / total.n * 100, 1), "%"),
      subtitle = "Suppressed Interactions",
      icon     = icon("arrow-down"),
      color    = "blue"
    )
  })
  output$diag.vbox.genomic = renderValueBox({
    pct_gen = if (exists("weight.layer.summary") && !is.null(weight.layer.summary$pct.genomic.weighted)) weight.layer.summary$pct.genomic.weighted else 0
    valueBox(
      value    = paste0(round(pct_gen, 1), "%"),
      subtitle = "Genomically Weighted",
      icon     = icon("dna"),
      color    = "purple"
    )
  })

  output$diag.weight.dir = renderPlotly({
    if (is.null(weight.direction.summary) || nrow(weight.direction.summary) == 0) {
      return(plotly::plot_ly() %>% plotly::layout(title = "No weight direction data available"))
    }
    df = weight.direction.summary %>%
      dplyr::mutate(receiver = gsub("scRNAseq ", "", receiver)) %>%
      tidyr::pivot_longer(tidyselect::any_of(c("Boosted", "Neutral", "Suppressed")),
                          names_to = "direction", values_to = "n") %>%
      dplyr::group_by(receiver) %>%
      dplyr::mutate(pct = n / sum(n) * 100) %>%
      dplyr::ungroup()

    plotly::plot_ly(df,
      x           = ~pct,
      y           = ~reorder(receiver, pct * (direction == "Boosted")),
      color       = ~direction,
      colors      = c(Boosted = "#E64B35FF", Neutral = "#B0B7BE", Suppressed = "#4DBBD5FF"),
      type        = "bar", orientation = "h",
      hovertemplate = "%{y}<br>%{data.name}: %{x:.1f}%<extra></extra>"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis   = list(title = "Proportion (%)", range = c(0, 100), titlefont = list(color = "#000000")),
        yaxis   = list(title = "", tickfont = list(size = 9, color = "#000000")),
        margin  = list(l = 200),
        legend  = list(title = list(text = "Direction", font = list(color = "#000000")), font = list(color = "#000000"))
      )
  })

  output$diag.rank.shift = renderPlotly({
    df = if (exists("rank.shift.summary") && is.data.frame(rank.shift.summary) && nrow(rank.shift.summary) > 0 &&
             !all(abs(rank.shift.summary$mean.rank.shift) < 1e-6)) {
      rank.shift.summary
    } else if (exists("full_interactions") && nrow(full_interactions) > 0 &&
               all(c("sample", "receiver", "base.score", "final.score") %in% colnames(full_interactions))) {
      full_interactions %>%
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
    } else if (exists("rank.shift.summary") && is.data.frame(rank.shift.summary) && nrow(rank.shift.summary) > 0) {
      rank.shift.summary
    } else {
      data.frame()
    }

    if (nrow(df) == 0) return(plotly::plot_ly() %>% plotly::layout(title = "No rank shift data available"))

    # If all mean.rank.shift are exactly 0, use mean absolute rank displacement
    use.abs = all(abs(df$mean.rank.shift) < 1e-6)
    plot.val = if (use.abs) df$mean.abs.rank.shift else df$mean.rank.shift
    x.title = if (use.abs) "Mean Absolute Rank Displacement" else "Mean Rank Shift (final vs base)"

    df = df %>%
      dplyr::mutate(
        receiver = gsub("scRNAseq ", "", receiver),
        metric.val = plot.val,
        bar.color = if (use.abs) "#E64B35FF" else ifelse(is.na(metric.val) | metric.val < 0, "#4DBBD5FF", "#E64B35FF")
      ) %>%
      dplyr::arrange(metric.val) %>%
      dplyr::mutate(receiver = factor(receiver, levels = unique(receiver)))

    plotly::plot_ly(df,
      x             = ~metric.val,
      y             = ~receiver,
      type          = "bar", orientation = "h",
      marker        = list(color = df$bar.color),
      text          = ~paste0(ifelse(metric.val > 0 & !use.abs, "+", ""), round(metric.val, 2), " (", round(pct.shifted, 1), "% shifted)"),
      textposition  = "outside",
      hovertemplate = "%{y}<br>Rank shift: %{x:.4f}<extra></extra>"
    ) %>%
      plotly::layout(
        xaxis  = list(title = x.title, titlefont = list(color = "#000000")),
        yaxis  = list(title = "", tickfont = list(size = 9, color = "#000000")),
        margin = list(l = 200, r = 100, t = 20, b = 50)
      )
  })

  output$diag.signal.health = renderPlotly({
    df = if (exists("mean.layer.weights.summary") && is.data.frame(mean.layer.weights.summary) && nrow(mean.layer.weights.summary) > 0) {
      mean.layer.weights.summary
    } else if (exists("full_interactions") && nrow(full_interactions) > 0 &&
               all(c("receiver", "weight.pathways", "weight.programs", "weight.genomic") %in% colnames(full_interactions))) {
      full_interactions %>%
        dplyr::group_by(receiver) %>%
        dplyr::summarize(
          Pathways = mean(weight.pathways, na.rm = TRUE),
          Programs = mean(weight.programs, na.rm = TRUE),
          Genomics = mean(weight.genomic, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      data.frame()
    }

    if (nrow(df) == 0) return(plotly::plot_ly() %>% plotly::layout(title = "No layer weight data available"))

    df = df %>%
      dplyr::mutate(receiver = gsub("scRNAseq ", "", receiver))

    df.long = df %>%
      tidyr::pivot_longer(cols = tidyselect::any_of(c("Pathways", "Programs", "Genomics")),
                          names_to = "Layer", values_to = "Weight") %>%
      dplyr::mutate(Layer = factor(Layer, levels = c("Pathways", "Programs", "Genomics")))

    layer.colors = c(Pathways = "#00A087FF", Programs = "#4DBBD5FF", Genomics = "#3C5488FF")

    plotly::plot_ly(df.long,
      x             = ~Weight,
      y             = ~receiver,
      color         = ~Layer,
      colors        = layer.colors,
      type          = "bar", orientation = "h",
      hovertemplate = "%{y}<br>%{data.name}: %{x:.4f}<extra></extra>"
    ) %>%
      plotly::layout(
        barmode = "stack",
        xaxis   = list(title = "Mean Weight Contribution", titlefont = list(color = "#000000")),
        yaxis   = list(title = "", tickfont = list(size = 9, color = "#000000")),
        margin  = list(l = 200, r = 50, t = 20, b = 50),
        legend  = list(title = list(text = "Functional Layer", font = list(color = "#000000")),
                       font = list(color = "#000000"))
      )
  })

  output$diag.receiver.source = renderPlotly({
    df = receiver.source.summary %>%
      dplyr::filter(!is.na(pct))
    
    plotly::plot_ly(df,
      x = ~pct,
      y = ~reorder(receiver.source, pct),
      type = "bar", orientation = "h",
      text = ~paste0(pct, "% (", scales::comma(n), ")"),
      textposition = "outside",
      marker = list(color = "#00A087FF"),
      hovertemplate = "%{y}<br>N: %{customdata}<br>%{x:.1f}%<extra></extra>",
      customdata = ~n
    ) %>%
      plotly::layout(
        xaxis = list(title = "Proportion (%)", range = c(0, max(df$pct) * 1.2)),
        yaxis = list(title = "", tickfont = list(size = 11)),
        margin = list(l = 150)
      )
  })



  output$diag.signal.dt = renderDT({
    df = input.signal.summary %>%
      dplyr::mutate(
        pct.zero = round(pct.zero, 1),
        pct.low = round(pct.low, 1),
        median.value = round(median.value, 3)
      )
    DT::datatable(df,
      rownames  = FALSE,
      selection = "none",
      colnames  = c("Modality", "Samples", "Features", "% Zero", "% Low/Zero", "Median Value"),
      options   = list(pageLength = 10, dom = "tp")
    )
  })





  # Observe biological scenario preset changes
  observeEvent(input$about.scenario, {
    req(input$about.scenario)
    sc = input$about.scenario
    if (sc == "checkpoint") {
      updateSliderInput(session, "eq.e.ligand", value = 4.5)
      updateSliderInput(session, "eq.e.receptor", value = 3.0)
      updateSliderInput(session, "eq.p.receiver", value = 8)
      updateSliderInput(session, "eq.g.prog", value = 600)
      updateSliderInput(session, "eq.cnv", value = 3.2)
      updateRadioButtons(session, "eq.snv", selected = "0")
    } else if (sc == "cytokine") {
      updateSliderInput(session, "eq.e.ligand", value = 5.2)
      updateSliderInput(session, "eq.e.receptor", value = 3.8)
      updateSliderInput(session, "eq.p.receiver", value = 18)
      updateSliderInput(session, "eq.g.prog", value = 400)
      updateSliderInput(session, "eq.cnv", value = 2.0)
      updateRadioButtons(session, "eq.snv", selected = "0")
    } else if (sc == "deletion") {
      updateSliderInput(session, "eq.e.ligand", value = 3.0)
      updateSliderInput(session, "eq.e.receptor", value = 2.5)
      updateSliderInput(session, "eq.p.receiver", value = 0)
      updateSliderInput(session, "eq.g.prog", value = 0)
      updateSliderInput(session, "eq.cnv", value = 1.0)
      updateRadioButtons(session, "eq.snv", selected = "0")
    } else if (sc == "baseline") {
      updateSliderInput(session, "eq.e.ligand", value = 3.5)
      updateSliderInput(session, "eq.e.receptor", value = 3.5)
      updateSliderInput(session, "eq.p.receiver", value = 0)
      updateSliderInput(session, "eq.g.prog", value = 0)
      updateSliderInput(session, "eq.cnv", value = 2.0)
      updateRadioButtons(session, "eq.snv", selected = "0")
    }
  }, ignoreInit = TRUE)

  about.score.calc = reactive({
    e.ligand = if (is.null(input$eq.e.ligand)) 4.5 else input$eq.e.ligand
    e.receptor = if (is.null(input$eq.e.receptor)) 3.0 else input$eq.e.receptor
    s.base = e.ligand + e.receptor

    p.receiver = if (is.null(input$eq.p.receiver)) 8 else input$eq.p.receiver
    c.pathway = if (is.null(input$eq.c.path) || input$eq.c.path <= 0) 10 else input$eq.c.path
    beta.pathway = if (is.null(input$eq.beta.path)) 5 else input$eq.beta.path
    w.pathway = beta.pathway * tanh(p.receiver / c.pathway)

    g.prog = if (is.null(input$eq.g.prog)) 600 else input$eq.g.prog
    c.program = if (is.null(input$eq.c.prog) || input$eq.c.prog <= 0) 1000 else input$eq.c.prog
    beta.program = if (is.null(input$eq.beta.prog)) 8 else input$eq.beta.prog
    w.program = beta.program * tanh(g.prog / c.program)

    cnv.val = if (is.null(input$eq.cnv)) 3.2 else input$eq.cnv
    snv.val = if (is.null(input$eq.snv)) 0 else as.numeric(input$eq.snv)
    w.genomic = snv.val + tanh(log2(max(cnv.val, 0.5) / 2))

    s.final = s.base + w.pathway + w.program + w.genomic

    list(
      e.ligand = e.ligand,
      e.receptor = e.receptor,
      s.base = s.base,
      w.pathway = w.pathway,
      w.program = w.program,
      w.genomic = w.genomic,
      s.final = s.final
    )
  })

  output$about.score.summary = renderUI({
    calc = about.score.calc()
    tags$div(style = "display: flex; gap: 10px; margin-bottom: 15px; flex-wrap: wrap;",
      tags$div(style = "flex: 1; min-width: 100px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 8px 12px; text-align: center;",
        tags$div(style = "font-size: 11px; color: #6c757d; text-transform: uppercase;", "S base"),
        tags$div(style = "font-size: 18px; color: #212529;", sprintf("%.2f", calc$s.base))
      ),
      tags$div(style = "flex: 1; min-width: 100px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 8px 12px; text-align: center;",
        tags$div(style = "font-size: 11px; color: #6c757d; text-transform: uppercase;", "W pathway"),
        tags$div(style = paste0("font-size: 18px; color: ", if (calc$w.pathway >= 0) "#28a745;" else "#dc3545;"),
          sprintf("%+.2f", calc$w.pathway)
        )
      ),
      tags$div(style = "flex: 1; min-width: 100px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 8px 12px; text-align: center;",
        tags$div(style = "font-size: 11px; color: #6c757d; text-transform: uppercase;", "W program"),
        tags$div(style = paste0("font-size: 18px; color: ", if (calc$w.program >= 0) "#28a745;" else "#dc3545;"),
          sprintf("%+.2f", calc$w.program)
        )
      ),
      tags$div(style = "flex: 1; min-width: 100px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px; padding: 8px 12px; text-align: center;",
        tags$div(style = "font-size: 11px; color: #6c757d; text-transform: uppercase;", "W genomic"),
        tags$div(style = paste0("font-size: 18px; color: ", if (calc$w.genomic >= 0) "#28a745;" else "#dc3545;"),
          sprintf("%+.2f", calc$w.genomic)
        )
      ),
      tags$div(style = "flex: 1; min-width: 100px; background: #e8f4fd; border: 1px solid #b8daff; border-radius: 6px; padding: 8px 12px; text-align: center;",
        tags$div(style = "font-size: 11px; color: #004085; text-transform: uppercase;", "S final"),
        tags$div(style = "font-size: 18px; color: #004085;", sprintf("%.2f", calc$s.final))
      )
    )
  })

  output$about.waterfall.plot = renderPlotly({
    calc = about.score.calc()

    x.labels = c(
      "Ligand E",
      "Receptor E",
      "S base",
      "W pathway",
      "W program",
      "W genomic",
      "S final"
    )

    y.vals = c(
      calc$e.ligand,
      calc$e.receptor,
      0,
      calc$w.pathway,
      calc$w.program,
      calc$w.genomic,
      0
    )

    text.labels = c(
      sprintf("+%.2f", calc$e.ligand),
      sprintf("+%.2f", calc$e.receptor),
      sprintf("= %.2f", calc$s.base),
      sprintf("%+.2f", calc$w.pathway),
      sprintf("%+.2f", calc$w.program),
      sprintf("%+.2f", calc$w.genomic),
      sprintf("= %.2f", calc$s.final)
    )

    plotly::plot_ly(
      type = "waterfall",
      orientation = "v",
      measure = c("relative", "relative", "total", "relative", "relative", "relative", "total"),
      x = x.labels,
      textposition = "outside",
      text = text.labels,
      y = y.vals,
      connector = list(line = list(color = "rgb(180, 180, 180)", width = 1, dash = "dot")),
      decreasing = list(marker = list(color = "#dc3545")),
      increasing = list(marker = list(color = "#28a745")),
      totals = list(marker = list(color = "#007bff"))
    ) %>%
      plotly::layout(
        title = list(text = "Decomposition of Final Interaction Score", font = list(color = "#212529", size = 15)),
        xaxis = list(title = "Score Component", tickfont = list(color = "#212529")),
        yaxis = list(title = "Interaction Score Contribution", tickfont = list(color = "#212529")),
        margin = list(t = 50, b = 50, l = 60, r = 30)
      )
  })

}

server
