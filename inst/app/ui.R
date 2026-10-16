# Integrated Immune-Tumor Signaling Viewer
# UI Definition

# Custom CSS ----
app.css = "
/* Header & sidebar */
.skin-blue .main-header .logo { background-color: #1a2940; font-weight: 700; letter-spacing: 0.5px; }
.skin-blue .main-header .logo:hover { background-color: #1a2940; }
.skin-blue .main-header .navbar { background-color: #1a2940; }
.skin-blue .main-sidebar { background-color: #22304a; }
.skin-blue .sidebar-menu > li.active > a,
.skin-blue .sidebar-menu > li > a:hover { background-color: #2c4a72; border-left-color: #3498db; }
.skin-blue .sidebar-menu > li > a { color: #c8d6e5; }

/* Instructional banner */
.banner-box {
  background: linear-gradient(135deg, #1a2940 0%, #2c6094 100%);
  color: white; padding: 20px 24px; border-radius: 8px;
  margin-bottom: 18px; box-shadow: 0 3px 8px rgba(0,0,0,0.2);
}
.banner-box h3 { margin-top: 0; font-size: 18px; }
.banner-box p  { margin-bottom: 6px; font-size: 13px; opacity: 0.92; }
.banner-box strong { opacity: 1; }

/* Tab description strip */
.tab-desc {
  background: #eaf4fb; border-left: 4px solid #3498db;
  padding: 9px 14px; border-radius: 0 4px 4px 0;
  margin-bottom: 14px; font-size: 13px; color: #2c3e50;
}

/* Metric definition cards */
.metric-def {
  background: #fdfefe; border: 1px solid #d5e8f5;
  border-radius: 6px; padding: 10px 14px; margin-bottom: 10px;
  font-size: 12.5px;
  color: #1a6ea8;
}
.metric-def .metric-title { font-weight: 700; color: #1a6ea8; }

/* Value boxes */
.small-box { border-radius: 6px; }

/* Clickable sample links in tables */
.sample-link { color: #1a6ea8; cursor: pointer; text-decoration: underline dotted; }

/* Role legend bullets */
.role-bullet { display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 5px; }

/* Circos grid */
.circos-panel { border: 1px solid #dce3ea; border-radius: 6px; padding: 6px; margin-bottom: 12px; background: white; }

/* Responsive plot panels */
.box-body { overflow: hidden; }
"

# UI ----
ui = dashboardPage(
  skin = "blue",

  # Header ----
  dashboardHeader(
    title = tags$span(icon("circle-nodes"), " INTERACT Viewer"),
    titleWidth = 240
  ),

  # Sidebar ----
  dashboardSidebar(
    width = 250,
    sidebarMenu(id = "tabs",
      menuItem("Home",                    tabName = "home",        icon = icon("home")),
      menuItem("Single Sample View",      tabName = "sample",     icon = icon("user")),
      menuItem("Multi-sample View",       tabName = "compare",     icon = icon("users")),
      menuItem("Clusters",                tabName = "clusters",    icon = icon("object-group")),
      menuItem("Pathways & Programs",     tabName = "pathway",     icon = icon("arrows-split-up-and-left")),
      menuItem("Model Diagnostics",       tabName = "diagnostics", icon = icon("square-root-variable")),
      menuItem("About Model",             tabName = "about",       icon = icon("info-circle"))
    )
  ),

  # Body ----
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #F8F9FA; }
        .tab-desc {
          margin-bottom: 15px; padding: 10px; background-color: #f1f7fd;
          border-left: 4px solid #007bff; border-radius: 4px; color: #333; font-size: 13.5px;
        }
        .metric-def { font-size: 13px; line-height: 1.5; color: #555; }
        .metric-title { font-weight: bold; color: #333; }
        .role-bullet { display: inline-block; width: 12px; height: 12px; border-radius: 50%; margin-right: 5px; }
      ")),
      tags$script(HTML("
        // Switch to a shinydashboard tab programmatically
        Shiny.addCustomMessageHandler('switch-to-tab', function(tabName) {
          $('.sidebar-menu a[data-value=\"' + tabName + '\"]').click();
        });
        // Set a selectize input value
        Shiny.addCustomMessageHandler('set-selectize', function(msg) {
          var sel = $('#' + msg.id).selectize();
          if (sel && sel[0] && sel[0].selectize) {
            sel[0].selectize.setValue(msg.value);
          }
        });
      "))
    ),

    tabItems(

      # ── Tab 1: Home ──────────────────────────────────────────────────────────
      tabItem(tabName = "home",

        div(class = "banner-box",
          h3(" INTERACT Signaling Detection Viewer"),
          p("This platform integrates gene expression across diverse cell types with tumor-specific somatic
            CNVs, SNVs, and transcriptional profiles to uncover immune-tumor interactions within the
            tumor microenvironment.")
        ),

        fluidRow(
          valueBoxOutput("vbox.samples",     width = 3),
          valueBoxOutput("vbox.interactions", width = 3),
          valueBoxOutput("vbox.celltypes",    width = 3),
          valueBoxOutput("vbox.genomic.pct",  width = 3)
        ),

        fluidRow(
          box(title = "Global Interaction Strength Heatmap", status = "primary",
              solidHeader = TRUE, width = 8, height = 510,
              div(class = "tab-desc",
                "Mean weighted interaction score between all sender-receiver cell type pairs,
                 aggregated across 241 samples. Darker cells indicate stronger signaling axes."
              ),
              withSpinner(plotlyOutput("home.heatmap", height = "390px"))
          ),
          box(title = "Sample Interaction Counts", status = "primary",
              solidHeader = TRUE, width = 4, height = 510,
              div(class = "tab-desc",
                "Top 10 and bottom 10 samples by total interaction count across the network."
              ),
              withSpinner(plotlyOutput("home.interaction.counts", height = "390px"))
          )
        ),

        fluidRow(
          box(title = "Top Cell-Cell interactions", status = "primary",
              solidHeader = TRUE, width = 4, height = 510,
              div(class = "tab-desc",
                "Highest ranked interactions involving tumor cells as either sender or receiver."
              ),
              withSpinner(plotlyOutput("home.bar.axes", height = "390px"))
          ),
          box(title = "Top Ligand-Receptors Flow", status = "primary",
              solidHeader = TRUE, width = 4, height = 510,
              div(class = "tab-desc",
                "Sankey diagram showing the flow of the top ligand-receptor interactions across the cohort."
              ),
              withSpinner(plotlyOutput("home.sankey", height = "390px"))
          ),
          box(title = "Top Programs", status = "primary",
              solidHeader = TRUE, width = 4, height = 510,
              div(class = "tab-desc",
                "Top biological programs enriched across the ecosystem. Click on a program name to view its source pathway database entry."
              ),
              withSpinner(DTOutput("home.programs.dt"))
          )
        ),

        fluidRow(
          box(title = "Top Cell-Cell Interactions", status = "primary",
              solidHeader = TRUE, width = 12,
              div(class = "tab-desc",
                "Highest-ranked ligand-receptor pairs by mean weighted score across all samples.
                 Click a gene name to view its GeneCards entry. Click a Ligand Axis to search OmnipathDB."
              ),
              fluidRow(
                column(3, numericInput("home.table.n", "Show top N interactions:",
                                       value = 50, min = 10, max = 500, step = 10)),
                column(3, selectInput("home.table.axis", "Filter by Ligand Axis:",
                                      choices = c("All" = "all"), selected = "all")),
                column(3, selectInput("home.table.sender", "Filter by Sender:",
                                      choices = c("All" = "all"), selected = "all")),
                column(3, selectInput("home.table.receiver", "Filter by Receiver:",
                                      choices = c("All" = "all"), selected = "all"))
              ),
              withSpinner(DTOutput("home.eco.table"))
          )
        )
      ),

      # ── Tab 2: Sample Clusters ───────────────────────────────────────────────
      tabItem(tabName = "clusters",

        div(class = "tab-desc",
          "Cluster samples based on cell-cell communication networks. Explore sample clusters in 2D projection space and examine microenvironment signaling patterns across clusters."
        ),

        fluidRow(
          column(width = 3,
            box(title = "Clustering & Display", status = "warning", solidHeader = TRUE, width = 12,
              sliderInput("cluster.res", "Louvain Resolution:",
                          min = 0.2, max = 2.0, value = 1.0, step = 0.2),
              uiOutput("cluster.proj.ui")
            )
          ),

          column(width = 9,
            fluidRow(
              column(width = 6,
                box(title = "Sample Cluster Map", status = "primary",
                    solidHeader = TRUE, width = 12, height = "520px",
                    uiOutput("cluster.map.subtitle"),
                    withSpinner(plotlyOutput("cluster.scatter", height = "430px"))
                )
              ),
              column(width = 6,
                box(title = "Sample Network Graph", status = "primary",
                    solidHeader = TRUE, width = 12, height = "520px",
                    uiOutput("cluster.network.subtitle"),
                    withSpinner(visNetworkOutput("cluster.network", height = "430px"))
                )
              )
            )
          )
        ),

        fluidRow(
          box(title = "Cluster-Specific Edge Attributions (GNN Integrated Gradients)", status = "primary",
              solidHeader = TRUE, width = 12,
              div(class = "tab-desc",
                "Feature attributions computed via Graph Neural Network (GNN) Path-Integrated Gradients.
                 Decomposes the predictive signaling edges that characterize patient microenvironment clusters
                 across five biological channels: RNA Expression, Pathway Activation, Phenotypic Programs, Somatic SNVs, and Copy Number (CNVs)."
              ),
              fluidRow(
                column(4,
                  uiOutput("cluster.ig.select.ui")
                ),
                column(4,
                  sliderInput("cluster.ig.top.n", "Top N Signaling Edges:",
                              min = 5, max = 30, value = 15, step = 5)
                ),
                column(4,
                  radioButtons("cluster.ig.barmode", "Display Mode:",
                               choices = c("Diverging Stacked" = "relative", "Grouped Channels" = "group"),
                               selected = "relative", inline = TRUE)
                )
              ),
              withSpinner(plotlyOutput("cluster.ig.plot", height = "520px"))
          )
        )
      ),

      # ── Tab 3: Sample Overview ───────────────────────────────────────────────
      tabItem(tabName = "sample",

        div(class = "tab-desc",
          " Deep-dive into a single sample's immune-tumor signaling network.
           Node color indicates functional role. Node size scales with hub + authority scores.
           Hover over nodes and edges for detailed metrics."
        ),

        fluidRow(
          column(width = 3,
            box(title = "Sample Selection", status = "warning", solidHeader = TRUE, width = 12,
              selectizeInput("sample.id", "Select Sample ID:",
                             choices = NULL,
                             options = list(
                               placeholder = "Search or select a sample...",
                               maxOptions = 300
                             )),
              sliderInput("sample.top.n", "Top N interactions to display:",
                          min = 20, max = 500, value = 150, step = 10),
              selectInput("sample.node.fill", "Node Fill Color:",
                          choices = c("Cell Type" = "cell_type", "Role" = "role", "Hub Score" = "hub", "Authority Score" = "authority"),
                          selected = "role"),
              selectInput("sample.node.outline", "Node Outline Color:",
                          choices = c("None" = "none", "Cell Type" = "cell_type", "Role" = "role", "Hub Score" = "hub", "Authority Score" = "authority"),
                          selected = "role"),
              selectInput("sample.edge.color", "Edge Color By:",
                          choices = c("Specific Pathway" = "pathway", "Weight" = "weight", "Genomic Hit" = "genomic"),
                          selected = "pathway"),
              selectInput("sample.network.layout", "Network Layout:",
                          choices = c("Nicely" = "layout_nicely",
                                      "Circle" = "layout_in_circle",
                                      "Force-Directed" = "physics",
                                      "Fruchterman-Reingold" = "layout_with_fr",
                                      "Kamada-Kawai" = "layout_with_kk"),
                          selected = "layout_nicely"),
              selectInput("sample.cell.highlight", "Highlight Cell Type:",
                          choices = c("None" = "None"), selected = "None"),
              selectizeInput("sample.ligand.search", "Search Ligand:",
                             choices = NULL, multiple = TRUE,
                             options = list(placeholder = "Type ligand gene name...")),
              selectizeInput("sample.receptor.search", "Search Receptor:",
                             choices = NULL, multiple = TRUE,
                             options = list(placeholder = "Type receptor gene name..."))
            ),

            box(title = "Network Legends", status = "info", solidHeader = TRUE, width = 12,
              h5(tags$b("Node Fill Color")),
              uiOutput("sample.node.legend"),
              hr(style = "margin: 10px 0;"),
              h5(tags$b("Node Outline Color")),
              uiOutput("sample.node.outline.legend"),
              hr(style = "margin: 10px 0;"),
              h5(tags$b("Edge Color")),
              uiOutput("sample.edge.legend"),
              hr(style = "margin: 10px 0;"),
              div(class = "metric-def",
                tags$b("Hub Score:"), " Kleinberg hub score — strength of a node as a signal broadcaster.", tags$br(),
                tags$b("Authority Score:"), " Kleinberg authority score — strength as a signal receiver/integrator."
              )
            ),
          ),

          column(width = 9,
            fluidRow(
              box(title = "Sample Signaling Network", status = "primary",
                  solidHeader = TRUE, width = 8, height = "580px",
                  withSpinner(visNetworkOutput("sample.network", height = "520px"))
              ),
              box(title = "Top Flow (Selected Cell)", status = "primary",
                  solidHeader = TRUE, width = 4, height = "580px",
                  uiOutput("sample.sankey.ui")
              )
            ),
            fluidRow(
              box(title = "Node Metrics & Roles", status = "primary",
                  solidHeader = TRUE, width = 6,
                  withSpinner(DTOutput("sample.metrics.dt"))
              ),
              box(title = "Top Ligand-Receptor Pairs", status = "primary",
                  solidHeader = TRUE, width = 6,
                  fluidRow(
                    column(4, numericInput("sample.lr.n", "Show top N:", value = 100, min = 10, max = 1000, step = 10)),
                    column(4, selectInput("sample.lr.sender", "Sender:", choices = c("All" = "all"), selected = "all")),
                    column(4, selectInput("sample.lr.receiver", "Receiver:", choices = c("All" = "all"), selected = "all"))
                  ),
                  withSpinner(DTOutput("sample.lr.dt"))
              )
            )
          )
        )
      ),

      # ── Tab 3b: Sample Comparison ────────────────────────────────────────────
      tabItem(tabName = "compare",

        div(class = "tab-desc",
          " Build a cohort of samples based on clinical or multi-omic profiles, and compare their immune-tumor signaling networks side-by-side.
            Select up to 12 samples from the cohort table below to generate their network views."
        ),

        fluidRow(
          column(width = 3,
            box(title = "Comparison Controls", status = "warning", solidHeader = TRUE, width = 12,
              sliderInput("compare.top.n", "Top N interactions:", min = 20, max = 500, value = 150, step = 10),
              selectInput("compare.node.fill", "Node Fill Color:",
                          choices = c("Cell Type" = "cell_type", "Role" = "role", "Hub Score" = "hub", "Authority Score" = "authority"), selected = "role"),
              selectInput("compare.node.outline", "Node Outline Color:",
                          choices = c("None" = "none", "Cell Type" = "cell_type", "Role" = "role", "Hub Score" = "hub", "Authority Score" = "authority"), selected = "role"),
              selectInput("compare.edge.color", "Edge Color By:",
                          choices = c("Specific Pathway" = "pathway", "Weight" = "weight", "Genomic Hit" = "genomic"), selected = "pathway"),
              selectInput("compare.network.layout", "Network Layout:",
                          choices = c("Nicely" = "layout_nicely",
                                      "Circle" = "layout_in_circle",
                                      "Force-Directed" = "physics",
                                      "Fruchterman-Reingold" = "layout_with_fr",
                                      "Kamada-Kawai" = "layout_with_kk"),
                          selected = "layout_nicely"),
              selectInput("compare.cell.highlight", "Highlight Cell Type:", choices = c("None" = "None"), selected = "None"),
              selectizeInput("compare.ligand.search", "Search Ligand:", choices = NULL, multiple = TRUE, options = list(placeholder = "Type ligand...")),
              selectizeInput("compare.receptor.search", "Search Receptor:", choices = NULL, multiple = TRUE, options = list(placeholder = "Type receptor..."))
            ),
            box(title = "Network Legend", status = "info", solidHeader = TRUE, width = 12,
              h5(tags$b("Node Fill Color")),
              uiOutput("compare.node.legend"),
              hr(style = "margin: 10px 0;"),
              h5(tags$b("Node Outline Color")),
              uiOutput("compare.node.outline.legend"),
              hr(style = "margin: 10px 0;"),
              h5(tags$b("Edge Color")),
              uiOutput("compare.legend"),
              hr(style = "margin: 10px 0;"),
              div(class = "metric-def",
                tags$b("Mean Edge Weight:"), " Average final functional score of all interactions.", tags$br(),
                tags$b("Kleinberg hub score:"), " Strength of a node as a signal broadcaster.", tags$br(),
                tags$b("Kleinberg authority score:"), " Strength as a signal receiver/integrator."
              )
            )
          ),
          column(width = 9,

            box(title = "Cohort Filters", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(2, uiOutput("cohort.filter.ligand.ui")),
                column(2, uiOutput("cohort.filter.receptor.ui")),
                column(2, uiOutput("cohort.filter.pathway.ui")),
                column(2, uiOutput("cohort.filter.program.ui")),
                column(4, uiOutput("cohort.filter.celltype.ui"))
              ),
              fluidRow(
                column(12, selectizeInput("compare.sample.ids", "Selected Samples (Max 12):",
                                          choices = NULL, multiple = TRUE,
                                          options = list(placeholder = "Search or select samples...", maxItems = 12)))
              )
            ),
            box(title = "Cohort Table", status = "primary", solidHeader = TRUE, width = 12, height = "800px",
                div(style = "overflow-x: auto;", withSpinner(DTOutput("compare.cohort.dt")))
            )
          )
        ),
        uiOutput("compare.network.ui")
      ),



      # ── Tab 5: Pathway Activity Map ───────────────────────────────────────────
      tabItem(tabName = "pathway",

        div(class = "tab-desc",
          " Pathway activity scores reflect the internal transcriptional signaling state of
           each cell type, inferred from gene expression footprints. These scores are computed
           ", tags$b("before"), " interaction weightings are applied and represent the baseline
           cellular context used to calibrate receiver weights."
        ),

        fluidRow(
          column(width = 3,
            box(title = "Display Options", status = "warning", solidHeader = TRUE, width = 12,
              radioButtons("pathway.score.type", "Score Type:",
                           choices = c("Pathways" = "pathway", "Biological Programs" = "program"),
                           selected = "pathway")
            )
          ),

          column(width = 9,
            box(title = "Pathway Activity by Cell Population",
                status = "primary", solidHeader = TRUE, width = 12,
                div(class = "tab-desc",
                  "Heatmap of mean pathway activity scores inferred from downstream transcriptional signatures for each cell population. Red indicates pathway activation, while blue indicates suppression."
                ),
                withSpinner(plotlyOutput("pathway.heatmap", height = "580px"))
            )
          )
        ),

        fluidRow(
          box(title = "Pathway Score Distribution Across Populations",
              status = "primary", solidHeader = TRUE, width = 12,
              div(class = "tab-desc",
                "Violin + box plot of pathway scores, showing variation across
                 all cell populations."
              ),
              withSpinner(plotlyOutput("pathway.violin", height = "380px"))
          )
        )
      ),



      # ── Tab 7: Signal Diagnostics ─────────────────────────────────────────────
      tabItem(tabName = "diagnostics",

        div(class = "tab-desc",
          " Quality-control metrics for the hybrid interaction scoring model. These plots
           reveal the relative contribution of ", tags$b("functional boosting"),
          " (Pathway + program weighting) and ", tags$b("genomic weighting"),
          " (SNV/CNV adjustments) to the final interaction scores across the cohort."
        ),

        fluidRow(
          valueBoxOutput("diag.vbox.boosted",    width = 4),
          valueBoxOutput("diag.vbox.suppressed",  width = 4),
          valueBoxOutput("diag.vbox.genomic",     width = 4)
        ),

        fluidRow(
          box(title = "Weighting Direction by Receiver Population",
              status = "primary", solidHeader = TRUE, width = 6,
              div(class = "tab-desc",
                "Proportion of interactions boosted, suppressed, or unweighted for each receiver.
                 Boosted = combined weight > 1; Suppressed = combined weight < 1."
              ),
              withSpinner(plotlyOutput("diag.weight.dir", height = "420px"))
          ),
          box(title = "Rank Shifting by Receiver Population",
              status = "primary", solidHeader = TRUE, width = 6,
              div(class = "tab-desc",
                "Mean displacement in interaction rank after applying functional and genomic weights.
                 Positive = promoted toward top; negative = demoted. % Shifted = fraction with
                 any rank change."
              ),
              withSpinner(plotlyOutput("diag.rank.shift", height = "420px"))
          )
        ),

        fluidRow(
          box(title = "Mean Layer Weight Contribution by Population",
              status = "primary", solidHeader = TRUE, width = 6,
              div(class = "tab-desc",
                "Average contribution of downstream pathway activation, phenotypic program enrichment, and somatic genomic alterations to the final interaction scores for each receiver population."
              ),
              withSpinner(plotlyOutput("diag.signal.health", height = "380px"))
          ),
          box(title = "Receiver Weight Source Distribution",
              status = "primary", solidHeader = TRUE, width = 6,
              div(class = "tab-desc",
                "Relative contribution of Genomics, Pathways, and Programs to the final interaction score for each receiver population."
              ),
              withSpinner(plotlyOutput("diag.receiver.source", height = "380px"))
          )
        ),

        fluidRow(
          box(title = "Input Signal Summary", status = "primary",
              solidHeader = TRUE, width = 12,
              withSpinner(DTOutput("diag.signal.dt"))
          )
        )
      ),

      # ── Tab 8: About Model ───────────────────────────────────────────────────
      tabItem(tabName = "about",
        tags$div(style = "background-color: #ffffff; padding: 30px; border-radius: 8px; font-size: 15px; box-shadow: 0 1px 3px rgba(0,0,0,0.1);",
          fluidRow(
            column(12,
              tags$img(src = "Logo.png", style = "float: right; width: 140px; margin-bottom: 15px; margin-left: 20px;"),
              tags$h2("INTERACT Framework", style = "margin-top: 0; color: #2c3e50; font-weight: normal;"),
              tags$p("INTERACT is a computational framework for inferring cell-cell communication networks in the tumor microenvironment. It integrates transcriptomic data across distinct cell populations (pseudobulk single-cell RNA-seq and bulk RNA-seq in log2 CPM) with optional tumor somatic SNVs and CNVs. Candidate ligand–receptor interactions are derived from a high-confidence OmniPath intercellular network supported by at least two literature references with compatible transmitter and receiver membrane topologies. Multimeric ligands and receptors are evaluated as complexes using the minimum expression of their required component subunits."),
              tags$hr(),
              tags$h3("Step-by-Step Scoring Methodology", style = "color: #2c3e50; margin-top: 20px; font-weight: normal;"),
              tags$div(style = "margin-bottom: 25px; line-height: 1.6;",
                tags$ol(style = "padding-left: 20px;",
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Baseline Communication Score (\\(S_{\\text{base}}\\)): "),
                    "Calculated for every observed sender–receiver cell population pair as the sum of sender ligand and receiver receptor expression values (\\(E\\)):",
                    tags$div(style = "margin: 8px 0;", "$$S_{\\text{base}} = E_{\\text{ligand}} + E_{\\text{receptor}}$$"),
                    "This represents expression-supported communication potential."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Downstream Pathway Activity (\\(W_{\\text{pathway}}\\)): "),
                    "Pathway activity in the receiver population is estimated using PROGENy (top 100 responsive genes). Each receptor is assigned the maximum activity score among pathways associated with that receptor (\\(P_{\\text{receiver}}\\)). Scores are scaled and bounded by a hyperbolic tangent (\\(\\tanh\\)) function:",
                    tags$div(style = "margin: 8px 0;", "$$W_{\\text{pathway}} = \\beta_{\\text{pathway}} \\tanh\\left(\\frac{P_{\\text{receiver}}}{c_{\\text{pathway}}}\\right)$$"),
                    "where \\(\\beta_{\\text{pathway}}\\) is the pathway weight and \\(c_{\\text{pathway}}\\) is the scaling parameter."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Phenotypic Programs (\\(W_{\\text{program}}\\)): "),
                    "Phenotypic states are estimated using single-sample gene-set enrichment analysis (ssGSEA) with gene sets from ImmPort, Reactome, and Gene Ontology under lineage constraints. Program scores linked to ligand components are assigned to sender populations (\\(G_{\\text{sender}}\\)) and scores linked to receptor components are assigned to receiver populations (\\(G_{\\text{receiver}}\\)). The contribution is calculated jointly:",
                    tags$div(style = "margin: 8px 0;", "$$W_{\\text{program}} = \\beta_{\\text{program}} \\tanh\\left(\\frac{G_{\\text{sender}} + G_{\\text{receiver}}}{c_{\\text{program}}}\\right)$$"),
                    "where \\(\\beta_{\\text{program}}\\) is the program weight and \\(c_{\\text{program}}\\) is the scaling parameter."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Genomic Modifiers (\\(W_{\\text{genomic}}\\)): "),
                    "When tumor somatic SNV and/or CNV data are provided, genomic alterations affecting ligand or receptor genes in the tumor compartment are incorporated. Somatic SNVs are encoded as +1 for activating alterations, -1 for repressive alterations, and 0 for wild type. CNVs use an absolute encoding (diploid = 2):",
                    tags$div(style = "margin: 8px 0;", "$$W_{\\text{genomic,sender}} = \\text{SNV}_{\\text{sender}} + \\tanh\\left[\\log_2\\left(\\frac{\\max(\\text{CNV}_{\\text{sender}}, 0.5)}{2}\\right)\\right]$$"),
                    tags$div(style = "margin: 8px 0;", "$$W_{\\text{genomic,receiver}} = \\text{SNV}_{\\text{receiver}} + \\tanh\\left[\\log_2\\left(\\frac{\\max(\\text{CNV}_{\\text{receiver}}, 0.5)}{2}\\right)\\right]$$"),
                    tags$div(style = "margin: 8px 0;", "$$W_{\\text{genomic}} = W_{\\text{genomic,sender}} + W_{\\text{genomic,receiver}}$$"),
                    "Interactions associated with homozygous deletion-level copy numbers are removed prior to scoring."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Integrated Score (\\(S_{\\text{final}}\\)) & Graph Construction: "),
                    "The final INTERACT score is calculated as the additive integration of expression, pathway, program, and genomic weights:",
                    tags$div(style = "margin: 8px 0;", "$$S_{\\text{final}} = S_{\\text{base}} + W_{\\text{pathway}} + W_{\\text{program}} + W_{\\text{genomic}}$$"),
                    "The resulting scores form sample-specific weighted, directed cell-cell communication graphs where nodes represent cell populations and edges denote active signaling circuits."
                  )
                )
              ),
              tags$hr(),
              tags$h3("Graph Neural Network (GNN) Architecture & Attribution", style = "color: #2c3e50; margin-top: 20px; font-weight: normal;"),
              tags$p("To capture higher-order microenvironment topology and discover predictive signaling circuits, INTERACT incorporates an optional Graph Neural Network (GNN) autoencoder with axiomatic edge attribution:"),
              tags$div(style = "margin-bottom: 25px; line-height: 1.6;",
                tags$ul(style = "padding-left: 20px;",
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Multi-Channel Graph Representation: "),
                    "Patient-specific cell-cell communication networks are encoded into multi-layer graph tensors where cell populations serve as nodes and scored interactions serve as directed edges with feature channels across RNA expression, pathway activation, phenotypic programs, SNVs, and CNVs."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Relational Graph Autoencoder (GINEConv): "),
                    "A multi-layer graph isomorphism autoencoder compresses each patient's ecosystem graph into a 32-dimensional latent embedding space, trained using link reconstruction and feature alignment objectives."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Unsupervised Ecosystem Subtyping: "),
                    "Multi-resolution Louvain community detection applied to the k-nearest neighbor patient graph in GNN latent space identifies recurrent patient microenvironment subtypes."
                  ),
                  tags$li(
                    tags$span(style = "font-weight: 600;", "Path-Integrated Gradients Edge Attribution: "),
                    "Cluster-defining signaling edges are discovered using Path-Integrated Gradients with respect to cluster centroids in latent space, attributing predictive importance to individual communication channels across all five biological layers."
                  )
                )
              ),
              tags$hr(),
              withMathJax(),
              tags$style(HTML(".MathJax_Display { text-align: left !important; margin-left: 0 !important; }")),
              tags$h3("Interactive Score Decomposition Explorer", style = "color: #2c3e50; margin-top: 20px; font-weight: normal;"),
              tags$p("This interactive tool visualizes how each biological layer additively builds the final communication score \\(S_{\\text{final}}\\). Select a biological scenario preset or adjust the parameters below to see the resulting waterfall decomposition."),
              fluidRow(
                box(title = "Scoring Parameters & Scenarios", status = "primary", solidHeader = TRUE, width = 4,
                  selectInput("about.scenario", "Biological Scenario Preset",
                    choices = c(
                      "Checkpoint Evasion (Tumor CD274 -> T cell PDCD1)" = "checkpoint",
                      "Pro-inflammatory Cytokine (Macrophage IL6 -> Tumor IL6R)" = "cytokine",
                      "Genomic Loss (Tumor FAS -> T cell FASLG)" = "deletion",
                      "Neutral Baseline (Pure Co-expression)" = "baseline",
                      "Custom Parameters" = "custom"
                    ),
                    selected = "checkpoint"
                  ),
                  tags$hr(style = "margin: 10px 0;"),
                  sliderInput("eq.e.ligand", "Ligand Expression (E ligand, log2 CPM)", min = 0, max = 10, value = 4.5, step = 0.1),
                  sliderInput("eq.e.receptor", "Receptor Expression (E receptor, log2 CPM)", min = 0, max = 10, value = 3.0, step = 0.1),
                  sliderInput("eq.p.receiver", "Receiver Pathway Score (P receiver)", min = -25, max = 25, value = 8, step = 1),
                  sliderInput("eq.g.prog", "Combined Program Score (G sender + G receiver)", min = -2000, max = 2000, value = 600, step = 50),
                  sliderInput("eq.cnv", "Tumor CNV (Diploid = 2.0)", min = 0.5, max = 5.0, value = 3.2, step = 0.1),
                  radioButtons("eq.snv", "Tumor Somatic SNV",
                    choices = c("Wild Type (0)" = "0", "Activating (+1)" = "1", "Repressive (-1)" = "-1"),
                    selected = "0", inline = TRUE
                  ),
                  tags$details(
                    tags$summary(style = "cursor: pointer; color: #007bff; margin-top: 10px;", "Advanced Weight Settings"),
                    tags$div(style = "padding-top: 10px;",
                      sliderInput("eq.beta.path", "Beta Pathway (max weight)", min = 1, max = 15, value = 5, step = 1),
                      sliderInput("eq.beta.prog", "Beta Program (max weight)", min = 1, max = 15, value = 8, step = 1),
                      numericInput("eq.c.path", "Scaling Parameter c pathway", value = 10, min = 1, max = 50),
                      numericInput("eq.c.prog", "Scaling Parameter c program", value = 1000, min = 100, max = 5000)
                    )
                  )
                ),
                box(title = "Waterfall Score Decomposition", status = "primary", solidHeader = TRUE, width = 8,
                  uiOutput("about.score.summary"),
                  plotlyOutput("about.waterfall.plot", height = "420px"),
                  tags$p(style = "color: #666; font-size: 13px; margin-top: 10px;",
                    "Bars illustrate the step-by-step assembly of the communication score: Baseline expression (S base) forms the foundation, followed by additive contributions from downstream pathway activation (W pathway), phenotypic program enrichment (W program), and somatic genomic alterations (W genomic) to produce the final edge score (S final)."
                  )
                )
              )
            )
          )
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage
