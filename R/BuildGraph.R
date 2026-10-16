#' Build an INTERACT iGraph for a Single Sample
#'
#' @description Converts the flat interaction data frame for a single sample into a
#' directed `igraph` object.
#'
#' @param sample.id String. The sample ID to convert to an iGraph.
#' @param interactions Data frame. The `interactions` output from `RunInteract()`.
#' @param snv.matrix Optional binary matrix of tumor somatic mutations (samples x features).
#' @param cnv.matrix Optional continuous matrix of tumor copy number alterations (samples x features).
#' @param tumor.id String. The exact name of the tumor modality in the dataset (default: "Tumor").
#' @param top.n Integer. The number of top scoring interactions to include in the graph (default: 100).
#'
#' @return A directed `igraph` object with annotated node and edge attributes.
#'
#' @importFrom dplyr %>% .data
#' @importFrom igraph graph_from_data_frame V
#' @importFrom utils head
#' @importFrom stats median
#' @export
BuildGraph = function(sample.id, interactions, snv.matrix = NULL, cnv.matrix = NULL, tumor.id = "Tumor", top.n = 100) {

  # 1. Filter and Subset Data ----

  if (is.null(interactions) || nrow(interactions) == 0) return(NULL)

  pt.data = interactions %>%
    dplyr::filter(.data$sample == sample.id) %>%
    dplyr::slice_max(order_by = .data$final.score, n = top.n, with_ties = FALSE)

  if (nrow(pt.data) == 0) {
    message("No interactions found for this sample.")
    return(NULL)
  }

  # 2. Aggregate Edges and Format Tooltips ----

  edges = pt.data %>%
    dplyr::group_by(.data$sender, .data$receiver) %>%
    dplyr::summarize(
      weight = sum(.data$final.score),
      n.pairs = dplyr::n(),
      top.interactions = paste(utils::head(paste0(
        "<b>", .data$source_genesymbol, " \u2192 ", .data$target_genesymbol, "</b><br>",
        "Score: ", round(.data$final.score, 2), " (Base: ", round(.data$base.score, 2), ")<br>",
        "Axis: ", .data$ligand.axis, " | ", .data$receptor.axis, "<br>",
        "Weights: Gen=", round(.data$weight.genomic, 2),
        " | Path=", round(.data$weight.pathways, 2),
        " | Prog=", round(.data$weight.programs, 2)
      ), 5), collapse = "<hr>"),
      .groups = "drop"
    )

  g = igraph::graph_from_data_frame(edges, directed = TRUE)

  # Store the sample ID directly in the graph object attributes
  g$sample.id = sample.id

  # 3. Inject Node Metadata ----

  for (v in igraph::V(g)$name) {
    v.idx = which(igraph::V(g)$name == v)

    # Track tumor status for downstream node sizing
    igraph::V(g)$is.tumor[v.idx] = (v == tumor.id)
    igraph::V(g)$genomic.info[v.idx] = "No genomic data available"

    if (v == tumor.id) {
      muts = "None detected"
      cnvs = "Diploid"

      if (!is.null(snv.matrix) && sample.id %in% rownames(snv.matrix)) {
        m = snv.matrix[sample.id, ]
        muts = paste(names(m[m == 1]), collapse = ", ")
        if (muts == "") muts = "None detected"
      }

      if (!is.null(cnv.matrix) && sample.id %in% rownames(cnv.matrix)) {
        cn = cnv.matrix[sample.id, ]
        pt.ploidy = stats::median(cn, na.rm = TRUE)
        gain.threshold = pt.ploidy + 0.5
        loss.threshold = pt.ploidy - 0.5

        gains = names(cn[cn > gain.threshold & !is.na(cn)])
        losses = names(cn[cn < loss.threshold & !is.na(cn)])

        cnvs = paste0("Gains: ", ifelse(length(gains) > 0, paste(gains, collapse=", "), "None"),
                      "<br>Losses: ", ifelse(length(losses) > 0, paste(losses, collapse=", "), "None"))
      }
      igraph::V(g)$genomic.info[v.idx] = paste0("<b>SNV:</b> ", muts, "<br><b>CNV:</b><br>", cnvs)
    }

    # Calculate degree statistics
    igraph::V(g)$sent.count[v.idx] = sum(edges$sender == v)
    igraph::V(g)$received.count[v.idx] = sum(edges$receiver == v)
  }

  return(g)
}
