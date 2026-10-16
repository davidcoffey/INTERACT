#' Run INTERACT Signaling Detection
#'
#' @description This function executes the core signaling pipeline for the INTERACT package,
#' integrating genomic alterations, pathway activity, and phenotypic states.
#'
#' @details
#' **Outputs a list of data frames:**
#' \itemize{
#'   \item \strong{interactions:} A data frame of all scored cell-cell interactions. Key columns include \code{source_genesymbol} (ligand), \code{target_genesymbol} (receptor), \code{sample}, \code{sender}, \code{receiver}, \code{base.score} (raw expression product), \code{weight.genomic} (DNA alteration modifier), \code{weight.pathways} (PROGENy modifier), \code{weight.programs} (phenotypic modifier), and \code{final.score} (fully integrated score).
#'   \item \strong{diagnostics:} A list of summary data frames for auditing the pipeline. Includes \code{input.signal} (raw matrix summaries), \code{weight.layers} (global impact of weights), \code{receiver.source} and \code{sender.source} (evidence breakdown), \code{weight.direction} (net functional multiplier direction), \code{rank.shift} (rank changes after weighting), and \code{genomic.hits} (genomically boosted interactions).
#'   \item \strong{topology.metrics:} A list of network metrics calculated per sample. Includes \code{node} (centrality and signaling roles per cell type), \code{edge} (communication weights between populations), and \code{graph} (ecosystem-wide density and transitivity).
#'   \item \strong{Additional Outputs:} \code{precompute} (data for Shiny app), \code{modality.progeny.scores} (raw pathway scores), \code{functional.program.scores} (raw program scores), and \code{program.gene.lists} (reference gene lists).
#' }
#'
#' @param nontumor.matrices A named list of expression matrices (samples x features) for microenvironment cells.
#' @param tumor.matrix Optional. A single expression matrix for the tumor cells.
#' @param snv.matrix Optional. A numeric matrix (samples x features) of tumor somatic mutations (e.g., 1 = activating, -1 = repressive, 0 = WT). Continuous probabilities are also supported.
#' @param cnv.matrix Optional. A matrix (samples x features) of tumor copy number alterations (absolute encoding, diploid = 2).
#' @param program.gene.lists Data frame containing gene sets (ImmPort/Reactome).
#' @param n.cores Integer. Number of cores for parallel processing.
#' @param beta.pathways Numeric. Weighting factor for pathway activity. Default is 5.
#' @param beta.programs Numeric. Weighting factor for phenotypic programs. Default is 8.
#' @param scale.pathways Numeric or "auto". Automatically tuning denominator to bring raw PROGENy scores into roughly `[-2, 2]` range before squashing.
#' @param scale.programs Numeric or "auto". Automatically tuning denominator to bring raw ssGSEA scores into roughly `[-2, 2]` range before squashing.
#' @param pathway.receptor.map Optional custom mapping of pathways to receptors.
#' @param program.constraints Optional data frame for lineage-specific program filtering.
#' @param omni.network Optional OmniPath network data frame to bypass API download.
#' @param min.expression Numeric. Minimum log-scale expression threshold for ligand or receptor subunits.
#' @param calibrate.n.samples Integer. Number of random samples to profile during auto-calibration. Default is 15.
#' @param run.autoencoder Logical. If TRUE, trains the multi-channel GNN Graph Autoencoder on the interaction network and extracts latent patient embeddings for downstream UMAP and Louvain clustering. Default is TRUE.
#' @param conda.env String. Conda environment containing PyTorch and PyTorch Geometric for the autoencoder. Default is "pytorch".
#' @param autoencoder.epochs Integer. Number of training epochs for the Graph Autoencoder. Default is 150.
#'
#' @return A list containing interactions, diagnostics, embeddings, precompute data, and topology metrics.
#'
#' @importFrom dplyr %>% .data
#' @importFrom stats setNames median var prcomp kmeans quantile
#' @importFrom rlang .data
#' @importFrom data.table as.data.table frank
#' @importFrom umap umap
#' @export
RunInteract <- function(nontumor.matrices,
                        tumor.matrix = NULL,
                        snv.matrix = NULL,
                        cnv.matrix = NULL,
                        program.gene.lists = NULL,
                        n.cores = max(1, parallelly::availableCores() - 2),
                        beta.pathways = 5,
                        beta.programs = 8,
                        scale.pathways = "auto",
                        scale.programs = "auto",
                        pathway.receptor.map = NULL,
                        program.constraints = NULL,
                        omni.network = NULL,
                        min.expression = 1,
                        calibrate.n.samples = 15,
                        run.autoencoder = TRUE,
                        conda.env = "pytorch",
                        autoencoder.epochs = 150L) {
  start.time <- Sys.time()

  old.pb <- pbapply::pboptions(type = "timer", char = "=", style = 3)
  on.exit(pbapply::pboptions(old.pb), add = TRUE)

  message("\n--- Running INTERACT Signaling Pipeline v1.0 ---")
  message(sprintf("Cores allocated: %d", n.cores))

  total.steps <- if (run.autoencoder) 9 else 8
  c.step <- 1

  # 1. Validation & Reference Import ----
  message(sprintf("[%d/%d] Validating inputs and generating OmniPath reference...", c.step, total.steps))
  c.step <- c.step + 1

  if (!is.list(nontumor.matrices)) stop("nontumor.matrices must be a named list of matrices.")

  if (is.null(program.gene.lists)) {
    gene.list.file = system.file("extdata", "ImmPort-gene-lists.txt", package = "INTERACT")
    if (gene.list.file == "") {
      gene.list.file = system.file("extdata", "ImmPort gene lists.txt", package = "INTERACT")
    }
    program.gene.lists = utils::read.delim(gene.list.file, header = FALSE, stringsAsFactors = FALSE)
    colnames(program.gene.lists) = c("program.name", "id", "source", "url", "count")
  }

  input <- nontumor.matrices
  tumor.modality <- "Tumor"
  if (!is.null(tumor.matrix)) input[[tumor.modality]] <- tumor.matrix

  genomic.modalities <- character()
  if (!is.null(snv.matrix)) genomic.modalities <- c(genomic.modalities, "SNV")
  if (!is.null(cnv.matrix)) genomic.modalities <- c(genomic.modalities, "CNV")

  apply.genomic <- (!is.null(tumor.matrix) && length(genomic.modalities) > 0)
  options(omnipathr.console_logging = FALSE)

  input.all <- input
  genomic.matrices <- list()
  if ("SNV" %in% genomic.modalities) {
    input.all[["SNV"]] <- snv.matrix
    genomic.matrices[["SNV"]] <- snv.matrix
  }
  if ("CNV" %in% genomic.modalities) {
    input.all[["CNV"]] <- cnv.matrix
    genomic.matrices[["CNV"]] <- cnv.matrix
  }

  all.features <- unique(unlist(lapply(input.all, colnames)))

  if (is.null(pathway.receptor.map)) pathway.receptor.map <- get.default.pathway.map()

  if (is.null(omni.network)) {
    lr.annotated <- OmnipathR::intercell_network(high_confidence = TRUE)
    lr.filtered <- OmnipathR::filter_intercell_network(
      lr.annotated,
      transmitter_topology = c("sec", "pmtm", "pmp"),
      receiver_topology = "pmtm",
      min_references = 2
    )
  } else {
    lr.filtered <- omni.network
  }

  lr.reference <- dplyr::filter(lr.filtered, is.mappable(.data$source_genesymbol, all.features), is.mappable(.data$target_genesymbol, all.features))
  lr.reference <- dplyr::group_by(lr.reference, .data$source_genesymbol, .data$target_genesymbol)
  lr.reference <- dplyr::summarize(lr.reference,
    ligand.axis = pick.axis(.data$category_intercell_source, .data$consensus_score_intercell_source),
    receptor.axis = pick.axis(.data$category_intercell_target, .data$consensus_score_intercell_target),
    .groups = "drop"
  )
  lr.reference <- dplyr::distinct(lr.reference, .data$source_genesymbol, .data$target_genesymbol, .keep_all = TRUE)

  unique.symbols <- unique(c(lr.reference$source_genesymbol, lr.reference$target_genesymbol))
  all.receptors <- unique(unlist(strsplit(gsub("COMPLEX:", "", lr.reference$target_genesymbol), "_", fixed = TRUE)))
  all.ligands <- unique(unlist(strsplit(gsub("COMPLEX:", "", lr.reference$source_genesymbol), "_", fixed = TRUE)))

  # 2. Functional Mapping & Indices ----
  expression.modalities <- setdiff(names(input), genomic.modalities)
  message(sprintf("[%d/%d] Precomputing optimized interaction mappings across %d populations...", c.step, total.steps, length(expression.modalities)))
  c.step <- c.step + 1

  receptor.pathways.lookup <- tibble::enframe(pathway.receptor.map, name = "pathway", value = "receptor")
  receptor.pathways.lookup <- tidyr::unnest(receptor.pathways.lookup, .data$receptor)
  receptor.pathways.lookup <- dplyr::filter(receptor.pathways.lookup, .data$receptor %in% all.receptors)

  valid.programs.df <- get.valid.programs(program.gene.lists, unique.symbols)
  programs.list <- split(valid.programs.df$gene.symbol, valid.programs.df$program.name)

  signature.receptor.lookup <- dplyr::filter(valid.programs.df, .data$gene.symbol %in% all.receptors) %>%
    dplyr::select(signature = .data$program.name, receptor = .data$gene.symbol) %>%
    dplyr::distinct()
  signature.ligand.lookup <- dplyr::filter(valid.programs.df, .data$gene.symbol %in% all.ligands) %>%
    dplyr::select(signature = .data$program.name, ligand = .data$gene.symbol) %>%
    dplyr::distinct()

  modality.to.lineage <- get.lineage.map(expression.modalities)
  if (is.null(program.constraints)) program.constraints <- get.program.constraints(names(programs.list))

  feature.universe <- toupper(unique(unlist(lapply(input.all, colnames))))

  symbol.parsing.list <- lapply(unique.symbols, function(s) {
    parts <- toupper(unlist(strsplit(gsub("COMPLEX:", "", s), "_", fixed = TRUE)))
    match(parts, feature.universe)
  })

  complex.lig.sigs <- lapply(unique.symbols, function(s) {
    parts <- unlist(strsplit(gsub("COMPLEX:", "", s), "_", fixed = TRUE))
    unique(signature.ligand.lookup$signature[signature.ligand.lookup$ligand %in% parts])
  })
  complex.rec.sigs <- lapply(unique.symbols, function(s) {
    parts <- unlist(strsplit(gsub("COMPLEX:", "", s), "_", fixed = TRUE))
    unique(signature.receptor.lookup$signature[signature.receptor.lookup$receptor %in% parts])
  })
  complex.rec.paths <- lapply(unique.symbols, function(s) {
    parts <- unlist(strsplit(gsub("COMPLEX:", "", s), "_", fixed = TRUE))
    unique(receptor.pathways.lookup$pathway[receptor.pathways.lookup$receptor %in% parts])
  })

  all.pops <- setdiff(names(input), genomic.modalities)

  pop.universe.idx <- lapply(input[all.pops], function(m) match(feature.universe, toupper(colnames(m))))
  genomic.universe.idx <- lapply(genomic.matrices, function(m) match(feature.universe, toupper(colnames(m))))

  # 3. State Inference ----
  message(sprintf("[%d/%d] Inferring modality-specific functional states (PROGENy & ssGSEA)...", c.step, total.steps))
  c.step <- c.step + 1

  modality.pathway.scores <- pbapply::pblapply(expression.modalities, function(pop) {
    tryCatch(
      {
        mat <- input[[pop]]
        if (is.null(mat) || nrow(mat) < 5) {
          return(NULL)
        }
        as.matrix(progeny::progeny(t(mat), scale = FALSE, top = 100))
      },
      error = function(e) NULL
    )
  }, cl = 1)
  names(modality.pathway.scores) <- expression.modalities

  functional.program.scores <- pbapply::pblapply(expression.modalities, function(pop) {
    tryCatch(
      {
        compute.modality.ssgsea(pop, input[[pop]], programs.list, modality.to.lineage, program.constraints)
      },
      error = function(e) NULL
    )
  }, cl = 1)
  names(functional.program.scores) <- expression.modalities

  if (!is.null(tumor.matrix)) {
    scoring.samples <- intersect(rownames(input[[tumor.modality]]), rownames(modality.pathway.scores[[tumor.modality]]))
    if (length(scoring.samples) == 0) stop("Sample ID mismatch across modalities.")
  } else {
    scoring.samples <- unique(unlist(lapply(nontumor.matrices, rownames)))
  }

  if (length(scoring.samples) == 0) stop("No samples found in input matrices. Please check rownames.")

  ctx <- list(
    lr.df = lr.reference,
    mats = input,
    gen.mats = genomic.matrices,
    all.path.scores = modality.pathway.scores,
    all.prog.scores = functional.program.scores,
    complex.lig.sigs = complex.lig.sigs,
    complex.rec.sigs = complex.rec.sigs,
    complex.rec.paths = complex.rec.paths,
    unique.symbols = unique.symbols,
    parsing.list = symbol.parsing.list,
    pop.idx = pop.universe.idx,
    gen.idx = genomic.universe.idx,
    lr.src.idx = match(lr.reference$source_genesymbol, unique.symbols),
    lr.tgt.idx = match(lr.reference$target_genesymbol, unique.symbols),
    apply.genomic = apply.genomic,
    tumor.pop = tumor.modality,
    gen.pops = genomic.modalities,
    beta.pathways = beta.pathways,
    beta.programs = beta.programs,
    min.expression = min.expression
  )

  # 4. Auto-Calibration mode ----
  calc.path <- identical(scale.pathways, "auto")
  calc.prog <- identical(scale.programs, "auto")

  if (calc.path || calc.prog) {
    message(sprintf("[%d/%d] Auto-calibrating functional score boundaries...", c.step, total.steps))
    c.step <- c.step + 1

    set.seed(42)
    calib.samples <- sample(scoring.samples, min(calibrate.n.samples, length(scoring.samples)))

    calib.df <- dplyr::bind_rows(lapply(calib.samples, function(pt) {
      bd <- tryCatch(compute.pt.base.state(pt, ctx), error = function(e) NULL)
      if (is.null(bd) || nrow(bd) == 0) {
        return(NULL)
      }
      bd <- bd[!((bd$s.cnv < 0.5) | (bd$r.cnv < 0.5)), , drop = FALSE]
      if (nrow(bd) == 0) {
        return(NULL)
      }
      data.frame(r.path = bd$r.path, prog.sum = bd$s.prog + bd$r.prog)
    }))

    if (is.null(calib.df) || nrow(calib.df) == 0) stop("Calibration failed: No valid edges found in sample subset.")

    if (calc.path) {
      q.path <- stats::quantile(abs(calib.df$r.path), probs = 0.95, na.rm = TRUE)
      ctx$scale.pathways <- max(1e-4, as.numeric(q.path) / 2)
      message(sprintf("      Auto-set scale.pathways to %.2f", ctx$scale.pathways))
    } else {
      ctx$scale.pathways <- as.numeric(scale.pathways)
    }

    if (calc.prog) {
      q.prog <- stats::quantile(abs(calib.df$prog.sum), probs = 0.95, na.rm = TRUE)
      ctx$scale.programs <- max(1e-4, as.numeric(q.prog) / 2)
      message(sprintf("      Auto-set scale.programs to %.2f", ctx$scale.programs))
    } else {
      ctx$scale.programs <- as.numeric(scale.programs)
    }
  } else {
    message(sprintf("[%d/%d] Using provided static scaling parameters...", c.step, total.steps))
    c.step <- c.step + 1
    ctx$scale.pathways <- as.numeric(scale.pathways)
    ctx$scale.programs <- as.numeric(scale.programs)
  }

  # 5. Final Assembly ----
  message(sprintf("[%d/%d] Computing final scored interactions...", c.step, total.steps))
  c.step = c.step + 1
  c.step <- c.step + 1

  ctx.size.gb = round(as.numeric(utils::object.size(ctx)) / 1e9, 2)
  safe.cores = max(1, min(n.cores, floor(32 / max(ctx.size.gb, 1))))
  ctx.size.gb <- round(as.numeric(utils::object.size(ctx)) / 1e9, 2)
  safe.cores <- max(1, min(n.cores, floor(32 / max(ctx.size.gb, 1))))

  use.fork = (safe.cores > 1) && can.support.fork()
  cluster.obj = NULL

  if (safe.cores > 1 && !use.fork) {
    cluster.obj = tryCatch(
      {
        cl.tmp = parallel::makeCluster(safe.cores, type = "PSOCK")
        parallel::clusterEvalQ(cl.tmp, {
          if (requireNamespace("INTERACT", quietly = TRUE)) {
            suppressPackageStartupMessages(library(INTERACT))
          }
        })
        cl.tmp
      },
      error = function(e) {
        warning("Could not initialize PSOCK cluster (", e$message, "). Falling back to single-core execution.")
        NULL
      }
    )
  }

  cl.arg = if (!is.null(cluster.obj)) {
    cluster.obj
  } else if (use.fork) {
    safe.cores
  } else {
    1
  }

  if (!is.null(cluster.obj)) {
    on.exit({
      if (!is.null(cluster.obj)) {
        try(parallel::stopCluster(cluster.obj), silent = TRUE)
      }
    }, add = TRUE)
  }

  score.patient.wrapper = function(pt) {
    base.df = tryCatch(
      compute.pt.base.state(pt, ctx),
      error = function(e) {
        warning(sprintf("[INTERACT sample %s error in base state]: %s", pt, e$message))
        NULL
      }
    )
    if (is.null(base.df) || nrow(base.df) == 0) {
      return(NULL)
    }
    res.weighted = tryCatch(
      apply.functional.weights(base.df, ctx),
      error = function(e) {
        warning(sprintf("[INTERACT sample %s error in weighting]: %s", pt, e$message))
        NULL
      }
    )
    return(res.weighted)
  }

  final.results = tryCatch(
    {
      pbapply::pblapply(scoring.samples, score.patient.wrapper, cl = cl.arg)
    },
    error = function(e) {
      if (grepl("fork", e$message, ignore.case = TRUE) && is.null(cluster.obj)) {
        message("Notice: Forked parallel execution is unsupported in this environment. Switching to PSOCK cluster...")
        psock.cl = tryCatch(parallel::makeCluster(safe.cores, type = "PSOCK"), error = function(e2) NULL)
        if (!is.null(psock.cl)) {
          parallel::clusterEvalQ(psock.cl, {
            if (requireNamespace("INTERACT", quietly = TRUE)) {
              suppressPackageStartupMessages(library(INTERACT))
            }
          })
          res = pbapply::pblapply(scoring.samples, score.patient.wrapper, cl = psock.cl)
          try(parallel::stopCluster(psock.cl), silent = TRUE)
          return(res)
        } else {
          return(pbapply::pblapply(scoring.samples, score.patient.wrapper, cl = 1))
        }
      } else {
        stop(e)
      }
    }
  )

  if (!is.null(cluster.obj)) {
    try(parallel::stopCluster(cluster.obj), silent = TRUE)
    cluster.obj = NULL
  }

  final.df = dplyr::bind_rows(Filter(function(x) is.data.frame(x) && nrow(x) > 0, final.results))

  if (is.null(final.df) || nrow(final.df) == 0) {
    stop("No interactions passed expression and mappability filters. Consider lowering min.expression.")
  }

  # Guardrail: Check if exact-match alignment failed silently for genomics
  if (apply.genomic && nrow(final.df) > 0) {
    tumor.edges <- final.df[final.df$sender == tumor.modality | final.df$receiver == tumor.modality, ]
    if (nrow(tumor.edges) > 0) {
      if ("SNV" %in% genomic.modalities && sum(abs(tumor.edges$s.snv) > 0 | abs(tumor.edges$r.snv) > 0, na.rm = TRUE) == 0) {
        warning("\n*** GENOMIC ALIGNMENT WARNING ***\nSNV matrix provided, but ALL resulting tumor interactions mapped to Wild-Type (0). Check if your SNV gene symbols perfectly match the expression matrix casing (e.g., HGNC vs Ensembl).\n")
      }
      if ("CNV" %in% genomic.modalities && sum(abs(tumor.edges$s.cnv - 2.0) > 0.1 | abs(tumor.edges$r.cnv - 2.0) > 0.1, na.rm = TRUE) == 0) {
        warning("\n*** GENOMIC ALIGNMENT WARNING ***\nCNV matrix provided, but ALL resulting tumor interactions mapped to Diploid (2.0). Check if your CNV gene symbols perfectly match the expression matrix casing.\n")
      }
    }
  }

  # 6. Diagnostics & Auditing ----
  message(sprintf("[%d/%d] Generating diagnostic and auditing summaries...", c.step, total.steps))
  c.step <- c.step + 1

  input.signal.list <- lapply(names(input.all), function(mod) {
    mat <- input.all[[mod]]
    if (is.null(mat)) {
      return(NULL)
    }
    mod.mean <- mean(mat, na.rm = TRUE)
    tibble::tibble(
      modality = mod, samples = nrow(mat), features = ncol(mat),
      pct.zero = mean(mat == 0, na.rm = TRUE) * 100,
      pct.low = mean(mat < mod.mean, na.rm = TRUE) * 100,
      median.value = stats::median(mat, na.rm = TRUE)
    )
  })
  input.signal.summary <- dplyr::bind_rows(input.signal.list)

  if (!is.null(final.df) && nrow(final.df) > 0) {
    weight.layer.summary <- dplyr::summarize(final.df,
      pct.pathway.weighted = mean(abs(.data$weight.pathways) > 1e-10) * 100,
      pct.program.weighted = mean(abs(.data$weight.programs) > 1e-10) * 100,
      pct.genomic.weighted = mean(abs(.data$weight.genomic) > 1e-10) * 100,
      pct.any.weighted = mean(abs(.data$weight.pathways) > 1e-10 | abs(.data$weight.programs) > 1e-10 | abs(.data$weight.genomic) > 1e-10) * 100
    )

    tmp.rec <- dplyr::mutate(final.df, receiver.source = dplyr::case_when(
      abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.pathways) > 1e-10 & abs(.data$weight.receiver.programs) > 1e-10 ~ "Genomic + Pathways + Programs",
      abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.pathways) > 1e-10 ~ "Genomic + Pathways",
      abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.receiver.programs) > 1e-10 ~ "Genomic + Programs",
      abs(.data$weight.receiver.genomic) > 1e-10 ~ "Genomic only",
      abs(.data$weight.pathways) > 1e-10 & abs(.data$weight.receiver.programs) > 1e-10 ~ "Pathways + Programs",
      abs(.data$weight.pathways) > 1e-10 ~ "Pathways only",
      abs(.data$weight.receiver.programs) > 1e-10 ~ "Programs only",
      TRUE ~ "No receiver functional weight"
    ))
    receiver.source.summary <- dplyr::mutate(dplyr::count(tmp.rec, .data$receiver.source), pct = round(.data$n / sum(.data$n) * 100, 1))

    tmp.snd <- dplyr::mutate(final.df, sender.source = dplyr::case_when(
      abs(.data$weight.sender.genomic) > 1e-10 & abs(.data$weight.sender.programs) > 1e-10 ~ "Genomic + Programs",
      abs(.data$weight.sender.genomic) > 1e-10 ~ "Genomic only",
      abs(.data$weight.sender.programs) > 1e-10 ~ "Programs only",
      TRUE ~ "No sender functional weight"
    ))
    sender.source.summary <- dplyr::mutate(dplyr::count(tmp.snd, .data$sender.source), pct = round(.data$n / sum(.data$n) * 100, 1))

    tmp.dir <- dplyr::mutate(final.df, net.modifier = .data$final.score - .data$base.score, weight.class = dplyr::case_when(.data$net.modifier > 0.01 ~ "Boosted", .data$net.modifier < -0.01 ~ "Suppressed", TRUE ~ "Neutral"))
    tmp.dir <- dplyr::group_by(tmp.dir, .data$receiver, .data$weight.class)
    weight.direction.summary <- tidyr::pivot_wider(dplyr::summarize(tmp.dir, n = dplyr::n(), .groups = "drop"), names_from = .data$weight.class, values_from = .data$n, values_fill = list(n = 0))

    tmp.rank <- dplyr::group_by(final.df, .data$sample) %>%
      dplyr::mutate(base.rank = dplyr::min_rank(dplyr::desc(.data$base.score)), final.rank = dplyr::min_rank(dplyr::desc(.data$final.score)), rank.shift = .data$base.rank - .data$final.rank) %>%
      dplyr::group_by(.data$receiver)
    rank.shift.summary <- dplyr::summarize(
      tmp.rank,
      mean.rank.shift = mean(.data$rank.shift, na.rm = TRUE),
      mean.abs.rank.shift = round(mean(abs(.data$rank.shift), na.rm = TRUE), 2),
      pct.shifted = round(mean(abs(.data$rank.shift) > 0, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )

    mean.layer.weights.summary <- final.df %>%
      dplyr::group_by(.data$receiver) %>%
      dplyr::summarize(
        Pathways = mean(.data$weight.pathways, na.rm = TRUE),
        Programs = mean(.data$weight.programs, na.rm = TRUE),
        Genomics = mean(.data$weight.genomic, na.rm = TRUE),
        .groups = "drop"
      )

    # Rebuilt dynamically using the raw directional states
    if (apply.genomic) {
      genomic.hits.summary <- final.df %>%
        dplyr::filter(abs(.data$s.snv) > 0 | abs(.data$r.snv) > 0 | abs(.data$s.cnv - 2.0) >= 0.5 | abs(.data$r.cnv - 2.0) >= 0.5) %>%
        dplyr::mutate(
          impacted = dplyr::case_when(
            .data$sender == tumor.modality & .data$receiver == tumor.modality ~ "Both",
            .data$sender == tumor.modality ~ "Sender",
            .data$receiver == tumor.modality ~ "Receiver",
            TRUE ~ "Unknown"
          ),
          alteration.type = dplyr::case_when(
            (abs(.data$s.snv) > 0 | abs(.data$r.snv) > 0) & (abs(.data$s.cnv - 2.0) >= 0.5 | abs(.data$r.cnv - 2.0) >= 0.5) ~ "SNV + CNV",
            .data$s.snv > 0 | .data$r.snv > 0 ~ "Activating SNV",
            .data$s.snv < 0 | .data$r.snv < 0 ~ "Repressive SNV",
            .data$s.cnv > 2.5 | .data$r.cnv > 2.5 ~ "CNV Gain",
            .data$s.cnv < 1.5 | .data$r.cnv < 1.5 ~ "CNV Loss",
            TRUE ~ "Other"
          )
        ) %>%
        dplyr::group_by(.data$source_genesymbol, .data$target_genesymbol, .data$sender, .data$receiver, .data$impacted, .data$alteration.type) %>%
        dplyr::summarize(
          n.samples = dplyr::n_distinct(.data$sample),
          mean.base.score = mean(.data$base.score, na.rm = TRUE),
          mean.final.score = mean(.data$final.score, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(.data$n.samples))
    } else {
      genomic.hits.summary <- tibble::tibble()
    }
  } else {
    weight.layer.summary <- receiver.source.summary <- weight.direction.summary <- rank.shift.summary <- sender.source.summary <- genomic.hits.summary <- mean.layer.weights.summary <- tibble::tibble()
  }

  # 7. Multi-Channel Graph Autoencoder ----
  gnn.embeddings <- NULL
  if (isTRUE(run.autoencoder)) {
    message(sprintf("[%d/%d] Training multi-channel Graph Autoencoder on interaction network...", c.step, total.steps))
    c.step <- c.step + 1

    gnn.embeddings <- tryCatch(
      {
        RunInteractAutoencoder(
          interact.results = final.df,
          conda.env = conda.env,
          epochs = as.integer(autoencoder.epochs),
          deterministic = TRUE,
          save.csv = FALSE
        )
      },
      error = function(e) {
        warning("\n*** GRAPH AUTOENCODER WARNING ***\nCould not train Graph Autoencoder: ", e$message,
                "\nFalling back to PCA for sample embeddings.\n")
        NULL
      }
    )
  }

  # 8. INTERACT Viewer Precompute ----
  message(sprintf("[%d/%d] Precomputing data for Shiny application...", c.step, total.steps))
  c.step <- c.step + 1

  sample.ids.vec <- sort(unique(final.df$sample))
  all.populations <- sort(unique(c(final.df$sender, final.df$receiver)))
  n.pops <- length(all.populations)

  pop.palette <- ggsci::pal_igv()(51)
  pop.colors <- pop.palette[seq_len(min(n.pops, 51))]
  if (n.pops > 51) pop.colors <- grDevices::colorRampPalette(pop.palette)(n.pops)
  names(pop.colors) <- all.populations

  all.gene.symbols <- sort(unique(c(final.df$source_genesymbol, final.df$target_genesymbol)))
  ligand.axes.available <- sort(unique(final.df$ligand.axis))

  sender.receiver.agg <- final.df %>%
    dplyr::group_by(sender, receiver) %>%
    dplyr::summarize(mean.score = mean(final.score, na.rm = TRUE), total.score = sum(final.score, na.rm = TRUE), n.pairs = dplyr::n(), pct.genomic = mean(as.numeric(abs(weight.genomic) > 1e-10), na.rm = TRUE) * 100, .groups = "drop")
  top.lr.global <- final.df %>%
    dplyr::group_by(source_genesymbol, target_genesymbol, ligand.axis, receptor.axis) %>%
    dplyr::summarize(mean.score = mean(final.score, na.rm = TRUE), n.samples = dplyr::n_distinct(sample), pct.genomic = mean(as.numeric(abs(weight.genomic) > 1e-10), na.rm = TRUE) * 100, .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean.score))
  top.ecosystem.events <- final.df %>%
    dplyr::group_by(source_genesymbol, target_genesymbol, sender, receiver, ligand.axis, receptor.axis) %>%
    dplyr::summarize(mean.final.score = mean(final.score, na.rm = TRUE), mean.base.score = mean(base.score, na.rm = TRUE), n.samples = dplyr::n_distinct(sample), genomic.rate = mean(as.numeric(abs(weight.genomic) > 1e-10), na.rm = TRUE) * 100, .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean.final.score)) %>%
    utils::head(500)

  pt.profile.long <- final.df %>%
    dplyr::group_by(sample, sender, receiver) %>%
    dplyr::summarize(total.score = sum(final.score, na.rm = TRUE), mean.score = mean(final.score, na.rm = TRUE), .groups = "drop")
  sample.interaction.counts <- final.df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(interaction.count = dplyr::n(), .groups = "drop")

  pt.profile.wide <- tidyr::pivot_wider(pt.profile.long, id_cols = sample, names_from = c(sender, receiver), names_sep = "__|__", values_from = mean.score, values_fill = list(mean.score = 0))
  mat.samples <- pt.profile.wide$sample
  pt.mat.raw <- as.matrix(pt.profile.wide[, -1])
  rownames(pt.mat.raw) <- mat.samples

  col.var <- apply(pt.mat.raw, 2, stats::var)
  pt.mat.clean <- pt.mat.raw[, col.var > 0, drop = FALSE]

  if (ncol(pt.mat.clean) > 1) {
    pt.pca.x <- stats::prcomp(scale(pt.mat.clean), center = FALSE, scale. = FALSE)$x
    n.pcs.use <- min(20, ncol(pt.pca.x))
  } else {
    pt.pca.x <- matrix(0, nrow = length(mat.samples), ncol = 2)
    rownames(pt.pca.x) <- mat.samples
    n.pcs.use <- 2
  }

  resolutions <- seq(0.2, 2.0, by = 0.2)
  all.louvain <- list()
  sample.knn.edges.df <- data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)

  if (!is.null(gnn.embeddings) && nrow(gnn.embeddings) > 0) {
    emb.match <- match(mat.samples, gnn.embeddings$sample)
    emb.df <- gnn.embeddings[emb.match, , drop = FALSE]
    emb.mat <- as.matrix(emb.df[, -1, drop = FALSE])
    rownames(emb.mat) <- mat.samples

    n.neighbors <- min(15, max(2, length(mat.samples) - 1))
    set.seed(42)
    umap.coords <- uwot::umap(emb.mat, n_neighbors = n.neighbors, min_dist = 0.3)
    embed.dim1 <- umap.coords[, 1]
    embed.dim2 <- umap.coords[, 2]
    embed.label <- c("GNN UMAP 1", "GNN UMAP 2")

    k_nn <- min(10, max(2, length(mat.samples) - 1))
    dist.mat <- as.matrix(stats::dist(emb.mat))
    edge.list <- lapply(seq_len(nrow(dist.mat)), function(i) {
      nn <- order(dist.mat[i, ])[2:(k_nn + 1)]
      data.frame(
        from = rownames(dist.mat)[i],
        to = rownames(dist.mat)[nn],
        weight = 1 / (1 + dist.mat[i, nn]),
        stringsAsFactors = FALSE
      )
    })
    sample.knn.edges.df <- do.call(rbind, edge.list)

    knn.g <- igraph::graph_from_data_frame(sample.knn.edges.df, directed = FALSE)
    knn.g <- igraph::simplify(knn.g, edge.attr.comb = "mean")

    for (res in resolutions) {
      res.key <- paste0("res_", sprintf("%.1f", res))
      cl <- tryCatch(igraph::cluster_louvain(knn.g, resolution = res), error = function(e) NULL)
      if (!is.null(cl)) {
        memb <- igraph::membership(cl)
        all.louvain[[res.key]] <- as.character(memb[mat.samples])
      } else {
        all.louvain[[res.key]] <- rep("1", length(mat.samples))
      }
    }
  } else {
    if (ncol(pt.pca.x) >= 2) {
      n.neighbors <- min(15, max(2, length(mat.samples) - 1))
      set.seed(42)
      umap.coords <- uwot::umap(pt.pca.x[, 1:n.pcs.use, drop = FALSE], n_neighbors = n.neighbors, min_dist = 0.3)
      embed.dim1 <- umap.coords[, 1]
      embed.dim2 <- umap.coords[, 2]
      embed.label <- c("UMAP 1", "UMAP 2")

      k_nn <- min(10, max(2, length(mat.samples) - 1))
      dist.mat <- as.matrix(stats::dist(pt.pca.x[, 1:n.pcs.use, drop = FALSE]))
      edge.list <- lapply(seq_len(nrow(dist.mat)), function(i) {
        nn <- order(dist.mat[i, ])[2:(k_nn + 1)]
        data.frame(
          from = rownames(dist.mat)[i],
          to = rownames(dist.mat)[nn],
          weight = 1 / (1 + dist.mat[i, nn]),
          stringsAsFactors = FALSE
        )
      })
      sample.knn.edges.df <- do.call(rbind, edge.list)

      knn.g <- igraph::graph_from_data_frame(sample.knn.edges.df, directed = FALSE)
      knn.g <- igraph::simplify(knn.g, edge.attr.comb = "mean")

      for (res in resolutions) {
        res.key <- paste0("res_", sprintf("%.1f", res))
        cl <- tryCatch(igraph::cluster_louvain(knn.g, resolution = res), error = function(e) NULL)
        if (!is.null(cl)) {
          memb <- igraph::membership(cl)
          all.louvain[[res.key]] <- as.character(memb[mat.samples])
        } else {
          all.louvain[[res.key]] <- rep("1", length(mat.samples))
        }
      }
    } else {
      embed.dim1 <- rep(0, length(mat.samples))
      embed.dim2 <- rep(0, length(mat.samples))
      embed.label <- c("Dim 1", "Dim 2")
      for (res in resolutions) {
        res.key <- paste0("res_", sprintf("%.1f", res))
        all.louvain[[res.key]] <- rep("1", length(mat.samples))
      }
    }
  }

  set.seed(42)
  all.kmeans = lapply(2:8, function(k) {
    if (nrow(pt.pca.x) >= k) {
      tryCatch(
        stats::kmeans(pt.pca.x[, 1:min(10, n.pcs.use)], centers = k, nstart = 25, iter.max = 100),
        error = function(e) list(cluster = rep(1, nrow(pt.pca.x)))
      )
    } else {
      list(cluster = rep(1, nrow(pt.pca.x)))
    }
  })
  names(all.kmeans) = paste0("k", 2:8)

  default.cluster <- if (!is.null(all.louvain$res_1.0)) paste0("Cluster ", all.louvain$res_1.0) else "Cluster 1"

  sample.coords.base <- data.frame(
    sample = mat.samples,
    Dim1 = embed.dim1,
    Dim2 = embed.dim2,
    PC1 = pt.pca.x[, 1],
    PC2 = pt.pca.x[, 2],
    cluster = default.cluster,
    stringsAsFactors = FALSE
  )

  modality.progeny.scores <- Filter(Negate(is.null), modality.pathway.scores)
  if (length(modality.progeny.scores) > 0) {
    progeny.mean.data <- dplyr::bind_rows(lapply(names(modality.progeny.scores), function(pop) tibble::tibble(population = pop, pathway = colnames(modality.progeny.scores[[pop]]), mean.score = colMeans(modality.progeny.scores[[pop]]))))
  } else {
    progeny.mean.data <- tibble::tibble()
  }

  functional.program.scores.filt <- Filter(Negate(is.null), functional.program.scores)
  if (length(functional.program.scores.filt) > 0) {
    program.mean.data <- dplyr::bind_rows(lapply(names(functional.program.scores.filt), function(pop) tibble::tibble(population = pop, program = colnames(functional.program.scores.filt[[pop]]), mean.score = colMeans(functional.program.scores.filt[[pop]]))))
  } else {
    program.mean.data <- tibble::tibble()
  }

  if (length(functional.program.scores.filt) == 0 || ncol(functional.program.scores.filt[[1]]) == 0) {
    pt.top.prog <- data.frame(sample = sample.ids.vec, `Top Program` = NA, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    pt.top.prog <- data.frame(sample = sample.ids.vec, stringsAsFactors = FALSE) %>%
      dplyr::left_join(dplyr::bind_rows(lapply(names(functional.program.scores.filt), function(pop) {
        mat <- functional.program.scores.filt[[pop]]
        if (is.null(mat) || nrow(mat) == 0) {
          return(NULL)
        }
        df <- as.data.frame(mat)
        df$sample <- rownames(mat)
        tidyr::pivot_longer(df, cols = -sample, names_to = "Program", values_to = "Score")
      })) %>% dplyr::group_by(sample, Program) %>% dplyr::summarize(TotalScore = sum(Score, na.rm = TRUE), .groups = "drop") %>% dplyr::group_by(sample) %>% dplyr::slice_max(TotalScore, n = 1, with_ties = FALSE) %>% dplyr::select(sample, `Top Program` = Program) %>% dplyr::ungroup(), by = "sample")
  }

  if (length(modality.progeny.scores) == 0 || ncol(modality.progeny.scores[[1]]) == 0) {
    pt.top.path <- data.frame(sample = sample.ids.vec, `Top Pathway` = NA, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    pt.top.path <- data.frame(sample = sample.ids.vec, stringsAsFactors = FALSE) %>%
      dplyr::left_join(dplyr::bind_rows(lapply(names(modality.progeny.scores), function(pop) {
        mat <- modality.progeny.scores[[pop]]
        if (is.null(mat) || nrow(mat) == 0) {
          return(NULL)
        }
        df <- as.data.frame(mat)
        df$sample <- rownames(mat)
        tidyr::pivot_longer(df, cols = -sample, names_to = "Pathway", values_to = "Score")
      })) %>% dplyr::group_by(sample, Pathway) %>% dplyr::summarize(TotalScore = sum(Score, na.rm = TRUE), .groups = "drop") %>% dplyr::group_by(sample) %>% dplyr::slice_max(TotalScore, n = 1, with_ties = FALSE) %>% dplyr::select(sample, `Top Pathway` = Pathway) %>% dplyr::ungroup(), by = "sample")
  }

  sample.stats.summary <- final.df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarize(`Number of interactions` = dplyr::n(), `Top Ligand` = source_genesymbol[which.max(final.score)], `Top Receptor` = target_genesymbol[which.max(final.score)], `Top Cell-Cell Interaction` = paste0(gsub("scRNAseq ", "", sender[which.max(final.score)]), " \u2192 ", gsub("scRNAseq ", "", receiver[which.max(final.score)])), `Mean Edge Weight` = round(mean(final.score, na.rm = TRUE), 4), .groups = "drop") %>%
    dplyr::left_join(pt.top.prog, by = "sample") %>%
    dplyr::left_join(pt.top.path, by = "sample")

  weight.dir.sum <- weight.direction.summary
  if (nrow(weight.dir.sum) == 0) weight.dir.sum <- data.frame(Boosted = numeric(), Suppressed = numeric(), Neutral = numeric())
  if (!"Boosted" %in% colnames(weight.dir.sum)) weight.dir.sum$Boosted <- 0
  if (!"Suppressed" %in% colnames(weight.dir.sum)) weight.dir.sum$Suppressed <- 0
  if (!"Neutral" %in% colnames(weight.dir.sum)) weight.dir.sum$Neutral <- 0

  tmp.rec.pre <- dplyr::mutate(final.df, receiver.source = dplyr::case_when(
    abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.pathways) > 1e-10 & abs(.data$weight.receiver.programs) > 1e-10 ~ "Genomic + Pathways + Programs",
    abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.pathways) > 1e-10 ~ "Genomic + Pathways",
    abs(.data$weight.receiver.genomic) > 1e-10 & abs(.data$weight.receiver.programs) > 1e-10 ~ "Genomic + Programs",
    abs(.data$weight.receiver.genomic) > 1e-10 ~ "Genomic only",
    abs(weight.pathways) > 1e-10 & abs(weight.receiver.programs) > 1e-10 ~ "Pathways + Programs",
    abs(weight.pathways) > 1e-10 ~ "Pathways only",
    abs(weight.receiver.programs) > 1e-10 ~ "Programs only",
    TRUE ~ "Baseline only"
  ))
  receiver.source.summary.pre <- dplyr::mutate(dplyr::count(tmp.rec.pre, receiver.source), pct = round(n / sum(n) * 100, 1))

  edge.gradient.summary = NULL
  if (!is.null(gnn.embeddings) && length(unique(default.cluster)) > 1) {
    message("Computing GNN Edge Integrated Gradients across patient clusters...")
    cluster.df = data.frame(sample = mat.samples, cluster = default.cluster, stringsAsFactors = FALSE)
    edge.gradient.summary = tryCatch(
      {
        RunInteractEdgeGradients(
          edge.tensors = final.df,
          clusters = cluster.df,
          conda.env = conda.env,
          epochs = 50L,
          steps = 25L,
          save.csv = FALSE
        )
      },
      error = function(e) {
        warning("Could not compute edge integrated gradients: ", e$message)
        NULL
      }
    )
  }

  precomputed_data <- list(
    sample.ids.vec = sample.ids.vec, all.populations = all.populations, pop.colors = pop.colors, all.gene.symbols = all.gene.symbols,
    sender.receiver.agg = as.data.frame(sender.receiver.agg), top.lr.global = as.data.frame(top.lr.global), top.ecosystem.events = as.data.frame(top.ecosystem.events),
    pt.profile.long = as.data.frame(pt.profile.long), sample.interaction.counts = as.data.frame(sample.interaction.counts),
    pt.pca.x = pt.pca.x, n.pcs.use = n.pcs.use, embed.dim1 = embed.dim1, embed.dim2 = embed.dim2, embed.label = embed.label,
    all.kmeans = all.kmeans, all.louvain = all.louvain, sample.knn.edges = sample.knn.edges.df, gnn.embeddings = gnn.embeddings,
    sample.coords.base = sample.coords.base, modality.progeny.scores = modality.progeny.scores, progeny.mean.data = progeny.mean.data,
    functional.program.scores = functional.program.scores.filt, program.mean.data = program.mean.data, ligand.axes.available = ligand.axes.available,
    input.signal.summary = input.signal.summary, weight.layer.summary = weight.layer.summary, weight.direction.summary = weight.dir.sum,
    rank.shift.summary = rank.shift.summary, receiver.source.summary = receiver.source.summary.pre, matched.program.summary = data.frame(), sample.stats.summary = as.data.frame(sample.stats.summary),
    mean.layer.weights.summary = mean.layer.weights.summary,
    edge.gradient.summary = edge.gradient.summary
  )

  # 9. Topology Metrics Calculation ----
  message(sprintf("[%d/%d] Calculating network topology metrics...", c.step, total.steps))
  c.step <- c.step + 1

  # Flooring final aggregated sum.score to strictly non-negative values to prevent Dijkstra/PageRank crashes
  cell.network.edges <- final.df |>
    dplyr::group_by(sample, sender, receiver) |>
    dplyr::summarize(sum.score = sum(final.score, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(sender != receiver) |>
    dplyr::mutate(sum.score = pmax(sum.score, 1e-6))

  cell.network.edges <- cell.network.edges |>
    dplyr::group_by(sample) |>
    dplyr::mutate(score.threshold = stats::quantile(sum.score, probs = 0.75)) |>
    dplyr::filter(sum.score >= score.threshold) |>
    dplyr::ungroup() |>
    dplyr::select(-score.threshold)

  calculate.patient.topology <- function(patient.edges) {
    current.sample <- unique(patient.edges$sample)[1]
    patient.graph <- igraph::graph_from_data_frame(d = patient.edges |> dplyr::select(sender, receiver, sum.score), directed = TRUE)

    # Distance weights safely computed (max distance 1,000,000)
    igraph::E(patient.graph)$weight <- igraph::E(patient.graph)$sum.score
    distance.weights <- 1 / igraph::E(patient.graph)$weight

    hits.results <- igraph::hits_scores(patient.graph)
    undirected.graph <- igraph::as_undirected(patient.graph, mode = "collapse", edge.attr.comb = list(weight = "mean", "ignore"))
    barrat.scores <- igraph::transitivity(undirected.graph, type = "barrat")

    in.str <- igraph::strength(patient.graph, mode = "in")
    out.str <- igraph::strength(patient.graph, mode = "out")
    raw.betweenness <- igraph::betweenness(patient.graph, directed = TRUE, weights = distance.weights)
    norm.betweenness <- if (max(raw.betweenness, na.rm = TRUE) > 0) raw.betweenness / max(raw.betweenness, na.rm = TRUE) else raw.betweenness

    category <- dplyr::case_when(
      norm.betweenness >= 0.3 ~ "Mediator",
      hits.results$hub > 0.5 & hits.results$auth > 0.5 ~ "Influencer",
      out.str > in.str * 1.5 ~ "Sender",
      in.str > out.str * 1.5 ~ "Receiver",
      TRUE ~ "Influencer"
    )

    node.metrics <- data.frame(
      cell.type = igraph::V(patient.graph)$name, degree = igraph::degree(patient.graph, mode = "all"), strength = igraph::strength(patient.graph, mode = "all"),
      betweenness = raw.betweenness, pagerank = igraph::page_rank(patient.graph, directed = TRUE)$vector, hub.score = hits.results$hub, authority.score = hits.results$auth,
      transitivity.barrat = barrat.scores, category = category, sample = current.sample, stringsAsFactors = FALSE
    )
    graph.metrics <- data.frame(sample = current.sample, density = igraph::edge_density(patient.graph), transitivity.global = igraph::transitivity(patient.graph, type = "global"), transitivity.barrat.mean = mean(barrat.scores, na.rm = TRUE), stringsAsFactors = FALSE)
    edge.metrics <- igraph::as_data_frame(patient.graph, what = "edges") |>
      dplyr::select(sender = from, receiver = to, sum.score) |>
      dplyr::mutate(edge.betweenness = igraph::edge_betweenness(patient.graph, directed = TRUE, weights = distance.weights), sample = current.sample)

    return(list(nodes = node.metrics, edges = edge.metrics, graph = graph.metrics))
  }

  patient.list <- split(cell.network.edges, cell.network.edges$sample)

  # TryCatch wrapper prevents a single problematic patient graph from failing the run
  topology.results <- lapply(patient.list, function(pe) tryCatch(calculate.patient.topology(pe), error = function(e) NULL))
  topology.results <- Filter(Negate(is.null), topology.results)

  if (length(topology.results) == 0) {
    node.metrics <- data.frame()
    edge.metrics <- data.frame()
    graph.metrics <- data.frame()
  } else {
    node.metrics <- dplyr::bind_rows(lapply(topology.results, function(x) x$nodes))
    edge.metrics <- dplyr::bind_rows(lapply(topology.results, function(x) x$edges))
    graph.metrics <- dplyr::bind_rows(lapply(topology.results, function(x) x$graph))
  }

  if (!is.null(node.metrics) && nrow(node.metrics) > 0) rownames(node.metrics) <- NULL
  if (!is.null(edge.metrics) && nrow(edge.metrics) > 0) rownames(edge.metrics) <- NULL
  if (!is.null(graph.metrics) && nrow(graph.metrics) > 0) rownames(graph.metrics) <- NULL

  topology.metrics <- list(node = node.metrics, edge = edge.metrics, graph = graph.metrics)

  diagnostics <- list(
    input.signal = input.signal.summary,
    weight.layers = weight.layer.summary,
    receiver.source = receiver.source.summary,
    sender.source = sender.source.summary,
    weight.direction = weight.direction.summary,
    rank.shift = rank.shift.summary,
    mean.layer.weights = mean.layer.weights.summary,
    genomic.hits = genomic.hits.summary,
    calibration.parameters = list(
      scale.pathways = ctx$scale.pathways,
      scale.programs = ctx$scale.programs
    )
  )

  elapsed <- round(difftime(Sys.time(), start.time, units = "mins"), 2)
  message(sprintf("\nAnalysis Complete in %.2f minutes. %d interactions found.", elapsed, if (!is.null(final.df)) nrow(final.df) else 0))

  return(list(
    interactions = final.df,
    diagnostics = diagnostics,
    embeddings = gnn.embeddings,
    precompute = precomputed_data,
    modality.progeny.scores = modality.pathway.scores,
    functional.program.scores = functional.program.scores,
    program.gene.lists = program.gene.lists,
    topology.metrics = topology.metrics,
    edge.gradients = edge.gradient.summary
  ))
}

# --- Core Logic Helpers (Internal) --- ----

#' @keywords internal
get.expr.optimized <- function(parsing.list, expr.vec) {
  vapply(parsing.list, function(idx) {
    vals <- expr.vec[idx]
    if (any(is.na(vals))) {
      return(0)
    }
    min(vals, na.rm = TRUE)
  }, numeric(1))
}

#' @keywords internal
get.snv.optimized <- function(parsing.list, snv.vec) {
  vapply(parsing.list, function(idx) {
    vals <- snv.vec[idx]
    if (any(is.na(vals))) {
      return(0)
    }
    # If any subunit is repressive, the complex is broken (-1)
    min_val <- min(vals, na.rm = TRUE)
    if (min_val < 0) {
      return(min_val)
    }
    max(vals, na.rm = TRUE)
  }, numeric(1))
}

#' @keywords internal
get.cnv.optimized <- function(parsing.list, cnv.vec) {
  vapply(parsing.list, function(idx) {
    vals <- cnv.vec[idx]
    if (any(is.na(vals))) {
      return(2.0)
    }
    min(vals, na.rm = TRUE)
  }, numeric(1))
}

#' @keywords internal
compute.pt.base.state <- function(pt.id, ctx) {
  present.pops <- names(ctx$pop.idx)[vapply(names(ctx$pop.idx), function(p) pt.id %in% rownames(ctx$mats[[p]]), logical(1))]
  if (length(present.pops) == 0) {
    return(NULL)
  }

  pt.pop.expr <- lapply(stats::setNames(present.pops, present.pops), function(pop) {
    get.expr.optimized(ctx$parsing.list, ctx$mats[[pop]][pt.id, ctx$pop.idx[[pop]]])
  })

  # Extract raw genomic state purely via the strictly separated gen.idx
  pt.snv.vec <- NULL
  pt.cnv.vec <- NULL
  if (ctx$apply.genomic && ctx$tumor.pop %in% present.pops) {
    if ("SNV" %in% names(ctx$gen.mats) && pt.id %in% rownames(ctx$gen.mats$SNV)) {
      pt.snv.vec <- get.snv.optimized(ctx$parsing.list, ctx$gen.mats$SNV[pt.id, ctx$gen.idx$SNV])
    }
    if ("CNV" %in% names(ctx$gen.mats) && pt.id %in% rownames(ctx$gen.mats$CNV)) {
      pt.cnv.vec <- get.cnv.optimized(ctx$parsing.list, ctx$gen.mats$CNV[pt.id, ctx$gen.idx$CNV])
    }
  }

  unit.weights <- stats::setNames(rep(0, length(ctx$unique.symbols)), ctx$unique.symbols)

  # Optimized O(1) Precomputed Signature Index Lookup
  pop.prog.scores <- lapply(stats::setNames(present.pops, present.pops), function(pop) {
    lig.vec <- unit.weights
    rec.vec <- unit.weights
    if (pop %in% names(ctx$all.prog.scores) && !is.null(ctx$all.prog.scores[[pop]]) && pt.id %in% rownames(ctx$all.prog.scores[[pop]])) {
      r.vals <- ctx$all.prog.scores[[pop]][pt.id, ]

      lig.vec[] <- vapply(ctx$complex.lig.sigs, function(sigs) {
        if (length(sigs) == 0) {
          return(0)
        }
        vals <- r.vals[sigs]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
          return(0)
        }
        max(vals)
      }, numeric(1))

      rec.vec[] <- vapply(ctx$complex.rec.sigs, function(sigs) {
        if (length(sigs) == 0) {
          return(0)
        }
        vals <- r.vals[sigs]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
          return(0)
        }
        max(vals)
      }, numeric(1))
    }
    list(ligand = lig.vec, receptor = rec.vec)
  })

  pop.path.scores <- lapply(stats::setNames(present.pops, present.pops), function(pop) {
    s.vec <- unit.weights
    if (pop %in% names(ctx$all.path.scores) && !is.null(ctx$all.path.scores[[pop]]) && pt.id %in% rownames(ctx$all.path.scores[[pop]])) {
      r.vals <- ctx$all.path.scores[[pop]][pt.id, ]

      s.vec[] <- vapply(ctx$complex.rec.paths, function(paths) {
        if (length(paths) == 0) {
          return(0)
        }
        vals <- r.vals[paths]
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
          return(0)
        }
        max(vals)
      }, numeric(1))
    }
    s.vec
  })

  l.idx <- ctx$lr.src.idx
  r.idx <- ctx$lr.tgt.idx

  pop.pairs <- expand.grid(sender = present.pops, receiver = present.pops, stringsAsFactors = FALSE)
  pop.pairs <- dplyr::distinct(pop.pairs)

  base.results <- list()

  for (i in seq_len(nrow(pop.pairs))) {
    s <- pop.pairs$sender[i]
    r <- pop.pairs$receiver[i]

    s.expr <- pt.pop.expr[[s]]
    r.expr <- pt.pop.expr[[r]]

    base.prod <- s.expr[l.idx] + r.expr[r.idx]
    invalid <- s.expr[l.idx] < ctx$min.expression | r.expr[r.idx] < ctx$min.expression
    base.prod[invalid] <- 0

    valid <- which(base.prod > 0 & !is.na(base.prod))
    if (length(valid) == 0) next

    res <- ctx$lr.df[valid, ]
    res$base.score <- base.prod[valid]
    res$sample <- pt.id
    res$sender <- s
    res$receiver <- r

    # Pass the raw unscaled genomic values into the edge frame seamlessly
    res$s.snv <- if (s == ctx$tumor.pop && !is.null(pt.snv.vec)) pt.snv.vec[l.idx[valid]] else 0
    res$s.cnv <- if (s == ctx$tumor.pop && !is.null(pt.cnv.vec)) pt.cnv.vec[l.idx[valid]] else 2.0

    res$r.snv <- if (r == ctx$tumor.pop && !is.null(pt.snv.vec)) pt.snv.vec[r.idx[valid]] else 0
    res$r.cnv <- if (r == ctx$tumor.pop && !is.null(pt.cnv.vec)) pt.cnv.vec[r.idx[valid]] else 2.0

    res$s.prog <- pop.prog.scores[[s]]$ligand[l.idx[valid]]
    res$r.prog <- pop.prog.scores[[r]]$receptor[r.idx[valid]]
    res$r.path <- pop.path.scores[[r]][r.idx[valid]]

    base.results[[paste(s, r, sep = "__")]] <- res
  }

  dplyr::bind_rows(base.results)
}

#' @keywords internal
apply.functional.weights <- function(base.df, ctx) {
  if (is.null(base.df) || nrow(base.df) == 0) {
    return(NULL)
  }

  # 1. Biological Knockout Rule: Hard removal of edges with homozygous deletions
  is.destroyed <- (base.df$s.cnv < 0.5) | (base.df$r.cnv < 0.5)
  base.df <- base.df[!is.destroyed, , drop = FALSE]
  if (nrow(base.df) == 0) {
    return(NULL)
  }

  # ============================================================================
  # The denominators `scale.pathways` and `scale.programs` shrink raw absolute
  # PROGENy and ssGSEA scores into the responsive [-2, +2] range of tanh().
  # Uses optimal bounds derived automatically in step 4 during runtime.
  # ============================================================================

  # 2. Scale and bound Pathway scores.
  weight.pathways <- tanh(base.df$r.path / ctx$scale.pathways) * ctx$beta.pathways

  # Programs bounded JOINTLY to prevent doubled caps
  prog.sum <- base.df$s.prog + base.df$r.prog
  weight.programs <- tanh(prog.sum / ctx$scale.programs) * ctx$beta.programs

  # Extract proportional split for pure diagnostic tracking (preserves contribution signs)
  frac_s <- ifelse(prog.sum == 0, 0.5, base.df$s.prog / prog.sum)
  frac_r <- ifelse(prog.sum == 0, 0.5, base.df$r.prog / prog.sum)
  weight.sender.programs <- weight.programs * frac_s
  weight.receiver.programs <- weight.programs * frac_r

  # 3. Structural Genomic Modifiers for the 1D visualization score
  # Activating variants add +1.0 in log space, Repressive variants subtract -1.0
  weight.sender.genomic <- (base.df$s.snv * 1.0) + tanh(log2(pmax(base.df$s.cnv, 0.5) / 2.0))
  weight.receiver.genomic <- (base.df$r.snv * 1.0) + tanh(log2(pmax(base.df$r.cnv, 0.5) / 2.0))
  weight.genomic <- weight.sender.genomic + weight.receiver.genomic

  # 4. Re-populate Diagnostic Columns (Ensures Steps 6 & 7 do not crash)
  base.df$weight.pathways <- weight.pathways
  base.df$weight.programs <- weight.programs
  base.df$weight.genomic <- weight.genomic
  base.df$weight.sender.programs <- weight.sender.programs
  base.df$weight.receiver.programs <- weight.receiver.programs
  base.df$weight.sender.genomic <- weight.sender.genomic
  base.df$weight.receiver.genomic <- weight.receiver.genomic

  # 5. Final Additive log-space integration
  # (Negative values here are mathematically valid for machine learning, but handled via pmax() for topology)
  base.df$final.score <- base.df$base.score + weight.pathways + weight.programs + weight.genomic

  # 6. Clean up temporary unscaled scores but KEEP raw genomic attributes for the PyG Tensor
  cols_to_keep <- setdiff(colnames(base.df), c("s.prog", "r.prog", "r.path"))
  base.df <- base.df[, cols_to_keep, drop = FALSE]

  return(base.df)
}

#' @keywords internal
compute.modality.ssgsea <- function(pop.name, mat, program.sigs, mod.map, prog.const) {
  if (is.null(mat) || nrow(mat) < 5) {
    return(NULL)
  }

  current.cat <- mod.map$lineage.type[mod.map$modality == pop.name]
  valid.sigs <- prog.const$signature[vapply(prog.const$required.category, function(rc) rc == "Broad" || rc == current.cat, logical(1))]
  sigs.filtered <- lapply(program.sigs[valid.sigs], function(g) intersect(g, colnames(mat)))
  sigs.filtered <- sigs.filtered[vapply(sigs.filtered, length, numeric(1)) >= 3]
  if (length(sigs.filtered) == 0) {
    return(NULL)
  }

  # normalize = FALSE explicitly added to prevent GSVA batch leakage
  as.matrix(t(GSVA::gsva(GSVA::ssgseaParam(
    exprData = t(mat),
    geneSets = sigs.filtered,
    minSize = 3,
    maxSize = 250,
    normalize = FALSE
  ), verbose = FALSE)))
}

#' @keywords internal
get.default.pathway.map <- function() {
  list(
    Androgen = c("AR"), EGFR = c("EGFR", "ERBB2", "ERBB3", "ERBB4"), Estrogen = c("ESR1", "ESR2", "GPER1"),
    "JAK-STAT" = c("IL2RA", "IL2RB", "IL2RG", "IL4R", "IL6R", "IL6ST", "IL7R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL21R", "IL23R", "IL27RA", "IL31RA", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "CSF1R", "CSF2RA", "CSF2RB", "CSF3R", "EPOR", "MPL", "LEPR", "OSMR", "LIFR", "CNTFR"),
    MAPK = c("EGFR", "ERBB2", "ERBB3", "ERBB4", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PDGFRA", "PDGFRB", "KIT", "MET", "RET", "NTRK1", "NTRK2", "NTRK3", "FLT3", "KDR", "FLT1", "FLT4", "IGF1R", "INSR"),
    NFkB = c("TNFRSF1A", "TNFRSF1B", "TNFRSF5", "TNFRSF11A", "TNFRSF13B", "TNFRSF17", "TNFRSF14", "LTBR", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFRSF12A", "TNFRSF13C", "TNFRSF18", "TNFRSF25", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TLR10", "IL1R1", "IL1RAP", "IL1RL1", "IL1RL2", "IL18R1", "IL18RAP", "IL17RA", "IL17RB", "IL17RC", "IL17RE"),
    PI3K = c("EGFR", "ERBB2", "ERBB3", "ERBB4", "IGF1R", "INSR", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "PDGFRA", "PDGFRB", "KIT", "MET", "KDR", "FLT1", "FLT4", "FLT3", "CSF1R", "IL2RA", "IL2RB", "IL2RG", "IL7R", "IL15RA", "IL21R", "CXCR4", "CCR1", "CCR5"),
    TGFb = c("TGFBR1", "TGFBR2", "TGFBR3", "ACVR1", "ACVR1B", "ACVR1C", "ACVRL1", "ACVR2A", "ACVR2B", "BMPR1A", "BMPR1B", "BMPR2", "ENG"),
    TNFa = c("TNFRSF1A", "TNFRSF1B"), Trail = c("TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D"),
    VEGF = c("KDR", "FLT1", "FLT4", "NRP1", "NRP2"), WNT = c("FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10", "LRP5", "LRP6", "ROR1", "ROR2", "RYK")
  )
}

#' @keywords internal
is.mappable <- function(symbols, feature.list) {
  clean <- gsub("COMPLEX:", "", symbols)
  feature.list.upper <- toupper(feature.list)
  vapply(clean, function(s) {
    parts <- unlist(strsplit(s, "_", fixed = TRUE))
    all(toupper(parts) %in% feature.list.upper)
  }, logical(1))
}

#' @keywords internal
pick.axis <- function(categories, consensus.scores) {
  if (length(categories) == 0) {
    return("Other")
  }
  axis.priority <- c("checkpoint", "cytokine", "chemokine", "interleukin", "growth_factor", "adhesion", "receptor", "ligand")
  options.df <- tibble::tibble(axis = categories, score = consensus.scores)
  options.df <- dplyr::mutate(options.df, priority = match(.data$axis, axis.priority, nomatch = 100))
  options.df <- dplyr::arrange(options.df, dplyr::desc(.data$score), .data$priority)
  options.df$axis[1]
}

#' @keywords internal
get.valid.programs <- function(program.gene.lists, unique.symbols) {
  reactome.map <- as.list(reactome.db::reactomePATHID2EXTID)

  reactome.long <- dplyr::filter(program.gene.lists, .data$source == "Reactome")
  reactome.long <- dplyr::mutate(reactome.long, entrez.id = lapply(.data$id, function(x) reactome.map[[x]]))
  reactome.long <- tidyr::unnest(reactome.long, .data$entrez.id)

  if (nrow(reactome.long) > 0) {
    reactome.long$gene.symbol <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(reactome.long$entrez.id), column = "SYMBOL", keytype = "ENTREZID"))
  } else {
    reactome.long$gene.symbol <- character()
  }

  go.ids <- program.gene.lists$id[program.gene.lists$source == "GO"]
  if (length(go.ids) > 0) {
    go.long <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = go.ids, columns = "SYMBOL", keytype = "GOALL"))
    go.long <- dplyr::rename(go.long, id = .data$GOALL, gene.symbol = .data$SYMBOL)
    go.long <- dplyr::left_join(go.long, dplyr::select(program.gene.lists, .data$id, .data$program.name), by = "id")
  } else {
    go.long <- data.frame(id = character(), gene.symbol = character(), program.name = character())
  }

  comb <- dplyr::bind_rows(dplyr::select(reactome.long, .data$program.name, .data$gene.symbol), dplyr::select(go.long, .data$program.name, .data$gene.symbol))
  comb <- dplyr::filter(comb, !is.na(.data$gene.symbol))
  comb <- dplyr::distinct(comb)
  comb <- dplyr::group_by(comb, .data$program.name)
  comb <- dplyr::filter(comb, any(.data$gene.symbol %in% unique.symbols), dplyr::n() >= 10, dplyr::n() <= 250)
  comb <- dplyr::ungroup(comb)
  comb
}

#' @keywords internal
get.lineage.map <- function(mods) {
  df <- tibble::tibble(modality = mods)
  dplyr::mutate(df, lineage.type = dplyr::case_when(
    grepl("Tumor|Plasma cells|Plasmablasts|B cells|pro-B|pre-B", .data$modality, ignore.case = TRUE) ~ "B",
    grepl("T cells|MAIT|Tcm|Tem|Regulatory T", .data$modality, ignore.case = TRUE) ~ "T",
    grepl("NK cells", .data$modality, ignore.case = TRUE) ~ "NK",
    grepl("monocytes|macrophages|Myeloid|Neutrophil|Promyelocytes|DC|pDC", .data$modality, ignore.case = TRUE) ~ "Myeloid",
    grepl("erythroid", .data$modality, ignore.case = TRUE) ~ "Erythroid",
    TRUE ~ "Other"
  ))
}

#' @keywords internal
get.program.constraints <- function(sigs) {
  df <- tibble::tibble(signature = sigs)
  dplyr::mutate(df, required.category = dplyr::case_when(
    grepl("TCR signaling|\\bT cell\\b|alpha-beta T", .data$signature, ignore.case = TRUE) ~ "T",
    grepl("BCR|\\bB cell\\b|Ig production|Germinal center|Plasma cell|Plasmablast", .data$signature, ignore.case = TRUE) ~ "B",
    grepl("Natural killer|NK", .data$signature, ignore.case = TRUE) ~ "NK",
    grepl("monocyte|macrophage|Myeloid|phagocytosis|Complement", .data$signature, ignore.case = TRUE) ~ "Myeloid",
    TRUE ~ "Broad"
  ))
}

#' @keywords internal
can.support.fork = function() {
  if (.Platform$OS.type != "unix") {
    return(FALSE)
  }
  # Positron IDE / Ark kernel detection
  if (identical(Sys.getenv("POSITRON"), "1") ||
      nzchar(Sys.getenv("POSITRON_VERSION")) ||
      nzchar(Sys.getenv("ARK_PORT")) ||
      "ark" %in% loadedNamespaces() ||
      identical(getOption("positron.package.version"), TRUE) ||
      !is.null(getOption("ark.version")) ||
      !is.null(getOption("positron.eval"))) {
    return(FALSE)
  }
  # Check if forking is disabled via R option
  if (isTRUE(getOption("parallel.disable.mcfork", FALSE))) {
    return(FALSE)
  }
  # Check parallelly multicore support
  if (!parallelly::supportsMulticore()) {
    return(FALSE)
  }
  TRUE
}
