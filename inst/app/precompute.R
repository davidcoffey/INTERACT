# Integrated Immune-Tumor Signaling Viewer
# Precomputation Script — Run Once Offline
# Generates Data/app-precomputed.Rda and Data/samples/*.rds
# David Coffey davidcoffey@miami.edu

t.total = proc.time()
message("=== TME Viewer Precompute ===")

# 0. Libraries ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(ggsci)
  library(uwot)
  library(igraph)
})

# 1. Load Source Data ----
message("[1/8] Loading RData...")
t0 = proc.time()
load("Data/Sample-functional-interactions-hybrid-diagnostics.RData")
message(sprintf("      Done (%.1fs)", (proc.time() - t0)[3]))

message("[1/8] Converting to data.table...")
final.dt = data.table::as.data.table(final.df)
if (!"genomic.weighted" %in% names(final.dt)) {
  final.dt[, genomic.weighted := sender.genomic.weighted | target.genomic.weighted]
}
data.table::setkey(final.dt, sample)

message("[1/8] Loading clinical data...")
clinical.raw = readr::read_csv("Data/Merged clinical data.csv", show_col_types = FALSE)
clinical.sample = clinical.raw %>%
  dplyr::group_by(SPECTRUM_SEQ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

sig.alt.gene = tryCatch(
  readr::read_csv("Data/Genomic alterations impacting gene expression.csv", show_col_types = FALSE),
  error = function(e) tibble::tibble(feature = character(), alteration = character(), med.fc = numeric())
)

# 2. Pre-Aggregate Interaction Data ----
message("[2/8] Pre-aggregating interaction data...")

sender.receiver.agg = final.dt[, .(
  mean.score  = mean(final.score, na.rm = TRUE),
  total.score = sum(final.score, na.rm = TRUE),
  n.pairs     = .N,
  pct.genomic = mean(as.numeric(genomic.weighted), na.rm = TRUE) * 100
), by = .(sender, receiver)]

top.lr.global = final.dt[, .(
  mean.score  = mean(final.score, na.rm = TRUE),
  n.samples  = data.table::uniqueN(sample),
  pct.genomic = mean(as.numeric(genomic.weighted), na.rm = TRUE) * 100
), by = .(source_genesymbol, target_genesymbol, ligand.axis, receptor.axis)][
  order(-mean.score)
]

top.ecosystem.events = final.dt[, .(
  mean.final.score = mean(final.score, na.rm = TRUE),
  mean.base.score  = mean(base.score, na.rm = TRUE),
  n.samples       = data.table::uniqueN(sample),
  genomic.rate     = mean(as.numeric(genomic.weighted), na.rm = TRUE) * 100
), by = .(source_genesymbol, target_genesymbol, sender, receiver, ligand.axis, receptor.axis)][
  order(-mean.final.score)
][1:min(.N, 500)]

# Sample sender-receiver profile (for cluster enrichment)
pt.profile.long = final.dt[, .(
  total.score = sum(final.score, na.rm = TRUE)
), by = .(sample, sender, receiver)]

# 3. Sample Clustering ----
message("[3/8] Building sample interaction matrix and running PCA...")

pt.profile.wide = data.table::dcast(
  pt.profile.long,
  sample ~ paste0(sender, "__|__", receiver),
  value.var = "total.score",
  fill = 0
)

sample.ids.vec = pt.profile.wide$sample
pt.mat.raw = as.matrix(pt.profile.wide[, -"sample"])
rownames(pt.mat.raw) = sample.ids.vec

col.var = apply(pt.mat.raw, 2, var)
pt.mat.clean  = pt.mat.raw[, col.var > 0, drop = FALSE]
pt.mat.scaled = scale(pt.mat.clean)

pt.pca  = prcomp(pt.mat.scaled, center = FALSE, scale. = FALSE)
n.pcs.use = min(20, ncol(pt.pca$x))
# Store only the PC scores (not the full prcomp object)
pt.pca.x = pt.pca$x

message("[3/8] Running UMAP...")
set.seed(42)
umap.coords = uwot::umap(pt.pca.x[, 1:n.pcs.use], n_neighbors = 15, min_dist = 0.3)
embed.dim1  = umap.coords[, 1]
embed.dim2  = umap.coords[, 2]
embed.label = c("UMAP1", "UMAP2")

message("[3/8] Running k-means for k = 2:8...")
set.seed(42)
all.kmeans = lapply(2:8, function(k) {
  stats::kmeans(pt.pca.x[, 1:min(10, n.pcs.use)],
                centers = k, nstart = 25, iter.max = 100)
})
names(all.kmeans) = paste0("k", 2:8)

# 4. Sample Coordinates ----
message("[4/8] Building sample coordinates with clinical data...")

sample.coords.base = data.frame(
  sample  = sample.ids.vec,
  Dim1    = embed.dim1,
  Dim2    = embed.dim2,
  PC1     = pt.pca.x[, 1],
  PC2     = pt.pca.x[, 2],
  cluster = paste0("Cluster ", all.kmeans$k4$cluster),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    clinical.sample %>%
      dplyr::select(
        SPECTRUM_SEQ, PUBLIC_ID,
        DEMOG_GENDER, DEMOG_SAMPLEAGE,
        chr17p_Loss, chr1q_Gain, t4_14_MMSET, t14_16_MAF,
        t14_20_MAFB, t11_14_CCND1, t_MYC, best.response,
        derived.pfs.day, censpfs, derived.os.day, censos
      ),
    by = c("sample" = "SPECTRUM_SEQ")
  )

survival.df = sample.coords.base %>%
  dplyr::filter(!is.na(derived.pfs.day), !is.na(censpfs)) %>%
  dplyr::mutate(
    pfs.months = derived.pfs.day / 30.44,
    os.months  = dplyr::if_else(!is.na(derived.os.day), derived.os.day / 30.44, NA_real_)
  )

# 5. PROGENy Data ----
message("[5/8] Aggregating PROGENy scores...")

progeny.mean.data = lapply(names(modality.progeny.scores), function(pop) {
  mat = modality.progeny.scores[[pop]]
  if (is.null(mat) || nrow(mat) == 0) return(NULL)
  means = colMeans(mat, na.rm = TRUE)
  tibble::tibble(population = pop, pathway = names(means), mean.score = as.numeric(means))
}) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!population %in% c("CNV", "SNV"))

message("[5/8] Aggregating Program scores...")
program.mean.data = lapply(names(functional.program.scores), function(pop) {
  mat = functional.program.scores[[pop]]
  if (is.null(mat) || nrow(mat) == 0) return(NULL)
  means = colMeans(mat, na.rm = TRUE)
  tibble::tibble(population = pop, program = names(means), mean.score = as.numeric(means))
}) %>%
  dplyr::bind_rows()

# 6. Genomic / Circos Data ----
message("[6/8] Preparing circos and genomic data...")

driver.genes = if (!is.null(sig.alt.gene) && nrow(sig.alt.gene) > 0) {
  unique(sig.alt.gene$feature)
} else character(0)

circos.data = final.dt[, .(
  total.score = sum(final.score, na.rm = TRUE),
  n.samples  = data.table::uniqueN(sample),
  pct.genomic = mean(as.numeric(genomic.weighted), na.rm = TRUE) * 100
), by = .(source_genesymbol, target_genesymbol, ligand.axis, receptor.axis)][
  order(-total.score)
]

ligand.axes.available = sort(unique(circos.data$ligand.axis))

# All unique gene symbols (for UI dropdowns)
all.gene.symbols = sort(unique(c(
  final.df$source_genesymbol,
  final.df$target_genesymbol
)))

# 7. Color Palettes ----
all.populations = sort(unique(c(final.df$sender, final.df$receiver)))
n.pops = length(all.populations)
pop.palette = ggsci::pal_igv()(51)
pop.colors  = pop.palette[seq_len(min(n.pops, 51))]
if (n.pops > 51) pop.colors = grDevices::colorRampPalette(pop.palette)(n.pops)
names(pop.colors) = all.populations

clinical.var.labels = c(
  DEMOG_GENDER       = "Sex",
  DEMOG_SAMPLEAGE   = "Age at Diagnosis",
  chr17p_Loss        = "del(17p)",
  chr1q_Gain         = "gain(1q)",
  t4_14_MMSET        = "t(4;14) MMSET",
  t14_16_MAF         = "t(14;16) MAF",
  t14_20_MAFB        = "t(14;20) MAFB",
  t11_14_CCND1       = "t(11;14) CCND1",
  t_MYC              = "MYC Rearrangement",
  best.response      = "Best Response",
  derived.pfs.day    = "PFS (days)",
  derived.os.day     = "OS (days)",
  line               = "Treatment Line",
  regimen            = "Regimen"
)

# 8. Save Per-Sample rds Files ----
message("[7/8] Saving per-sample interaction files...")

samples.dir = "Data/samples"
if (!dir.exists(samples.dir)) dir.create(samples.dir)

sample.ids.all = unique(final.dt$sample)
n.samples = length(sample.ids.all)

pb = txtProgressBar(min = 0, max = n.samples, style = 3)
for (i in seq_along(sample.ids.all)) {
  pt.id  = sample.ids.all[i]
  pt.dt  = final.dt[pt.id]
  saveRDS(pt.dt, file = file.path(samples.dir, paste0(pt.id, ".rds")), compress = TRUE)
  setTxtProgressBar(pb, i)
}
close(pb)
message(sprintf("      Saved %d sample files to %s/", n.samples, samples.dir))

# Verify file sizes
pt.sizes = file.size(file.path(samples.dir, paste0(sample.ids.all, ".rds")))
message(sprintf("      Median file size: %.1f KB", median(pt.sizes) / 1024))

# 9. Calculate Interaction Counts ----
sample.interaction.counts = final.dt[, .N, by = sample]
data.table::setnames(sample.interaction.counts, "N", "interaction.count")

# 10. Save Precomputed App Data ----
message("[8/8] Saving app-precomputed.Rda...")

# Diagnostic objects are already loaded from the original RData
save(
  # Global app metadata
  sample.ids.vec,
  all.populations,
  pop.colors,
  all.gene.symbols,
  clinical.var.labels,

  # Clinical data
  clinical.sample,
  sig.alt.gene,

  # Interaction aggregations (Tabs 1, 2, 4)
  sender.receiver.agg,
  top.lr.global,
  top.ecosystem.events,
  pt.profile.long,
  functional.program.scores,
  program.mean.data,
  sample.interaction.counts,

  # Sample clustering (Tabs 2, 3, 6)
  pt.pca.x,
  n.pcs.use,
  embed.dim1,
  embed.dim2,
  embed.label,
  all.kmeans,
  sample.coords.base,

  # PROGENy (Tab 5)
  modality.progeny.scores,
  progeny.mean.data,

  # Survival (Tab 6)
  survival.df,

  # Genomic / circos (Tab 4)
  circos.data,
  ligand.axes.available,
  driver.genes,

  # Diagnostics — from original RData (Tab 7)
  input.signal.summary,
  weight.layer.summary,
  weight.direction.summary,
  rank.shift.summary,
  receiver.source.summary,

  # Biological validation (Tab 6)
  driver.report,
  driver.category.report,

  file = "Data/app-precomputed.Rda"
)

rda.size = file.size("Data/app-precomputed.Rda") / 1e6
message(sprintf("      Saved (%.1f MB)", rda.size))
message(sprintf("\nPrecompute complete in %.1f seconds.", (proc.time() - t.total)[3]))
message("Run shiny::runApp('.') to launch the optimized app.")
