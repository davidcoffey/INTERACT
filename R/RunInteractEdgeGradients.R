#' Compute GNN Edge Integrated Gradients
#'
#' @description Trains the multi-channel Graph Autoencoder (GINEConv) on patient
#' microenvironment networks and computes path-integrated gradients (attribution scores)
#' for each cell-cell interaction with respect to cluster centroids in latent embedding space.
#' This quantifies the exact contribution of each biological communication channel (RNA,
#' Reactome pathways, PROGENy programs, somatic SNVs, and CNVs) toward shaping each
#' patient cluster.
#'
#' @param edge.tensors Edge tensors data.table, data.frame, CSV file path, or an INTERACT results object.
#' @param node.features Node features data.table, data.frame, or CSV file path. Optional if edge.tensors is an INTERACT results object.
#' @param edge.tensors.val Optional validation edge tensors data.table, data.frame, or CSV file path (for cell mapping). Default: NULL.
#' @param clusters Patient cluster assignments data.table, data.frame, named vector, or CSV file path. Optional if edge.tensors contains clusters.
#' @param conda.env String. Conda environment containing PyTorch and PyTorch Geometric (default: "pytorch").
#' @param epochs Integer. Number of training epochs (default: 150L).
#' @param batch.size Integer. Mini-batch size (default: 32L).
#' @param steps Integer. Number of interpolation steps for Riemann path integral (default: 50L).
#' @param learning.rate Numeric. Optimizer learning rate (default: 0.005).
#' @param weight.decay Numeric. L2 weight decay regularization (default: 1e-4).
#' @param hidden.dim Integer. Hidden dimension for GINEConv MLPs (default: 32L).
#' @param emb.dim Integer. Latent embedding dimension per cell type (default: 4L).
#' @param channel.weights Numeric vector of length 5. Loss weights for `[RNA, Pathway, Program, SNV, CNV]` (default: c(1.0, 3.0, 3.0, 5.0, 5.0)).
#' @param seed Integer. Random seed for reproducibility (default: 42L).
#' @param deterministic Logical. If TRUE, enforces single-threaded deterministic operations for 100% reproducibility (default: TRUE).
#' @param auto.install Logical. If TRUE, automatically checks and installs required Python packages if missing (default: FALSE).
#' @param save.csv Logical. If TRUE, saves summary attributions to output.file (default: FALSE).
#' @param output.file String. Path to output CSV file (default: "Data/INTERACT Edge Integrated Gradients.csv").
#'
#' @return A data.frame of mean interaction attributions per cluster across all 5 biological channels.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run Integrated Gradients using INTERACT results
#' ig.results = RunInteractEdgeGradients(
#'   edge.tensors = results,
#'   conda.env = "pytorch",
#'   epochs = 150L
#' )
#' }
RunInteractEdgeGradients = function(
  edge.tensors = "Data/2 INTERACT Discovery Tensors.csv",
  node.features = "Data/2 INTERACT Discovery Node Features.csv",
  edge.tensors.val = NULL,
  clusters = "Data/6 Louvain Clusters.csv",
  conda.env = "pytorch",
  epochs = 150L,
  batch.size = 32L,
  steps = 50L,
  learning.rate = 0.005,
  weight.decay = 1e-4,
  hidden.dim = 32L,
  emb.dim = 4L,
  channel.weights = c(1.0, 3.0, 3.0, 5.0, 5.0),
  seed = 42L,
  deterministic = TRUE,
  auto.install = FALSE,
  save.csv = FALSE,
  output.file = "Data/INTERACT Edge Integrated Gradients.csv"
) {
  # Helper to load data from path, object, or data.table
  load.data = function(x) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) {
      if (!file.exists(x)) stop("File does not exist: ", x)
      return(as.data.frame(data.table::fread(x)))
    }
    return(as.data.frame(x))
  }

  extracted.clusters = NULL

  # Check if edge.tensors is an INTERACT results object
  if (is.list(edge.tensors) && !is.data.frame(edge.tensors) &&
      ("interactions" %in% names(edge.tensors) || "precompute" %in% names(edge.tensors))) {
    message("Extracting edge tensors and node features from INTERACT results object...")
    tensors.list = format.interact.tensors(edge.tensors)
    df.edge.disc = as.data.frame(tensors.list$edge.tensors)
    df.node.disc = as.data.frame(tensors.list$node.features)

    # Extract clusters if available and not explicitly provided
    if (is.null(clusters) || (is.character(clusters) && !file.exists(clusters))) {
      if (!is.null(edge.tensors$precompute$sample.coords.base$cluster)) {
        extracted.clusters = data.frame(
          sample = edge.tensors$precompute$sample.coords.base$sample,
          cluster = edge.tensors$precompute$sample.coords.base$cluster,
          stringsAsFactors = FALSE
        )
      } else if (!is.null(edge.tensors$precompute$all.louvain$res_1.0)) {
        extracted.clusters = data.frame(
          sample = edge.tensors$precompute$sample.ids.vec,
          cluster = paste0("Cluster ", edge.tensors$precompute$all.louvain$res_1.0),
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    df.edge.disc = load.data(edge.tensors)
    # If raw interactions data frame without tensor columns, format on the fly
    req.tensor.cols = c("edge_rna", "edge_pathway", "edge_program", "edge_snv", "edge_cnv")
    if (!all(req.tensor.cols %in% names(df.edge.disc))) {
      message("Formatting raw interaction table into multi-channel edge and node tensors...")
      tensors.list = format.interact.tensors(df.edge.disc)
      df.edge.disc = as.data.frame(tensors.list$edge.tensors)
      df.node.disc = as.data.frame(tensors.list$node.features)
    } else {
      df.node.disc = load.data(node.features)
    }
  }

  df.edge.val = load.data(edge.tensors.val)

  # Resolve cluster assignments
  if (!is.null(extracted.clusters)) {
    df.clusters = extracted.clusters
  } else if (!is.null(clusters)) {
    if (is.character(clusters) && file.exists(clusters)) {
      df.clusters = as.data.frame(data.table::fread(clusters))
    } else if (is.data.frame(clusters) || data.table::is.data.table(clusters)) {
      df.clusters = as.data.frame(clusters)
    } else if (is.vector(clusters) && !is.null(names(clusters))) {
      df.clusters = data.frame(
        sample = names(clusters),
        cluster = as.character(clusters),
        stringsAsFactors = FALSE
      )
    } else {
      stop("Could not parse cluster assignments. Please provide a data.frame with 'sample' and 'cluster' columns, a named vector, or a valid CSV file path.")
    }
  } else {
    stop("Cluster assignments must be provided via 'clusters' argument or within the INTERACT results object.")
  }

  # Ensure sample and cluster column names
  if (!"sample" %in% names(df.clusters) || !"cluster" %in% names(df.clusters)) {
    if (ncol(df.clusters) >= 2) {
      names(df.clusters)[1:2] = c("sample", "cluster")
    } else {
      stop("Clusters table must have 'sample' and 'cluster' columns.")
    }
  }

  # Python Environment Setup
  if (isTRUE(auto.install)) {
    InstallInteractEnv(conda.env = conda.env, method = "auto")
  }

  # Initialize Python environment
  custom.python = Sys.getenv("RETICULATE_PYTHON")
  conda.envs = tryCatch(reticulate::conda_list()$name, error = function(e) character())
  venvs = tryCatch(reticulate::virtualenv_list(), error = function(e) character())

  if (nzchar(custom.python) && file.exists(custom.python)) {
    message("Using Python specified by RETICULATE_PYTHON: ", custom.python)
    reticulate::use_python(custom.python, required = TRUE)
  } else if (!is.null(conda.env) && conda.env %in% conda.envs) {
    message("Initializing Conda environment: ", conda.env)
    reticulate::use_condaenv(conda.env, required = TRUE)
  } else if (!is.null(conda.env) && conda.env %in% venvs) {
    message("Initializing Virtualenv: ", conda.env)
    reticulate::use_virtualenv(conda.env, required = TRUE)
  }

  req.pkgs = c("torch", "torch_geometric", "pandas", "numpy", "sklearn")
  missing.pkgs = character(0)
  for (pkg in req.pkgs) {
    if (!reticulate::py_module_available(pkg)) {
      missing.pkgs = c(missing.pkgs, pkg)
    }
  }

  if (length(missing.pkgs) > 0) {
    stop(paste0(
      "\nMissing required Python package(s): ", paste(missing.pkgs, collapse = ", "),
      "\nPlease install them using:\n",
      "  INTERACT::InstallInteractEnv(conda.env = '", conda.env, "')\n",
      "or run RunInteractEdgeGradients(..., auto.install = TRUE)\n"
    ))
  }

  # Configure Python session and pass data
  reticulate::py_run_string("
import os
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GINEConv
from sklearn.preprocessing import StandardScaler
")

  py = reticulate::py
  py$r_edge_disc = df.edge.disc
  py$r_node_disc = df.node.disc
  py$r_edge_val = df.edge.val
  py$r_clusters = df.clusters
  py$r_epochs = as.integer(epochs)
  py$r_batch_size = as.integer(batch.size)
  py$r_steps = as.integer(steps)
  py$r_lr = as.numeric(learning.rate)
  py$r_wd = as.numeric(weight.decay)
  py$r_hidden_dim = as.integer(hidden.dim)
  py$r_emb_dim = as.integer(emb.dim)
  py$r_weights = as.numeric(channel.weights)
  py$r_seed = as.integer(seed)
  py$r_deterministic = as.logical(deterministic)

  message("Executing GNN Edge Integrated Gradients in PyTorch Geometric...")
  reticulate::py_run_string("
if r_deterministic:
    torch.set_num_threads(1)

torch.manual_seed(r_seed)
np.random.seed(r_seed)

class TME_Autoencoder(nn.Module):
    def __init__(self, num_node_features, hidden_dim, emb_dim, num_edge_features, num_cells):
        super(TME_Autoencoder, self).__init__()
        self.num_cells = num_cells
        self.emb_dim = emb_dim
        mlp1 = nn.Sequential(nn.Linear(num_node_features, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, hidden_dim))
        self.conv1 = GINEConv(nn=mlp1, edge_dim=num_edge_features)
        mlp2 = nn.Sequential(nn.Linear(hidden_dim, emb_dim), nn.ReLU(), nn.Linear(emb_dim, emb_dim))
        self.conv2 = GINEConv(nn=mlp2, edge_dim=num_edge_features)
        self.edge_decoder = nn.Sequential(nn.Linear(emb_dim * 2, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, num_edge_features))

    def encode(self, x, edge_index, edge_attr):
        h = F.elu(self.conv1(x, edge_index, edge_attr=edge_attr))
        return self.conv2(h, edge_index, edge_attr=edge_attr)

    def decode(self, z, edge_index):
        return self.edge_decoder(torch.cat([z[edge_index[0]], z[edge_index[1]]], dim=-1))

    def forward(self, x, edge_index, edge_attr, batch):
        z = self.encode(x, edge_index, edge_attr)
        return self.decode(z, edge_index), z.view(-1, self.num_cells * self.emb_dim)

def prep_dataset(edge_df, nodes_dict, feature_cols, cell_to_idx, tumor_idx, num_cells):
    dataset = []
    for pt, pt_edges in edge_df.groupby('sample'):
        senders = [cell_to_idx[s] for s in pt_edges['sender']]
        receivers = [cell_to_idx[r] for r in pt_edges['receiver']]
        edge_index = torch.tensor([senders, receivers], dtype=torch.long)
        attr_tensor = torch.tensor(pt_edges[feature_cols].values, dtype=torch.float32)
        x = np.eye(num_cells, dtype=np.float32)
        extra_genomic = np.zeros((num_cells, 2), dtype=np.float32)
        if tumor_idx is not None and pt in nodes_dict:
            extra_genomic[tumor_idx, 0] = nodes_dict[pt]['Tumor_Mutational_Burden']
            extra_genomic[tumor_idx, 1] = nodes_dict[pt]['Tumor_Aneuploidy_Var']
        x = np.hstack([x, extra_genomic])
        x_tensor = torch.tensor(x, dtype=torch.float32)
        dataset.append(Data(x=x_tensor, edge_index=edge_index, edge_attr=attr_tensor, sample_id=pt))
    return dataset

edge_df_disc = r_edge_disc.copy()
nodes_df_disc = r_node_disc.copy()

feature_cols = ['edge_rna', 'edge_pathway', 'edge_program', 'edge_snv', 'edge_cnv']
scaler = StandardScaler()
edge_df_disc[feature_cols] = scaler.fit_transform(edge_df_disc[feature_cols])

unique_cells = set(edge_df_disc['sender']).union(set(edge_df_disc['receiver']))
if r_edge_val is not None and len(r_edge_val) > 0:
    edge_df_val = r_edge_val.copy()
    unique_cells = unique_cells.union(set(edge_df_val['sender'])).union(set(edge_df_val['receiver']))

unique_cells = sorted(list(unique_cells))
num_cells = len(unique_cells)
cell_to_idx = {cell: i for i, cell in enumerate(unique_cells)}
tumor_idx = cell_to_idx.get('Tumor')

nodes_dict_disc = nodes_df_disc.set_index('sample').to_dict('index')
disc_dataset = prep_dataset(edge_df_disc, nodes_dict_disc, feature_cols, cell_to_idx, tumor_idx, num_cells)

loader = DataLoader(disc_dataset, batch_size=r_batch_size, shuffle=False)

node_dim = num_cells + 2
device = torch.device('cpu')
model = TME_Autoencoder(node_dim, r_hidden_dim, r_emb_dim, 5, num_cells).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=r_lr, weight_decay=r_wd)

# Train Autoencoder
model.train()
weights_t = torch.tensor(r_weights, dtype=torch.float32, device=device)
for epoch in range(r_epochs):
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred_edges, _ = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
        loss = torch.sum(torch.mean((pred_edges - batch.edge_attr)**2, dim=0) * weights_t)
        loss.backward()
        optimizer.step()

# Compute baseline embeddings
model.eval()
baseline_embs = []
patient_ids = []
with torch.no_grad():
    for batch in loader:
        batch = batch.to(device)
        _, z = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
        baseline_embs.append(z)
        patient_ids.extend(batch.sample_id)
baseline_embs = torch.cat(baseline_embs, dim=0)

clusters_df = r_clusters.copy()
cluster_dict = dict(zip(clusters_df['sample'], clusters_df['cluster']))

# Calculate cluster centroids
centroids = {}
for c_name in set(cluster_dict.values()):
    c_idx = [idx for idx, pid in enumerate(patient_ids) if cluster_dict.get(pid) == c_name]
    if len(c_idx) > 0:
        centroids[c_name] = baseline_embs[c_idx].mean(dim=0)

# Compute Integrated Gradients
importances = []
idx_to_cell = {v: k for k, v in cell_to_idx.items()}

for batch in loader:
    batch = batch.to(device)
    baseline_attr = torch.zeros_like(batch.edge_attr)
    edge_attr_target = batch.edge_attr.clone()
    total_grads = torch.zeros_like(batch.edge_attr)

    batch_centroids = [
        centroids[cluster_dict[pid]] if pid in cluster_dict and cluster_dict[pid] in centroids
        else torch.zeros_like(list(centroids.values())[0])
        for pid in batch.sample_id
    ]
    batch_centroids = torch.stack(batch_centroids).to(device)

    for alpha in np.linspace(1.0 / r_steps, 1.0, r_steps):
        attr_step = baseline_attr + alpha * (edge_attr_target - baseline_attr)
        attr_step.requires_grad_(True)
        _, z = model(batch.x, batch.edge_index, attr_step, batch.batch)
        scores = -torch.sum((z - batch_centroids)**2, dim=1)
        model.zero_grad()
        scores.sum().backward()
        total_grads += attr_step.grad.detach()

    avg_grads = total_grads / r_steps
    ig = (edge_attr_target - baseline_attr) * avg_grads

    edge_importance = ig.sum(dim=1).cpu().numpy()
    ig_rna = ig[:, 0].cpu().numpy()
    ig_pathway = ig[:, 1].cpu().numpy()
    ig_program = ig[:, 2].cpu().numpy()
    ig_snv = ig[:, 3].cpu().numpy()
    ig_cnv = ig[:, 4].cpu().numpy()

    graph_indices = batch.batch[batch.edge_index[0]].cpu().numpy()
    senders = batch.edge_index[0].cpu().numpy()
    receivers = batch.edge_index[1].cpu().numpy()
    senders_orig = senders % num_cells
    receivers_orig = receivers % num_cells

    for i in range(len(edge_importance)):
        g_idx = graph_indices[i]
        pid = batch.sample_id[g_idx]
        c_name = cluster_dict.get(pid)
        if c_name is not None:
            sender_name = idx_to_cell[senders_orig[i]]
            receiver_name = idx_to_cell[receivers_orig[i]]
            edge_name = f'{sender_name} -> {receiver_name}'
            importances.append({
                'cluster': c_name,
                'edge_name': edge_name,
                'sender': sender_name,
                'receiver': receiver_name,
                'ig_score': edge_importance[i],
                'ig_rna': ig_rna[i],
                'ig_pathway': ig_pathway[i],
                'ig_program': ig_program[i],
                'ig_snv': ig_snv[i],
                'ig_cnv': ig_cnv[i]
            })

df_imp = pd.DataFrame(importances)
df_summary = df_imp.groupby(['cluster', 'edge_name', 'sender', 'receiver'])[
    ['ig_score', 'ig_rna', 'ig_pathway', 'ig_program', 'ig_snv', 'ig_cnv']
].mean().reset_index()

df_summary = df_summary.sort_values(['cluster', 'ig_score'], ascending=[True, False])
py_ig_summary = df_summary
")

  ig.summary = as.data.frame(py$py_ig_summary)

  if (save.csv) {
    out.dir = dirname(output.file)
    if (!dir.exists(out.dir)) {
      dir.create(out.dir, recursive = TRUE)
    }
    data.table::fwrite(ig.summary, output.file)
    message("Saved Integrated Gradients attributions to: ", output.file)
  }

  return(ig.summary)
}
