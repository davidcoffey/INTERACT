#' Format INTERACT Raw Interactions into Multi-Channel Edge Tensors and Node Features
#'
#' @description Formats INTERACT interaction tables and tumor genomic features into
#' multi-channel edge tensors and node features suitable for input into the Graph Autoencoder.
#'
#' @param interact.data INTERACT results object, .Rda file path, or interactions data.frame/data.table.
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{edge.tensors:} A data.table containing edge features (RNA, pathway, program, SNV, CNV) across cell pairs.
#'   \item \strong{node.features:} A data.table containing patient-level tumor genomic node features.
#' }
#' @importFrom data.table as.data.table is.data.table := copy
#' @importFrom stats var
#' @rawNamespace export(format.interact.tensors)
format.interact.tensors = function(interact.data) {
  # If input is a character path, load the .Rda file
  if (is.character(interact.data)) {
    message("Loading INTERACT results from: ", interact.data)
    loaded.env = new.env()
    load(interact.data, envir = loaded.env)
    var.names = ls(loaded.env)
    # Find results object
    if ("results" %in% var.names) {
      target.obj = loaded.env$results
    } else if ("results.bms" %in% var.names) {
      target.obj = loaded.env$results.bms
    } else {
      target.obj = loaded.env[[var.names[1]]]
    }
  } else {
    target.obj = interact.data
  }
  
  # Extract interaction table
  if (is.list(target.obj) && "interactions" %in% names(target.obj)) {
    dt = data.table::copy(data.table::as.data.table(target.obj$interactions))
  } else if (data.table::is.data.table(target.obj) || is.data.frame(target.obj)) {
    dt = data.table::copy(data.table::as.data.table(target.obj))
  } else {
    stop("Input must be an INTERACT results object with an $interactions element, an .Rda file path, or a data.frame.")
  }
  
  message("Aggregating multi-channel edge tensors across cell pairs...")
  dt[, linear_base := 2^base.score]
  dt[, s.cnv.dev := s.cnv - 2.0]
  dt[, r.cnv.dev := r.cnv - 2.0]
  
  sample.agg = dt[, .(
    sum_linear = sum(linear_base, na.rm = TRUE) + 1e-8,
    edge_pathway = sum(linear_base * weight.pathways, na.rm = TRUE),
    edge_program = sum(linear_base * weight.programs, na.rm = TRUE),
    edge_snv_mean = sum(linear_base * (s.snv + r.snv), na.rm = TRUE),
    edge_cnv_mean = sum(linear_base * (s.cnv.dev + r.cnv.dev), na.rm = TRUE)
  ), by = .(sample, sender, receiver)]
  
  sample.agg = sample.agg[sender != receiver]
  
  sample.agg[, `:=`(
    edge_rna = log2(sum_linear),
    edge_pathway = edge_pathway / sum_linear,
    edge_program = edge_program / sum_linear,
    edge_snv = edge_snv_mean / sum_linear,
    edge_cnv = edge_cnv_mean / sum_linear
  )]
  
  edge.tensors = sample.agg[, .(sample, sender, receiver, edge_rna, edge_pathway, edge_program, edge_snv, edge_cnv)]
  
  message("Extracting patient-level tumor genomic node features...")
  tumor.edges = dt[sender == "Tumor" | receiver == "Tumor"]
  if (nrow(tumor.edges) > 0) {
    node.features = tumor.edges[, .(
      Tumor_Mutational_Burden = mean(abs(s.snv) + abs(r.snv), na.rm = TRUE),
      Tumor_Aneuploidy_Var = stats::var(c(s.cnv, r.cnv), na.rm = TRUE)
    ), by = .(sample)]
  } else {
    node.features = data.table::data.table(
      sample = unique(dt$sample),
      Tumor_Mutational_Burden = 0,
      Tumor_Aneuploidy_Var = 0
    )
  }
  
  node.features[is.na(Tumor_Aneuploidy_Var), Tumor_Aneuploidy_Var := 0]
  node.features[is.na(Tumor_Mutational_Burden), Tumor_Mutational_Burden := 0]

  all.samples = unique(edge.tensors$sample)
  missing.samples = setdiff(all.samples, node.features$sample)
  if (length(missing.samples) > 0) {
    extra.rows = data.table::data.table(
      sample = missing.samples,
      Tumor_Mutational_Burden = 0,
      Tumor_Aneuploidy_Var = 0
    )
    node.features = rbind(node.features, extra.rows)
  }
  
  return(list(edge.tensors = edge.tensors, node.features = node.features))
}

# 2 Environment Installation Helper ----

#' Install Python Dependencies for INTERACT Graph Autoencoder
#'
#' @description Sets up a Conda or Virtualenv environment with the necessary Python
#' dependencies (PyTorch, PyTorch Geometric, pandas, numpy, scikit-learn) required
#' to run the INTERACT multi-channel Graph Autoencoder.
#'
#' @param conda.env String. Name of Conda or Virtualenv environment to create or install into (default: "pytorch").
#' @param method String. Package management method: "auto", "conda", or "virtualenv" (default: "auto").
#' @param python.version String. Python version to use when creating a new environment (default: "3.10").
#' @param pip Logical. If TRUE, packages are installed via pip (default: TRUE).
#' @param ... Additional arguments passed to \code{reticulate::conda_install()} or \code{reticulate::virtualenv_install()}.
#'
#' @return Invisible TRUE upon successful installation and verification.
#'
#' @importFrom reticulate conda_binary conda_create conda_install conda_list py_run_string virtualenv_create virtualenv_install virtualenv_list use_condaenv use_virtualenv
#' @export
InstallInteractEnv = function(
  conda.env = "pytorch",
  method = c("auto", "conda", "virtualenv"),
  python.version = "3.10",
  pip = TRUE,
  ...
) {
  method = match.arg(method)
  req.pkgs = c("torch", "torch_geometric", "pandas", "numpy", "scikit-learn")

  message("--- INTERACT Python Environment Setup ---")
  message("Target environment: ", conda.env)
  message("Required packages: ", paste(req.pkgs, collapse = ", "))

  if (method == "auto") {
    conda.bin = tryCatch(reticulate::conda_binary(), error = function(e) NULL)
    if (!is.null(conda.bin) && nzchar(conda.bin)) {
      method = "conda"
    } else {
      method = "virtualenv"
    }
  }

  if (method == "conda") {
    message("Using Conda package manager...")
    existing.envs = tryCatch(reticulate::conda_list()$name, error = function(e) character())
    if (!conda.env %in% existing.envs) {
      message("Creating Conda environment '", conda.env, "' with Python ", python.version, "...")
      reticulate::conda_create(envname = conda.env, python_version = python.version)
    } else {
      message("Found existing Conda environment '", conda.env, "'.")
    }

    message("Installing required packages into Conda environment '", conda.env, "'...")
    reticulate::conda_install(
      envname = conda.env,
      packages = req.pkgs,
      pip = pip,
      ...
    )
    reticulate::use_condaenv(conda.env, required = TRUE)
  } else {
    message("Using Virtualenv package manager...")
    existing.venv = tryCatch(reticulate::virtualenv_list(), error = function(e) character())
    if (!conda.env %in% existing.venv) {
      message("Creating Virtualenv '", conda.env, "' with Python ", python.version, "...")
      reticulate::virtualenv_create(envname = conda.env, python_version = python.version)
    } else {
      message("Found existing Virtualenv '", conda.env, "'.")
    }

    message("Installing required packages into Virtualenv '", conda.env, "'...")
    reticulate::virtualenv_install(
      envname = conda.env,
      packages = req.pkgs,
      ...
    )
    reticulate::use_virtualenv(conda.env, required = TRUE)
  }

  message("\nVerifying Python dependencies...")
  verification = tryCatch(
    {
      reticulate::py_run_string("
import torch
import torch_geometric
import pandas
import numpy
import sklearn
print('INTERACT Python dependencies successfully verified!')
")
      TRUE
    },
    error = function(e) {
      warning("Verification failed: ", e$message)
      FALSE
    }
  )

  if (verification) {
    message("Setup complete! You can now run RunInteractAutoencoder(..., conda.env = '", conda.env, "').")
  }

  invisible(verification)
}

# 3 Main Autoencoder Function ----

#' Run INTERACT Graph Autoencoder
#'
#' Trains a multi-channel Graph Autoencoder (GINEConv) on patient tumor microenvironment
#' communication networks and extracts low-dimensional latent embeddings.
#'
#' @param interact.results INTERACT results object, .Rda file path, or pre-computed edge tensor data.table for the reference cohort.
#' @param interact.node.features Optional node features data.table (if pre-computed tensors are passed).
#' @param interact.results.query Optional query cohort INTERACT results, .Rda path, or edge tensor data.table to project into the reference embedding space.
#' @param interact.node.features.query Optional query cohort node features data.table.
#' @param conda.env Conda or virtualenv environment containing PyTorch and PyG (default: "pytorch").
#' @param auto.install Logical; if TRUE, automatically attempts to create and install the environment if dependencies are missing.
#' @param epochs Number of training epochs (default: 150).
#' @param batch.size Mini-batch size for training (default: 32).
#' @param learning.rate Adam optimizer learning rate (default: 0.005).
#' @param weight.decay Weight decay L2 regularization (default: 1e-4).
#' @param hidden.dim Hidden dimension for GINEConv MLPs (default: 32).
#' @param emb.dim Latent embedding dimension per cell type (default: 4).
#' @param channel.weights Channel loss weights for `[RNA, Pathway, Program, SNV, CNV]` (default: c(1.0, 3.0, 3.0, 5.0, 5.0)).
#' @param seed Random seed for reproducibility (default: 42).
#' @param deterministic Logical; if TRUE, enforces single-threaded deterministic operations for 100% reproducibility.
#' @param save.csv Logical; if TRUE, saves embeddings as CSV files in output.dir.
#' @param output.dir Directory to save CSV embeddings (default: "Data").
#'
#' @return A data.frame of patient embeddings (or a list with `$reference` and `$query` data.frames when `interact.results.query` is provided).
#'
#' @importFrom data.table as.data.table data.table fwrite
#' @importFrom reticulate use_condaenv use_virtualenv use_python conda_list virtualenv_list py_module_available py_run_string py
#' @export
RunInteractAutoencoder = function(
  interact.results = NULL,
  interact.node.features = NULL,
  interact.results.query = NULL,
  interact.node.features.query = NULL,
  conda.env = "pytorch",
  auto.install = FALSE,
  epochs = 150L,
  batch.size = 32L,
  learning.rate = 0.005,
  weight.decay = 1e-4,
  hidden.dim = 32L,
  emb.dim = 4L,
  channel.weights = c(1.0, 3.0, 3.0, 5.0, 5.0),
  seed = 42L,
  deterministic = TRUE,
  save.csv = FALSE,
  output.dir = "Data"
) {
  if (is.null(interact.results)) {
    stop("interact.results must be provided.")
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
  } else if (isTRUE(auto.install)) {
    message("Environment '", conda.env, "' not found. Automatically creating and installing required packages...")
    InstallInteractEnv(conda.env = conda.env)
  } else {
    tryCatch(reticulate::use_condaenv(conda.env, required = FALSE), error = function(e) NULL)
  }

  # Verify required Python packages
  req.pkgs = c("torch", "torch_geometric", "pandas", "numpy", "sklearn")
  missing.pkgs = req.pkgs[!vapply(req.pkgs, reticulate::py_module_available, logical(1))]

  if (length(missing.pkgs) > 0) {
    if (isTRUE(auto.install)) {
      message("Missing Python package(s): ", paste(missing.pkgs, collapse = ", "), ". Installing automatically...")
      InstallInteractEnv(conda.env = conda.env)
    } else {
      stop(
        "\nMissing required Python package(s): ", paste(missing.pkgs, collapse = ", "), ".\n",
        "To install them automatically from R, run:\n",
        "  INTERACT::InstallInteractEnv(conda.env = '", conda.env, "')\n",
        "Or pass `auto.install = TRUE`:\n",
        "  RunInteractAutoencoder(..., conda.env = '", conda.env, "', auto.install = TRUE)\n",
        "Alternatively, configure a Conda environment in terminal:\n",
        "  conda create -n ", conda.env, " python=3.10 -y\n",
        "  conda activate ", conda.env, "\n",
        "  pip install torch torch_geometric pandas numpy scikit-learn\n"
      )
    }
  }
  
  # Check if discovery input needs tensor formatting
  req.tensor.cols = c("edge_rna", "edge_pathway", "edge_program", "edge_snv", "edge_cnv")
  if (is.character(interact.results) ||
      (is.list(interact.results) && "interactions" %in% names(interact.results)) ||
      (is.data.frame(interact.results) && !all(req.tensor.cols %in% names(interact.results)))) {
    disc.prepped = format.interact.tensors(interact.results)
    edge.disc = disc.prepped$edge.tensors
    node.disc = disc.prepped$node.features
  } else {
    edge.disc = data.table::as.data.table(interact.results)
    if (is.null(interact.node.features)) {
      stop("When passing pre-computed edge tensors, interact.node.features must also be provided.")
    }
    node.disc = data.table::as.data.table(interact.node.features)
  }

  # Check if query input is provided
  has.val = !is.null(interact.results.query)
  if (has.val) {
    if (is.character(interact.results.query) ||
        (is.list(interact.results.query) && "interactions" %in% names(interact.results.query)) ||
        (is.data.frame(interact.results.query) && !all(req.tensor.cols %in% names(interact.results.query)))) {
      val.prepped = format.interact.tensors(interact.results.query)
      edge.val = val.prepped$edge.tensors
      node.val = val.prepped$node.features
    } else {
      edge.val = data.table::as.data.table(interact.results.query)
      if (is.null(interact.node.features.query)) {
        stop("When passing pre-computed query edge tensors, interact.node.features.query must also be provided.")
      }
      node.val = data.table::as.data.table(interact.node.features.query)
    }
  } else {
    edge.val = data.table::data.table()
    node.val = data.table::data.table()
  }
  
  # Pass data into Python
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
  
  # Transfer tables to Python
  py = reticulate::py
  py$r_edge_disc = as.data.frame(edge.disc)
  py$r_node_disc = as.data.frame(node.disc)
  py$r_has_val = has.val
  py$r_edge_val = if (has.val) as.data.frame(edge.val) else data.frame()
  py$r_node_val = if (has.val) as.data.frame(node.val) else data.frame()
  py$r_epochs = as.integer(epochs)
  py$r_batch_size = as.integer(batch.size)
  py$r_lr = as.numeric(learning.rate)
  py$r_wd = as.numeric(weight.decay)
  py$r_hidden_dim = as.integer(hidden.dim)
  py$r_emb_dim = as.integer(emb.dim)
  py$r_weights = as.numeric(channel.weights)
  py$r_seed = as.integer(seed)
  py$r_deterministic = as.logical(deterministic)
  
  message("Configuring PyTorch Geometric model and training...")
  reticulate::py_run_string("
if r_deterministic:
    torch.set_num_threads(1)

torch.manual_seed(r_seed)
np.random.seed(r_seed)
g = torch.Generator()
g.manual_seed(r_seed)

class TME_Autoencoder(nn.Module):
    def __init__(self, num_node_features, hidden_dim, emb_dim, num_edge_features, num_cells):
        super(TME_Autoencoder, self).__init__()
        self.num_cells = num_cells
        self.emb_dim = emb_dim
        
        mlp1 = nn.Sequential(
            nn.Linear(num_node_features, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )
        self.conv1 = GINEConv(nn=mlp1, edge_dim=num_edge_features)
        
        mlp2 = nn.Sequential(
            nn.Linear(hidden_dim, emb_dim),
            nn.ReLU(),
            nn.Linear(emb_dim, emb_dim)
        )
        self.conv2 = GINEConv(nn=mlp2, edge_dim=num_edge_features)
        
        self.edge_decoder = nn.Sequential(
            nn.Linear(emb_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, num_edge_features)
        )

    def encode(self, x, edge_index, edge_attr):
        h = F.elu(self.conv1(x, edge_index, edge_attr=edge_attr))
        h = self.conv2(h, edge_index, edge_attr=edge_attr)
        return h

    def decode(self, z, edge_index):
        edge_features = torch.cat([z[edge_index[0]], z[edge_index[1]]], dim=-1)
        return self.edge_decoder(edge_features)

    def forward(self, x, edge_index, edge_attr, batch):
        z = self.encode(x, edge_index, edge_attr)
        predicted_edges = self.decode(z, edge_index)
        patient_emb = z.view(-1, self.num_cells * self.emb_dim)
        return predicted_edges, patient_emb

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
        
        data = Data(x=x_tensor, edge_index=edge_index, edge_attr=attr_tensor, sample_id=pt)
        dataset.append(data)
    return dataset

edge_df_disc = r_edge_disc.copy()
nodes_df_disc = r_node_disc.copy()

feature_cols = ['edge_rna', 'edge_pathway', 'edge_program', 'edge_snv', 'edge_cnv']
scaler = StandardScaler()
edge_df_disc[feature_cols] = scaler.fit_transform(edge_df_disc[feature_cols])

if r_has_val:
    edge_df_val = r_edge_val.copy()
    nodes_df_val = r_node_val.copy()
    edge_df_val[feature_cols] = scaler.transform(edge_df_val[feature_cols])
    unique_cells = sorted(list(set(edge_df_disc['sender']).union(set(edge_df_disc['receiver']))
                               .union(set(edge_df_val['sender'])).union(set(edge_df_val['receiver']))))
else:
    unique_cells = sorted(list(set(edge_df_disc['sender']).union(set(edge_df_disc['receiver']))))

num_cells = len(unique_cells)
cell_to_idx = {cell: i for i, cell in enumerate(unique_cells)}
tumor_idx = cell_to_idx.get('Tumor')

nodes_dict_disc = nodes_df_disc.set_index('sample').to_dict('index')
disc_dataset = prep_dataset(edge_df_disc, nodes_dict_disc, feature_cols, cell_to_idx, tumor_idx, num_cells)

if r_deterministic:
    loader = DataLoader(disc_dataset, batch_size=r_batch_size, shuffle=True, generator=g)
else:
    loader = DataLoader(disc_dataset, batch_size=r_batch_size, shuffle=True)

infer_loader_disc = DataLoader(disc_dataset, batch_size=1, shuffle=False)

if r_has_val:
    nodes_dict_val = nodes_df_val.set_index('sample').to_dict('index')
    val_dataset = prep_dataset(edge_df_val, nodes_dict_val, feature_cols, cell_to_idx, tumor_idx, num_cells)
    infer_loader_val = DataLoader(val_dataset, batch_size=1, shuffle=False)

node_dim = num_cells + 2
device = torch.device('cpu')

model = TME_Autoencoder(node_dim, r_hidden_dim, r_emb_dim, 5, num_cells).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=r_lr, weight_decay=r_wd)

model.train()
weights_t = torch.tensor(r_weights, dtype=torch.float32, device=device)
for epoch in range(r_epochs):
    total_loss = 0.0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred_edges, _ = model(batch.x, batch.edge_index, batch.edge_attr, batch.batch)
        channel_mse = torch.mean((pred_edges - batch.edge_attr)**2, dim=0)
        loss = torch.sum(channel_mse * weights_t)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch.num_graphs

model.eval()

def extract_emb(data_loader):
    embeddings = []
    with torch.no_grad():
        for data in data_loader:
            data = data.to(device)
            _, patient_vector = model(data.x, data.edge_index, data.edge_attr, torch.zeros(data.x.size(0), dtype=torch.long).to(device))
            sid = data.sample_id[0]
            if isinstance(sid, torch.Tensor):
                sid = sid.item()
            embeddings.append([sid] + patient_vector.cpu().numpy().flatten().tolist())
    flat_dim = num_cells * r_emb_dim
    col_names = ['sample'] + [f'GNN.Dim.{i+1}' for i in range(flat_dim)]
    return pd.DataFrame(embeddings, columns=col_names)

py_disc_emb = extract_emb(infer_loader_disc)
py_val_emb = extract_emb(infer_loader_val) if r_has_val else pd.DataFrame()
")
  
  ref.embeddings = as.data.frame(py$py_disc_emb)
  
  if (save.csv) {
    if (!dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE)
    }
    ref.out = file.path(output.dir, "4 Autoencoder Reference Embeddings.csv")
    data.table::fwrite(ref.embeddings, ref.out)
    message("Saved Reference embeddings to: ", ref.out)
  }
  
  if (has.val) {
    query.embeddings = as.data.frame(py$py_val_emb)
    if (save.csv) {
      query.out = file.path(output.dir, "4 Autoencoder Query Embeddings.csv")
      data.table::fwrite(query.embeddings, query.out)
      message("Saved Query embeddings to: ", query.out)
    }
    return(list(reference = ref.embeddings, query = query.embeddings))
  }
  
  return(ref.embeddings)
}
