# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "monocle3", "purrr", "mclust", "flexclust")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
NCORES <- 8
BiocParallel::register(MulticoreParam(NCORES))

# Load data ----
source(here("Simulations", "Scripts", "01-create_data.R"))
source(here("Simulations", "Scripts", "02-clusterings-helper.R"))
set.seed(19038)
nCells <- 5000
sce <- create_simple_balanced_data(nCells = nCells, nClus = 30, 
                                   DE =.05, seed = runif(1, 0, 100))

# Run clustering
run_clusterings(sce, id = 3)

# Do the consensus
rsec <- readRDS(here("Simulations", "Data", paste0("Merger_", 3, ".rds")))
run_merging_methods(rsec, sce, 3)

# Do the measures
evaluate_clustering_methods(sce, 3)