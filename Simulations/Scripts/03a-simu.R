# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "monocle3", "purrr", "mclust", "flexclust")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
NCORES <- 2
BiocParallel::register(MulticoreParam(NCORES))

# Load data ----
source(here("Simulations", "Scripts", "01-create_data.R"))
source(here("Simulations", "Scripts", "02-clusterings-helper.R"))
set.seed(118617)
nCells <- 5000
sce <- create_simple_balanced_data(nCells = nCells, nClus = 10, seed = 197)

# Run clustering
run_clusterings(sce, id = 1)