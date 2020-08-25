# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "monocle3", "purrr", "mclust", "flexclust")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
NCORES <- 10
BiocParallel::register(MulticoreParam(NCORES))

# Load data ----
set.seed(118617)
nCells <- 5000
sce <- create_hard_balanced_data(nCells = nCells, nClus = 10, seed = 77865)

# Run clustering
run_clusterings(sce, id = 3)