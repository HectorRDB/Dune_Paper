# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "BiocParallel",
          "ggplot2", "purrr", "purrr", "mclust")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
NCORES <- 32
BiocParallel::register(MulticoreParam(NCORES))

# Create data ----
source(here("Simulations", "Scripts", "01-create_data.R"))
source(here("Simulations", "Scripts", "04-size-helper.R"))
set.seed(118617)
nCells <- c(100, 200, 500, 1000, 2000, 5000, 10000)
sces <- lapply(nCells, function(nCell) {
  sce <- create_simple_balanced_data(nCells = nCell, nClus = 30, 
                                     DE =.1, seed = runif(1, 0, 100))
})

# With 3 inputs
clusterings <- purrr::map(sces, run_clusterings)

# Do the consensus
Dunes <- purrr::map(clusterings, run_Dune)

# Do the measures
ARIs <- purrr::map2(sces, Dunes, evaluate_clustering_methods)
names(ARIs) <- nCells
ARIs <- bind_rows(ARIs, .id = "nCells")
write.csv(x = ARIs, file = here("Simulations", "Data", "Sizes.csv"))