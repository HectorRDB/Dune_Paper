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
print("Creating datasets")
sces <- lapply(nCells, function(nCell) {
  print(paste0(".. ", nCell))
  sce <- create_simple_balanced_data(nCells = nCell, nClus = 30, 
                                     DE =.1, seed = runif(1, 0, 100))
})

print("Running clustering methods")
clusterings <- purrr::map(sces, run_clusterings)

# Do the consensus
print("Running Dune")
Dunes <- purrr::map(clusterings, function(clustering) {
  df <- data.frame("sc3" = clustering$sc3[, "40"],
                   "UMAP_KMEANS" = clustering$UMAP_KMEANS[, "40"],
                   "TSNE_KMEANS" = clustering$TSNE_KMEANS[, "40"]
                  )
  return(run_Dune(df))
})

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map2(sces, Dunes, evaluate_clustering_methods)
names(ARIs) <- nCells
ARIs <- bind_rows(ARIs, .id = "nCells")
write.csv(x = ARIs, file = here("Simulations", "Data", "Sizes.csv"))