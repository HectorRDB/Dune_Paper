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
source(here("Simulations", "Scripts", "04-size-helper.R"))
set.seed(118616)
nCells <- c(100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900,
            1000, 1500, 2000, 2500, 3000, 5000)
print("Creating datasets")
sces <- lapply(nCells, function(nCell) {
  print(paste0(".. ", nCell))
  return(create_and_clean(nCell))
})

print("Running clustering methods")
clusterings <- purrr::map(sces, run_clusterings)

# 3 methods, different sizes and parameters ----
# Do the consensus
print("Running and evaluating Dunes")
ks <- as.character(seq(30, 50, 5))
names(ks) <- ks
ARIs <- purrr::map(seq_along(nCells), function(i) {
  clustering <- clusterings[[i]]
  Dunes <- purrr::map(ks, function(k) {
    df <- data.frame(cells = clustering$sc3$cells,
                     "sc3" = clustering$sc3[, k],
                     "UMAP_KMEANS" = clustering$UMAP_KMEANS[, k],
                     "TSNE_KMEANS" = clustering$tSNE_KMEANS[, k]
    )
    return(run_Dune(df))
  })
  ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[i]])
  names(ARIs) <- ks
  ARIs <- bind_rows(ARIs, .id = "param")
  return(ARIs)
})

# Do the measures
names(ARIs) <- nCells
ARIs <- bind_rows(ARIs, .id = "nCells")
write.csv(x = ARIs, file = here("Simulations", "Data", "Imp.csv"))