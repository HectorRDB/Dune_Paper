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

# 3 methods, different sizes ----
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

# 3 methods, changing the parameters ----
# Do the consensus
print("Running Dune")
ks <- as.character(seq(25, 50, 5))
names(ks) <- ks
Dunes <- purrr::map(seqks, function(k) {
  df <- data.frame("sc3" = clusterings[[1]]$sc3[, k],
                   "UMAP_KMEANS" = clusterings[[1]]$UMAP_KMEANS[, k],
                   "TSNE_KMEANS" = clusterings[[1]]$TSNE_KMEANS[, k]
  )
  return(run_Dune(df))
})

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[1]])
names(ARIs) <- ks
ARIs <- bind_rows(ARIs, .id = "param")
write.csv(x = ARIs, file = here("Simulations", "Data", "Param.csv"))

# 3 methods, changing the parameters ----
# Do the consensus
print("Running Dune")
ks <- as.character(seq(30, 50, 5))
names(ks) <- ks
Dunes <- purrr::map(seqks, function(k) {
  df <- data.frame("sc3" = clusterings[[1]]$sc3[, k],
                   "UMAP_KMEANS" = clusterings[[1]]$UMAP_KMEANS[, k],
                   "TSNE_KMEANS" = clusterings[[1]]$TSNE_KMEANS[, k]
  )
  return(run_Dune(df))
})

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[1]])
names(ARIs) <- ks
ARIs <- bind_rows(ARIs, .id = "param")
write.csv(x = ARIs, file = here("Simulations", "Data", "Param.csv"))

# Changing the number of methods ----
# Do the consensus
print("Running Dune")
df <- bind_cols(clusterings[[1]]$sc3[, c("35", "45")],
                clusterings[[1]]$UMAP_KMEANS[, c("35", "45")],
                clusterings[[1]]$TSNE_KMEANS[, c("35", "45")])
colnames(df) <- paste0(rep(c("sc3_", "UMAP_KMEANS_", "TSNE_KMEANS_"), each = 3),
                       c("35", "45"))
Dunes <- list()
clusMat <- data.frame(cells = clusterings[[1]]$sc3$cells,
                      "sc3_40" = clusterings[[1]]$sc3[, "40"],
                      "UMAP_KMEANS_40" = clusterings[[1]]$UMAP_KMEANS[, "40"],
                      "TSNE_KMEANS_40" = clusterings[[1]]$TSNE_KMEANS[, "40"]
                      )
for (i in 1:7) {
  Dunes[[as.character(i + 2)]] <- run_Dune(clusMat)
  k <- sample(ncol(df), 1)
  clusMat[, colnames(df)[k]] <- df[, k]
  df <- df[, -k]
}

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[1]])
names(ARIs) <- 3:9
ARIs <- bind_rows(ARIs, .id = "ns")
write.csv(x = ARIs, file = here("Simulations", "Data", "Ns.csv"))