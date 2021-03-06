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
nCells <- c(100, 200, 500, 1000, 2000, 5000)
print("Creating datasets")
sces <- lapply(nCells, function(nCell) {
  print(paste0(".. ", nCell))
  return(create_and_clean(nCell))
})

print("Running clustering methods")
# clusterings <- purrr::map(sces, run_clusterings)
# saveRDS(clusterings, here("Simulations", "Data", "clusterings.rds"))
clusterings  <- readRDS(here("Simulations", "Data", "clusterings.rds"))

# 3 methods, different sizes ----
# Do the consensus
print("Running Dune")
Dunes <- purrr::map(clusterings, function(clustering) {
  df <- data.frame(cells = clustering$sc3$cells,
                   "sc3" = clustering$sc3[, "40"],
                   "UMAP_KMEANS" = clustering$UMAP_KMEANS[, "40"],
                   "TSNE_KMEANS" = clustering$tSNE_KMEANS[, "40"]
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
ks <- as.character(seq(30, 50, 5))
names(ks) <- ks
Dunes <- purrr::map(ks, function(k) {
  df <- data.frame(cells = clusterings[[6]]$sc3$cells,
                   "sc3" = clusterings[[6]]$sc3[, k],
                   "UMAP_KMEANS" = clusterings[[6]]$UMAP_KMEANS[, k],
                   "TSNE_KMEANS" = clusterings[[6]]$tSNE_KMEANS[, k]
  )
  return(run_Dune(df))
})

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[6]])
names(ARIs) <- ks
ARIs <- bind_rows(ARIs, .id = "param")
write.csv(x = ARIs, file = here("Simulations", "Data", "Param.csv"))

# Changing the number of methods ----
# Do the consensus
print("Running Dune")
df <- cbind(clusterings[[6]]$sc3$cells,
            clusterings[[6]]$sc3[, c("35", "40", "45")],
            clusterings[[6]]$UMAP_KMEANS[, c("35", "40", "45")],
            clusterings[[6]]$tSNE_KMEANS[, c("35", "40", "45")])
colnames(df) <- c("cells",
                  paste0(rep(c("sc3_", "UMAP_KMEANS_", "tSNE_KMEANS_"), each = 3),
                         c("35", "40", "45")))
df <- as.data.frame(df)
Dunes <- list()
for (i in 2:9) {
  groups <- combn(9, i)
  for (j in sample(seq_len(ncol(groups)), min(5, ncol(groups)))) {
    clusMat <- df[, c(1, groups[, j] + 1)]
    print(ncol(clusMat))
    Dunes[[paste0(i, "_", j)]] <- run_Dune(clusMat)
  }
}

# Do the measures
print("Evaluating Dune")
ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[6]])
names(ARIs) <- names(Dunes)
ARIs <- bind_rows(ARIs, .id = "ns")
write.csv(x = ARIs, file = here("Simulations", "Data", "Ns.csv"))