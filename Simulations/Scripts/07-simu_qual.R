# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "BiocParallel",
          "ggplot2", "purrr", "mclust", "Dune")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
NCORES <- 32
BiocParallel::register(MulticoreParam(NCORES))

# Create data ----
source(here("Simulations", "Scripts", "04-size-helper.R"))
set.seed(118616)
nCells <- 5000
print("Creating datasets")
sces <- lapply(1:2, function(i) {
  print(paste0(".. ", i))
  return(create_and_clean(nCells))
})
names(sces) <- paste0("Dataset_", 1:2)

print("Running clustering methods")
clusterings <- purrr::map(sces, run_clusterings)

# 3 methods, different sizes and parameters ----
# Do the consensus
print("Running and evaluating Dunes")
n_randoms <- 1:3
size_random <- 1:5/5
set.seed(2)
params <- expand.grid(n_randoms, size_random)
rownames(params) <- paste0(params[, 1], "_", params[,2])
ARIs <- purrr::map(1:2, function(i) {
  clustering <- clusterings[[i]]
  source <- cbind(clustering$sc3[, c("35", "45")],
                  clustering$UMAP_KMEANS[, c("35", "45")],
                  clustering$tSNE_KMEANS[, c("35", "45")])
  df <- data.frame(cells = clustering$sc3$cells,
                   "sc3" = clustering$sc3[, "40"],
                   "UMAP_KMEANS" = clustering$UMAP_KMEANS[, "40"],
                   "TSNE_KMEANS" = clustering$tSNE_KMEANS[, "40"])
  Dunes <- purrr::map(seq_len(nrow(params)), function(k) {
    n <- params[k, 1]
    r <- params[k, 2]
    random <- source[, sample(1:6, n, replace = F)] %>% as.matrix()
    colnames(random) <- paste0("random_", 1:n)
    random_rows <- sample(nrow(random), nrow(random) * r)
    random[random_rows, ] <- random[sample(random_rows), ]
    df2 <- bind_cols(df, as.data.frame(random))
    return(run_Dune(df2))
  })
  ARIs <- purrr::map(Dunes, evaluate_clustering_methods, sce = sces[[i]])
  names(ARIs) <- rownames(params)
  ARIs[["ref"]] <- evaluate_clustering_methods(sce = sces[[i]], merger = run_Dune(df))
  ARIs <- bind_rows(ARIs, .id = "param")
  return(ARIs)
})

# Do the measures
names(ARIs) <- paste0("Dataset_", 1:2)
ARIs <- bind_rows(ARIs, .id = "sim")
write.csv(x = ARIs, file = here("Simulations", "Data", "Random.csv"))