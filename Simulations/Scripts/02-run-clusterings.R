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
source(here("Simulations", "Scripts", "01-create_data.R"))
set.seed(101)
nCells <- 5000
sce1 <- create_simple_balanced_data(nCells = nCells, nClus = 10, seed = 197)
sce2 <- create_simple_balanced_data(nCells = nCells, nClus = 20, seed = 197)
sce3 <- create_hard_balanced_data(nCells = nCells, nClus = 10, seed = 77865)
sce4 <- create_unbalanced_data(nCells = nCells, nClus = 10, nBatches = 1,
                               DE = .2, seed = 45678)

# Run clustering
Seurats <- list()
SC3s <- list()
Monocles <- list()
for (dataset in c(paste0("sce", 1:4))) {
  print(dataset)
  # Create data ----
  # We follow the workflow from Duo et al 2018
  # https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
  sce <- get(dataset)
  clusters <- sce$Group
  keep_features <- rowSums(counts(sce) > 0) > 0
  sce <- sce[keep_features, ]
  df <- perCellQCMetrics(sce)
  df$libsize.drop <- isOutlier(df$total, nmads = 3,
                                         type = "lower", log = TRUE)
  df$feature.drop <- isOutlier(df$detected, nmads = 3,
                                         type = "lower", log = TRUE)
  df <- as.data.frame(df)
  sce <- sce[, !(df$libsize.drop | df$feature.drop)]
  clusters <- clusters[!(df$libsize.drop | df$feature.drop)]
  sce <- computeSumFactors(sce, sizes = pmin(ncol(sce), seq(20, 120, 20)),
                           min.mean = 0.1)
  logcounts(sce) <- scater::normalizeCounts(sce)
  # Running Seurat ----
  print("... Running Seurat")
  sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K')
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
  sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                  dispersion.function = LogVMR, do.plot = F)
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = c("nCount_RNA", "human"))
  sce <- as.SingleCellExperiment(sSeurat)
  
  sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)
  clusterMatrix <- NULL
  for (RESOLUTION in seq(from = 0.3, to = 2.5, by = .1)) {
    print(paste0("...... ", RESOLUTION))
    for (K.PARAM in c(30, 50, 100)) {
      print(paste0("......... ", K.PARAM))
      sSeurat_star <- FindNeighbors(sSeurat, dims = 1:K.PARAM)
      sSeurat_star <- FindClusters(sSeurat_star, resolution = RESOLUTION)
      clusterMatrix <- cbind(clusterMatrix, Idents(sSeurat_star))
      colnames(clusterMatrix)[ncol(clusterMatrix)] <- paste(RESOLUTION, K.PARAM,
                                                            sep = ",")
    }
  }
  
  clusterMatrix <- as.data.frame(clusterMatrix)
  clusterMatrix$cells <- colnames(sce)
  Seurats[[dataset]] <- clusterMatrix
  
  # Running SC3 ----
  print("... Running SC3")
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  sce <- sc3_estimate_k(sce)
  K <- metadata(sce)$sc3$k_estimation
  
  print(paste0("...... The optimal number of clusters defined by sc3 is ", K))
  ks <- 10:(K + 10)
  names(ks) <- ks - K
  
  sc3 <- map_df(ks, function(k){
    SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE, 
               gene_filter = F, n_cores = NCORES, rand_seed = 786907)
    SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
    return(SC3)
  })
  
  sc3$cells <- colnames(sce)
  SC3s[[dataset]] <- sc3
  
  # Running ZinbWave ----
  vars <- matrixStats::rowVars(logcounts(sce))
  print("... Running Zinbwave")
  cat("...... Running with K = 0 on the full data\n")
  cat("...... Time to run zinbwave (seconds):\n")
  
  ind <- vars > sort(vars,decreasing = TRUE)[1000]
  whichGenes <- rownames(sce)[ind]
  sceVar <- sce[ind,]
  
  zinbWs <- lapply( 1:5 * 10, function(zinbDim) {
    cat("Running with K = ", zinbDim, " on the filtered data\n")
    cat("Time to run zinbwave (seconds):\n")
    print(system.time(zinb <- zinbwave(sceVar, K = zinbDim)))
    return(zinb)
  })
  
  
  for (i in 1:5) {
    type <- paste0("zinb-K-", i * 10)
    reducedDim(sce, type = type) <- zinbW <- reducedDim(zinbWs[[i]])
    TNSE <- Rtsne(zinbW, initial_dims = i  * 10)
    df <- data.frame(x = TNSE$Y[, 1], y = TNSE$Y[, 2], col = clusters)
    p <- ggplot(df, aes(x = x, y = y, col = col)) +
      geom_point(size = .4, alpha = .3) +
      theme_classic() +
      labs(x = "dim1", y = "dim2")
    ggsave(here("Simulations", "Figures",
                paste0(dataset, "_K_", i * 10, ".png"), p))
  }
  
  # Running Monocle ----
  pd <- as.data.frame(sce@colData)
  fd <- data.frame(gene_short_name = rownames(assays(sce)$counts))
  zinbW <- reducedDim(sce, type = reducedDimNames(sce)[3])
  rownames(fd) <- rownames(assays(sce)$counts)
  sce <- new_cell_data_set(assays(sce)$counts,
                           cell_metadata = pd,
                           gene_metadata = fd)
  sce@reducedDims <- SimpleList("PCA" = zinbW)
  sce <- reduce_dimension(sce)
  ks <- c(2:9, seq(from = 10, to = 100, by = 5))
  names(ks) <- paste0("k_", ks)
  clusterMatrix <- map_df(ks, function(k){
    sce2 <- cluster_cells(sce,
                          k = k,
                          louvain_iter = 2,
                          verbose = F)
    return(sce2@clusters$UMAP$clusters %>% as.numeric())
  })
  
  clusterMatrix$cells <- colnames(sce)
  Monocles[[dataset]] <- clusterMatrix
}

Seurats <- bind_rows(Seurats, .id = dataset)
write.csv(Seurats, here("Simulations", "Data", "Seurats.csv"))
SC3s <- bind_rows(SC3s, .id = dataset)
write.csv(Seurats, here("Simulations", "Data", "SC3s.csv"))
Monocles <- bind_rows(Monocles, .id = dataset)
write.csv(Seurats, here("Simulations", "Data", "Monocles.csv"))
