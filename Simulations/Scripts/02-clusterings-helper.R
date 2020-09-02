# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "purrr", "mclust", "flexclust", "monocle3")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)

run_clusterings <- function(sce, id) {
  # Create data ----
  # We follow the workflow from Duo et al 2018
  # https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
  
  ## Pre-processing ----
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
  sce_og <- sce
  
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
  Seurats <- clusterMatrix
  
  # Running SC3 ----
  print("... Running SC3")
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  sce <- sc3_estimate_k(sce)
  K <- metadata(sce)$sc3$k_estimation
  
  print(paste0("...... The optimal number of clusters defined by sc3 is ", K))
  ks <- min(K, 10):(K + 10)
  names(ks) <- ks - K
  
  sc3 <- map_df(ks, function(k){
    SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE, 
               gene_filter = F, n_cores = NCORES, rand_seed = 786907)
    SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
    return(SC3)
  })
  
  sc3$cells <- colnames(sce)
  SC3s <- sc3
  
  # Running ZinbWave ----
  sce <- sce_og
  vars <- matrixStats::rowVars(logcounts(sce))
  ind <- vars > sort(vars,decreasing = TRUE)[1000]
  whichGenes <- rownames(sce)[ind]
  sceVar <- sce[ind,]
  
  cat("Running with K =  on the filtered data\n")
  cat("Time to run zinbwave (seconds):\n")
  sceVar$Batch <- as.factor(sceVar$Batch)
  if (nlevels(sceVar$Batch) > 1) {
    print(system.time(zinbW <- zinbwave(Y = sceVar, X = "~Batch", K = )))  
  } else {
    print(system.time(zinbW <- zinbwave(Y = sceVar, K = )))
  }
  
  
  reducedDim(sce, type = "zinb-K-") <- reducedDim(zinbW)
  TNSE <- Rtsne(reducedDim(zinbW), initial_dims = )
  df <- data.frame(x = TNSE$Y[, 1], y = TNSE$Y[, 2], col = clusters)
  p <- ggplot(df, aes(x = x, y = y, col = col)) +
    geom_point(size = .4, alpha = .3) +
    theme_classic() +
    labs(x = "dim1", y = "dim2")
  ggsave(here("Simulations", "Figures", paste0(id, "_K_.png")), p)
  
  # Running Monocle ----
  pd <- as.data.frame(sce@colData)
  fd <- data.frame(gene_short_name = rownames(assays(sce)$counts))
  zinbW <- reducedDim(sce, type = "zinb-K-")
  rownames(fd) <- rownames(assays(sce)$counts)
  sce <- new_cell_data_set(assays(sce)$counts,
                           cell_metadata = pd,
                           gene_metadata = fd)
  reducedDim(sce, "Aligned") <- zinbW
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
  Monocles <- clusterMatrix
  
  # Save results ----
  write.csv(Seurats, here("Simulations", "Data", paste0("Seurat", id, ".csv")))
  write.csv(Seurats, here("Simulations", "Data", paste0("SC3", id, ".csv")))
  write.csv(Seurats, here("Simulations", "Data", paste0("Monocle", id, ".csv")))
  
  # RSEC ----
  sequential <- FALSE
  subsample <- T
  clusterFunction <- "pam"
  sce <- sce_og
  reduceMeth <- reducedDimNames(sce)
  
  print(system.time(
    sce <- RSEC(sce, k0s = seq(10, 50, by = 5), alphas = c(0.1, 0.3),
                reduceMethod = reduceMeth, sequential = sequential,
                subsample = subsample, minSizes = 1, betas = c(0.8), 
                clusterFunction = clusterFunction, ncores = NCORES, run = TRUE,
                isCount = FALSE, dendroReduce = reduceMeth[length(reduceMeth)],
                dendroNDims = 50, consensusProportion = 0.7, verbose = TRUE,
                random.seed = 23578, mergeMethod = "adjP", mergeCutoff = 0.05,
                subsampleArgs = list(resamp.num = 50, clusterFunction = "kmeans"),
                mergeLogFCcutoff = 1, consensusMinSize = 10)
  ))
  
  # Saving objects
  saveRDS(sce, here("Simulations", "Data", paste0("Merger_", id, ".rds")))
  return()
}