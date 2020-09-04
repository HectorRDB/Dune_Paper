# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "purrr", "mclust", "flexclust", "monocle3", "scater")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

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
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = c("nCount_RNA", "Batch"))
  sce <- as.SingleCellExperiment(sSeurat)
  sce_Seurat <- sce
  
  sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)
  sSeurat <- RunUMAP(sSeurat, verbose = FALSE, dims = 2)
  p <- UMAPPlot(sSeurat, group.by = clusters)
  ggsave(p, here("Simulations", "Figures", paste0("UMAP_", id, ".png")))
  
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
  write.csv(Seurats, here("Simulations", "Data", paste0("Seurat", id, ".csv")))
  
  # Running SC3 ----
  print("... Running SC3")
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  sce <- sc3_estimate_k(sce)
  K <- metadata(sce)$sc3$k_estimation
  
  print(paste0("...... The optimal number of clusters defined by sc3 is ", K))
  ks <- seq(from = 20, to = 50, by = 5)
  names(ks) <- ks - K
  
  sc3 <- map_df(ks, function(k){
    SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE, 
               gene_filter = F, n_cores = NCORES, rand_seed = 786907)
    SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
    return(SC3)
  })
  
  sc3$cells <- colnames(sce)
  SC3s <- sc3
  write.csv(SC3s, here("Simulations", "Data", paste0("SC3", id, ".csv")))
  
  # RSEC ----
  sequential <- FALSE
  subsample <- T
  clusterFunction <- "pam"
  sce <- sce_og
  sce <- scater::runPCA(sce)
  sce <- scater::runUMAP(sce)
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
  
  # K-Means ----
  sce <- sce_Seurat
  sce <- scater::runUMAP(sce, ncomponents = 3)
  UMAP <- reducedDim(sce, "UMAP")
  ks <- seq(5, 50, 5)
  names(ks) <- ks
  K_MEANS <- map_dfc(ks, function(k){
    return(kmeans(UMAP, centers = k)$cluster)
  })
  K_MEANS$cells <- rownames(UMAP)
  
  write.csv(K_MEANS, here("Simulations", "Data", paste0("KMEANS", id, ".csv")))
  return()
  
}

run_merging_methods <- function(rsec, id, sce) {
  # Input clustering results -----
  Seurat <- read.csv(here("Simulations", "Data", paste0("Seurat", id, ".csv")))
  SC3 <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")))
  KMEANS <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")))

  # # Running Dune ----
  # clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
  # rownames(clusMat) <- Names
  # BPPARAM <- BiocParallel::MulticoreParam(32)
  # merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE)
  # Names <- as.character(Names)
  # chars <- c("sc3", "Monocle", "Seurat")
  # levels <- seq(from = 0, to = 1, by = .05)
  # stopMatrix <- lapply(levels, function(p){
  #   print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
  #   mat <- intermediateMat(merger = merger, p = p)
  #   suppressWarnings(rownames(mat) <- mat$cells)
  #   mat <- mat[Names, ]
  #   mat <- mat %>%
  #     select(-cells) %>%
  #     as.matrix()
  #   return(mat)
  # }) %>%
  #   do.call('cbind', args = .)
  # colnames(stopMatrix) <- lapply(levels, function(p){
  #   i <- as.character(round(100 * p))
  #   if (nchar(i) == 1) {
  #     i <- paste0("0", i)
  #   }
  #   return(paste(chars, i, sep = "-"))
  # }) %>% unlist()
  # print("...Full matrix")
  # mat <- cbind(as.character(Names), stopMatrix)
  # colnames(mat)[1] <- "cells"
  # 
  # write_csv(x = as.data.frame(mat), path = here("Simulations", "Data", paste0(id, "_Dune.csv")))
  # 
  # # Running Dune NMI ----
  # clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
  # rownames(clusMat) <- Names
  # BPPARAM <- BiocParallel::MulticoreParam(32)
  # merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE, metric = "NMI")
  # Names <- as.character(Names)
  # chars <- c("sc3", "Monocle", "Seurat")
  # levels <- seq(from = 0, to = 1, by = .05)
  # stopMatrix <- lapply(levels, function(p){
  #   print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
  #   mat <- intermediateMat(merger = merger, p = p)
  #   suppressWarnings(rownames(mat) <- mat$cells)
  #   mat <- mat[Names, ]
  #   mat <- mat %>%
  #     select(-cells) %>%
  #     as.matrix()
  #   return(mat)
  # }) %>%
  #   do.call('cbind', args = .)
  # colnames(stopMatrix) <- lapply(levels, function(p){
  #   i <- as.character(round(100 * p))
  #   if (nchar(i) == 1) {
  #     i <- paste0("0", i)
  #   }
  #   return(paste(chars, i, sep = "-"))
  # }) %>% unlist()
  # print("...Full matrix")
  # mat <- cbind(as.character(Names), stopMatrix)
  # colnames(mat)[1] <- "cells"
  # 
  # write_csv(x = as.data.frame(mat),
  #           path = here("Simulations", "Data", paste0(id, "_Dune_NMI.csv")))
  # 
  # # Do hierarchical merging with fraction of DE----
  # for (clustering in c("sc3", "Monocle", "Seurat")) {
  #   Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
  # }
  # cutoffs <- seq(from = 0, to = .01, by = .01)
  # res <- list()
  # for (clustering in c("sc3", "Monocle", "Seurat")) {
  #   print(clustering)
  #   Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  #   names(cutoffs) <- paste(clustering, cutoffs, sep = "_")
  #   res[[clustering]] <- map_dfc(cutoffs,
  #                                function(i){
  #                                  print(paste0("...", i))
  #                                  Rsec3 <- mergeClusters(Rsec2,
  #                                                         mergeMethod = "adjP",
  #                                                         plotInfo = "adjP",
  #                                                         cutoff = i,
  #                                                         clusterLabel = "Clusters",
  #                                                         plot = F,
  #                                                         DEMethod = "limma")
  #                                  return(Rsec3@clusterMatrix[,"Clusters"])
  #                                })
  # }
  # 
  # res <- do.call('cbind', res) %>% as.data.frame()
  # res$cells <- colnames(Rsec)
  # write_csv(res, path = here("Simulations", "Data", paste0(id, "_hierarchical_DE.csv")))
  # 
  # # Do hierarchical merging with cutting the tree ----
  # res <- list()
  # for (clustering in c("sc3", "Monocle", "Seurat")) {
  #   print(clustering)
  #   n <- n_distinct(get(clustering))
  #   cutoffs <- 5:n
  #   Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  #   Tree <- as.hclust(convertToDendrogram(Rsec2))
  #   names(cutoffs) <- paste(clustering, n - cutoffs, sep = "_")
  #   res[[clustering]] <- map_dfc(cutoffs,
  #                                function(cutoff){
  #                                  print(paste0("...", cutoff))
  #                                  return(cutree(Tree, k = cutoff))
  #                                })
  # }
  # res <- do.call('cbind', res) %>% as.data.frame()
  # res$cells <- colnames(Rsec)
  # write_csv(res, path = here("Simulations", "Data", paste0(id, "_hierarchical_Dist.csv")))
}