# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "purrr", "mclust", "flexclust", "monocle3", "scater", "readr")
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
  df <- colData(sce) %>% as.data.frame() %>%
    mutate(Batch = as.factor(Batch))
  rownames(df) <- df$Cell
  sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K',
                                meta.data = df)
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
  sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                  dispersion.function = LogVMR, do.plot = F)
  if (nlevels(df$Batch) > 1) {
    sSeurat <- ScaleData(object = sSeurat, vars.to.regress = c("nCount_RNA", "Batch"))
  } else {
    sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nCount_RNA")
  }
  sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)
  sSeurat <- RunUMAP(sSeurat, verbose = FALSE, dims = 1:100)
  p <- UMAPPlot(sSeurat, group.by = "Group")
  sce <- as.SingleCellExperiment(sSeurat)
  sce_Seurat <- sce
  ggsave(here("Simulations", "Figures", paste0("UMAP_", id, ".png")), p)
  
  # TSNE K-Means ----
  print("... Running RtsneKmeans")
  sce_Seurat <- scater::runPCA(sce_Seurat, ntop = 2000)
  sce_Seurat <- scater::runTSNE(sce_Seurat, ntop = 2000, ncomponents = 3, perplexity = 30)
  TSNE <- reducedDim(sce_Seurat, "TSNE")
  ks <- seq(20, 50, 5)
  names(ks) <- ks
  K_MEANS <- map_dfc(ks, function(k){
    return(kmeans(TSNE, centers = k)$cluster)
  })
  K_MEANS$cells <- rownames(TSNE)
  write.csv(K_MEANS, here("Simulations", "Data", paste0("tSNE-KMEANS", id, ".csv")),
            col.names = TRUE, row.names = FALSE)
  
  # UMAP K-Means ----
  sce_Seurat <- scater::runUMAP(sce_Seurat, ntop = 2000, ncomponents = 3)
  UMAP <- reducedDim(sce_Seurat, "UMAP")
  ks <- seq(20, 50, 5)
  names(ks) <- ks
  K_MEANS <- map_dfc(ks, function(k){
    return(kmeans(UMAP, centers = k)$cluster)
  })
  K_MEANS$cells <- rownames(UMAP)
  write.csv(K_MEANS, here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")),
            col.names = TRUE, row.names = FALSE)
  
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
  write.csv(SC3s, here("Simulations", "Data", paste0("SC3", id, ".csv")),
            col.names = TRUE, row.names = FALSE)
  
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
  return()
  
}

run_merging_methods <- function(Rsec, id, sce) {
  # Input clustering results -----
  SC3 <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")),
                  stringsAsFactors = FALSE) %>%
    arrange(cells)
  UMAP_KMEANS <- read.csv(here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  TSNE_KMEANS <- read.csv(here("Simulations", "Data", paste0("TSNE-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)

  # Running Dune ----
  Names <- SC3$cells
  clusMat <- data.frame("SC3" = SC3$`40`,
                        "UMAP_KMEANS" = UMAP_KMEANS$`40`,
                        "TSNE_KMEANS" = TSNE_KMEANS$`40`)
  rownames(clusMat) <- Names
  BPPARAM <- BiocParallel::MulticoreParam(32)
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE)
  Names <- as.character(Names)
  chars <- c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")
  levels <- seq(from = 0, to = 1, by = .05)
  stopMatrix <- lapply(levels, function(p){
    print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
    mat <- intermediateMat(merger = merger, p = p)
    suppressWarnings(rownames(mat) <- mat$cells)
    mat <- mat[Names, ]
    mat <- mat %>%
      select(-cells) %>%
      as.matrix()
    return(mat)
  }) %>%
    do.call('cbind', args = .)
  colnames(stopMatrix) <- lapply(levels, function(p){
    i <- as.character(round(100 * p))
    if (nchar(i) == 1) {
      i <- paste0("0", i)
    }
    return(paste(chars, i, sep = "-"))
  }) %>% unlist()
  print("...Full matrix")
  mat <- cbind(as.character(Names), stopMatrix)
  colnames(mat)[1] <- "cells"

  write_csv(x = as.data.frame(mat),
            path = here("Simulations", "Data", paste0("Dune_", id, ".csv")), 
            col_names = TRUE)

  # Running Dune NMI ----
  BPPARAM <- BiocParallel::MulticoreParam(32)
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE, metric = "NMI")
  Names <- as.character(Names)
  chars <- c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")
  levels <- seq(from = 0, to = 1, by = .05)
  stopMatrix <- lapply(levels, function(p){
    print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
    mat <- intermediateMat(merger = merger, p = p)
    suppressWarnings(rownames(mat) <- mat$cells)
    mat <- mat[Names, ]
    mat <- mat %>%
      select(-cells) %>%
      as.matrix()
    return(mat)
  }) %>%
    do.call('cbind', args = .)
  colnames(stopMatrix) <- lapply(levels, function(p){
    i <- as.character(round(100 * p))
    if (nchar(i) == 1) {
      i <- paste0("0", i)
    }
    return(paste(chars, i, sep = "-"))
  }) %>% unlist()
  print("...Full matrix")
  mat <- cbind(as.character(Names), stopMatrix)
  colnames(mat)[1] <- "cells"

  write_csv(x = as.data.frame(mat),
            path = here("Simulations", "Data", paste0("Dune_NMI_", id, ".csv")), 
            col_names = TRUE)

  # Do hierarchical merging with fraction of DE----
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
  }
  cutoffs <- seq(from = 0, to = .01, by = .01)
  res <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    print(clustering)
    Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
    names(cutoffs) <- paste(clustering, cutoffs, sep = "_")
    res[[clustering]] <- map_dfc(cutoffs,
                                 function(i){
                                   print(paste0("...", i))
                                   Rsec3 <- mergeClusters(Rsec2,
                                                          mergeMethod = "adjP",
                                                          plotInfo = "adjP",
                                                          cutoff = i,
                                                          clusterLabel = "Clusters",
                                                          plot = F,
                                                          DEMethod = "limma")
                                   return(Rsec3@clusterMatrix[,"Clusters"])
                                 })
  }

  res <- do.call('cbind', res) %>% as.data.frame()
  res$cells <- colnames(Rsec)
  write_csv(res, col_names = TRUE,
            path = here("Simulations", "Data", paste0("DE_", id, ".csv")))

  # Do hierarchical merging with cutting the tree ----
  res <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    print(clustering)
    n <- n_distinct(get(clustering))
    cutoffs <- 5:n
    Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
    Tree <- as.hclust(convertToDendrogram(Rsec2))
    names(cutoffs) <- paste(clustering, n - cutoffs, sep = "_")
    res[[clustering]] <- map_dfc(cutoffs,
                                 function(cutoff){
                                   print(paste0("...", cutoff))
                                   return(cutree(Tree, k = cutoff))
                                 })
  }
  res <- do.call('cbind', res) %>% as.data.frame()
  res$cells <- colnames(Rsec)
  write_csv(res, col_names = TRUE,
            path = here("Simulations", "Data", paste0("Dist_", id, ".csv")))
}

evaluate_clustering_methods <- function(id, sce) {
  # Input clustering results -----
  SC3 <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")),
                  stringsAsFactors = FALSE) %>%
    arrange(cells)
  UMAP_KMEANS <- read.csv(here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  TSNE_KMEANS <- read.csv(here("Simulations", "Data", paste0("TSNE-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  
  ref <- data.frame(groups = sce$Group, 
                    cells = colnames(sce)) %>%
    arrange(cells) %>%
    filter(cells %in% SC3$cells)
  
  res <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    df <- get(clustering) %>% select(-cells)
    res[[clustering]] <- data.frame(
      "ARI" = lapply(df, adjustedRandIndex, y = ref$groups) %>% unlist(),
      "n_clus" = lapply(df, n_distinct) %>% unlist())
  }
  
  res <- bind_rows(res, .id = "method")
  write_csv(res, here("Simulations", "Data", paste0("ARI_Param_", id, ".csv")))
  
  # ARI with ref for the methods
  
}