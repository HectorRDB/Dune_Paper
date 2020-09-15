# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne", "clusterExperiment",
          "BiocParallel", "zinbwave", "matrixStats", "ggplot2", "reticulate",
          "purrr", "mclust", "flexclust", "monocle3", "scater", "readr", "aricode",
          "Dune", "cluster")
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
  df <- colData(sce) %>% as.data.frame()
  rownames(df) <- df$Cell
  sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K',
                                meta.data = df)
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
  sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                  dispersion.function = LogVMR, do.plot = F)
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nCount_RNA")
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
  write.csv(x = K_MEANS, row.names = FALSE,
            file = here("Simulations", "Data", paste0("tSNE-KMEANS", id, ".csv")))

  # UMAP K-Means ----
  sce_Seurat <- scater::runUMAP(sce_Seurat, ntop = 2000, ncomponents = 3)
  UMAP <- reducedDim(sce_Seurat, "UMAP")
  ks <- seq(20, 50, 5)
  names(ks) <- ks
  K_MEANS <- map_dfc(ks, function(k){
    return(kmeans(UMAP, centers = k)$cluster)
  })
  K_MEANS$cells <- rownames(UMAP)
  write.csv(x = K_MEANS, row.names = FALSE,
            file = here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")))

  # RSEC ----
  sequential <- FALSE
  subsample <- T
  clusterFunction <- "pam"
  # sce <- sce_og
  # sce <- scater::runPCA(sce)
  # sce <- scater::runUMAP(sce)
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
                mergeLogFCcutoff = 1, consensusMinSize = 10
                )
  ))

  # Saving objects
  saveRDS(sce, here("Simulations", "Data", paste0("Merger_", id, ".rds")))

  # Running SC3 ----
  print("... Running SC3")
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  sce <- sc3_estimate_k(sce)
  K <- metadata(sce)$sc3$k_estimation

  print(paste0("...... The optimal number of clusters defined by sc3 is ", K))
  ks <- seq(from = 20, to = 50, by = 5)
  names(ks) <- ks

  sc3 <- map_df(ks, function(k){
    SC3 <- sc3(sce, ks = k, svm_max = ncol(sce) + 1, biology = FALSE,
               gene_filter = F, n_cores = NCORES, rand_seed = 786907)
    SC3 <- colData(SC3)[, paste0("sc3_", k, "_clusters")] %>% as.numeric()
    return(SC3)
  })
  sc3$cells <- colnames(sce)
  SC3s <- sc3
  write.csv(x = SC3s, row.names = FALSE,
            file = here("Simulations", "Data", paste0("SC3", id, ".csv")))

  return()
  
}

run_merging_methods <- function(Rsec, sce, id) {
  # Input clustering results -----
  SC3 <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")),
                  stringsAsFactors = FALSE) %>%
    arrange(cells)
  UMAP_KMEANS <- read.csv(here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  TSNE_KMEANS <- read.csv(here("Simulations", "Data", paste0("tSNE-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)

  # Running Dune ----
  Names <- SC3$cells
  clusMat <- data.frame("SC3" = SC3$X40,
                        "UMAP_KMEANS" = UMAP_KMEANS$X40,
                        "TSNE_KMEANS" = TSNE_KMEANS$X40)
  rownames(clusMat) <- Names
  BPPARAM <- BiocParallel::MulticoreParam(8)
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
      dplyr::select(-cells) %>%
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
  BPPARAM <- BiocParallel::MulticoreParam(8)
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
      dplyr::select(-cells) %>%
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
  Rsec <- Rsec[, Names]
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    Rsec <- addClusterings(Rsec, clusMat[,clustering], clusterLabels = clustering)
  }
  cutoffs <- seq(from = 0, to = .15, by = .005)
  res <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    print(clustering)
    Rsec2 <- Rsec
    counts(Rsec2) <- as.matrix(counts(Rsec2))
    Rsec2 <- makeDendrogram(Rsec2, whichCluster = clustering)
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
    n <- n_distinct(clusMat[, clustering])
    cutoffs <- 5:n
    Rsec2 <- Rsec
    counts(Rsec2) <- as.matrix(counts(Rsec2))
    Rsec2 <- makeDendrogram(Rsec2, whichCluster = clustering)
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

evaluate_clustering_methods <- function(sce, id) {
  # Input clustering results -----
  SC3 <- read.csv(here("Simulations", "Data", paste0("SC3", id, ".csv")),
                  stringsAsFactors = FALSE) %>%
    arrange(cells)
  UMAP_KMEANS <- read.csv(here("Simulations", "Data", paste0("UMAP-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  TSNE_KMEANS <- read.csv(here("Simulations", "Data", paste0("tSNE-KMEANS", id, ".csv")),
                          stringsAsFactors = FALSE) %>%
    arrange(cells)
  
  ref <- data.frame(groups = sce$Group, 
                    cells = colnames(sce)) %>%
    arrange(cells) %>%
    filter(cells %in% SC3$cells)
  
  # ARI with ref for the methods ----
  params <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    df <- get(clustering) %>% dplyr::select(-cells)
    params[[clustering]] <- data.frame(
      "Value" = lapply(df, adjustedRandIndex, y = ref$groups) %>% unlist(),
      "n_clus" = lapply(df, n_distinct) %>% unlist())
  }
  params <- bind_rows(params, .id = "clustering") %>%
    mutate(method = "param")
  ARI <- list()
  for (method in c("Dune", "Dune_NMI", "DE", "Dist")) {
    df <- read.csv(here("Simulations", "Data", paste0(method, "_", id, ".csv")),
                    stringsAsFactors = FALSE) %>%
      arrange(cells)
      
    ARI[[method]] <- data.frame(
      "Value" = lapply(df %>% dplyr::select(-cells), 
                       adjustedRandIndex, y = ref$groups) %>% unlist(),
      "Name" = colnames(df %>% dplyr::select(-cells)),
      "n_clus" = lapply(df %>% dplyr::select(-cells), n_distinct) %>% unlist(),
      "clustering" = word(colnames(df %>% dplyr::select(-cells)), 1, sep = "\\."),
      stringsAsFactors = FALSE
      ) %>%
      mutate(level = str_remove_all(Name, clustering),
             level = str_remove(level, "^\\."),
             level = str_remove(level, "^\\_") %>% as.numeric())
  }
  ARI <- bind_rows(ARI, .id = "method")
  ARI <- bind_rows(ARI, params)
  
  # NMI with ref for the methods ----
  params <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    df <- get(clustering) %>% dplyr::select(-cells)
    params[[clustering]] <- data.frame(
      "Value" = lapply(df, NMI, c2 = ref$groups, variant = "sum") %>% unlist(),
      "n_clus" = lapply(df, n_distinct) %>% unlist())
  }
  params <- bind_rows(params, .id = "clustering") %>%
    mutate(method = "param")
  NMI_ <- list()
  for (method in c("Dune", "Dune_NMI", "DE", "Dist")) {
    df <- read.csv(here("Simulations", "Data", paste0(method, "_", id, ".csv")),
                   stringsAsFactors = FALSE) %>%
      arrange(cells)
    NMI_[[method]] <- data.frame(
      "Value" = lapply(df %>% dplyr::select(-cells), NMI, c2 = ref$groups, variant = "sum") %>% unlist(),
      "Name" = colnames(df %>% dplyr::select(-cells)),
      "n_clus" = lapply(df %>% dplyr::select(-cells), n_distinct) %>% unlist(),
      "clustering" = word(colnames(df %>% dplyr::select(-cells)), 1, sep = "\\."),
      stringsAsFactors = FALSE
    ) %>%
      mutate(level = str_remove_all(Name, clustering),
             level = str_remove(level, "^\\."),
             level = str_remove(level, "^\\_") %>% as.numeric())
  }
  NMI_ <- bind_rows(NMI_, .id = "method")
  NMI_ <- bind_rows(NMI_, params)
  
  # Silhouettes ----
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
  df <- colData(sce) %>% as.data.frame()
  rownames(df) <- df$Cell
  sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K',
                                meta.data = df)
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize")
  sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                  dispersion.function = LogVMR, do.plot = F)
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nCount_RNA")
  sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 100)
  sSeurat <- RunUMAP(sSeurat, verbose = FALSE, dims = 1:100)
  UMAP <- Embeddings(sSeurat, reduction = "umap")
  UMAP <- UMAP[SC3$cells, ]
  dist_mat <- dist(UMAP)
  params <- list()
  for (clustering in c("SC3", "UMAP_KMEANS", "TSNE_KMEANS")) {
    df <- get(clustering) %>% dplyr::select(-cells)
    params[[clustering]] <- data.frame(
      "Value" = lapply(df, function(i) {
        silhouette(i, dist = dist_mat)[,3] %>% 
          mean() %>%
          return()
      }) %>% unlist(),
      "n_clus" = lapply(df, n_distinct) %>% unlist())
  }
  params <- bind_rows(params, .id = "clustering") %>%
    mutate(method = "param")
  SL <- list()
  for (method in c("Dune", "Dune_NMI", "DE", "Dist")) {
    df <- read.csv(here("Simulations", "Data", paste0(method, "_", id, ".csv")),
                   stringsAsFactors = FALSE) %>%
      arrange(cells)
    SL[[method]] <- data.frame(
      "Value" = lapply(df %>% select(-cells), function(i) {
        if(n_distinct(i) == 1) return(0)
        silhouette(i, dist = dist_mat)[,3] %>% 
          mean() %>%
          return()
      }) %>% unlist(),
      "Name" = colnames(df %>% dplyr::select(-cells)),
      "n_clus" = lapply(df %>% dplyr::select(-cells), n_distinct) %>% unlist(),
      "clustering" = word(colnames(df %>% dplyr::select(-cells)), 1, sep = "\\."),
      stringsAsFactors = FALSE
    ) %>%
      mutate(level = str_remove_all(Name, clustering),
             level = str_remove(level, "^\\."),
             level = str_remove(level, "^\\_") %>% as.numeric())
  }
  SL <- bind_rows(SL, .id = "method")
  SL <- bind_rows(SL, params)
  
  # Save results ----
  res <- bind_rows("ARI" = ARI, "NMI" = NMI_, "SL" = SL, .id = "Metric")
  write_csv(res, col_names = TRUE, 
            path = here("Simulations", "Data", paste0("Metric_Ref_", id, ".csv"))
  )
}
