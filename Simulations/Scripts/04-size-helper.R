# Packages ----
libs <- c("splatter", "here", "scater", "scran", "Seurat", "dplyr",
          "stringr", "SingleCellExperiment", "SC3", "Rtsne",
          "BiocParallel", "matrixStats", "ggplot2", "reticulate",
          "purrr", "mclust", "flexclust", "monocle3", "scater", "readr", "aricode",
          "Dune", "cluster")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
source(here("Simulations", "Scripts", "01-create_data.R"))
create_and_clean <- function(nCells) {
  sce <- create_simple_balanced_data(nCells = nCells, nClus = 30, 
                                     DE =.1, seed = runif(1, 0, 100))
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
  sce <- computeSumFactors(sce, min.mean = 0.1,
                           sizes = pmin(ncol(sce), seq(20, 120, 20) )%>% unique())
  logcounts(sce) <- scater::normalizeCounts(sce)
  df <- colData(sce) %>% as.data.frame()
  rownames(df) <- df$Cell
  sSeurat <- CreateSeuratObject(counts = assays(sce)$counts, project = 'allen40K',
                                meta.data = df)
  sSeurat <- NormalizeData(object = sSeurat, normalization.method = "LogNormalize",
                           verbose = FALSE)
  sSeurat <- FindVariableFeatures(object = sSeurat, mean.function = ExpMean,
                                  dispersion.function = LogVMR, verbose = FALSE)
  sSeurat <- ScaleData(object = sSeurat, vars.to.regress = "nCount_RNA",
                       verbose = FALSE)
  sSeurat <- RunPCA(object = sSeurat, ndims.print = 1, npcs = 98, verbose = FALSE)
  sSeurat <- RunUMAP(sSeurat, verbose = FALSE, dims = 1:98)
  sce <- as.SingleCellExperiment(sSeurat)

}

run_clusterings <- function(sce) {
  print(paste0(".. ", ncol(sce)))
  NCORES <- 32
  # Create data ----
  # We follow the workflow from Duo et al 2018
  # https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
  # TSNE K-Means ----
  print("... Running RtsneKmeans")
  sce <- scater::runPCA(sce, ntop = 2000)
  sce <- scater::runTSNE(sce, ntop = 2000, ncomponents = 3, perplexity = 30)
  TSNE <- reducedDim(sce, "TSNE")
  ks <- seq(30, 50, 5)
  names(ks) <- ks
  tSNE_KMEANS <- map_dfc(ks, function(k){
    return(kmeans(TSNE, centers = k)$cluster)
  })
  tSNE_KMEANS$cells <- rownames(TSNE)
  
  # UMAP K-Means ----
  sce <- scater::runUMAP(sce, ntop = 2000, ncomponents = 3)
  UMAP <- reducedDim(sce, "UMAP")
  ks <- seq(30, 50, 5)
  names(ks) <- ks
  UMAP_KMEANS <- map_dfc(ks, function(k){
    return(kmeans(UMAP, centers = k)$cluster)
  })
  UMAP_KMEANS$cells <- rownames(UMAP)
  
  # Running SC3 ----
  print("... Running SC3")
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))
  
  ks <- seq(from = 30, to = 50, by = 5)
  names(ks) <- ks
  
  SC3 <- sc3(sce, ks = ks, svm_max = ncol(sce) + 1, biology = FALSE,
        gene_filter = F, n_cores = NCORES, rand_seed = 786907)
  sc3 <- colData(SC3)[, paste0("sc3_", ks, "_clusters")]
  colnames(sc3) <- names(ks)
  sc3$cells <- colnames(sce)
  return(list("sc3" = sc3, "UMAP_KMEANS" = UMAP_KMEANS, "tSNE_KMEANS" = tSNE_KMEANS))
  
}

run_Dune <- function(clusMat) {
  print(paste0(".. ", nrow(clusMat)))
  Names <- clusMat$cells
  # Running Dune NMI ----
  BPPARAM <- BiocParallel::MulticoreParam(32)
  clusMat <- clusMat %>% select(-cells) %>% map_dfc(as.numeric) %>%
    as.matrix()
  rownames(clusMat) <- Names
  colnames(clusMat) <- c("sc3", "UMAP_KMEANS", "tSNE_KMEANS")
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE,
                  metric = "NMI")
  Names <- as.character(Names)
  chars <- c("sc3", "UMAP_KMEANS", "tSNE_KMEANS")
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
  mat <- cbind(as.character(Names), stopMatrix)
  colnames(mat)[1] <- "cells"
  return(as.data.frame(mat))
}

evaluate_clustering_methods <- function(sce, merger) {
  print(paste0(".. ", ncol(sce)))
  ref <- data.frame(groups = sce$Group,
                    cells = colnames(sce)) %>%
    arrange(cells) %>%
    filter(cells %in% merger$cells)
  merger <- merger %>% arrange(cells)
  ARI <- data.frame(
    "Value" = lapply(merger %>% dplyr::select(-cells),
                     adjustedRandIndex, y = ref$groups) %>% unlist(),
    "Name" = colnames(merger %>% dplyr::select(-cells)),
    "n_clus" = lapply(merger %>% dplyr::select(-cells), n_distinct) %>% unlist(),
    "clustering" = word(colnames(merger %>% dplyr::select(-cells)), 1, sep = "-"),
    stringsAsFactors = FALSE
    ) %>%
    mutate(level = str_remove_all(Name, clustering),
           level = str_remove(level, "-") %>% as.numeric())
  return(ARI)
}
