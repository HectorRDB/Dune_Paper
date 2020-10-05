# Plotting ----
plot_UMAP <- function(df, type) {
  p <- ggplot(bind_cols(df, data.frame("type" = type)),
              aes(x = V1, y = V2, col = type)) +
    geom_point(size = .5, alpha = .5) +
    theme_bw() +
    labs(x = "UMAP_1", y = "UMAP_2")
  return(p)
}

# Seurat ----
map_seurat_proba <- function(ref, target, labels) {
  ## Normalize ----
  dtlist <- list()
  dtlist$ref <- CreateSeuratObject(counts = counts(ref),
                                   min.cells = 0,
                                   min.features = 0,
                                   project = "mapping")
  dtlist$target <- CreateSeuratObject(counts = counts(target),
                                      min.cells = 0,
                                      min.features = 0,
                                      project = "mapping")
  dtlist$ref <- AddMetaData(dtlist$ref, as.data.frame(colData(ref)))
  dtlist <- lapply(dtlist, function(sce){
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
    if ("human" %in% sce@meta.data) {
      sce <- ScaleData(object = sce, vars.to.regress = c("nCount_RNA", "human"))  
    } else {
      sce <- ScaleData(object = sce, vars.to.regress = "nCount_RNA")
    }
    return(sce)
  })
  ## Mapping
  sim.anchors <- FindTransferAnchors(reference = dtlist$ref,
                                     query = dtlist$target,
                                     dims = 1:30,
                                     project.query = FALSE)
  prediction <- TransferData(anchorset = sim.anchors,
                             refdata = as.character(labels),
                             dims = 1:30)
  return(prediction$prediction.score.max)
}


# SingleR----
map_singleR_proba <- function(ref, target, labels) {
  preds <- SingleR(test = logNormCounts(target), ref = logNormCounts(ref),
                   labels = labels)
  return(preds$tuning.scores$first)
}

# Together
start_finish_clustering <- function(ref, target, merger, clustering) {
  print(paste0(".... ", clustering))
  labels <- data.frame(init = merger[, paste0(clustering, ".00")],
                       final = merger[, paste0(clustering, ".100")])
  print("...... SingleR")
  proba_singleR <- map_dfc(labels, map_singleR_proba, ref = ref, target = target)
  imp_singleR <- mean(proba_singleR$final - proba_singleR$init)
  print("...... Seurat")
  proba_seurat <- map_dfc(labels, map_seurat_proba, ref = ref, target = target)
  imp_seurat <- mean(proba_seurat$final - proba_seurat$init)
  return(data.frame('singleR' = imp_singleR, "Seurat" = imp_seurat))
}

start_finish <- function(ref, target, merger) {
  clusterings <- c("sc3", "Seurat", "Monocle")
  merger <- merger[colnames(ref), ]
  df <- lapply(clusterings, start_finish_clustering,
               ref = ref, target = target, merger = merger)
  df <- bind_rows(df, .id = "clustering")
}

switch_refs <- function(comp, sce1, sce2, n1, n2, m_loc) {
  df1 <- read.csv(file = paste0(m_loc, opt$n1, "_", comp, "_Dune_NMI.csv"))
  df2 <- read.csv(file = paste0(m_loc, opt$n2, "_", comp, "_Dune_NMI.csv"))
  print(comp)
  print(paste0(".. ", n1))
  N1 <- start_finish(ref = sce1, target = sce2, merger = df1)
  print(paste0(".. ", n2))
  N2 <- start_finish(ref = sce2, target = sce1, merger = df2)
  probas <- list(N1, N2)
  names(probas) <- c(n1, n2)
  df <- bind_rows(probas, .id = "Ref")
  return(df)
}

all_comps <- function(sce1, sce2, n1, n2, m_loc, comps) {
  names(comps) <- comps
  df <- lapply(comps, switch_refs, 
               sce1 = sce1, sce2 = sce2, n1 = n1, n2 = n2, m_loc = m_loc)
  df <- bind_rows(df, .id = "comps")
}