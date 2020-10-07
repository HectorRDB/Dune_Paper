# Packages to load ----
libs <- c("here", "tidyverse", "Seurat", "Dune")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper ----
compute_silhouette <- function(ref, label) {
  sce <- CreateSeuratObject(counts = counts(ref),
                            min.cells = 0,
                            min.features = 0,
                            project = "sl")
  sce <- AddMetaData(sce, as.data.frame(colData(ref)))
  sce <- NormalizeData(sce, verbose = FALSE)
  sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000,
                              verbose = FALSE)
  if ("human" %in% sce@meta.data) {
    sce <- ScaleData(object = sce, vars.to.regress = c("nCount_RNA", "human"),
                     verbose = FALSE)  
  } else {
    sce <- ScaleData(object = sce, vars.to.regress = "nCount_RNA",
                     verbose = FALSE)
  }
  sce <- RunPCA(object = sce, ndims.print = 1, npcs = 100,
                verbose = FALSE)
  sce <- RunUMAP(sce, verbose = FALSE, dims = 1:100)
  UMAP <- Embeddings(sce, reduction = "umap")
  dist_mat <- dist(UMAP)
  SL <- cluster::silhouette(label, dist = dist_mat)[,3] %>%  mean()
  return(SL)
}


# Functions ----
meanMethod_comp <- function(dataset, comp = "") {
  if (str_detect(dataset, "SMART")) {
    where <- "/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/"
    if (comp == "") {
      df <- readRDS(paste0(where, "Dune/", dataset, "_NMI_mergers.rds"))
    } else {
      df <- readRDS(paste0(where, "singleTree/", dataset, comp, "_NMI_merger.rds"))
    }
  } else {
    df <- readRDS(
      paste0("/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/",
             dataset, "_", comp, "_NMI_merger.rds"))
  }
  NMI <- NMIs(as.matrix(df$currentMat)) %>% colMeans()
  return(data.frame(clustering = names(NMI),
                    Value = NMI))
}

sihouette_comp <- function(loc, dataset, comp = "") {
  if (str_detect(dataset, "SMART")) {
    where <- "/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/"
    if (comp == "") {
      df <- readRDS(paste0(where, "Dune/", dataset, "_NMI_mergers.rds"))
    } else {
      df <- readRDS(paste0(where, "singleTree/", dataset, comp, "_NMI_merger.rds"))
    }
    sce <- readRDS(file = paste0(loc, dataset, "_filt.rds"))
    sce <- sce[, rownames(df$initialMat)]
  } else {
    df <- readRDS(
      paste0("/accounts/projects/epurdom/singlecell/Pancreas/Data/Dune/",
             dataset, "_", comp, "_NMI_merger.rds"))
    sce <- readRDS(file = paste0(loc, dataset, "_filt.rds"))
    sce <- sce[, rownames(df$initialMat)]
  }
  clusterings <- colnames(df$initialMat)
  names(clusterings) <- clusterings
  SLs <- lapply(clusterings, function(clus) {
   return(compute_silhouette(sce, df$currentMat[, clus])) 
  }) %>% unlist()
  return(data.frame(clustering = clusterings,
                    Value = SLs))
}

run_both <- function(loc, dataset, comp = "") {
  print(paste0("...", comp))
  df <- bind_rows(
    "MeanNMI" = meanMethod_comp(dataset, comp),
    "SL" = sihouette_comp(loc, dataset, comp),
    .id = "Metric"
  )
  return(df)
}
# Run Pancreas ----
run_pancreas_comps <- function(dataset) {
  print(dataset)
  loc <- "/scratch/users/singlecell/Pancreas/ProcessedData/"
  df <- bind_rows(
    "comp1" = run_both(loc, dataset, comp = "comp1"),
    "comp2" = run_both(loc, dataset, comp = "comp2"),
    "comp3" = run_both(loc, dataset, comp = "comp3"),
    .id = "comp"
  )
}

df <- bind_rows(
  "baron" = run_pancreas_comps(dataset = "baron"),
  "segerstolpe" = run_pancreas_comps(dataset = "segerstolpe"),
  .id = "dataset"
)

write.table(df, here("BestMethod", "Data", "Pancreas.txt"), row.names = FALSE)
# Run Brain ----
run_brain_comps <- function(dataset) {
  loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/"
  df <- bind_rows(
    "comp1" = run_both(loc, dataset, comp = ""),
    "comp2" = run_both(loc, dataset, comp = "_large2"),
    "comp3" = run_both(loc, dataset, comp = "_large3"),
    .id = "comp"
  )
}

df <- bind_rows(
  "SMARTer_cells_MOp" = run_brain_comps(dataset = "SMARTer_cells_MOp"),
  "SMARTer_nuclei_MOp" = run_brain_comps(dataset = "SMARTer_nuclei_MOp"),
  .id = "dataset"
)

write.table(df, here("BestMethod", "Data", "Brain.txt"), row.names = FALSE)
