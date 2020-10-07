# Packages to load ----
libs <- c("here", "tidyverse", "Seurat", "Dune", "SingleCellExperiment", "parallel", "mclust")
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

load_Dune <- function(dataset, comp) {
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
  return(df)
}

load_sce <- function(dataset, loc, comp) {
  if (str_detect(dataset, "SMART")) {
    where <- "/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/"
    sce <- readRDS(file = paste0(loc, dataset, "_filt.rds"))
  } else {
    sce <- readRDS(file = paste0(loc, dataset, "_filt.rds"))
  }
  return(sce)
}

load_gold_standard <- function(dataset) {
  if (str_detect(dataset, "SMART")) {
    where <- "/accounts/projects/epurdom/singlecell/allen/allen40K/Pipeline_Brain/data/"
    allen_clusters <- read.csv(paste0(where, "Smart-Seq/", dataset, "_cluster.membership.csv"),
                               col.names = c("cells", "cluster_id"))
    clusters <- read.csv(paste0(where, "Smart-Seq/", dataset, "_cluster.annotation.csv"),
                         header = T)
    allen_clusters <- full_join(allen_clusters, clusters) %>%
      select(cells, subclass_label) %>%
      dplyr::rename(ref = subclass_label)
    return(allen_clusters)
  } else {
    gold_clusters <- read.csv(
      paste0("/accounts/projects/epurdom/singlecell/Pancreas/Data/", 
             str_to_title(dataset), "/", dataset, "_meta.csv")) %>% 
      select(X, cell_type1) %>%
      dplyr::rename(cells = X, ref = cell_type1)
    return(gold_clusters)
  }
}

# Functions ----
meanMethod_comp <- function(dataset, comp = "") {
  df <- load_Dune(dataset, comp)
  NMI <- NMIs(as.matrix(df$currentMat)) %>% colMeans()
  return(data.frame(clustering = names(NMI),
                    Value = NMI))
}

sihouette_comp <- function(loc, dataset, comp = "") {
  df <- load_Dune(dataset, comp)
  sce <- load_sce(dataset, loc, comp)
  sce <- sce[, rownames(df$initialMat)]
  clusterings <- colnames(df$initialMat)
  names(clusterings) <- clusterings
  SLs <- lapply(clusterings, function(clus) {
   return(compute_silhouette(sce, df$currentMat[, clus])) 
  }) %>% unlist()
  return(data.frame(clustering = clusterings, Value = SLs))
}

ARI_ref <- function(dataset, comp = "") {
  df <- load_Dune(dataset, comp)
  ref <- load_gold_standard(dataset)
  rownames(ref) <- ref$cells
  ref <- ref[rownames(df$currentMat), ]
  clusterings <- colnames(df$initialMat)
  ARI <- lapply(df$currentMat, adjustedRandIndex, y = ref$ref) %>%
    unlist()
  return(data.frame(clustering = clusterings, Value = ARI))
}

rep_comp <- function(dataset, comp) {
  if (str_detect(dataset, "SMART")) {
    where <- paste0("/accounts/projects/epurdom/singlecell/allen/allen40K/",
                    "Pipeline_Brain/data/Replicability/Dune_Smart/")
    if (comp == "") {
      comp <- "Normal"
    } else {
      comp <- str_remove(comp, "^_") %>% str_to_title()
    }
    df <- read.table(
      paste0(where, comp, "_NMI/consensus_cluster_replicability.txt"),
      header = TRUE)
  } else {
    where <- "/accounts/projects/epurdom/singlecell/Pancreas/Data/Replicability/Dune_NMI/"
    df <- read.table(
      paste0(where, comp, "/consensus_cluster_replicability.txt"),
      header = TRUE)
  }
  df <- df  %>%
    filter(level == 100) %>%
    select(clustering_name, fraction_replicable_cells) %>%
    dplyr::rename(clustering = clustering_name,
                  Value = fraction_replicable_cells)
  return(df)
}

run_all <- function(loc, dataset, comp = "") {
  print(paste0("...", comp))
  df <- bind_rows(
    "MeanNMI" = meanMethod_comp(dataset, comp),
    "SL" = sihouette_comp(loc, dataset, comp),
    "ARI" = ARI_ref(dataset, comp),
    "Rep" = rep_comp(dataset, comp),
    .id = "Metric"
  )
  return(df)
}
# Run Pancreas ----
run_pancreas_comps <- function(dataset) {
  print(dataset)
  loc <- "/scratch/users/singlecell/Pancreas/ProcessedData/"
  df <- bind_rows(
    "comp1" = run_all(loc, dataset, comp = "comp1"),
    "comp1" = run_all(loc, dataset, comp = "comp2"),
    "comp1" = run_all(loc, dataset, comp = "comp3"),
    .id = "comp"
  )
}

# df <- bind_rows(
#   "baron" = run_pancreas_comps(dataset = "baron"),
#   "segerstolpe" = run_pancreas_comps(dataset = "segerstolpe"),
#   .id = "dataset"
# )

# write.table(df, here("BestMethod", "Data", "Pancreas.txt"), row.names = FALSE)
# Run Brain ----
run_brain_comps <- function(dataset) {
  loc <- "/scratch/users/singlecell/MiniAtlas/data/rds/"
  df <- bind_rows(
    "comp1" = run_all(loc, dataset, comp = ""),
    "comp1" = run_all(loc, dataset, comp = "_large2"),
    "comp1" = run_all(loc, dataset, comp = "_large3"),
    .id = "comp"
  )
}

df <- bind_rows(
  "SMARTer_cells_MOp" = run_brain_comps(dataset = "SMARTer_cells_MOp"),
  "SMARTer_nuclei_MOp" = run_brain_comps(dataset = "SMARTer_nuclei_MOp"),
  .id = "dataset"
)

write.table(df, here("BestMethod", "Data", "Brain.txt"), row.names = FALSE)
