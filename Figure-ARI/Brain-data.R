# Packages ----
library(knitr)
libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "pracma", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper functions ----
setwd(here("Brain"))
reload(inst("here"))
interpolate <- function(df, ns) {
  if (any(df$n_clus == ns)) {
    return(df %>% filter(n_clus >= ns))
  } else {
    df <- df %>% arrange(desc(n_clus))
    cutoff <- which(df$n_clus < ns)[1]
    slope <- (df$ARI[cutoff] - df$ARI[cutoff - 1]) / 
      (df$n_clus[cutoff] - df$n_clus[cutoff - 1])
    intercept <- df$ARI[cutoff]
    filt <- df %>% filter(n_clus >= ns) %>%
      add_row(n_clus = ns,
              ARI = slope * (ns - df$n_clus[cutoff]) + intercept)
    return(filt %>% arrange(n_clus))
  }
}
## function to get the ARI ----
comp_merger_with_ref <- function(clusMat, clus, ref) {
  Mat <- clusMat %>%
    arrange(cells) %>%
    as.matrix()
  return(adjustedRandIndex(Mat[,clus], ref))
}
n_clus  <- function(clusMat, clus) {
  return(n_distinct(as.matrix(clusMat)[,clus]))
}
comp_tree_with_ref <- function(x, ref) {
  return(c("n_clus" = n_distinct(x), "ARI" = adjustedRandIndex(x, ref)))
}

## Load Dune ----
comp_dune_ref <- function(dataset, comp = "", ref, metric = "") {
  if (comp == "") {
    df <- readRDS(here("data", "Dune",
                       paste0(dataset, metric,  "_mergers.rds")))
  } else {
    df <- readRDS(here("data", "singleTree",
                       paste0(dataset, comp, metric, "_merger.rds")))
  }
  if (is.null(df$metric)) df$metric <- "ARI"
  if (is.null(df$ImpMetric)) df$ImpMetric <- df$ImpARI
  
  ARI_ref_sc3 <- data.frame(
    "n_clus" = functionTracking(df, n_clus, clus = "sc3"),
    "ARI" = functionTracking(df, comp_merger_with_ref, clus = "sc3",
                             ref = ref)) %>%
    arrange(n_clus) %>%
    distinct()
  ARI_ref_sc3 <- trapz(x = ARI_ref_sc3$n_clus, y = ARI_ref_sc3$ARI)
  
  ARI_ref_seurat <- data.frame(
    "n_clus" = functionTracking(df, n_clus, clus = "Seurat"),
    "ARI" = functionTracking(df, comp_merger_with_ref, clus = "Seurat",
                             ref = ref)) %>%
    arrange(n_clus) %>%
    distinct()
  ARI_ref_seurat <- trapz(x = ARI_ref_seurat$n_clus, y = ARI_ref_seurat$ARI)
  
  ARI_ref_monocle <- data.frame(
    "n_clus" = functionTracking(df, n_clus, clus = "Monocle"),
    "ARI" = functionTracking(df, comp_merger_with_ref, clus = "Monocle",
                             ref = ref)) %>%
    arrange(n_clus) %>%
    distinct()
  ARI_ref_monocle <- trapz(x = ARI_ref_monocle$n_clus, y = ARI_ref_monocle$ARI)
  
  df <- data.frame("comp" = comp,
                   "method" = c("SC3", "Seurat", "Monocle"),
                   "AUARIC" = c(ARI_ref_sc3, ARI_ref_seurat, ARI_ref_monocle))
  return(df)
}

## Load Tree data ----
comp_tree_ref <- function(dataset, comp, ref, type) {
  if (comp == "") {
    merger <- readRDS(here("data", "Dune",
                           paste0(dataset, "_mergers.rds")))
  } else {
    merger <- readRDS(here("data", "singleTree",
                           paste0(dataset, comp, "_merger.rds")))
  }
  ns <- lapply(merger$currentMat, n_distinct)
  
  ARI_ref_sc3 <- read.csv(here("Data", "singleTree",
                               paste0(dataset, comp, "_hierarchical_",
                                      type, ".csv"))) %>%
    arrange(cells) %>%
    dplyr::select(starts_with("sc3")) %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    distinct() %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$sc3)
  ARI_ref_sc3 <- trapz(x = ARI_ref_sc3$n_clus, y = ARI_ref_sc3$ARI)
  
  ARI_ref_seurat <- read.csv(here("Data", "singleTree",
                                  paste0(dataset, comp, "_hierarchical_",
                                         type, ".csv"))) %>%
    arrange(cells) %>%
    dplyr::select(starts_with("seurat")) %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    distinct() %>%
    arrange(n_clus) %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$Seurat)
  ARI_ref_seurat <- trapz(x = ARI_ref_seurat$n_clus, y = ARI_ref_seurat$ARI)
  
  ARI_ref_monocle <- read.csv(here("Data", "singleTree",
                                   paste0(dataset, comp, "_hierarchical_",
                                          type, ".csv"))) %>%
    arrange(cells) %>%
    dplyr::select(starts_with("monocle")) %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    distinct() %>%
    arrange(n_clus) %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$Monocle)
  ARI_ref_monocle <- trapz(x = ARI_ref_monocle$n_clus, y = ARI_ref_monocle$ARI)
  
  df <- data.frame("comp" = comp,
                   "method" = c("SC3", "Seurat", "Monocle"),
                   "AUARIC" = c(ARI_ref_sc3, ARI_ref_seurat, ARI_ref_monocle))
  
  return(df)
}

comp_DE_tree <- function(dataset, comp, ref) {
  return(comp_tree_ref(dataset, comp, ref, type = "DE"))
}

comp_Dist_tree <- function(dataset, comp, ref) {
  return(comp_tree_ref(dataset, comp, ref, type = "Dist"))
}

# Load param data (if we end up using it) ----
load_sc3 <- function(dataset, dune) {
  sc3_init <- dune[, "sc3.00"]
  SC3 <- read.csv(here("Data", "singleMethod",
                       paste0(dataset, "_SC3.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(SC3) <- str_remove(colnames(SC3), "X") %>%
    str_replace("\\.", "-") %>% unlist
  sc3_param <- colnames(SC3)[which.max(colSums(SC3 == sc3_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(SC3)) <= sc3_param
  SC3 <- SC3[, keep]
  colnames(SC3) <- paste("sc3", colnames(SC3), sep = "_")
  return(SC3)
}

load_monocle <- function(dataset, dune) {
  Monocle_init <- dune[, "Monocle.00"]
  Monocle <- read.csv(here("Data", "singleMethod",
                           paste0(dataset, "_Monocle.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(Monocle) <- str_remove(colnames(Monocle), "k_") %>% unlist()
  Monocle_param <- colnames(Monocle)[which.max(colSums(Monocle == Monocle_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(Monocle)) >= Monocle_param
  Monocle <- Monocle[, keep]
  colnames(Monocle) <- paste("Monocle", colnames(Monocle), sep = "_")
  return(Monocle)
}

load_seurat <- function(dataset, dune) {
  Seurat_init <- dune[, "Seurat.00"]
  Seurat <- read.csv(here("Data", "singleMethod",
                          paste0(dataset, "_Seurat.csv"))) %>%
    select(-X) %>%
    arrange(cells) %>%
    select(-cells)
  colnames(Seurat) <- str_remove(colnames(Seurat), "X") %>% unlist()
  Seurat_param <- colnames(Seurat)[which.max(colSums(Seurat == Seurat_init))]
  k_seurat <- word(Seurat_param, 3, sep = "\\.")
  keep <- str_detect(colnames(Seurat), k_seurat) %>% unlist()
  Seurat <- Seurat[, keep]
  colnames(Seurat) <- str_remove(colnames(Seurat), paste0(".", k_seurat))
  Seurat_param <- colnames(Seurat)[which.max(colSums(Seurat == Seurat_init))] %>%
    as.numeric()
  keep <- as.numeric(colnames(Seurat)) <= Seurat_param
  Seurat <- Seurat[, keep]
  colnames(Seurat) <- paste("Seurat", colnames(Seurat), sep = "_")
  return(Seurat)
}

comp_single_ref <- function(dataset, comp, ref) {
  if (comp == "") {
    merger <- readRDS(here("data", "Dune",
                           paste0(dataset, "_mergers.rds")))
    dune <- read.csv(here("data", "Dune",
                          paste0(dataset, ".csv")))
  } else {
    dune <- read.csv(here("data", "singleTree",
                          paste0(dataset, comp, "_Dune.csv")))
    merger <- readRDS(here("data", "singleTree",
                           paste0(dataset, comp, "_merger.rds")))
  }
  
  ns <- lapply(merger$currentMat, n_distinct)
  
  dune <- dune %>% arrange(cells)
  Seurat <- load_seurat(dataset, dune)
  SC3 <- load_sc3(dataset, dune)
  Monocle <- load_monocle(dataset, dune)
  
  ARI_ref_sc3 <- SC3 %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
    distinct() %>%
    arrange(n_clus) %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$sc3)
  ARI_ref_sc3 <- trapz(x = ARI_ref_sc3$n_clus, y = ARI_ref_sc3$ARI)
  
  ARI_ref_seurat <- Seurat %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
    distinct() %>%
    arrange(n_clus) %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$Seurat)
  ARI_ref_seurat <- trapz(x = ARI_ref_seurat$n_clus, y = ARI_ref_seurat$ARI)
  
  ARI_ref_monocle <- Monocle %>%
    map_df(., comp_tree_with_ref, ref = ref) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
    distinct() %>%
    arrange(n_clus)  %>%
    group_by(n_clus) %>%
    summarise(ARI = mean(ARI)) %>%
    interpolate(df = ., ns = ns$Monocle)
  ARI_ref_monocle <- trapz(x = ARI_ref_monocle$n_clus, y = ARI_ref_monocle$ARI)
  
  df <- data.frame("comp" = comp,
                   "method" = c("SC3", "Seurat", "Monocle"),
                   "AUARIC" = c(ARI_ref_sc3, ARI_ref_seurat, ARI_ref_monocle))
  
  return(df)
}

# Create the data used for the table ----
comp_all <- function(dataset, comp, ref){
  df <- bind_rows(
    "Dune" = comp_dune_ref(dataset, comp, ref),
    "Dune_NMI" = comp_dune_ref(dataset, comp, ref, metric = "_NMI"),
    "DE" = comp_DE_tree(dataset, comp, ref),
    "Dist" = comp_Dist_tree(dataset, comp, ref),
    .id = "Merge_method"
  )
  return(df)
}

comp_dataset <- function(dataset){
  allen_clusters <- read.csv(here("data", "Smart-Seq",
                                  paste0(dataset, "_cluster.membership.csv")),
                             col.names = c("cells", "cluster_id"))
  clusters <- read.csv(here("data", "Smart-Seq",
                            paste0(dataset, "_cluster.annotation.csv")),
                       header = T)
  allen_clusters <- full_join(allen_clusters, clusters) %>%
    arrange(cells) %>%
    `$`("subclass_label") %>%
    as.character()
  df <- bind_rows(
    comp_all(dataset, comp = "", ref = allen_clusters),
    comp_all(dataset, comp = "_large2", ref = allen_clusters),
    comp_all(dataset, comp = "_large3", ref = allen_clusters)
  )  %>%
    mutate(comp = case_when(comp == "" ~ "x1",
                            comp == "_large2" ~ "x2",
                            comp == "_large3" ~ "x3"))
  return(df)
}

df <- bind_rows(
    "SMARTer_cells_MOp" = comp_dataset(dataset = "SMARTer_cells_MOp"),
    "SMARTer_nuclei_MOp" = comp_dataset(dataset = "SMARTer_nuclei_MOp"),
    .id = "dataset"
)

setwd("..")
reload(inst("here"))
write.table(df, here("Figure-ARI", "data", "Brain.txt"), row.names = FALSE)
