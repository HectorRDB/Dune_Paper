# Packages ----
library(knitr)
libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "pracma", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)
# Load ref ----
# monocle comp1 smarter-nuclei
allen_clusters <- read.csv(here("Brain", "data", "Smart-Seq",
                                "SMARTer_cells_MOp_cluster.membership.csv"),
                           col.names = c("cells", "cluster_id"))
clusters <- read.csv(here("Brain", "data", "Smart-Seq",
                          "SMARTer_cells_MOp_cluster.annotation.csv"),
                     header = T)
allen_clusters <- full_join(allen_clusters, clusters) %>%
  arrange(cells) %>%
  `$`("subclass_label") %>%
  as.character()

# Helper function ----
load_monocle <- function(dataset, dune) {
  Monocle_init <- dune[, "Monocle.Initial"]
  Monocle <- read.csv(here("Brain", "Data", "singleMethod",
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
comp_merger_with_ref <- function(clusMat, clus, ref) {
  return(adjustedRandIndex(as.matrix(clusMat)[,clus], ref))
}
n_clus  <- function(clusMat, clus) {
  return(n_distinct(as.matrix(clusMat)[,clus]))
}
comp_tree_with_ref <- function(x, ref) {
  return(c("n_clus" = n_distinct(x), "ARI" = adjustedRandIndex(x, ref)))
}

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
# Load Monocle Param ---- 
dune <- read.csv(here("Brain", "data", "Dune", "SMARTer_cells_MOp.csv"))
merger <- readRDS(here("Brain", "data", "Dune",
                       "SMARTer_cells_MOp_mergers.rds"))
ns <- lapply(merger$currentMat, n_distinct)

dune <- dune %>% arrange(cells)
Monocle <- load_monocle("SMARTer_cells_MOp", dune)

Params <- Monocle %>%
  map_df(., comp_tree_with_ref, ref = allen_clusters) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
  distinct() %>%
  arrange(n_clus) %>%
  group_by(n_clus) %>%
  summarise(ARI = mean(ARI)) %>%
  interpolate(df = ., ns = ns$Monocle)

# Load Dune ----
merger$initialMat <- merger$initialMat[sort(rownames(merger$initialMat)), ]
dune <- data.frame(
  "n_clus" = functionTracking(merger, n_clus, clus = "Monocle"),
  "ARI" = functionTracking(merger, comp_merger_with_ref, clus = "Monocle",
                           ref = allen_clusters)) %>%
  arrange(n_clus) %>%
  distinct()

# Load Trees ----
DE_Tree <- read.csv(here("Brain", "Data", "singleTree",
                                 "SMARTer_cells_MOp_hierarchical_DE.csv")) %>%
  arrange(cells) %>%
  dplyr::select(starts_with("monocle")) %>%
  map_df(., comp_tree_with_ref, ref = allen_clusters) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
  distinct() %>%
  arrange(desc(n_clus)) %>%
  group_by(n_clus) %>%
  mutate(ARI = mean(ARI)) %>%
  distinct() %>%
  ungroup() %>%
  interpolate(df = ., ns = ns$Monocle)

Dist_Tree <- read.csv(here("Brain", "Data", "singleTree",
                         "SMARTer_cells_MOp_hierarchical_Dist.csv")) %>%
  arrange(cells) %>%
  dplyr::select(starts_with("monocle")) %>%
  map_df(., comp_tree_with_ref, ref = allen_clusters) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename("n_clus" = V1, "ARI" = V2) %>%
  distinct() %>%
  arrange(n_clus) %>%
  filter(n_clus >= ns$Monocle)

# data for plot ----
df <- bind_rows(
  "Params" = Params,
  "DE" = DE_Tree,
  "Dist" = Dist_Tree,
  "Dune" = dune,
  .id = "mergeMethod"
)
write.table(df, here("Figure-ARI", "data", "sc3-comp1.txt"), row.names = F)
