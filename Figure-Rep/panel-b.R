# Packages ----
library(knitr)
libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "pracma", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper ----
toRank <- function(i) {
  case_when(
    i == "Initial" ~ 0,
    i == "Final" ~ 100,
    TRUE ~ as.numeric(i)
  )
}

interpolate <- function(df, ns) {
  if (any(df$nb_cluster == ns)) {
    return(df %>% filter(nb_clusters >= ns))
  } else {
    df <- df %>% arrange(desc(nb_clusters))
    cutoff <- which(df$nb_clusters < ns)[1]
    slope <- (df$fraction_replicable_cells[cutoff] - 
                df$fraction_replicable_cells[cutoff - 1]) / 
      (df$nb_clusters[cutoff] - df$nb_clusters[cutoff - 1])
    intercept <- df$fraction_replicable_cells[cutoff]
    filt <- df %>% filter(nb_clusters >= ns) %>%
      add_row(nb_clusters = ns,
              fraction_replicable_cells = slope * 
                (ns - df$nb_clusters[cutoff]) + intercept)
    return(filt %>% arrange(nb_clusters))
  }
}

n_clus  <- function(clusMat, clus) {
  return(n_distinct(as.matrix(clusMat)[,clus]))
}

## Load Dune ----
comp_dune <- function() {
  df <- read.table(here("Brain", "Data", "Replicability", "Dune_Smart",
                        "Large2", "consensus_cluster_replicability.txt"))
  df <- df %>% filter(!str_detect(clustering_method, "Consensus")) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters) %>%
    select(nb_clusters, fraction_replicable_cells, clustering_name)
    
  return(df)
}
## Load Tree data ----
comp_tree <- function(comp, type) {
  dune <- read.table(here("Brain", "Data", "Replicability", "Dune_Smart",
                          "Large2", "consensus_cluster_replicability.txt")) %>%
    filter(!str_detect(clustering_method, "Consensus"),
           level == 100) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters)
  ns <- dune$nb_clusters
  names(ns) <- dune$clustering_name
  
  sep <- if_else(type == "DE", "\\.", "_")
  df <- read.table(here("Brain", "Data", "Replicability", "SingleTree",
                        paste0("Large2_", type),
                        "consensus_cluster_replicability.txt")) %>%
    mutate(level = word(clustering_method, 2, sep = sep) %>%
             as.numeric(),
           nb_clusters = replicable_clusters + non_replicable_clusters,
           clustering_name = word(clustering_method, 1, sep = sep))
  
  Sc3 <- df %>% 
    filter(clustering_name == "sc3") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["sc3"])
  
  Seurat <- df %>% 
    filter(clustering_name == "Seurat") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["Seurat"])
  
  Monocle <- df %>% 
    filter(clustering_name == "Monocle") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["Monocle"])
  
  df <- bind_rows("sc3" = Sc3,
                  "Monocle" = Monocle,
                  "Seurat" = Seurat,
                  .id = "clustering_name")
  
  return(df)
}

comp_DE_tree <- function(comp) {
  return(comp_tree(comp, type = "DE"))
}

comp_Dist_tree <- function(comp) {
  return(comp_tree(comp, type = "Dist"))
}
