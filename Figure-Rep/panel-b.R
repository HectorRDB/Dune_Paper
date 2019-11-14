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
load_Dune <- function() {
  df <- read.table(here("Brain", "data", "Replicability", "Dune_Smart",
                        "Large2", "consensus_cluster_replicability.txt"))
  df <- df %>% filter(!str_detect(clustering_method, "Consensus")) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters)
  df$level <- lapply(df$level, toRank) %>% unlist()
  df$level <- as.numeric(df$level)
  return(df)
}

load_DE_tree <- function() {
  df <- read.table(here("Brain", "data", "Replicability", "SingleTree",
                        "Large2_DE", "consensus_cluster_replicability.txt"))
  df <- df %>% filter(!str_detect(clustering_method, "Consensus")) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters,
           level = str_remove(clustering_method, "^.*\\."),
           level = if_else(nchar(level) == 1, paste0(level, "0"), level),
           level = as.numeric(level))
  return(df)
}

load_Dist_tree <- function() {
  df <- read.table(here("Brain", "data", "Replicability", "SingleTree",
                        "Large2_Dist", "consensus_cluster_replicability.txt"))
  df <- df %>% filter(!str_detect(clustering_method, "Consensus")) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters,
           level = str_remove(clustering_method, "^.*_"),
           level = as.numeric(level),
           clustering_method = str_remove(clustering_method, "_.*$"))
  return(df)
}

load_single_method <- function() {
  df <- bind_rows("Dune" = load_Dune(),
                  "Hierarchical_DE" = load_DE_tree(),
                  "Hierarchical_Dist" = load_Dist_tree(),
                  .id = "Method") %>%
    arrange(level) %>%
    mutate(clustering_method = str_remove_all(clustering_method, "\\..*$"))
  return(df)
}