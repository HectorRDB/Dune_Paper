# Packages to load ---
libs <- c("here", "tidyverse", "DailyHRB", "devtools", "pkgload", "pracma")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper functions ----
setwd(here("Brain"))
reload(inst("here"))
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
comp_dune <- function(comp = "comp1", metric = "") {
  df <- read.table(here("Data", "Replicability", "Dune_Smart", 
                        paste0(comp, metric),
                        "consensus_cluster_replicability.txt"))
  df <- df %>% filter(!str_detect(clustering_method, "Consensus")) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters)
  
  
  ARI_ref_sc3 <- df %>%
    filter(clustering_name == "sc3") %>%
    arrange(nb_clusters)
  ARI_ref_sc3 <- trapz(x = ARI_ref_sc3$nb_clusters,
                       y = ARI_ref_sc3$fraction_replicable_cells)
  
  ARI_ref_seurat <- df %>%
    filter(clustering_name == "Seurat") %>%
    arrange(nb_clusters)
  ARI_ref_seurat <- trapz(x = ARI_ref_seurat$nb_clusters,
                          y = ARI_ref_seurat$fraction_replicable_cells)
  
  ARI_ref_monocle <- df %>%
    filter(clustering_name == "Monocle") %>%
    arrange(nb_clusters)
  ARI_ref_monocle <- trapz(x = ARI_ref_monocle$nb_clusters,
                           y = ARI_ref_monocle$fraction_replicable_cells)
  
  df <- data.frame("comp" = comp,
                   "method" = c("SC3", "Seurat", "Monocle"),
                   "AUARIC" = c(ARI_ref_sc3, ARI_ref_seurat, ARI_ref_monocle))
  return(df)
}
## Load Tree data ----
comp_tree <- function(comp, type) {
  dune <- read.table(here("Data", "Replicability", "Dune_Smart", comp,
                          "consensus_cluster_replicability.txt")) %>%
    filter(!str_detect(clustering_method, "Consensus"),
           level == 100) %>%
    mutate(nb_clusters = replicable_clusters + non_replicable_clusters)
  ns <- dune$nb_clusters
  names(ns) <- dune$clustering_name
  
  sep <- if_else(type == "DE", "\\.", "_")
  df <- read.table(here("Data", "Replicability", "SingleTree",
                        paste0(comp, "_", type),
                        "consensus_cluster_replicability.txt")) %>%
    mutate(level = word(clustering_method, 2, sep = sep) %>%
             as.numeric(),
           nb_clusters = replicable_clusters + non_replicable_clusters,
           clustering_name = word(clustering_method, 1, sep = sep))
  
  ARI_ref_sc3 <- df %>% 
    filter(clustering_name == "sc3") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["sc3"])
  ARI_ref_sc3 <- trapz(x = ARI_ref_sc3$nb_clusters,
                       y = ARI_ref_sc3$fraction_replicable_cells)
  
  ARI_ref_seurat <- df %>% 
    filter(clustering_name == "Seurat") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["Seurat"])
  ARI_ref_seurat <- trapz(x = ARI_ref_seurat$nb_clusters,
                          y = ARI_ref_seurat$fraction_replicable_cells)
  
  ARI_ref_monocle <- df %>% 
    filter(clustering_name == "Monocle") %>%
    select(nb_clusters, fraction_replicable_cells) %>%
    interpolate(df = ., ns = ns["Monocle"])
  ARI_ref_monocle <- trapz(x = ARI_ref_monocle$nb_clusters,
                           y = ARI_ref_monocle$fraction_replicable_cells)
  
  df <- data.frame("comp" = comp,
                   "method" = c("SC3", "Seurat", "Monocle"),
                   "AUARIC" = c(ARI_ref_sc3, ARI_ref_seurat, ARI_ref_monocle))
  
  return(df)
}

comp_DE_tree <- function(comp) {
  return(comp_tree(comp, type = "DE"))
}

comp_Dist_tree <- function(comp) {
  return(comp_tree(comp, type = "Dist"))
}
# Create the data used for the table ----
comp_all <- function(comp){
  df <- bind_rows(
    "Dune" = comp_dune(comp),
    "Dune_NMI" = comp_dune(comp, metric = "_NMI"),
    "DE" = comp_DE_tree(comp),
    "Dist" = comp_Dist_tree(comp),
    .id = "Merge_method"
  )
  return(df)
}

df <- bind_rows(
  comp_all(comp = "Normal"),
  comp_all(comp = "large2"),
  comp_all(comp = "large3")
) %>%
  mutate(comp = case_when(comp == "Normal" ~ "x1",
                          comp == "large2" ~ "x2",
                          comp == "large3" ~ "x3"))

setwd("..")
reload(inst("here"))
write.table(df, here("Figure-Rep", "data", "Brain.txt"), row.names = FALSE)
