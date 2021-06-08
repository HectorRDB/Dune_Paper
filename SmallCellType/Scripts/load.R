load_pancreas <- function(dataset = "Baron") {
  meta <- read.csv(here("Pancreas", "Data", dataset, 
                        paste0(dataset, "_meta.csv"))) %>%
    dplyr::select(X, cell_type1) %>%
    rename("cells" = "X",
           "celltype" = cell_type1)
  df <- read.csv(here("Pancreas", "Data", "tSNE", paste0(dataset, ".csv"))) %>%
    dplyr::select(-X)
  df <- full_join(df, meta)
  return(df)
}

load_res_pancreas <- function(dataset = "Baron") {
  comps <- c("comp1", "comp2", "comp3")
  df <- load_pancreas(dataset)
  small_clus <- small_distinct(df)
  mergers <- 
    lapply(comps, function(comp){
      merger <- readRDS(here("Pancreas", "Data", "Dune", 
                             paste(dataset, comp, "NMI_merger.rds", sep = "_")))
      return(merger)
    })
  names(mergers) <- comps
  res <- lapply(mergers, jaccard_imp, df = df, small_clus = small_clus)
  return(bind_rows(res, .id = "comp"))
}

load_brain <- function(dataset = "SMARTer_cells") {
  info <- read.csv(here("Brain", "data", "Smart-Seq", 
                        paste0(dataset, "_MOp_cluster.annotation.csv")))
  meta <- read.csv(here("Brain", "data", "Smart-Seq", 
                        paste0(dataset, "_MOp_cluster.membership.csv")),
                   col.names = c("cells", "cluster_id"))
  meta <- full_join(meta, info) %>%
    dplyr::select(cells, subclass_label) %>%
    rename("celltype" = "subclass_label")
    
  df <- read.csv(here("Brain", "data", "tSNE", paste0(dataset, "_MOp_tnse.csv"))) %>%
    dplyr::select(-X)
  df <- full_join(df, meta)
  return(df)
}

load_res_brain <- function(dataset = "SMARTer_cells") {
  comps <- c("", "large2", "large3")
  df <- load_brain(dataset)
  small_clus <- small_distinct(df)
  mergers <- list(
    readRDS(here("Brain", "data", "Dune", paste0(dataset, "_MOp_mergers.rds"))),
    readRDS(here("Brain", "data", "singleTree", paste0(dataset, "_MOp_large2_merger.rds"))),
    readRDS(here("Brain", "data", "singleTree", paste0(dataset, "_MOp_large3_merger.rds")))
  )
  mergers <- lapply(mergers, function(merger) {
    merger$metric <- "NMI"
    if (is.numeric(merger$merges$clusteringLabel)) {
      merger$merges$clusteringLabel <- convert(vec = merger$merges$clusteringLabel, 
                                               df = merger$initialMat)
    }
    return(merger)
  })
  names(mergers) <- c("comp1", "comp2", "comp3")
  res <- lapply(mergers, jaccard_imp, df = df, small_clus = small_clus)
  return(bind_rows(res, .id = "comp"))
}

convert <- function(vec, df) {
  if(length(vec) > 1) return(lapply(vec, convert, df = df) %>% unlist())
  return(colnames(df)[vec])
}
