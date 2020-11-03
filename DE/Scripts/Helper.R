# Real scripts ----
Prep <- function(sce) {
  se <- CreateSeuratObject(counts = counts(sce),
                           min.cells = 0,
                           min.features = 0,
                           project = "de")
  se <- AddMetaData(se, as.data.frame(colData(sce)))
  se <- NormalizeData(se, verbose = FALSE)
  se <- FindVariableFeatures(se, selection.method = 'vst', nfeatures = 10000,
                             verbose = FALSE)
  if ("human" %in% se@meta.data) {
    se <- ScaleData(object = se, vars.to.regress = c("nCount_RNA", "human"),
    verbose = FALSE)  
  } else {
    se <- ScaleData(object = se, vars.to.regress = "nCount_RNA",
                     verbose = FALSE)
  }
  return(as.SingleCellExperiment(se))
}

de_label <- function(label, sce) {
  markers <- scran::findMarkers(
    sce, groups = label, pval.type = "any", direction = "up"
  )
  markers <- lapply(markers, function(df){
     df %>% 
      as.data.frame() %>%
      dplyr::filter(FDR < .05) %>%
      rownames() %>%
      data.frame(markers = .) %>%
    return(.)
  })
  markers <- do.call('bind_rows', markers) %>%
    distinct()
  return(markers)
}

# All together----
compute_de_comp <- function(i, sce, comps, m_locs, f) {
  df <- read.csv(file = paste0(m_locs[i], f, "_", comps[i]),
                 stringsAsFactors = FALSE)
  print(paste0("...", names(comps)[i]))
  df <- df %>% select("cells", paste0(c("Monocle", "sc3", "Seurat"), 
                                      rep(c(".00", ".100"), each = 3)))
  sce <- sce[, df$cells]
  de_genes <- lapply(df %>% select(-cells), de_label, sce = sce)
  de_genes <- bind_rows(de_genes, .id = "label")
  return(de_genes)
}

compute_de_sce <- function(i, sces, comps, m_locs) {
  names(comps) <- word(comps, 1, sep = "_")
  f <- names(sces)[i]
  print(f)
  sce <- Prep(sces[[i]])
  de_genes <- lapply(1:3, compute_de_comp, sce = sce,
                     comps = comps, m_locs = m_locs, f = f)
  names(de_genes) <- names(comps)
  de_genes <- bind_rows(de_genes, .id = "comp")
  return(de_genes)
}

all_de <- function(sce1, sce2, f, s, m_locs, comps) {
  sces <- list(sce1, sce2)
  names(sces) <- c(f, s)
  de_genes <- lapply(1:2, compute_de_sce, sces = sces, 
                    comps = comps, m_locs = m_locs)
  names(de_genes) <- names(sces)
  de_genes <- bind_rows(de_genes, .id = "dataset")
  de_genes <- de_genes %>%
    mutate(method = paste0(dataset, "_", comp, "_", label),
           Clustering = word(label, 1, sep = "\\."),
           Level = if_else(str_detect(label, "100"), "Final", "Initial"),
           Group = paste0(dataset, "_", comp, "_", Level))
  clusterings <- unique(de_genes$Clustering)
  res <- lapply(unique(de_genes$Group), function(group){
    de_genes_group <- de_genes %>% filter(Group == group)
    dist_mat <- matrix(0, nrow = 3, ncol = 3)
    for (i in 1:2) {
      for (j in i:3) {
        genes_i <- de_genes_group %>% filter(Clustering == clusterings[i]) 
        genes_i <- genes_i$markers 
        genes_j <- de_genes_group %>% filter(Clustering == clusterings[j]) 
        genes_j <- genes_j$markers
        dist_mat[i, j] <- sum(genes_i %in% genes_j) / 
          (length(genes_i) + sum(!genes_j %in% genes_i))
      }
    }
    return(sum(dist_mat[upper.tri(dist_mat)]))
  })
  res <- data.frame(unique(de_genes$Group),
                    Concordance = unlist(res))
  return(res)
}