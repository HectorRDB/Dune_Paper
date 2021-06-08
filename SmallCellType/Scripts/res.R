small_distinct <- function(df) {
  clus <- df %>% 
    group_by(celltype) %>%
    summarise(size = n()) %>%
    ungroup() %>%
    arrange(celltype) %>%
    mutate(number = row_number())
  sil <- silhouette(df$celltype %>% as.factor() %>% as.numeric,
                    dist = dist(as.matrix(df[, c("x", "y")]))) %>%
    summary() %>%
    `$`(., "clus.avg.widths")
  sil <- data.frame(number = seq_along(sil), sil = sil)
  clus <- full_join(clus, sil) %>% 
    filter(size < nrow(df) / 20, 
           size > 5,
           sil > 0) %>%
    dplyr::select(-number)
  return(clus)
}

check_if_there <- function(df, mat, small_clus) {
  mat %>%
    mutate(cells = rownames(.)) %>%
    full_join(df) %>%
    pivot_longer(c("sc3", "Monocle", "Seurat"), names_to = "method", values_to = "cluster") %>%
    dplyr::select(method, celltype, cluster) %>%
    group_by(method, cluster, celltype) %>%
    summarise(intersect = n()) %>%
    ungroup() %>%
    group_by(method, cluster) %>%
    mutate(A = sum(intersect)) %>%
    ungroup() %>%
    group_by(celltype, method) %>%
    mutate(B = sum(intersect)) %>%
    ungroup() %>%
    mutate(jaccard = intersect / (A + B - intersect)) %>%
    inner_join(small_clus) %>%
    group_by(celltype, method) %>%
    filter(jaccard == max(jaccard)) %>%
    dplyr::select(method, cluster, celltype, jaccard, sil) %>%
    return()
}

jaccard_imp <- function(df, merger, small_clus) {
  init <- check_if_there(df, merger$initialMat, small_clus)
  final <- check_if_there(df, merger$currentMat, small_clus)
  conversion <- clusterConversion(merger) %>%
    bind_rows(.id = "method")
  init <- inner_join(init, conversion,
                    by = c("cluster"="old", "method"="method")) %>%
    dplyr::select(-cluster) %>%
    rename("cluster" = "new")
  res <- left_join(init, final, by = c('method', "cluster","celltype", "sil")) %>%
    rename("Og_jaccard" = "jaccard.x",
           "New_jaccard" = "jaccard.y") %>%
    replace_na(list(New_jaccard = 0)) %>%
    mutate(Delta = New_jaccard - Og_jaccard)
}