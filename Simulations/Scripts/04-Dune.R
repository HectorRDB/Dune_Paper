suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the output"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-n", "--nCores"),
              action = "store", default = 1, type = "integer",
              help = "Number of cores to use"
  ),
  make_option(c("-S", "--SeuratParam"),
              action = "store", default = NA, type = "character",
              help = "Parameter to use for Seurat"
  ),
  make_option(c("-C", "--C3"),
              action = "store", default = NA, type = "character",
              help = "SC3 parameter"
  ),
  make_option(c("-m", "--monocle"),
              action = "store", default = NA, type = "character",
              help = "Monocle parameter"
  ),
  make_option(c("-r", "--rsec"),
              action = "store", default = NA, type = "character",
              help = "Location of the Rsec object"
  ),
  make_option(c("-p", "--plot"),
              action = "store", default = "None", type = "character",
              help = "Where and if to store plot"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$l)) {
  loc <- opt$l
  cat("The selected dataset is located at ", loc, "\n")
} else {
  stop("Missing l argument")
}

if (!is.na(opt$o)) {
  output <- opt$o
} else {
  stop("Missing o argument")
  cat("The output will be stored at ", output, "\n")
}

if (!is.na(opt$C)) {
  sc3_p <- opt$C
} else {
  stop("Missing C argument")
}

if (!is.na(opt$C)) {
  seurat_p <- opt$S
} else {
  stop("Missing S argument")
}

if (!is.na(opt$m)) {
  monocle_p <- opt$m
} else {
  stop("Missing m argument")
}

library(SummarizedExperiment)
library(parallel)
library(matrixStats)
library(clusterExperiment)
library(tidyverse)
library(Dune)
library(mclust)

# Load Data ----
# Load sc3 clustering results
sc3 <- read.csv(paste0(loc, "_SC3.csv"))[, -1]
colnames(sc3) <- str_remove(colnames(sc3), "^X") %>% str_replace("\\.", "-")
if (opt$p != "None") {
  ggsave(filename = paste0(opt$p, "_SC3_ARI.png"),
         plot = plotARIs(sc3 %>% select(-cells), values = FALSE,
                         numericalLabels = TRUE))  
}
Names <- sc3$cells
sc3 <- sc3[, sc3_p] %>% as.numeric()

# Load Seurat clustering results
Seurat <- read.csv(paste0(loc, "_Seurat.csv"))[, -1]
colnames(Seurat) <- str_remove(colnames(Seurat), "^X")
if (opt$p != "None") {
  ggsave(filename = paste0(opt$p, "_Seurat_ARI.png"),
         plot = plotARIs(Seurat %>% select(-cells), values = FALSE))
}
Seurat <- Seurat[, seurat_p] %>% as.numeric()

# Load Monocle clustering results
Monocle <- read.csv(paste0(loc, "_Monocle.csv"))[, -1]
if (opt$p != "None") {
  ggsave(filename = paste0(opt$p, "_Monocle_ARI.png"),
         plot = plotARIs(Monocle %>% select(-cells), values = FALSE))
}
Monocle <- as.data.frame(Monocle)[, monocle_p] %>% as.numeric()

# Get the final clustering labels
clusMat <- data.frame("sc3" = sc3, "Monocle" = Monocle, "Seurat" = Seurat)
rownames(clusMat) <- Names  

# Do the consensus clustering with ARI ----
BPPARAM <- BiocParallel::MulticoreParam(opt$n)
print(system.time(
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE)
))
saveRDS(merger,  paste0(output, "_merger.rds"))
cat("Finished Consensus Merge\n")

# Save the matrix with all the consensus steps ----
Names <- as.character(Names)
chars <- c("sc3", "Monocle", "Seurat")
levels <- seq(from = 0, to = 1, by = .05)
stopMatrix <- lapply(levels, function(p){
  print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
  mat <- intermediateMat(merger = merger, p = p)
  suppressWarnings(rownames(mat) <- mat$cells)
  mat <- mat[Names, ]
  mat <- mat %>%
    select(-cells) %>%
    as.matrix()
  return(mat)
}) %>%
  do.call('cbind', args = .)
colnames(stopMatrix) <- lapply(levels, function(p){
  i <- as.character(round(100 * p))
  if (nchar(i) == 1) {
    i <- paste0("0", i)
  }
  return(paste(chars, i, sep = "-"))
}) %>% unlist()
print("...Full matrix")
mat <- cbind(as.character(Names), stopMatrix)
colnames(mat)[1] <- "cells"

write_csv(x = as.data.frame(mat), path = paste0(output, "_Dune.csv"))

# Do the consensus clustering with NMI ----
BPPARAM <- BiocParallel::MulticoreParam(opt$n)
print(system.time(
  merger <- Dune(clusMat = clusMat, BPPARAM = BPPARAM, parallel = TRUE,
                 metric = "NMI")
))
saveRDS(merger,  paste0(output, "_NMI_merger.rds"))
cat("Finished Consensus Merge\n")

# Save the matrix with all the consensus steps for NMI----
Names <- as.character(Names)
chars <- c("sc3", "Monocle", "Seurat")
levels <- seq(from = 0, to = 1, by = .05)
stopMatrix <- lapply(levels, function(p){
  print(paste0("...Intermediary consensus at ", round(100 * p), "%"))
  mat <- intermediateMat(merger = merger, p = p)
  suppressWarnings(rownames(mat) <- mat$cells)
  mat <- mat[Names, ]
  mat <- mat %>%
    select(-cells) %>%
    as.matrix()
  return(mat)
}) %>%
  do.call('cbind', args = .)
colnames(stopMatrix) <- lapply(levels, function(p){
  i <- as.character(round(100 * p))
  if (nchar(i) == 1) {
    i <- paste0("0", i)
  }
  return(paste(chars, i, sep = "-"))
}) %>% unlist()
print("...Full matrix")
mat <- cbind(as.character(Names), stopMatrix)
colnames(mat)[1] <- "cells"

write_csv(x = as.data.frame(mat), path = paste0(output, "_Dune_NMI.csv"))

# Do hierarchical merging with fraction of DE----
Rsec <- readRDS(opt$r)

for (clustering in c("sc3", "Monocle", "Seurat")) {
  Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
}

# Doing the merges
cutoffs <- seq(from = 0, to = .01, by = .01)
res <- list()
for (clustering in c("sc3", "Monocle", "Seurat")) {
  print(clustering)
  Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  names(cutoffs) <- paste(clustering, cutoffs, sep = "_")
  res[[clustering]] <- map_dfc(cutoffs,
                              function(i){
                                print(paste0("...", i))
                                Rsec3 <- mergeClusters(Rsec2,
                                                       mergeMethod = "adjP",
                                                       plotInfo = "adjP",
                                                       cutoff = i,
                                                       clusterLabel = "Clusters",
                                                       plot = F,
                                                       DEMethod = "limma")
                                return(Rsec3@clusterMatrix[,"Clusters"])
                              })
}

res <- do.call('cbind', res) %>% as.data.frame()
res$cells <- colnames(Rsec)
write_csv(res, path = paste0(output, "_hierarchical_DE.csv"))

# Do hierarchical merging with cutting the tree ----
Rsec <- readRDS(opt$r)

for (clustering in c("sc3", "Monocle", "Seurat")) {
  Rsec <- addClusterings(Rsec, get(clustering), clusterLabels = clustering)
}

# Doing the merges
res <- list()
for (clustering in c("sc3", "Monocle", "Seurat")) {
  print(clustering)
  n <- n_distinct(get(clustering))
  cutoffs <- 5:n
  Rsec2 <- makeDendrogram(Rsec, whichCluster = clustering)
  Tree <- as.hclust(convertToDendrogram(Rsec2))
  names(cutoffs) <- paste(clustering, n - cutoffs, sep = "_")
  res[[clustering]] <- map_dfc(cutoffs,
                              function(cutoff){
                                print(paste0("...", cutoff))
                                return(cutree(Tree, k = cutoff))
                              })
}

res <- do.call('cbind', res) %>% as.data.frame()
res$cells <- colnames(Rsec)
write_csv(res, path = paste0(output, "_hierarchical_Dist.csv"))