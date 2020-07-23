# Packages ----
libs <- c("splatter", "here", "scater", "scran")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Load data ----
source(here("Simulations", "01-create_data.R"))
set.seed(101)
nCells <- 5000
sce1 <- create_simple_balanced_data(nCells = nCells, nClus = 10, seed = 197)
sce2 <- create_medium_balanced_data(nCells = nCells, nClus = 10, seed = 971)
sce3 <- create_hard_balanced_data(nCells = nCells, nClus = 10, seed = 77865)
sce4 <- create_unbalanced_data(nCells = nCells, nClus = 10, nBatches = 2,
                               DE = .2, seed = 45678)

# Run clustering
for (dataset in c(paste0("sce", 1:4))) {
  # We follow the workflow from Duo et al 2018
  # https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison
  sce <- get(dataset)
  keep_features <- rowSums(counts(sce) > 0) > 0
  sce <- sce[keep_features, ]
  df <- perCellQCMetrics(sce)
  df$libsize.drop <- isOutlier(df$total, nmads = 3,
                                         type = "lower", log = TRUE)
  df$feature.drop <- isOutlier(df$detected, nmads = 3,
                                         type = "lower", log = TRUE)
  df <- as.data.frame(df)
  sce <- sce[, !(df$libsize.drop | df$feature.drop)]
  sce <- computeSumFactors(sce, sizes = pmin(ncol(sce), seq(20, 120, 20)),
                           min.mean = 0.1)
  sce <- normalizeCounts(sce)
  
}
