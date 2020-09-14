# Packages ----
libs <- c("splatter")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Create data functions ----
create_simple_balanced_data <- function(nCells, nClus, DE =.1, seed = 19070) {
  sce <- splatSimulate(group.prob = rep(1, nClus) / nClus,
                       method = "groups",
                       verbose = FALSE,
                       nGenes = 10^4,
                       de.prob = DE,
                       seed = seed)
  return(sce)
}

# Balanced data functions ----
create_unbalanced_data <- function(nCells, nClus, DE, seed = 19070) {
  set.seed(seed)
  groupProb <- table(sample(seq_len(nClus), nCells, replace = TRUE)) / nCells
  groupProb <- as.vector(groupProb)
  deProb <- runif(nClus, .75 * DE, 1.25 * DE)
  sce <- splatSimulate(group.prob = groupProb,
                       method = "groups",
                       verbose = FALSE,
                       nGenes = 10^4,
                       de.prob = deProb,
                       seed = seed)
  return(sce)
}
