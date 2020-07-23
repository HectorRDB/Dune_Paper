# Packages ----
libs <- c("splatter")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Create data functions ----
.create_balanced_data <- function(nCells, nClus, nBatches, DE, seed = 19070) {
  sce <- splatSimulate(group.prob = rep(1, nClus) / nClus,
                       method = "groups",
                       verbose = FALSE,
                       nGenes = 10^4,
                       batchCells = rep(nCells / nBatches, nBatches),
                       de.prob = DE,
                       seed = seed)
  return(sce)
}

create_simple_balanced_data <- function(nCells, nClus, seed = 19070) {
  return(.create_balanced_data(nCells = nCells, nClus = nClus, nBatches = 1,
                               DE = .25, seed = seed))
}


create_hard_balanced_data <- function(nCells, nClus, seed = 19070) {
  return(.create_balanced_data(nCells = nCells, nClus = nClus, nBatches = 4,
                               DE = .1, seed = seed))
}

# Balanced data functions ----
create_unbalanced_data <- function(nCells, nClus, nBatches, DE, seed = 19070) {
  set.seed(seed)
  groupProb <- table(sample(seq_len(nClus), nCells, replace = TRUE)) / nCells
  groupProb <- as.vector(groupProb)
  deProb <- rnorm(nClus, DE, sd = .1)
  deProb[deProb <= 0] <- DE
  sce <- splatSimulate(group.prob = groupProb,
                       method = "groups",
                       verbose = FALSE,
                       nGenes = 10^4,
                       batchCells = rep(nCells / nBatches, nBatches),
                       de.prob = deProb,
                       seed = seed)
  return(sce)
}
