#' Compute the distribution of the ARI between the partition and a uniformly
#' distributed partition
#'
#' @param clus the partition
#' @param rep Number of repetitions for the Monte Carlo simulations
#' @importFrom magrittr %>%
#' @importFrom dplyr n_distinct
#' @importFrom mclust adjustedRandIndex
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed.
#'  @examples 
#'  clus <- sample(x = 20, size = 300, replace = T)
#'  meanARI(clus)
meanARI <- function(clus, rep = 1000) {
  # How many clusters do we have ?
  ns <- n_distinct(clus)
  # Simulate rep ARI values
  qs <- lapply(seq_len(rep), function(i){
    rd <- sample(x = ns, size = length(clus), replace = T)
    return(adjustedRandIndex(clus, rd))
  })
  return(mean(unlist(qs)))
}
