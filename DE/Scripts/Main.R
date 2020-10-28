suppressWarnings(library(optparse))

# Arguments for R Script ----
option_list <- list(
  make_option(c("-o", "--output"),
              action = "store", default = NA, type = "character",
              help = "Where to store the object after running"
  ),
  make_option(c("-l", "--location"),
              action = "store", default = NA, type = "character",
              help = "The location of the data"
  ),
  make_option(c("-f", "--first_name"),
              action = "store", default = NA, type = "character",
              help = "The name of the first dataset"
  ),
  make_option(c("-s", "--second_name"),
              action = "store", default = NA, type = "character",
          help = "The name of the second dataset"
  ),
  make_option(c("-m", "--merge"),
              action = "store", default = NA, type = "character",
              help = "The location of the merged files "
  ),
  make_option(c("-n", "--number_two_merge"),
              action = "store", default = NA, type = "character",
              help = "The location of the second merged files "
  )
  
)

opt <- parse_args(OptionParser(option_list = option_list))
# Packages ----
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)
library(here)
library(tidyr)
library(purrr)
library(scran)
library(Seurat)
source(here("DE", "Scripts", "Helper.R"))
# Load data ----
sce1 <- readRDS(file = paste0(opt$l, opt$f, "_filt.rds"))
sce2 <- readRDS(file = paste0(opt$l, opt$s, "_filt.rds"))

if (is.na(opt$n)) {
  comps <- list.files(opt$m) %>%
    str_subset("Dune_NMI.csv")  %>%
    str_remove("_Dune_NMI.csv") %>%
    unlist() %>%
    word(-1, sep = "_") %>%
    unlist() %>%
    unique() %>%
    paste0("_Dune_NMI.csv")
  m_locs <- rep(opt$m, length(comps))
} else {
  comps <- list.files(opt$m) %>%
    str_subset("NMI_Dune.csv")  %>%
    str_remove("_NMI_Dune.csv") %>%
    unlist() %>%
    word(-1, sep = "_") %>%
    unlist() %>%
    unique() %>%
    paste0("_NMI_Dune.csv")
  comps <- c(comps, "NMI.csv")
  m_locs <- c(rep(opt$m, length(comps) -1 ), opt$n)
}

# Make predictions
dist_mat <- all_de(sce1, sce2, opt$f, opt$s, m_locs, comps)
write.csv(dist_mat, file = opt$o, col.names = TRUE, row.names = TRUE)