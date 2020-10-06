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
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
# Packages ----
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)
library(SingleR)
library(Seurat)
library(ggplot2)
library(here)
library(tidyr)
library(purrr)
source(here("Mapping", "Scripts", "Helper.R"))
# Load data ----
sce1 <- readRDS(file = paste0(opt$l, opt$f, "_filt.rds"))
sce2 <- readRDS(file = paste0(opt$l, opt$s, "_filt.rds"))

comps <- list.files(opt$m) %>%
  str_subset("Dune_NMI.csv")  %>%
  str_remove("_Dune_NMI.csv") %>%
  unlist() %>%
  word(-1, sep = "_") %>%
  unlist() %>%
  unique()

# Make predictions
probas <- all_comps(sce1, sce2, opt$f, opt$s, opt$m, comps)
write.csv(probas, file = opt$o, col.names = TRUE, row.names = FALSE)