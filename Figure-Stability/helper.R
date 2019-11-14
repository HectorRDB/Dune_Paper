libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "DailyHRB", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Load rep data -----
setwd(here("Pancreas"))
reload(inst("here"))
no_garb <- read.table(here("Data", "Replicability", "Dune", "Comp1",
                           "consensus_cluster_replicability.txt"))
for (i in 1:3) {
  df <- read.table(
    here("Data", "Replicability", "Garbage", paste0("bad_", i),
         "consensus_cluster_replicability.txt")
    )
  assign(paste0("garb_", i), df)
}
df <- bind_rows(
  "0" = no_garb,
  "1" = garb_1,
  "2" = garb_2,
  "3" = garb_3,
  .id = "nb_garb"
) %>%
  filter(!str_detect(clustering_name, "garbage"))

rm(no_garb, garb_1, garb_2, garb_3, i)

setwd(here(".."))
reload(inst("here"))

write.table(df, here("Figure-Stability", "data", "garbage.txt"))

# Load ari imp data -----
setwd(here("Pancreas"))
reload(inst("here"))
ari_imp <- function(clusMat) {
  return(mean(adjustedRandIndex(clusMat$sc3, clusMat$Monocle),
              adjustedRandIndex(clusMat$sc3, clusMat$Seurat),
              adjustedRandIndex(clusMat$Seurat, clusMat$Monocle)))
}

for (dataset in c("baron", "segerstolpe")) {
  print(dataset)
  merger <- readRDS(here("Data", "Dune", paste0(dataset, "_comp1_merger.rds")))
  no_garb <- data.frame(ARI = ARIImp(merger),
                        x = 1:length(ARIImp(merger)) / length(ARIImp(merger)))
  
  lapply(1:3, function(i){
    print(paste0("...", i))
    merger <- readRDS(here("data", "Garbage",
                           paste0(dataset, "_", i, "-bad_merger.rds")))
    ari_with <- functionTracking(merger, ari_imp, p = 1)
    df <- data.frame(ARI = ari_with,
                     x = 1:length(ari_with) / length(ari_with))
    assign(paste0("garb_", i), df, envir = parent.env(environment()))
  })
  df <- bind_rows(
    "0" = no_garb,
    "1" = garb_1,
    "2" = garb_2,
    "3" = garb_3,
    .id = "nb_garb"
  )
  assign(dataset, df)
}

df <- bind_rows(
  "segerstolpe" = segerstolpe,
  "baron" = baron,
  .id = "dataset"
)

setwd(here(".."))
reload(inst("here"))

write.table(df, here("Figure-Stability", "data", "garbage_ari.txt"))
