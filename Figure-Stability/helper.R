libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "DailyHRB", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Load data -----
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