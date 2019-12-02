# Packages ----
libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "pracma", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper functions ----
setwd(here("Brain"))
reload(inst("here"))

# Load Dune ----
ARI_methods <- function(dataset) {
  comps <- c("Normal", "Large2", "Large3")
  names(comps) <- comps
  dunes <- map(comps, function(comp) {
    if (comp == "Normal") {
      df <- read.csv(here("data", "Dune", paste0(dataset, ".csv")))
    } else {
      df <- read.csv(here("data", "singleTree",
                          paste(dataset, comp, "Dune.csv", sep = "_")))
    }
    df <- df %>% arrange(cells)
    return(df)
  })
  
  chars <- c("sc3", "Monocle", "Seurat")
  names(chars) <- chars
  df <- map_df(chars, function(method){
    ARI <- rep(0, 21)
    for (i in 1:21) {
      clusMat <- cbind(
        dunes$Normal[, str_detect(colnames(dunes$Normal), method)][, i],
        dunes$Large2[, str_detect(colnames(dunes$Large2), method)][, i],
        dunes$Large3[, str_detect(colnames(dunes$Large3), method)][, i])
      ARI[i] <- mean(ARIs(clusMat = clusMat)[upper.tri(diag(1, 3, 3))])
    }
    return(ARI)
  })
  return(df %>% mutate(steps = seq(from = 0, to = 1, by = .05)))
}

df <- bind_rows(
  "SMARTer_cells" = ARI_methods("SMARTer_cells_MOp"),
  "SMARTer_nuclei" = ARI_methods("SMARTer_nuclei_MOp"),
  .id = "dataset"
)

setwd("..")
reload(inst("here"))
write.table(df, here("Figure-Stability", "data", "Pancreas-methods.txt"), row.names = FALSE)
