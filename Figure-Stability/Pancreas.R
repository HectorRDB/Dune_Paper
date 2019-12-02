# Packages ----
libs <- c("here", "dplyr", "ggplot2", "tidyr", "stringr", "readr", "cowplot",
          "mclust", "RColorBrewer", "purrr", "Dune", "pracma", "devtools",
          "pkgload")
suppressMessages(
  suppressWarnings(sapply(libs, require, character.only = TRUE))
)
rm(libs)

# Helper functions ----
setwd(here("Pancreas"))
reload(inst("here"))

# Load Dune ----
ARI_methods <- function(dataset) {
  comps <- paste0("comp", 1:3)
  names(comps) <- comps
  dunes <- map(comps, function(comp) {
    df <- read.csv(here("data", "Dune",
                        paste(dataset, comp, "Dune.csv", sep = "_")))
    df <- df %>% arrange(cells)
    return(df)
  })
  chars <- c("sc3", "Monocle", "Seurat")
  names(chars) <- chars
  df <- map_df(chars, function(method){
    ARI <- rep(0, 21)
    for (i in 1:21) {
      if (method != "sc3") {
        clusMat <- cbind(
          dunes$comp1[, str_detect(colnames(dunes$comp1), method)][, i],
          dunes$comp2[, str_detect(colnames(dunes$comp2), method)][, i],
          dunes$comp3[, str_detect(colnames(dunes$comp3), method)][, i])
        ARI[i] <- mean(ARIs(clusMat = clusMat)[upper.tri(diag(1, 3, 3))])  
      } else {
        ARI[i] <- adjustedRandIndex(
          x = dunes$comp1[, str_detect(colnames(dunes$comp1), method)][, i],
          y = dunes$comp3[, str_detect(colnames(dunes$comp3), method)][, i]
        )
      }
    }
    return(ARI)
  })
  return(df %>% mutate(steps = seq(from = 0, to = 1, by = .05)))
}

df <- bind_rows(
  "baron" = ARI_methods("baron"),
  "segerstolpe" = ARI_methods("segerstolpe"),
  .id = "dataset"
)

setwd("..")
reload(inst("here"))
write.table(df, here("Figure-Stability", "data", "Pancreas-methods.txt"), row.names = FALSE)
