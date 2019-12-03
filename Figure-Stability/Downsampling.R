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

comps <- 1:3
names(comps) <- 1:3
df <- map(comps, function(comp) {
  print(paste0("Comp ", comp))
  fractions <- c(.01, .05, 1:10 / 10) * 100
  names(fractions) <- fractions
  df <- map(fractions, function(fraction){
    return(read.table(
      here("Data", "Replicability", "Downsampling",
           paste("Fraction", comp, fraction, sep = "-"),
           "consensus_cluster_replicability.txt")))
  })
  df <- bind_rows(df, .id = "fraction") %>%
    mutate(fraction = as.numeric(fraction))
  return(df)
})

df <- bind_rows(df, .id = "Comp")

setwd(here(".."))
reload(inst("here"))

write.table(df, here("Figure-Stability", "data", "downsampling.txt"))
