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
names(comps)
df <- map(comps, function(comp) {
  print(paste0("Comp ", comp))
  no_garb <- read.table(here("Data", "Replicability", "Dune",
                             paste0("comp", comp), 
                             "consensus_cluster_replicability.txt")) %>%
    mutate(rep = 0, nb_garb = 0)
  reps <- 1:10
  names(reps) <- as.character(1:10)
  df <- map(reps, function(rep){
    sizes <- 1:3
    names(sizes) <- sizes
    sizes <- map(sizes, function(size){
      return(read.table(
        here("Data", "Replicability", "Garbage",
             paste("Bad", comp, size, rep, sep = "-"),
             "consensus_cluster_replicability.txt")))
    })
    sizes <- bind_rows(sizes, .id = "size") %>%
      rename("nb_garb" = size) %>%
      mutate(nb_garb = as.numeric(nb_garb))
    return(sizes)
  })
  df <- bind_rows(df, .id = "rep") %>%
    mutate(rep = as.numeric(rep))
  df <- bind_rows(df, no_garb)
})

df <- bind_rows(df, .id = "Comp") %>%
  filter(!str_detect(clustering_name, "Garbage"))

setwd(here(".."))
reload(inst("here"))

write.table(df, here("Figure-Stability", "data", "garbage.txt"))
