#' Script to re-randomise positions of F3s.
#' 
#' Previously we grew F3 and S3 plants in November 2020.
#' We need to repeat this using the same families. I hav manually updated ID 
#' names (e.g. F3.01.005).
#' This script randomises positions within each tray
library(tidyverse)

f3s <- read_csv("03_randomisation/output/F3_randomisation_november_2023.csv")

for(t in 1:4){
  f3s[f3s$tray == t,] <- f3s %>% 
    filter(tray == t) %>% 
    mutate(coord = coord[sample(1:48, replace = F)])
}

write_csv(f3s, file = "03_randomisation/output/F3_randomisation_november_2023_rerandomised.csv")
