#' To fill four trays of 48 pots = 2 plates of 96 = 192 pots in total
#' 
#' 7 F2 crosses, 5 families per cross, 4 individuals each = 140 plants
#' Each cross is replicated 20 times
#' Total of 35 families
#' Suggest taking one family of each cross from the same block, for five blocks.
#' 
#' 4 S2 lines, 4 families per line, 3 individuals each = 48 plants
#' Each parent replicated 12 times.
#' This is complicated because there are sometimes two S2s per parents, so do 
#' 2 families from each subline.
#' This is also complicated because for 1158 there are no S2s; use S1s or parents instead.
#' 
#' That gives 188 plants, with four blank spaces
#' 
#' Grow one big plant and plate tissue to the blank wells, to control for between
#' plate differences
#' 
#' Each tray of 48 should include 1 plant of each F2 family = 5 plants for each cross, so 35 plants
#' Plus 12 S2s; 3 trays should have one plant, and a fourth zero 
#' Plus 1 control well

library(data.table)
library('tidyverse')
set.seed(420)


f2 <- 7 * 5 * 4 
s2 <- 4 * 3 * 4
blanks <- 2
positive <- 2

f2 + s2 + blanks + positive

getPlateCoord <- function(x){
  xx <- (x - 1)%%48
  nrow <- 8
  ncol <- 6
  letters <- strsplit(intToUtf8(c(65:72)),"")[[1]]
  row <- as.integer(xx/ncol) + 1
  col <- xx%%ncol + 1 + 6 * (x > 48)
  return(paste(letters[row], col, sep = ""))
}

#tray = A-H / 1-12

nReps <- 4
dat <- fread("/groups/nordborg/projects/epiclines/002.pedigree/03_randomisation/F3_lines_to_grow.csv")
dat[, i := 1:.N, by = F2_family]
dat <- dat[rep(1:.N, nReps)][, rep := 1:.N, by = c("F2_family", "i")]
dat[, plant_id := paste(F3_family, ".00", rep, sep = "")]
dat[, tray := rep]
dat[, sequencing_plate := ceiling(tray/2)]

dat[, order := 1:.N, by = rep]
dat[, random := sample(order, .N, replace = F), by = rep]
dat[, coord := getPlateCoord(random + 48*((tray - 1)%%2))]

dat <- dat[order(F3_family, rep)]

write_csv(dat, "/groups/nordborg/projects/epiclines/002.pedigree/03_randomisation/output/F3_randomisation.csv")


dat %>% 
  group_by(tray, F3_family) %>% 
  summarise(
    n= n()
  ) %>%  View()

dat %>% 
  arrange(coord) %>%  View()
