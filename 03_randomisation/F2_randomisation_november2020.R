#' Tom Ellis, October 2020
#' 
#' Script to perform randomisation of F2s from 7 crosses, plus selfed offspring
#' of the parents (S2).
#' 
#' We consider a block as two trays of 6x8 pots, so that a pair of trays
#' corresponds to a 12x8 sequencing plate.
#' In those 96 pots we can have 13 F2s for each of seven crosses; 7x13=91
#' Doing 15 blocks (30 trays) gives 13x15=195 F2s per cross.
#' 
#' That leaves 75 spots for S2s (5 per block), but we have six S2 families to 
#' 

library('tidyverse')

# List of the F2 families used
f2s <- c(
  'F2.05',
  'F2.32',
  'F2.25',
  'F2.27',
  #F2.33 no seeds for this cross,
  'F2.16',
  'F2.35',
  'F2.03',
  'S2' # For S2s insert a placeholder we will change later
)
f2s <- rep(f2s, times=c(13,13,13,13,13,13,13,5))

output <- vector('list', 15)

# For each block, set up row and column labels and randomise genotypes
for(block in 1:15){
  ix <- sample(1:96, replace = F)
  output[[block]] <-  tibble(
    block  = block,
    tray   = c(rep(2*block-1, 48), rep(2*block, 48)),
    row    = rep(LETTERS[8:1], 12), # reversed to match a sequencing plate
    column = rep(1:12, each=8),
    family = as.character(f2s[ix])
  )
}
# concatenate list to a data frame
output <- do.call('rbind', output)
# Add a variable to recapitulate the current order
output$original_order <- 1:1440

# There are 75 spots for S2s, but 6 S2 families
# Take 13 of three families and 12 of the other three.
s2s <- c('S2.12', 'S2.14',  'S2.01',  'S2.06', 'S2.15', 'S2.13')
s2s <- rep(s2s, 13)[sample(1:75, 75, replace = F)]
output$family[output$family == "S2"] <- s2s

# Add Unique IDs for each plant
output <- split(output, output$family)
output <- lapply(output, function(x){
  three_digit_index <-sprintf("%03d", 1:nrow(x))
  cbind(x,
        id = paste(x$family, three_digit_index, sep = ".")
  )
})
output <- do.call("rbind", output)

# Add colours for the labels
output <- output %>% 
  mutate(label_colour = "white") %>% 
  mutate(label_colour = ifelse(family == "F2.16", "green", label_colour)) %>% 
  mutate(label_colour = ifelse(family %in% c("F2.32", "F2.05"), "blue", label_colour)) %>% 
  mutate(label_colour = ifelse(family %in% c("F2.25", "F2.27"), "pink", label_colour)) %>% 
  mutate(label_colour = ifelse(family %in% c("F2.35", "F2.03"), "orange", label_colour))

write_csv(output, "data_raw/F2_randomisation_november2020.csv")
