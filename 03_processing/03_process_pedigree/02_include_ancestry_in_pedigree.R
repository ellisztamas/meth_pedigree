#' Update the pedigree file with ancestry information
#' 
#' The pedigree file in 01_data/01_cross_design/pedigree.xlsx gives the parents
#' of each individual. This script uses that to piece together the full lineage
#' for each individual, back to the original cross.
#' 
#' Input:
#'   CSV file with a row for each individual listing the mother and father, and 
#'   the generation.
#' Output:
#'   The same CSV with two additional columns:
#'     - crossID: Accession IDs of the original parents crosses (e.g. 8249x6024).
#'         If the individual is selfed, only one parent is given.
#'     - ancestry: String giving all the individuals in the lineage separated by
#'         semicolons (e.g. 1158_2;6184_1;H1.2 x H2.1 - 2;F2.03.013) in alphabetical
#'         order.
#'         
#' Tom Ellis, 18th March 2024

library('readxl')
library('tidyverse')

source("03_processing/03_process_pedigree/01_ancestry_functions.R")

# Pedigree file giving the parents of every individual in the study
pedigree <- read_excel("01_data/01_cross_design/pedigree.xlsx", sheet = "Individuals") %>% 
  dplyr::select(plantID, mother, father, type, generation)

# Empty columns in pedigree to track cross type and the full list of ancestors
pedigree$crossID     <- character(nrow(pedigree))
pedigree$ancestors <- character(nrow(pedigree))
# Get the cross ID (e.g. '8249x6024') and the full list of ancestors for each individual
pb <- txtProgressBar(min = 1, max = nrow(pedigree), style = 3)
for(i in 1:nrow(pedigree)){
  lineage <- track_ancestors(pedigree, row_number = i)
  pedigree$crossID[i]     <- cross_ancestry(lineage)
  pedigree$ancestors[i] <- paste(unique(lineage$plantID), collapse = ";")
  setTxtProgressBar(pb, i)
}
close(pb)

# Write the updated pedigree file to disk
outfile="03_processing/03_process_pedigree/output/pedigree_with_ancestry.csv"
dir.create(dirname(outfile), showWarnings = FALSE)
write_csv(pedigree, file = outfile)

