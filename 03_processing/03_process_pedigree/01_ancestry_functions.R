#' List ancestors through a pedigree
#' 
#' List ancestors of an individual in a pedigree file.
#' 
#' @param pedigree Dataframe giving the full pedigree of the data. This should 
#' contain a row for each individual with the header 'plantID', with columns
#' giving the mother and father of that individual, and another giving the 
#' generation (0 for parents, 1 for F1s, etc). Additional columns are ignored. 
#' All ancestors should be included in the pedigree file.
#' @param row_number Integer row number of the individual in the pedigree whose
#' ancestry is to be recovered.
#' 
#' @return Dataframe giving the rows of the original dataframe that are the 
#' ancestors of the focal individual, including the focal individual
#' 
#' @author Tom Ellis
track_ancestors <- function(pedigree, row_number){
  # Pull out the row in the pedigree for the focal individual
  this_row <- pedigree[row_number,]
  final_generation <- this_row$generation
  # Empty list to store subsequent generations
  lineage <- vector('list', length = final_generation+1)
  lineage[[ length(lineage) ]] <- this_row
  # For F1s onward, iteratively retreive the rows for the parents of each generation
  if(final_generation >= 1){
    for(i in (length(lineage)-1):1){
      parents <- unique(c(this_row$mother, this_row$father))
      parent_rows <- pedigree %>% 
        filter(plantID %in% parents)
      lineage[[ i ]] <- parent_rows
      this_row <- parent_rows
    }
  }
  # Return a single dataframe.
  lineage <- do.call(what = "rbind", lineage)
  lineage
}  

#' Identify the original parents of the cross.
#'
#' Retreive and format the genotype of the parent(s) at generation zero. This 
#' also checks that the output is one of twelve valid values.
#' 
#' @param lineage Dataframe giving pedigree information backwards to generation
#' zero for a single individual. This should be the output of track_ancestors().
#' 
#' @return String giving the geneotypes of the parents at generation zero.
#' If the offspring is the result of a cross this will be the two accession IDs
#' separated by a 'x' (e.g. '8249x6024'). If the offspring is the selfed progeny
#' of a single parent, only a single ID will be returned (e.g. '8249').
#' 
#' @author Tom Ellis
cross_ancestry <- function(lineage){
  if( ! any( lineage$generation %in% c(0,1) ) ){
    stop("This family tree appears to be incomplete. No individals are at geneartion 0 or 1.")
  }
  # Parents (plants at generation zero.)
  if(nrow(lineage) == 1){
    if(lineage$generation != 0){
      stop("This family tree appears to be incomplete. There is only one row, but it is not generation 0.")
    } else {
        cross <- lineage %>% 
          mutate(
            plantID = str_split_i(plantID, "_", 1)
          ) %>% 
          pull(plantID)
      }
  # Offspring (plants at generation 1 or more)
  } else {
    cross <- lineage %>%
      filter(generation == 1) %>% 
      mutate(
        mother = str_split_i(mother, "_", 1),
        father = str_split_i(father, "_", 1),
        cross = ifelse(
          .$mother == .$father, 
          mother,             
          paste(mother, father, sep="x")
        )
      ) %>%
      pull(cross)
  }
  # # Check that the output is one of the possible values.
  # valid_cross_names <- c(
  #   '6024', '8249', '1158', '6184',
  #   '6024x8249', '8249x6024',
  #   '6024x1158', '1158x6024',
  #   '8249x6184', '6184x8249',
  #   '1158x6184', '6184x1158'
  # )
  # if( ! cross %in% valid_cross_names ){
  #   stop(paste0("The apparent parent(s) (", cross, ") of the cross are not valid"))
  # }
  # 
  
  return(cross)
}
