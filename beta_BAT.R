beta_BAT <- function(comm_data,trait_data) {
    abund <- comm_data
  abund[abund>0] <- 1
  
  colnames(trait_data) <- c("trait1", "trait2")
  abundance <- as.matrix(abund)
  trait_data <- as.matrix(trait_data)
  
  comm <- hull.build(comm=abundance, trait = trait_data) ##build convex hull
  beta <- hull.beta(comm, func = "jaccard", comp = FALSE) ##compute jaccard family diversity indices
  return(beta)
  
}