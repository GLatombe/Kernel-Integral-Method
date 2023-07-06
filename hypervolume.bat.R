hypervolume.bat <- function(comm_data,trait_data, site1 = "Site1", site2 = "Site2",kde.bandwidth.type = "default"){
  
  first <- comm_data[site1, comm_data[site1,] > 0]
  second <- comm_data[site2, comm_data[site2,] > 0]
  trait1 <- trait_data[names(first),]
  trait2 <- trait_data[names(second),]
  
  if(kde.bandwidth.type=="uniform"){
    trait12 <- rbind(trait1,trait2) ##put all species together
    kde.bandwidth1 <- kde.bandwidth2 <- estimate_bandwidth(trait12)/2 ##compute common bandwidth
    hyp <- kernel.build(comm = comm_data, trait = trait_data, distance = "gower", method = "gaussian", abund = TRUE,kde.bandwidth=kde.bandwidth1) ##compute hypervolume using the bandwidth
  }else{
    hyp <- kernel.build(comm = comm_data, trait = trait_data, distance = "gower", method = "gaussian", abund = TRUE)  ##compute hypervolume - default
  }
  
  beta <- kernel.beta(comm= hyp) ##compute beta diversity indices
  
  final <- list()
  final$hypervol <- hyp
  final$indices <- data.frame(Site1 = combn(colnames(as.matrix(beta$Brich)),2)[1,],
                              Site2 = combn(colnames(as.matrix(beta$Brich)),2)[2,],
                              Btotal = as.numeric(beta$Btotal),
                              Bturn = as.numeric(beta$Brepl), 
                              Bnest = as.numeric(beta$Brich))
  
  return(final)
}