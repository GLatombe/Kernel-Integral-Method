all_FT_function <- function(config,sizediff="same",npoints){
  
  ##Generate the different configurations, as in the paper
  if(config=="overlap"){
    ##overlapping
    Site1AX <- c(1,1,5,5)
    Site1AY <- c(4,8,4,8)
    
    if(sizediff=="same"){
      Site2AX <- c(4,4,8,8)
      Site2AY <- c(1,5,1,5)
      
      ##same size
      par(mfrow=c(2,2))
      comms_similar <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site2AX, Y= Site2AY),mode = "Similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_uniform <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site2AX, Y= Site2AY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_diff <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site2AX, Y= Site2AY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
    }else{
      Site3AX <- c(4,4,6,6)
      Site3AY <- c(3,5,3,5)
      
      ##different sizes
      par(mfrow=c(2,2))
      comms_similar <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_uniform <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_diff <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
    }
  }else if(config=="disjoint"){
    
    ##disjoint
    Site1BX <- c(1,1,5,5)
    Site1BY <- c(5,9,5,9)
    
    if(sizediff=="same"){
      
      Site2BX <- c(6,6,10,10)
      Site2BY <- c(0,4,0,4)
      
      ##same size
      par(mfrow=c(2,2))
      comms_similar <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Disjoint similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_uniform <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_diff <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
    }else{
      Site3BX <- c(6,6,8,8)
      Site3BY <- c(2,4,2,4)
    ##different sizes
      par(mfrow=c(2,2))
      comms_similar <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Disjoint similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_uniform <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_diff <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
    }
  
  }else if(config=="nested"){
    ##similar
    Site1BX <- c(1,1,5,5)
    Site1BY <- c(5,9,5,9)
    
    if(sizediff=="same"){
      Site2BX <- c(1,1,5,5)
      Site2BY <- c(5,9,5,9)
      
      ##same size
      par(mfrow=c(2,2))
      comms_similar <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_uniform <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      comms_diff <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site2BX, Y= Site2BY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
      }else{
        Site3BX <- c(2.5,2.5,4.5,4.5)
        Site3BY <- c(5.5,7.5,5.5,7.5)
        
        ##different sizes
        par(mfrow=c(2,2))
        comms_similar <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
        comms_uniform <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
        comms_diff <- creation_data(comm1 =  data.frame(X = Site1BX, Y = Site1BY),comm2 = data.frame(X = Site3BX, Y= Site3BY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
    }
  }
  
  
  ###########################
  ##compute dissimilarities##
  ###########################
  
  ##CONVEX HULL
  #convex hull results -- Uses the hull.beta function from the BAT packages
  cvxhull.res <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  beta_convex <- beta_BAT(comms_diff[[1]],comms_diff[[2]])
  cvxhull.res[1,2:4] <- c(beta_convex$Btotal[1],beta_convex$Brepl[1],beta_convex$Brich[1])
  beta_convex <- beta_BAT(comms_uniform[[1]],comms_uniform[[2]])
  cvxhull.res[2,2:4] <- c(beta_convex$Btotal[1],beta_convex$Brepl[1],beta_convex$Brich[1])
  beta_convex <- beta_BAT(comms_similar[[1]],comms_similar[[2]])
  cvxhull.res[3,2:4] <- c(beta_convex$Btotal[1],beta_convex$Brepl[1],beta_convex$Brich[1])
  
  
  
  
  
  
  ##HYPERVOLUME FROM BAT
  #hypervolume default -- This is the hypervolume computation from the BAT package
  hypervol.res.default <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  hypervol.res.all.default <- hypervolume.bat(comms_diff[[1]],comms_diff[[2]]) 
  hypervol.res.default[1,2:4] <- hypervol.res.all.default$indices[3:5]
  hypervol.res.all.diff.default <- hypervol.res.all.default
  hypervol.res.all.default <- hypervolume.bat(comms_uniform[[1]],comms_uniform[[2]]) 
  hypervol.res.default[2,2:4] <- hypervol.res.all.default$indices[3:5]
  hypervol.res.all.unif.default <- hypervol.res.all.default
  hypervol.res.all.default <- hypervolume.bat(comms_similar[[1]],comms_similar[[2]]) 
  hypervol.res.default[3,2:4] <- hypervol.res.all.default$indices[3:5]
  hypervol.res.all.simil.default <- hypervol.res.all.default
  
  
  #hypervolume adjusted bandwidth -- Here we force the KDEs to have the same bandwidth for the two communities
  hypervol.res <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  hypervol.res.all <- hypervolume.bat(comms_diff[[1]],comms_diff[[2]],kde.bandwidth.type="uniform")
  hypervol.res[1,2:4] <- hypervol.res.all$indices[3:5]
  hypervol.res.all.diff <- hypervol.res.all
  hypervol.res.all <- hypervolume.bat(comms_uniform[[1]],comms_uniform[[2]],kde.bandwidth.type="uniform")
  hypervol.res[2,2:4] <- hypervol.res.all$indices[3:5]
  hypervol.res.all.unif <- hypervol.res.all
  hypervol.res.all <- hypervolume.bat(comms_similar[[1]],comms_similar[[2]],kde.bandwidth.type="uniform")
  hypervol.res[3,2:4] <- hypervol.res.all$indices[3:5]
  hypervol.res.all.simil <- hypervol.res.all
  
  
  
  
  
  
  ##NEW KERNEL-BASED METHODS
  #Kernel-based function using the same parameters as the original hypervolume method from BAT: each community KDE has a different bandwidth and resampling of random points to uniform distribution
  bandwith.fac <- 1
  kernel.res <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  kernel.res.diff <- for.multi(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac)
  kernel.res.unif <- for.multi(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac)
  kernel.res.simil <- for.multi(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac)
  
  kernel.res[1,2:4] <- kernel.res.diff$indices
  kernel.res[2,2:4] <- kernel.res.unif$indices
  kernel.res[3,2:4] <- kernel.res.simil$indices
  
  
  #Kernel-based function adjusted for bandwidth: same bandwidth for the 2 communities and resampling of random points to uniform distribution
  kernel.res.adjusted <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  kernel.res.diff.adjusted <- for.multi(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac,kde.bandwidth.type="uniform")
  kernel.res.unif.adjusted <- for.multi(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac,kde.bandwidth.type="uniform")
  kernel.res.simil.adjusted <- for.multi(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac,kde.bandwidth.type="uniform")
  
  kernel.res.adjusted[1,2:4] <- kernel.res.diff.adjusted$indices
  kernel.res.adjusted[2,2:4] <- kernel.res.unif.adjusted$indices
  kernel.res.adjusted[3,2:4] <- kernel.res.simil.adjusted$indices
  
  
  #Kernel-based function fully adjusted: same bandwidth for the 2 communities and no resampling of random points to uniform distribution (1000 points for each species, leading to changes in density in the functional space)
  kernel.res.nonunif <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  kernel.res.diff.nonunif <- for.multi.nonunif(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac, scales=1)
  kernel.res.unif.nonunif <- for.multi.nonunif(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac, scales=1)
  kernel.res.simil.nonunif <- for.multi.nonunif(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac, scales=1)
  
  kernel.res.nonunif[1,2:4] <- kernel.res.diff.nonunif$indices
  kernel.res.nonunif[2,2:4] <- kernel.res.unif.nonunif$indices
  kernel.res.nonunif[3,2:4] <- kernel.res.simil.nonunif$indices
  
  
  
  
  
  
  ##TREE-BASED FT
  Tree.FT.diff <- treeFT(comm_data=comms_diff[[1]], trait_data = comms_diff[[2]])
  Tree.FT.unif <- treeFT(comm_data=comms_uniform[[1]], trait_data = comms_uniform[[2]])
  Tree.FT.simil <- treeFT(comm_data=comms_similar[[1]], trait_data = comms_similar[[2]])
  
  Tree.FT<- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
  Tree.FT[1,2:4] <- Tree.FT.diff[c(2,5,6)]
  Tree.FT[2,2:4] <- Tree.FT.unif[c(2,5,6)]
  Tree.FT[3,2:4] <- Tree.FT.simil[c(2,5,6)]
  
  
  
  
  
  
  ##RATIO OF THE WILLIAMS AND JACCARD INDICES
  cvxhull.res$Nest.ratio <- cvxhull.res$Nest/cvxhull.res$Jac
  hypervol.res.default$Nest.ratio <- hypervol.res.default$Nest/hypervol.res.default$Jac
  hypervol.res$Nest.ratio <- hypervol.res$Nest/hypervol.res$Jac
  kernel.res$Nest.ratio <- kernel.res$Nest/kernel.res$Jac
  kernel.res.adjusted$Nest.ratio <- kernel.res.adjusted$Nest/kernel.res.adjusted$Jac
  kernel.res.nonunif$Nest.ratio <- kernel.res.nonunif$Nest/kernel.res.nonunif$Jac
  Tree.FT$Nest.ratio <- Tree.FT$Nest/Tree.FT$Jac
  
  
  ##Reformat
  cvxhull.res.mat <- as.matrix(cvxhull.res[,2:4])
  hypervol.res.default.mat <- as.matrix(hypervol.res.default[,2:4])
  hypervol.res.mat <- as.matrix(hypervol.res[,2:4])
  kernel.res.mat <- as.matrix(kernel.res[,2:4])
  kernel.res.adjusted.mat <- as.matrix(kernel.res.adjusted[,2:4])
  kernel.res.nonunif.mat <- as.matrix(kernel.res.nonunif[,2:4])
  Tree.FT.mat <- as.matrix(Tree.FT[,2:4])
  
  row.names(cvxhull.res) <- cvxhull.res[,1]
  row.names(hypervol.res.default) <- cvxhull.res[,1]
  row.names(hypervol.res) <- cvxhull.res[,1]
  row.names(kernel.res) <- cvxhull.res[,1]
  row.names(kernel.res.adjusted) <- cvxhull.res[,1]
  row.names(kernel.res.nonunif) <- cvxhull.res[,1]
  row.names(Tree.FT) <- cvxhull.res[,1]
  
  #results
  r <- abind(cvxhull.res.mat,hypervol.res.default.mat,hypervol.res.mat,kernel.res.mat,kernel.res.adjusted.mat,kernel.res.nonunif.mat,Tree.FT.mat,along=3)
  return(r)
  
  
}