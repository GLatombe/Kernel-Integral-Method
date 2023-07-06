rm(list=ls(all.names = TRUE))
gc()
graphics.off()
#quartz()


##I admit I lost track of all the functions that are called, and not all packages may not be necessary, but better safe than sorry
library("BAT") #version 2.0.1
library("FD") #version 1.1-12
library("ggplot2")
library("gridExtra")
library("labdsv")
library("hypervolume") #version 2.0.11
library("psych")
library("StatMatch")
library("TPD") #version 1.1.0
library("geometry")
library("mFD")
library("MASS")
library("rgl")
library("reshape2")
library("TreeTools")
library("StatMatch")
library("picante")
library("alphahull")
library("data.table")
library("progress")
library("abind")

##function to generate the theoretical data
source("creation_data.R") 

##Functions for the tree-based method
source("treeFT.R") ##wrapper function calling all the others
  source("drop.tip.custom.R") ##adapt drop.tip so that the id of the nodes and tips is kept the same, to be able to compute the intersection of the trees. The resulting trees are not valid, but this will be fixed by reorder.phy.
  source("keep.tip.custom.R") ##like keep.tip(), but calls drop.tip.custom() instead
  source("reorder.phy.R") ##sort the id of the nodes and tips to make a valid phylo object
  source("merge.trees.R") ##compute the intersection between i trees


##Functions for the functional-space methods
source("all_FT_function.R") ##This is the function generating the theoretical point distributions and computing the different indices of functional turnover
  source("beta_BAT.R") ##wrapper function to compute convex hull indices
  source("hypervolume.bat.R") ##wrapper function to compute hypervolume-based indices
  source("for.multi.R") ##function to compute KIM V1 and V2 indices 
  source("for.multi.nonunif.R") ##function to compute KIM V3 indices (no resampling of random points)
    source("hypervolume_gaussian.custom.R") ##compute point distributions using a modification of the hypervolume_gaussian() function from package hypervolume, to stop everything before resampling to uniform density
      source("sample_model_ellipsoid.custom.R") ##modification of sample_model_ellipsoid() to get the correct points
        source("predict_function_gaussian.custom.R") ##modification of predict_function_gaussian() needed for sample_model_ellipsoid.custom()
        source("estimate_bandwidth_silent.R") ##just to avoid some annoying message that is not an error



#creation of different square configurations
configs <- c("overlap","disjoint","nested") ##relative position of the squares
sizediffs <- c("same","different") ##size of the two squares
npointss <- c(10,40,70,100) ##number of species-ponts

is.res4 <-0
for(config in configs){
  is.res3 <-0
  for(sizediff in sizediffs){
    is.res2 <-0
    for(npoints in npointss){
      is.res1 <-0
      for(r in 1:50){ ##replicates
        if(is.res1==0){
          results1 <- all_FT_function(config=config,sizediff=sizediff,npoints=npoints)
          is.res1 <- 1
        }else{
          res <- all_FT_function(config=config,sizediff=sizediff,npoints=npoints)
          results1 <- abind(results1,res,along=4)
        }
      }
      if(is.res2==0){
        results2 <- results1
        is.res2 <- 1
      }else{
        results2 <- abind(results2,results1,along=5)
      }
    }
    if(is.res3==0){
      results3 <- results2
      is.res3 <- 1
    }else{
      results3 <- abind(results3,results2,along=6)
    }
  }
  if(is.res4==0){
    results4 <- results3
    is.res4 <- 1
  }else{
    results4 <- abind(results4,results3,along=7)
  }
}


# save(results4,file="wkspace.RData")





