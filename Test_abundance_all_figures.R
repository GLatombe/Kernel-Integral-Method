rm(list=ls(all.names = TRUE))
gc()
graphics.off()
quartz()


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





configs <- c("overlap") ##relative position of the squares
sizediffs <- c("different") ##size of the two squares
#npointss <- c(10) ##number of species-ponts
npoints <- 2
bandwith.fac <- 1

Site1AX <- c(1,1,5,5)
Site1AY <- c(4,8,4,8)

Site3AX <- c(4,4,6,6)
Site3AY <- c(3,5,3,5)

##different sizes
par(mfrow=c(2,2))
comms_similar <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Similar",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
comms_uniform <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Uniform",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)
comms_diff <- creation_data(comm1 =  data.frame(X = Site1AX, Y = Site1AY),comm2 = data.frame(X = Site3AX, Y= Site3AY),mode = "Different",points_number = npoints,tol.abs = 1,tol.rel = NULL, plot=TRUE)



##If you play with the comments, you can test different combinations: same abundance and trait variability for all species, or make on or both vary for each species. In practice, you will want to play with the values to get something that works for your data

abundance <- 1000 ##same abundance for all species (1000 random points generated)
trait.var <- 1 ##same trait variability for all species and for the two traits (that will make circles with a radius of 1)

#Kernel-based function fully adjusted: same bandwidth for the 2 communities and no resampling of random points to uniform distribution (1000 points for each species, leading to changes in density in the functional space)
kernel.res.nonunif <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
kernel.res.diff.nonunif <- for.multi.nonunif(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.unif.nonunif <- for.multi.nonunif(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.simil.nonunif <- for.multi.nonunif(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)

kernel.res.nonunif[1,2:4] <- kernel.res.diff.nonunif$indices
kernel.res.nonunif[2,2:4] <- kernel.res.unif.nonunif$indices
kernel.res.nonunif[3,2:4] <- kernel.res.simil.nonunif$indices

png(filename="All_figures.png",width=600,height=1000)
par(mfrow=c(5,3))
plot(comms_uniform[[2]],col=c(rep("black",npoints+4),rep("red",npoints+4)),xlim=c(0,7),ylim=c(2,9),xlab="Trait 1",ylab="Trait 2",main="Initial species points")
plot.new()
plot.new()
plot(kernel.res.unif.nonunif$points1,pch=20,cex=0.1,xlim=c(0,7),ylim=c(2,9),xlab="Trait 1",ylab="Trait 2",main="Fixed abundance and radii")
points(kernel.res.unif.nonunif$points2,pch=20,cex=0.1,col="red")
image(x=kernel.res.unif.nonunif$kernel1$eval.points[[1]],y=kernel.res.unif.nonunif$kernel1$eval.points[[2]],z=kernel.res.unif.nonunif$kernel1$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 1",xlab="Trait 1",ylab="Trait 2")
image(x=kernel.res.unif.nonunif$kernel2$eval.points[[1]],y=kernel.res.unif.nonunif$kernel2$eval.points[[2]],z=kernel.res.unif.nonunif$kernel2$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 2",xlab="Trait 1",ylab="Trait 2")


#########

abundance <- matrix(ceiling(1000*runif(npoints+4)),2,npoints+4) ## different abundance for all species (within [1,1000])
trait.var <- 1 ##same trait variability for all species and for the two traits (that will make circles with a radius of 1)

#Kernel-based function fully adjusted: same bandwidth for the 2 communities and no resampling of random points to uniform distribution (1000 points for each species, leading to changes in density in the functional space)
kernel.res.nonunif <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
kernel.res.diff.nonunif <- for.multi.nonunif(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.unif.nonunif <- for.multi.nonunif(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.simil.nonunif <- for.multi.nonunif(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)

kernel.res.nonunif[1,2:4] <- kernel.res.diff.nonunif$indices
kernel.res.nonunif[2,2:4] <- kernel.res.unif.nonunif$indices
kernel.res.nonunif[3,2:4] <- kernel.res.simil.nonunif$indices

plot(kernel.res.unif.nonunif$points1,pch=20,cex=0.1,xlim=c(0,7),ylim=c(2,9),xlab="Trait 1",ylab="Trait 2",main="Varying abundance and fixed radii")
points(kernel.res.unif.nonunif$points2,pch=20,cex=0.1,col="red")
image(x=kernel.res.unif.nonunif$kernel1$eval.points[[1]],y=kernel.res.unif.nonunif$kernel1$eval.points[[2]],z=kernel.res.unif.nonunif$kernel1$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 1",xlab="Trait 1",ylab="Trait 2")
image(x=kernel.res.unif.nonunif$kernel2$eval.points[[1]],y=kernel.res.unif.nonunif$kernel2$eval.points[[2]],z=kernel.res.unif.nonunif$kernel2$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 2",xlab="Trait 1",ylab="Trait 2")

#########

abundance <- 1000 ##same abundance for all species (1000 random points generated)
trait.var <- list() ## different trait variability for all species and for each trait (within [0,1])
trait.var[[1]] <- matrix(runif((npoints+4)*2),npoints+4,2)
trait.var[[2]] <- matrix(runif((npoints+4)*2),npoints+4,2)

#Kernel-based function fully adjusted: same bandwidth for the 2 communities and no resampling of random points to uniform distribution (1000 points for each species, leading to changes in density in the functional space)
kernel.res.nonunif <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
kernel.res.diff.nonunif <- for.multi.nonunif(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.unif.nonunif <- for.multi.nonunif(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.simil.nonunif <- for.multi.nonunif(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)

kernel.res.nonunif[1,2:4] <- kernel.res.diff.nonunif$indices
kernel.res.nonunif[2,2:4] <- kernel.res.unif.nonunif$indices
kernel.res.nonunif[3,2:4] <- kernel.res.simil.nonunif$indices

plot(kernel.res.unif.nonunif$points1,pch=20,cex=0.1,xlim=c(0,7),ylim=c(2,9),xlab="Trait 1",ylab="Trait 2",main="Fixed abundance and varying radii")
points(kernel.res.unif.nonunif$points2,pch=20,cex=0.1,col="red")
image(x=kernel.res.unif.nonunif$kernel1$eval.points[[1]],y=kernel.res.unif.nonunif$kernel1$eval.points[[2]],z=kernel.res.unif.nonunif$kernel1$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 1",xlab="Trait 1",ylab="Trait 2")
image(x=kernel.res.unif.nonunif$kernel2$eval.points[[1]],y=kernel.res.unif.nonunif$kernel2$eval.points[[2]],z=kernel.res.unif.nonunif$kernel2$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 2",xlab="Trait 1",ylab="Trait 2")

##########


#abundance <- 1000 ##same abundance for all species (1000 random points generated)
abundance <- matrix(ceiling(1000*runif(npoints+4)),2,npoints+4) ## different abundance for all species (within [1,1000])
# trait.var <- 1 ##same trait variability for all species and for the two traits (that will make circles with a radius of 1)
trait.var <- list() ## different trait variability for all species and for each trait (within [0,1])
trait.var[[1]] <- matrix(runif((npoints+4)*2),npoints+4,2)
trait.var[[2]] <- matrix(runif((npoints+4)*2),npoints+4,2)

#Kernel-based function fully adjusted: same bandwidth for the 2 communities and no resampling of random points to uniform distribution (1000 points for each species, leading to changes in density in the functional space)
kernel.res.nonunif <- data.frame(Config=c("Joint_different","Joint_uniform","Joint_similar"),Jac=numeric(3),Turn=numeric(3),Nest=numeric(3))
kernel.res.diff.nonunif <- for.multi.nonunif(comms_diff[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_diff[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.unif.nonunif <- for.multi.nonunif(comms_uniform[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_uniform[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)
kernel.res.simil.nonunif <- for.multi.nonunif(comms_similar[[1]], site1 = "Site1", site2 = "Site2", trait_data = comms_similar[[2]],res = 0.05, bandwith.fac = bandwith.fac,chunk.size=abundance,scales=trait.var)

kernel.res.nonunif[1,2:4] <- kernel.res.diff.nonunif$indices
kernel.res.nonunif[2,2:4] <- kernel.res.unif.nonunif$indices
kernel.res.nonunif[3,2:4] <- kernel.res.simil.nonunif$indices

plot(kernel.res.unif.nonunif$points1,pch=20,cex=0.1,xlim=c(0,7),ylim=c(2,9),xlab="Trait 1",ylab="Trait 2",main="Varying abundance and varying radii")
points(kernel.res.unif.nonunif$points2,pch=20,cex=0.1,col="red")
image(x=kernel.res.unif.nonunif$kernel1$eval.points[[1]],y=kernel.res.unif.nonunif$kernel1$eval.points[[2]],z=kernel.res.unif.nonunif$kernel1$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 1",xlab="Trait 1",ylab="Trait 2")
image(x=kernel.res.unif.nonunif$kernel2$eval.points[[1]],y=kernel.res.unif.nonunif$kernel2$eval.points[[2]],z=kernel.res.unif.nonunif$kernel2$estimate,col = hcl.colors(24, "YlOrRd", rev = TRUE),main="KDE community 2",xlab="Trait 1",ylab="Trait 2")

dev.off()









