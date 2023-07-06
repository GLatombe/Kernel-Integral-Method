#treeFT2 <- function(comm_data,site1,site2,trait_data){
treeFT <- function(comm_data,trait_data){

  du <- comm_data[1,]
  du <- du[du > 0] 
  da <- comm_data[2,]
  da <- da[da > 0]
 
  
  
  ##compute the dendrogram from the clustering algorithm (used the same as in the original paper) and convert it to a tree of class "phylo"
  sp.names <- row.names(trait_data)
  gower.mat <- gower.dist(trait_data)
  gower.mat <- as.dist(gower.mat) ##needs to be in a dist format
  trait.tree <- as.phylo(hclust(d = gower.mat, method = "average"))
  trait.tree$tip.label <- sp.names
  
  
  
  ##prepare the 2 trees for combination i
  tips <- names(du)
  trait.tree.du <- keep.tip.custom(phy = trait.tree, tip = tips)
  tips <- names(da)
  trait.tree.da <- keep.tip.custom(phy = trait.tree, tip = tips)
  
  trait.tree.du.2 <- reorder.phy(trait.tree.du,trait.tree)
  trait.tree.da.2 <- reorder.phy(trait.tree.da,trait.tree)
  
  trait.tree.subs <- list()
  trait.tree.subs[[1]] <- trait.tree.du
  trait.tree.subs[[2]] <- trait.tree.da
  
  
  ##Compute the intersection of the 2 trees
  trait.tree.intersect <- merge.trees(trees = trait.tree.subs, tree.all = trait.tree)
  if(length(trait.tree.intersect$tip.label)==0){
    tree.pd.Will <- tree.pd.Jac <- tree.pd.Sor <- tree.pd.Sim <- tree.pd <- matrix(c(1,1),1,2)
    tree.pd.Nest <- matrix(c(0,0),1,2)
  }else{
    tree.sample <- matrix(1,1,length(trait.tree.intersect$tip.label))
    colnames(tree.sample) <- trait.tree.intersect$tip.label
    tree.pd <- pd(samp = tree.sample, tree = trait.tree.intersect) ##phylogenetic/functional diversity of the resulting tree
    
    ##compute the PD/FD of the 2 trees - needed for the Sorensen and Simpson versions
    tree.sample.indiv <- comm_data
    colnames(tree.sample.indiv) <- trait.tree$tip.label
    tree.pd.indiv <- pd(samp = tree.sample.indiv, tree = trait.tree)
    
    ##compute the PD/FD of the union of the 2 trees (i.e. we keep all the branches leading to the species present in at least one of the ord sites) - needed for the Jaccard version
    tree.sample.union <- matrix(as.numeric(colSums(comm_data)>0),1,ncol(comm_data))
    colnames(tree.sample.union) <- trait.tree$tip.label
    tree.pd.union <- pd(samp = tree.sample.union, tree = trait.tree)
    
    tree.pd.Jac <- tree.pd.Sor <- tree.pd.Sim <- tree.pd.Nest  <- tree.pd.Will <- tree.pd 
    tree.pd.Jac[,1] <- 1-tree.pd[,1]/tree.pd.union[1,1]
    tree.pd.Sor[,1] <- 1-tree.pd[,1]/mean(tree.pd.indiv[,1])
    tree.pd.Sim[,1] <- 1-tree.pd[,1]/min(tree.pd.indiv[,1])
    # tree.pd.Will[,1] <- 2*(min(tree.pd.indiv[,1])-tree.pd[,1])/tree.pd.union[1,1]
    tree.pd.Will[,1] <- 2*(min(tree.pd.indiv[,1]) - tree.pd[,1])/tree.pd.union[1,1]
    
    tree.pd.Nest[,1] <- tree.pd.Jac[,1]-tree.pd.Will[,1]
  }
  
  return(data.frame(Raw=tree.pd[1,1],Jac=tree.pd.Jac[1,1],Sor=tree.pd.Sor[1,1],Sim=tree.pd.Sim[1,1],Will=tree.pd.Will[1,1],Nest=tree.pd.Nest[1,1]))
  
}
