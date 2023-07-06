merge.trees <- function(trees,tree.all){
  
  tree.m <- trees[[1]]
  for(i in 2:length(trees)){
    tree.m$edge <- as.matrix(merge(tree.m$edge,trees[[i]]$edge,sort=FALSE))
    if(!is.null(tree.m$edge.length)){
      tree.m$edge.length <- trees[[i]]$edge.length[which(trees[[i]]$edge[,1] %in% tree.m$edge[,1] & trees[[i]]$edge[,2] %in% tree.m$edge[,2])]
    }
  }
  tree.m$tip.label <- as.character(tree.m$edge[which(!(tree.m$edge[,2] %in% tree.m$edge[,1])),2])
  tree.m$Nnode <- length(unique(c(tree.m$edge))) - length(which(!(tree.m$edge[,2] %in% tree.m$edge[,1])))
  
  tree.m <- reorder.phy(tree.m, tree.all)
  
  tree.m
}