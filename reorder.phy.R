reorder.phy <- function(phy,phy1){
  
  
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  Nnode <- phy$Nnode
  Ntip <- length(phy$tip.label)
  NEWROOT <- length(phy1$tip.label) + 1
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[which(newNb > 0) - 
                                       Ntip]
  phy
}