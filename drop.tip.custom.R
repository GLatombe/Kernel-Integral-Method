drop.tip.custom <- function (phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, 
          rooted = is.rooted(phy), collapse.singles = TRUE, interactive = FALSE) 
{
  
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)
  
    if (is.character(tip)) 
      tip <- which(phy$tip.label %in% tip)
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) 
    return(phy)
  if (length(tip) == Ntip) {
    if (Nnode(phy) < 3 || trim.internal) {
      warning("drop all tips of the tree: returning NULL")
      return(NULL)
    }
  }
  wbl <- !is.null(phy$edge.length)
  if (!rooted && subtree) {
    phy <- root(phy, (1:Ntip)[-tip][1])
    root.edge <- 0
  }
  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  if (trim.internal) {
    ints <- edge2 > Ntip
    repeat {
      sel <- !(edge2 %in% edge1[keep]) & ints & keep
      if (!sum(sel)) 
        break
      keep[sel] <- FALSE
    }
    if (subtree) {
      subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
      keep[subt] <- TRUE
    }
    if (root.edge && wbl) {
      degree <- tabulate(edge1[keep])
      if (degree[ROOT] == 1) {
        j <- integer(0)
        repeat {
          i <- which(edge1 == NEWROOT & keep)
          j <- c(i, j)
          NEWROOT <- edge2[i]
          if (degree[NEWROOT] > 1) 
            break
        }
        keep[j] <- FALSE
        j <- j[1:root.edge]
        NewRootEdge <- sum(phy$edge.length[j])
        if (length(j) < root.edge && !is.null(phy$root.edge)) 
          NewRootEdge <- NewRootEdge + phy$root.edge
        phy$root.edge <- NewRootEdge
      }
    }
  }
  if (!root.edge) 
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  
  phy
}
