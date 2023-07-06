keep.tip.custom <- function (phy, tip) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) {
    idx <- match(tip, phy$tip.label)
    if (anyNA(idx)) {
      um <- c("umatched tip labels:\n", paste(tip[is.na(idx)], 
                                              collapse = " "))
      stop(um)
    }
    tip <- idx
  }else {
    out.of.range <- tip > Ntip
    if (any(out.of.range)) {
      warning("some tip numbers were larger than the number of tips: they were ignored")
      tip <- tip[!out.of.range]
    }
  }
  toDrop <- setdiff(1:Ntip, tip)
  drop.tip.custom(phy, toDrop)
}
