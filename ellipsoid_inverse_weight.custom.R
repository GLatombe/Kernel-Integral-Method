ellipsoid_inverse_weight.custom <- function (samples, centers, scales, scales.long, verbose) 
{
  
  # samples = full_samples
  # centers = data
  # scales = scales
  # scales.long = scales.long
  
  
  for (i in 1:ncol(centers)) {
    samples[, i] = samples[, i]/scales.long[, i]
    centers[, i] = centers[, i]/scales[, i]
  }
  tree <- hypervolume:::kdtree_build(centers, verbose = verbose)
  query <- hypervolume:::kdtree_ball_query_multiple(tree, t(samples), nrow(samples), 
                                      ncol(samples), r = 1, verb = verbose)
  rm(tree)
  return(query)
}
