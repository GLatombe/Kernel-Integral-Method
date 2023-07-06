sample_model_ellipsoid.custom <- function (predict_function = NULL, data, scales, min.value, 
          samples.per.point, chunk.size = chunk.size, verbose = TRUE, return.full = FALSE, kde.bandwidth=estimate_bandwidth_silent(data)) 
{
  
  data = na.omit(as.matrix(data))
  d = ncol(data)
  N.samples <- ceiling(samples.per.point * nrow(data))
  if (is.null(dimnames(data)[[2]])) {
    dimnames(data) <- list(NULL, paste("X", 1:d, sep = ""))
  }
  if (verbose == TRUE) {
    pb <- progress_bar$new(total = N.samples)
    pb$tick(0)
  }
  samples = list()
  volume_sampling_extent_all = list()
  total_accepted <- 0
  total_tried <- 0
  #while (total_accepted < N.samples) {
    if (verbose == TRUE) {
      if (!pb$finished == TRUE) {
        pb$update(total_accepted/N.samples)
      }
    }
  
    if(length(chunk.size)==1){
      full_samples = lapply(1:nrow(data), function(i) {
        se = hypervolume:::sample_ellipsoid(data[i, ], chunk.size, scales = scales[i,])
        return(data.frame(se))
      })
    }else{
      full_samples = lapply(1:nrow(data), function(i) {
        se = hypervolume:::sample_ellipsoid(data[i, ], chunk.size[i], scales = scales[i,])
        return(data.frame(se))
      })
    }
   
  
  full_samples = as.matrix(rbindlist(full_samples))
  
  
  if (verbose == TRUE) {
    pb$terminate()
  }
  
  return(list(full_samples = full_samples))
}
