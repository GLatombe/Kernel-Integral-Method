sample_model_ellipsoid.custom.abundance <- function (predict_function = NULL, data, scales, min.value, 
          samples.per.point, chunk.size = 1000, verbose = TRUE, return.full = FALSE, kde.bandwidth=estimate_bandwidth_silent(data)) 
{
  
  data = na.omit(as.matrix(data))
  d = ncol(data)
  N.samples <- samples.per.point#ceiling(samples.per.point * nrow(data))
  if (is.null(dimnames(data)[[2]])) {
    dimnames(data) <- list(NULL, paste("X", 1:d, sep = ""))
  }
  # if (verbose == TRUE) {
  #   pb <- progress_bar$new(total = N.samples)
  #   pb$tick(0)
  # }
  samples = list()
  volume_sampling_extent_all = list()
  total_accepted <- 0
  total_tried <- 0
  #while (total_accepted < N.samples) {
    # if (verbose == TRUE) {
    #   if (!pb$finished == TRUE) {
    #     pb$update(total_accepted/N.samples)
    #   }
    # }
    full_samples = lapply(1:nrow(data), function(i) {
      se = hypervolume:::sample_ellipsoid(data[i, ], chunk.size[i], scales = scales)
      return(data.frame(se))
    })
    full_samples = as.matrix(rbindlist(full_samples))
    iweight = hypervolume:::ellipsoid_inverse_weight(full_samples, centers = data, 
                                       scales = scales, verbose = verbose)
    weight = 1/(iweight/min(iweight))
    mean_weight = sum(1/iweight)/nrow(full_samples)
    volume_sampling_extent = hypervolume:::ellipsoid_volume(scales) * 
      nrow(data) * mean_weight
    included = as.logical(rbinom(length(iweight), size = 1, 
                                 prob = weight))
    samples_retained = full_samples[included, , drop = FALSE]
    #predicted_values <- predict_function(samples_retained)
    predicted_values <- hypervolume:::calculate_density(s = samples_retained, means = data, kernel_sd = kde.bandwidth, weight = weight)
    included_thresholded = (as.numeric(predicted_values) > 
                              min.value)
    samples_retained_thresholded = samples_retained[included_thresholded, 
                                                    , drop = FALSE]
    predicted_values_thresholded = predicted_values[included_thresholded]
    samples_final_this = cbind(samples_retained_thresholded, 
                               predicted_values_thresholded)
    dimnames(samples_final_this) <- list(NULL, c(dimnames(data)[[2]], 
                                                 "value"))
    total_tried <- total_tried + nrow(samples_retained)
    total_accepted <- total_accepted + nrow(samples_retained_thresholded)
    samples <- c(samples, list(data.frame(samples_final_this)))
    volume_sampling_extent_all <- c(volume_sampling_extent_all, 
                                    volume_sampling_extent)
  #}
  samples <- as.matrix(rbindlist(samples))
  #samples <- samples[sample(1:nrow(samples), N.samples), , drop = FALSE]
  volume_sampling_extent_all_mean = mean(unlist(volume_sampling_extent_all), 
                                         na.rm = T)
  volume = volume_sampling_extent_all_mean * total_accepted/total_tried
  if (verbose == TRUE) {
    pb$terminate()
  }
  if (return.full == TRUE) {
    return(list(samples = samples, full_samples = full_samples, 
                volume = volume))
  }
  else {
    return(list(samples = samples, volume = volume, full_samples = full_samples))
  }
}
