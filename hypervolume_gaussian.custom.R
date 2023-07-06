hypervolume_gaussian.custom <- function (data, name = NULL, weight = NULL, samples.per.point = ceiling((10^(3+sqrt(ncol(data))))/nrow(data)), kde.bandwidth = estimate_bandwidth(data), scales = kde.bandwidth * (sd.count), sd.count = 3, quantile.requested = 0.95, quantile.requested.type = "probability", chunk.size = chunk.size, verbose = TRUE, ...) 
{
  
  if(length(scales)==1)
    scales = matrix(scales,nrow(data),ncol(data))
  
  
  data = as.matrix(data)
  d = ncol(data)
  np = nrow(data)
  if (is.null(weight)) {
    weight <- rep(1/nrow(data), nrow(data))
  }
  if (is.null(dimnames(data)[[2]])) {
    dimnames(data)[[2]] <- paste("X", 1:ncol(data))
  }
  if (sd.count < 3) {
    warning(sprintf("Values of sd.count (%d) is low.\nRecommended minimum value is 3, with higher values giving better performance.\nBoundaries and volumes may be inaccurate.", 
                    sd.count))
  }
  if (ncol(data) != length(kde.bandwidth)) {
    stop("data and kde.bandwidth must have same dimensionality")
  }
  if (any(kde.bandwidth == 0)) {
    stop("Bandwidth must be non-zero.")
  }
  names(kde.bandwidth) <- dimnames(data)[[2]]
  if (length(weight) != nrow(data)) {
    stop("The length of the weights must be equal to the number of observations.")
  }
  if (abs(sum(weight) - 1) > 10^-10) {
    warning("The sum of the weights must be equal to 1. Normalizing the weights.")
    weight <- weight/sum(weight)
  }
  kde.method = attr(kde.bandwidth, "method")
  if (is.null(kde.method)) {
    stop("KDE bandwidth must be set using estimate_bandwidth(), not entered manually.")
  }
  predict_function_gaussian <- function(x) {
    return(calculate_density(s = x, means = data, kernel_sd = kde.bandwidth, 
                             weight = weight))
  }
  
  samples_all = sample_model_ellipsoid.custom(predict_function = predict_function_gaussian.custom, data = data, scales = scales, min.value = 0, samples.per.point = samples.per.point, chunk.size = chunk.size, verbose = verbose)
  
  random_points_full = samples_all$full_samples
  return(list(random_points_full))
}
