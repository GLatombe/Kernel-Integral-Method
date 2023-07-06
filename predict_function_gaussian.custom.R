predict_function_gaussian.custom <- function(x) {
  return(hypervolume:::calculate_density(s = x, means = data, kernel_sd = kde.bandwidth, 
                           weight = weight))
}
