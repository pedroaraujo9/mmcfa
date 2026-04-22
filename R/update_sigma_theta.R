update_sigma_theta = function(epsilon, mu, z, model_data, smooth, add_cluster = TRUE) {

  G = model_data$dims$G
  H = model_data$dims$K
  n = model_data$dims$n

  if(add_cluster == TRUE) {

    ess2 = colSums((epsilon - mu[z, ])^2)

  }else{

    ess2 = colSums((epsilon)^2)

  }

  a = 2 + 0.5 * n
  b = 2 + 0.5 * ess2
  sigma_theta = sqrt(1/rgamma(H, shape = a, rate = b))

  return(sigma_theta)

}
