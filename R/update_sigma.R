update_sigma = function(epsilon,
                        mu,
                        z,
                        model_data,
                        sigma_a = 2,
                        sigma_b = 2) {

  G = model_data$dims$G
  H = ncol(epsilon)

  res = (epsilon - mu[z, ])^2

  sigma = numeric(G)

  for(g in 1:G) {

    if(sum(z == g) == 0) {

      sigma[g, ] = 1/sqrt(rgamma(n = 1, shape = sigma_a, rate = sigma_b))

    }else{

      ss = colSums(rbind(res[z == g, ]))
      ng = sum(z == g)
      sigma[g, ] = 1/sqrt(rgamma(n = H, shape = rep(sigma_a + ng/2, H), rate = sigma_b + ss/2))

    }

  }

  return(sigma)

}

# update_sigma = function(mu, model_data, sigma_a = 1, sigma_b = 1) {
#
#   G = model_data$dims$G
#   H = ncol(mu)
#
#   ess = rowSums(mu^2)
#   sigma = 1/sqrt(rgamma(n = G, shape = sigma_a + H/2, rate = sigma_b + ess/2))
#   return(sigma)
#
# }
