update_mu = function(epsilon, z, sigma, model_data) {

  G = model_data$dims$G
  H = ncol(epsilon)

  mu = matrix(0, nrow = G, ncol = H)

  for(g in 1:G) {

    ng = sum(z == g)

    if(ng == 0) {

      mu[g, ] = rnorm(H, 0, 1)

    } else {

      epsilon_bar = colMeans(epsilon[z == g, , drop = FALSE])

      data_prec  = ng / (sigma[g]^2)
      prior_prec = 1

      post_prec = data_prec + prior_prec
      post_var  = 1 / post_prec
      post_mean = post_var * data_prec * epsilon_bar

      mu[g, ] = rnorm(H, post_mean, sqrt(post_var))
    }
  }

  return(mu)

}
