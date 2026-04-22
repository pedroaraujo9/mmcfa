update_mu = function(epsilon, z, sigma, model_data) {

  G = model_data$dims$G
  H = model_data$dims$K

  mu = matrix(0, nrow = G, ncol = H)

  for(g in 1:G) {

    ng = sum(z == g)

    if(ng == 0) {

      mu[g, ] = rnorm(H, mean = 0, sd = 1)

    } else {

      epsilon_bar = colMeans(epsilon[z == g, , drop = FALSE])
      prec = 1/(sigma[g, ]^2)

      nu = 1 / (1 + ng/prec)
      m = nu * epsilon_bar * (ng/prec)
      mu[g, ] = rnorm(H, m, sqrt(nu))

    }
  }

  return(mu)

}
