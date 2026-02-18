update_z = function(epsilon, mu, sigma, w, beta, model_data, add_prob_spline) {

  H = ncol(epsilon)
  G = model_data$dims$G
  M = model_data$dims$M
  n = model_data$dims$n
  id = model_data$data$id
  B = model_data$theta_spline$B_theta

  if(G > 1) {

    ll = matrix(0, nrow = n, ncol = G)

    for(g in 1:G) {

      avg = matrix(mu[g, ], nrow = n, ncol = H, byrow = T)
      std = matrix(sigma[g], nrow = n, ncol = H)

      ll[, g] = dnorm(epsilon, mean = avg, sd = std, log = T) |> rowSums()

    }

    if(add_prob_spline == TRUE) {

      W = create_dummy(w, M)
      X = kronecker(W, B)
      prob = mclust::softmax(X %*% beta)
      #prob = predict_prob_cpp(M, w, B, beta)
      ll = ll + log(prob)

    }else{

      ll = ll - log(G)

    }

    ll = ll - matrix(
      mclust::logsumexp(ll), nrow = n, ncol = G, byrow = F
    )

    z = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)

  }else{

    z = rep(1L, n)

  }


  return(z)

}
