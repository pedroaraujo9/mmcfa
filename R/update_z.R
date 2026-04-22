update_z = function(epsilon,
                    mu,
                    sigma,
                    alpha,
                    psi,
                    pz,
                    w,
                    beta,
                    model_data,
                    mixscat_prior,
                    integrate = FALSE) {

  H = ncol(epsilon)
  G = model_data$dims$G
  M = model_data$dims$M
  n = model_data$dims$n
  id = model_data$data$id
  B = model_data$theta_spline$B_theta
  #R = model_data$theta_spline$R
  #S = (R %*% t(R)) + diag(1, n_time)
  # s = matrix(sqrt(rep(diag(S), model_data$dims$n_id)), nrow = n, ncol = H, byrow = FALSE)

  if(G > 1) {

    ll = matrix(0, nrow = n, ncol = G)

    for(g in 1:G) {

      if(integrate == FALSE) {

        avg = matrix(mu[g, ], nrow = n, ncol = H, byrow = T)
        # ll[, g] = dnorm(epsilon, mean = avg, sd = sigma[g], log = T) |> rowSums()
        ll[, g] = dnorm(epsilon, mean = avg, sd = 1, log = T) |> rowSums()

      }else{

        ll[, g] = mvtnorm::dmvnorm(
          x = model_data$data$y,
          mean = mu[g, ] %*% t(alpha),
          sigma = (sigma[g]^2) * tcrossprod(alpha) + diag(psi),
          log = TRUE
        )


      }



    }

    if(mixscat_prior == TRUE) {

      prob = compute_probs(w = w, M = M, B = B, beta = beta)

    }else{

      prob = matrix(pz, nrow = n, ncol = G, byrow = T)

    }


    ll = ll + log(prob)

    ll = ll - matrix(
      mclust::logsumexp(ll), nrow = n, ncol = G, byrow = F
    )

    z = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)

  }else{

    z = rep(1L, n)

  }


  return(z)

}
