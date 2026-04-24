update_z = function(H,
                    w,
                    theta,
                    mu,
                    sigma,
                    pz,
                    beta,
                    model_data,
                    mixscat_prior) {

  G = model_data$dims$G
  M = model_data$dims$M
  n = model_data$dims$n

  id = model_data$data$id
  B = model_data$clustering$B
  idx = as.integer(id)

  if(mixscat_prior == TRUE) {

    prob = compute_probs(w = w, M = M, B = B, beta = beta)

  }else{

    prob = matrix(as.numeric(pz), nrow = n, ncol = G, byrow = T)

  }

  logP = lapply(1:G, function(g){

    mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
    s = sigma[cbind(idx, g)]

    logp = dnorm(theta, mean = mu_g, sd = s, log = TRUE)
    rowSums(logp)

  }) |> (\(x) do.call(cbind, x))()

  logP = norm_mat(logP + log(prob))
  z_post_prob = exp(logP)
  z = sample_cat(logP)

  out = list(z = z, z_post_prob = z_post_prob)

  return(out)

}
