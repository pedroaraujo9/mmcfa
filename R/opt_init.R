opt_init = function(iters_opt = 100,
                    alpha,
                    init_list,
                    sample_list,
                    model_data) {
  w = init_list$w
  z = init_list$z
  alpha = init_list$alpha
  mu = init_list$mu
  H = model_data$dims$K

  beta = sample_list$beta[1,,]
  pw = sample_list$pw[1,]
  sigma = sample_list$sigma[1,]
  epsilon = sample_list$epsilon[1,,]
  alpha_precision = sample_list$alpha_precision[1,,]
  psi = sample_list$psi[1,]

  for(iter in 1:iters_opt) {

    cat(iter, "\r")

    update = update_chain(
      H = H,
      epsilon = epsilon,
      alpha = alpha,
      alpha_precision = alpha_precision,
      psi = psi,
      z = z,
      mu = mu,
      sigma = sigma,
      beta = beta,
      w = w,
      pw = pw,
      model_data = model_data,
      add_prob_spline = TRUE,
      add_cluster = TRUE,
      update_w_params = TRUE,
      update_z = FALSE,
      update_w = FALSE
    )

    epsilon = update$epsilon
    alpha = update$alpha
    psi = update$psi
    mu = update$mu
    beta = update$beta
    pw = update$pw

  }

  out = list(
    epsilon = epsilon,
    alpha = alpha,
    psi = psi,
    w = w,
    beta = beta,
    pw = pw,
    mu = mu
  )

  return(out)

}
