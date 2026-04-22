single_run_zw = function(model_data,
                         iters,
                         burn_in,
                         thin,
                         tune = 0,
                         alpha_prior,
                         z = NULL,
                         w = NULL,
                         adapt_H,
                         alpha0,
                         alpha1,
                         init_list = NULL,
                         seed) {

  init_time = Sys.time()

  if(is.null(init_list)) {

    init_list = create_init_values2(
      y = model_data$data$y,
      id = model_data$data$id,
      time = model_data$data$time,
      K_max = model_data$dims$K,
      G_max = model_data$dims$G,
      M_max = model_data$dims$M,
      n_basis = model_data$theta_spline$n_basis
    )

    sample_list = create_sample_list(
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      model_data = model_data,
      alpha_prior = alpha_prior,
      init_list = NULL,
      mixscat_prior = TRUE,
      seed = seed
    )

    alpha = init_list$alpha

  }else{

    sample_list = create_sample_list(
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      model_data = model_data,
      alpha_prior = alpha_prior,
      init_list = init_list,
      mixscat_prior = TRUE,
      seed = seed
    )

    alpha = sample_list$alpha[1,,]

  }

  G = model_data$dims$G
  K = model_data$dims$K
  M = model_data$dims$M

  if(!is.null(z)) {
    z_fixed = TRUE
  }else{
    z_fixed = FALSE
    z = sample_list$z[1,]
  }

  if(!is.null(w)) {
    w_fixed = TRUE
  }else{
    w_fixed = FALSE
    w = sample_list$w[1,]

  }


  sigma = sample_list$sigma[1,]
  mu = sample_list$mu[1,,]
  beta = sample_list$beta[1,,]
  pw = sample_list$pw[1,]
  epsilon = sample_list$epsilon[1,,]
  alpha_precision = sample_list$alpha_precision[1,,]
  psi = sample_list$psi[1,]
  pz = sample_list$pz[1, ]
  omega = sample_list$omega[1, ]
  ind = sample_list$ind[1, ]
  sigma_theta = sample_list$sigma_theta[1, ]
  v = sample_list$v[1, ]
  z_w = sample_list$z_w[1,,]

  H_max = model_data$dims$K

  if(adapt_H == TRUE) {

    H = H_max - 1

  }else{

    H = H_max

  }

  i = 1

  for(iter in (-tune + 1):iters) {

    cat(iter, "\r")

    update = update_chain_zw(
      H = H,
      H_max = H_max,
      epsilon = epsilon,
      alpha = alpha,
      alpha_precision = alpha_precision,
      psi = psi,
      z = z,
      mu = mu,
      omega = omega,
      sigma = sigma,
      beta = beta,
      w = w,
      pw = pw,
      pz = pz,
      model_data = model_data,
      z_fixed = z_fixed,
      w_fixed = w_fixed,
      alpha_fixed = FALSE,
      alpha_prior = alpha_prior,
      adapt_H = (runif(1) < exp(-(alpha0 + alpha1 * iter))) & (adapt_H == TRUE) & (iter > 0)
    )

    theta = update$theta
    epsilon = update$epsilon
    alpha = update$alpha
    psi = update$psi
    z = update$z
    mu = update$mu
    sigma = update$sigma
    beta = update$beta
    w = update$w
    pw = update$pw
    pz = update$pz
    omega = update$omega
    alpha_precision = update$alpha_precision
    H_active = update$H_active
    H = update$H
    w_post_prob = update$w_post_prob
    z_w = update$z_w

    if(iter %in% sample_list$iters_vec) {

      sample_list$alpha[i,,] = alpha
      sample_list$theta[i,,] = theta
      sample_list$epsilon[i,,] = epsilon
      sample_list$alpha_precision[i,,] = alpha_precision
      sample_list$psi[i, ] = psi
      sample_list$z[i, ] = z
      sample_list$mu[i,,] = mu
      sample_list$sigma[i,] = sigma
      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = as.numeric(pw)
      sample_list$pz[i, ] = pz
      sample_list$H_active[i] = H_active
      sample_list$H[i] = H
      sample_list$omega[i, ] = omega
      sample_list$z_w[i,,] = z_w
      sample_list$w_post_prob[i,,] = w_post_prob


      i = i + 1

    }

  }

  rotation_list = get_rotation_list(
    sample_list$alpha, reference_matrix = init_list$alpha
  )

  for(param in c("alpha", "theta", "mu", "epsilon")){
    sample_list[[param]] = apply_rotation(
      sample_list[[param]], rotation_list = rotation_list
    )
  }

  out = list(
    sample_list = sample_list,
    init = init_list
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)


}
