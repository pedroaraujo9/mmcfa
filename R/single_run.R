single_run = function(model_data,
                      iters,
                      burn_in,
                      thin,
                      alpha_prior,
                      add_prob_spline,
                      seed) {

  init_time = Sys.time()

  init_list = create_init_values(model_data = model_data)

  sample_list = create_sample_list(
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    model_data = model_data,
    alpha_prior = alpha_prior,
    init_list = NULL,
    seed = seed
  )

  G = model_data$dims$G
  K = model_data$dims$K
  M = model_data$dims$M

  alpha = init_list$alpha

  #z = init_list$z
  #w = init_list$w
  #mu = init_list$mu
  z = sample_list$z[1,] / sample_list$z[1,]
  w = sample_list$w[1,] / sample_list$w[1,]
  mu = sample_list$mu[1,,]

  beta = sample_list$beta[1,,]
  pw = sample_list$w[1,]
  epsilon = sample_list$epsilon[1,,]
  sigma = sample_list$sigma[1,]
  alpha_precision = sample_list$alpha_precision[1,,]
  psi = sample_list$psi[1,]

  H = model_data$dims$K
  i = 1

  for(iter in 1:iters) {

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
      add_prob_spline = add_prob_spline,
      add_cluster = TRUE,
      update_w_params = TRUE,
      update_z = TRUE,
      update_w = TRUE
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

    if(iter %in% sample_list$iters_vec) {

      sample_list$alpha[i,,] = alpha
      sample_list$theta[i,,] = theta
      sample_list$epsilon[i,,] = epsilon
      sample_list$alpha_precision[i,,] = alpha_precision
      sample_list$psi[i, ] = psi
      sample_list$z[i, ] = z
      sample_list$mu[i,,] = mu
      sample_list$sigma[i, ] = sigma
      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = pw

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
