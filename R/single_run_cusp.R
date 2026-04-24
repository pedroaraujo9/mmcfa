single_run_cusp = function(model_data,
                            iters,
                            burn_in,
                            thin,
                            cusp_adapt_H,
                            cusp_alpha0,
                            cusp_alpha1,
                            est_sigma = TRUE,
                            init_list = NULL,
                            seed = NULL,
                            verbose = TRUE) {

  init_time = Sys.time()

  if(is.null(init_list)) {

    init_list = find_init_FA(
      y = model_data$data$y,
      time = model_data$data$time,
      id = model_data$data$id,
      K_max = model_data$dims$K
    )

    sample_list = create_sample_list(
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      model_data = model_data,
      alpha_prior = "cusp",
      init_list = NULL,
      seed = seed
    )

    alpha = init_list$alpha

  }else{

    sample_list = create_sample_list(
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      model_data = model_data,
      alpha_prior = "cusp",
      init_list = init_list,
      seed = seed
    )

    alpha = sample_list$alpha[1,,]

  }

  sigma = sample_list$sigma[1,,]

  if(est_sigma == FALSE) {
    sigma = matrix(1, nrow = model_data$dims$n_id, ncol = model_data$dims$G)
  }

  alpha_precision = sample_list$alpha_precision[1,]
  psi = sample_list$psi[1,]
  omega = sample_list$omega[1, ]
  v = sample_list$v[1, ]

  H_max = model_data$dims$K

  if(cusp_adapt_H == TRUE) {
    H = H_max - 1
  }else{
    H = H_max
  }

  i = 1

  for(iter in 1:iters) {

    if(verbose == TRUE) cat(iter, "\r")

    update = update_chain_cusp(
      H = H,
      H_max = H_max,
      alpha = alpha,
      alpha_precision = alpha_precision,
      psi = psi,
      sigma = sigma,
      omega = omega,
      v = v,
      model_data = model_data,
      adapt_H = (runif(1) < exp(-(cusp_alpha0 + cusp_alpha1 * iter))) & (cusp_adapt_H == TRUE),
      est_sigma = est_sigma
    )

    theta = update$theta
    alpha = update$alpha
    psi = update$psi
    omega = update$omega
    alpha_precision = update$alpha_precision
    H_active = update$H_active
    H = update$H
    sigma = update$sigma
    v = update$v

    if(iter %in% sample_list$iters_vec) {

      sample_list$alpha[i,,] = alpha
      sample_list$theta[i,,] = theta
      sample_list$alpha_precision[i,] = alpha_precision
      sample_list$psi[i, ] = psi
      sample_list$sigma[i,,] = sigma
      sample_list$H_active[i] = H_active
      sample_list$H[i] = H
      sample_list$omega[i, ] = omega
      sample_list$v[i, ] = v

      i = i + 1

    }

  }

  if(cusp_adapt_H == FALSE) {

    rotation_list = get_rotation_list(
      sample_list$alpha,
      reference_matrix = init_list$alpha
    )

    for(param in c("alpha", "theta")) {
      sample_list[[param]] = apply_rotation(
        sample_list[[param]],
        rotation_list = rotation_list
      )
    }

  }

  out = list(
    sample_list = sample_list,
    init = init_list
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)

}
