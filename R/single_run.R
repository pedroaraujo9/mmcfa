single_run = function(model_data,
                      iters,
                      burn_in,
                      thin,
                      alpha_prior,
                      mixscat_prior,
                      cusp_adapt_H,
                      cusp_alpha0,
                      cusp_alpha1,
                      add_cluster = TRUE,
                      init_list = NULL,
                      z = NULL,
                      w = NULL,
                      w_prior = NULL,
                      seed = NULL,
                      verbose = TRUE) {

  init_time = Sys.time()

  if(is.null(init_list)) {

    init_list = find_init_FA(y = model_data$data$y,
                             time = model_data$data$time,
                             id = model_data$data$id,
                             K = model_data$dims$K)

    sample_list = create_sample_list(
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      model_data = model_data,
      alpha_prior = alpha_prior,
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
      alpha_prior = alpha_prior,
      init_list = init_list,
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


  mu = sample_list$mu[1,,]
  sigma = sample_list$sigma[1,,]
  beta = sample_list$beta[1,,]
  pw = sample_list$pw[1,]
  alpha_precision = sample_list$alpha_precision[1,]
  psi = sample_list$psi[1,]
  pz = sample_list$pz[1, ]
  omega = sample_list$omega[1, ]
  ind = sample_list$ind[1, ]
  v = sample_list$v[1, ]
  z_post_prob = sample_list$z_post_prob[1,,]

  H_max = model_data$dims$K

  if(cusp_adapt_H == TRUE) {

    H = H_max - 1

  }else{

    H = H_max

  }

  i = 1

  for(iter in 1:iters) {

    if(verbose == TRUE) cat(iter, "\r")

    update = update_chain(
      H = H,
      H_max = H_max,
      alpha = alpha,
      alpha_precision = alpha_precision,
      psi = psi,
      z = z,
      mu = mu,
      sigma,
      omega = omega,
      beta = beta,
      w = w,
      pw = pw,
      pz = pz,
      model_data = model_data,
      z_fixed = z_fixed,
      w_fixed = w_fixed,
      mixscat_prior = mixscat_prior,
      alpha_prior = alpha_prior,
      w_prior = w_prior,
      adapt_H = (alpha_prior == "cusp") & (runif(1) < exp(-(cusp_alpha0 + cusp_alpha1 * iter))) & (cusp_adapt_H == TRUE),
      add_cluster = add_cluster
    )

    theta = update$theta
    alpha = update$alpha
    psi = update$psi
    z = update$z
    mu = update$mu
    beta = update$beta
    w = update$w
    pw = update$pw
    pz = update$pz
    omega = update$omega
    alpha_precision = update$alpha_precision
    H_active = update$H_active
    H = update$H
    sigma = update$sigma
    w_post_prob = update$w_post_prob
    z_post_prob = update$z_post_prob

    if(iter %in% sample_list$iters_vec) {

      sample_list$alpha[i,,] = alpha
      sample_list$theta[i,,] = theta
      sample_list$alpha_precision[i,] = alpha_precision
      sample_list$psi[i, ] = psi
      sample_list$z[i, ] = z
      sample_list$mu[i,,] = mu
      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = as.numeric(pw)
      sample_list$sigma[i,,] = sigma
      sample_list$pz[i, ] = pz
      sample_list$w_post_prob[i,,] = w_post_prob
      sample_list$H_active[i] = H_active
      sample_list$H[i] = H
      sample_list$omega[i, ] = omega
      sample_list$z_post_prob[i,,] = z_post_prob

      i = i + 1

    }

  }

  if(cusp_adapt_H == FALSE) {

    rotation_list = get_rotation_list(
      sample_list$alpha, reference_matrix = init_list$alpha
    )

    for(param in c("alpha", "theta", "mu")){
      sample_list[[param]] = apply_rotation(
        sample_list[[param]], rotation_list = rotation_list
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
