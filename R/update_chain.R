update_chain = function(H,
                        H_max,
                        epsilon,
                        alpha,
                        alpha_precision,
                        psi,
                        sigma_theta,
                        z,
                        mu,
                        sigma,
                        omega,
                        beta,
                        w,
                        pw,
                        pz,
                        model_data,
                        alpha_fixed = FALSE,
                        z_fixed = FALSE,
                        w_fixed = FALSE,
                        mixscat_prior = TRUE,
                        alpha_prior = "cusp",
                        add_cluster = TRUE,
                        logP_proposal = NULL,
                        smooth,
                        global_update,
                        add_sigma_w,
                        est_epsilon = FALSE,
                        time_pz,
                        sigma_state,
                        sigma_mu,
                        clust_var,
                        w_prior,
                        adapt_H) {


  theta = update_theta_new(
    mu = mu[, 1:H],
    sigma = sigma,
    alpha = alpha[, 1:H],
    psi = psi,
    z = z,
    H = H,
    w = w,
    sigma_state = sigma_state,
    model_data = model_data,
    smooth = smooth,
    clust_var = clust_var
  ) %>% complete_dim(dimension = H_max)

  theta_sd = apply(theta, 2, sd)
  theta = scale(theta, center = FALSE, scale = theta_sd)

  mu = update_mu_new(
    H = H,
    z = z,
    w = w,
    sigma = sigma,
    theta = theta[, 1:H],
    model_data = model_data,
    sigma_state = sigma_state,
    clust_var = clust_var,
    sigma_mu = sigma_mu,
    smooth = smooth
  ) %>% complete_dim(dimension = H_max)

  # sigma_mu = update_sigma_mu(
  #   mu = mu[, 1:H], model_data = model_data, H = H
  # ) %>% complete_dim(dimension = H_max)

  if(clust_var == TRUE) {

    sigma_state = update_sigma_clust(
      H = H,
      theta = theta[, 1:H],
      mu = mu[, 1:H],
      z = z,
      model_data = model_data,
      smooth = smooth,
      w = w
    )

  }else{

    sigma = update_sigma_new(
      theta = theta[, 1:H],
      H = H,
      mu = mu[, 1:H],
      w = w,
      z = z,
      model_data = model_data,
      smooth = smooth
    )

  }

  if(z_fixed == FALSE) {

    z_out = update_z_new(
      H = H,
      theta = theta[, 1:H],
      mu = mu[, 1:H],
      z = z,
      sigma = sigma,
      pz = pz,
      model_data = model_data,
      beta = beta,
      mixscat_prior = mixscat_prior,
      w = w,
      clust_var = clust_var,
      sigma_state = sigma_state,
      logP_proposal = logP_proposal,
      global_update = global_update
    )

    z = z_out$z
    accept = z_out$accept
    log_accept_prob = z_out$log_accept_prob
    logP = z_out$logP


  }else{
    accept = NA
    logP = NULL
    accept = 1
    log_accept_prob = 0
  }

  if(alpha_fixed == FALSE) {

    alpha = update_alpha(
      theta = cbind(theta[, 1:H]),
      psi = psi,
      prior_precision = cbind(alpha_precision[, 1:H]),
      model_data = model_data
    ) %>% complete_dim(dimension = H_max)

    alpha = scale(alpha, center = FALSE, scale = 1/theta_sd)

    if(alpha_prior == "normal") {

      alpha_precision = matrix(
        1, nrow = model_data$dims$J, ncol = model_data$dims$K
      )

      H_active = model_data$dims$K

    }else if(alpha_prior == "cusp") {

      out_cusp = update_cusp(
        alpha = cbind(alpha[, 1:H]),
        omega = omega[1:H],
        nu = model_data$cusp$nu,
        a = model_data$cusp$a,
        b = model_data$cusp$b,
        min_var = model_data$cusp$min_var,
        H = H
      )

      ind = out_cusp$ind
      H_active = sum(ind[1:H] > 1:H, na.rm = TRUE)
      omega = out_cusp$omega[1:H] %>% complete_dim(dimension = H_max)
      v = out_cusp$v[1:H] %>% complete_dim(dimension = H_max)

      alpha_precision = matrix(
        out_cusp$prec,
        nrow = model_data$dims$J,
        ncol = H,
        byrow = TRUE
      ) %>% complete_dim(dimension = H_max)

      if(adapt_H == TRUE) {

        out_adapt = adapt_cusp(
          H_active = H_active,
          H = H,
          H_max = H_max,
          ind = ind,
          alpha = alpha,
          mu = mu,
          theta = theta,
          omega = omega,
          sigma = sigma,
          alpha_precision = alpha_precision,
          v = v,
          z = z,
          nu = model_data$cusp$nu,
          min_var = model_data$cusp$min_var,
          model_data = model_data
        )

        alpha = out_adapt$alpha
        mu = out_adapt$mu
        theta = out_adapt$theta
        omega = out_adapt$omega
        alpha_precision = out_adapt$alpha_precision
        v = out_adapt$v
        H = out_adapt$H

      }

    }

  }

  psi = update_psi(
    theta = cbind(theta[, 1:H]),
    alpha = cbind(alpha[, 1:H]),
    model_data = model_data,
    a = 1,
    b = 1
  )

  if(mixscat_prior == TRUE) {

    beta = update_beta(
      beta = beta,
      z = z,
      w = w,
      model_data = model_data
    )

    if(is.null(w_prior)) {
      pw = update_pw(
        w = w,
        model_data = model_data, epsilon = model_data$theta_spline$w_dirichlet
      )
    }

    if(w_fixed == FALSE) {

      w_out = update_w(
        beta = beta,
        z = z,
        theta = theta,
        mu = mu,
        sigma = sigma,
        pw = rbind(pw),
        w_prior = w_prior,
        model_data = model_data
      )

      w = w_out$w
      w_post_prob = w_out$w_post_prob

    }

  }else{

    pz = update_pz(
      z = z, z_dir = 1, model_data = model_data, time_pz = FALSE
    )

  }

  out = list(
    theta = theta,
    alpha = alpha,
    psi = psi,
    z = z,
    w = w,
    mu = mu,
    sigma = sigma,
    beta = beta,
    sigma_mu = sigma_mu,
    sigma_state = sigma_state,
    log_accept_prob = log_accept_prob,
    pw = pw,
    pz = pz,
    alpha_precision = alpha_precision,
    omega = omega,
    H_active = H_active,
    H = H,
    accept = accept,
    logP = logP
  )

  return(out)

}


