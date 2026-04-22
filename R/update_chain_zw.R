update_chain_zw = function(H,
                           H_max,
                           epsilon,
                           alpha,
                           alpha_precision,
                           psi,
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
                           alpha_prior = "cusp",
                           add_cluster = TRUE,
                           adapt_H) {

  theta = update_theta_new(
    mu = mu[, 1:H],
    sigma = sigma,
    alpha = alpha[, 1:H],
    psi = psi,
    z = z,
    H = H,
    w = w,
    mixscat_prior = TRUE,
    model_data = model_data,
    add_cluster = add_cluster,
    smooth = TRUE
  ) %>% complete_dim(dimension = H_max)

  theta_sd = apply(theta, 2, sd)
  theta = scale(theta, center = FALSE, scale = theta_sd)

  mu = update_mu_new(
    H = H,
    z = z,
    w = w,
    mixscat_prior = TRUE,
    sigma = sigma,
    theta = theta[, 1:H],
    model_data = model_data,
    smooth = TRUE
  ) %>% complete_dim(dimension = H_max)

  sigma = update_sigma_new(
    theta = theta[, 1:H],
    H = H,
    mu = mu[, 1:H],
    mixscat_prior = TRUE,
    w = w,
    z = z,
    model_data = model_data,
    smooth = TRUE
  )

  if(z_fixed == FALSE) {

    z_out = update_zw(
      w = w,
      theta = theta,
      mu = mu,
      sigma = sigma,
      H = H,
      pz = pz,
      beta = beta,
      model_data = model_data
    )

    z = z_out$z
    z_w = z_out$z_w


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

  # beta = update_beta(
  #   beta = beta,
  #   z = z,
  #   w = w,
  #   model_data = model_data
  # )

  pw = update_pw(
    w = w,
    model_data = model_data, epsilon = model_data$theta_spline$w_dirichlet
  )

  if(w_fixed == FALSE) {

    w_out = update_wz(
      beta = beta,
      z_w = z_w,
      pw = pw,
      theta = theta,
      mu = mu,
      sigma = sigma,
      model_data = model_data,
      smooth = TRUE
    )

    w = w_out$w
    w_post_prob = w_out$w_post_prob

  }


  pz = update_pz(
    z = as.numeric(z_w), z_dir = 0.01, model_data = model_data
  )


  out = list(
    theta = theta,
    epsilon = epsilon,
    alpha = alpha,
    psi = psi,
    z = z,
    mu = mu,
    sigma = sigma,
    beta = beta,
    w = w,
    pw = pw,
    pz = pz,
    alpha_precision = alpha_precision,
    omega = omega,
    H_active = H_active,
    H = H,
    z_w = z_w,
    w_post_prob = w_post_prob
  )

  return(out)

}


