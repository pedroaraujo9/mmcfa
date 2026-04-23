update_chain = function(H,
                        H_max,
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
                        w_prior,
                        model_data,
                        z_fixed,
                        w_fixed,
                        mixscat_prior,
                        alpha_prior,
                        adapt_H,
                        add_cluster = TRUE
                        ) {


  theta = update_theta(
    H = H,
    mu = mu[, 1:H],
    sigma = sigma,
    alpha = alpha[, 1:H],
    psi = psi,
    z = z,
    model_data = model_data,
    add_cluster = add_cluster
  ) %>% complete_dim(dimension = H_max)

  theta_sd = apply(theta, 2, sd)
  theta = scale(theta, center = FALSE, scale = theta_sd)

  z_post_prob = NA
  w_post_prob = NA

  if(add_cluster == TRUE) {

      mu = update_mu(
      H = H,
      z = z,
      sigma = sigma,
      theta = theta[, 1:H],
      model_data = model_data
    ) %>% complete_dim(dimension = H_max)

    sigma = update_sigma(
      H = H,
      theta = theta[, 1:H],
      mu = mu[, 1:H],
      z = z,
      model_data = model_data,
      add_cluster = add_cluster
    )

    if(z_fixed == FALSE) {

      z_out = update_z(
        H = H,
        w = w,
        theta = theta[, 1:H],
        mu = mu[, 1:H],
        sigma = sigma,
        pz = pz,
        beta = beta,
        model_data = model_data,
        mixscat_prior = mixscat_prior
      )

      z = z_out$z
      z_post_prob = z_out$z_post_prob

    }else{
      z_post_prob = NA
    }

  }else{

    sigma = update_sigma(
      H = H,
      theta = theta[, 1:H],
      mu = mu[, 1:H],
      z = z,
      model_data = model_data,
      add_cluster = add_cluster
    )

  }

  alpha = update_alpha(
    H = H,
    theta = cbind(theta[, 1:H]),
    psi = psi,
    prior_precision = alpha_precision[1:H],
    model_data = model_data
  ) %>% complete_dim(dimension = H_max)

  alpha = scale(alpha, center = FALSE, scale = 1/theta_sd)

  if(alpha_prior == "normal") {

    alpha_precision = rep(1, model_data$dims$K)
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
    alpha_precision = out_cusp$prec %>% complete_dim(dimension = H_max)

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
        model_data = model_data,
        add_cluster = add_cluster
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

  psi = update_psi(
    theta = cbind(theta[, 1:H]),
    alpha = cbind(alpha[, 1:H]),
    model_data = model_data,
    a = 2,
    b = 2
  )

  if(add_cluster == TRUE) {

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
          model_data = model_data,
          epsilon = model_data$clustering$w_dirichlet
        )
      }

      if(w_fixed == FALSE) {

        w_out = update_w(
          beta = beta,
          z = z,
          pw = pw,
          model_data = model_data,
          w_prior = w_prior
        )

        w = w_out$w
        w_post_prob = w_out$w_post_prob

      }else{
        w_post_prob = NA
      }

      }else{

      pz = update_pz(
        z = z, z_dir = model_data$clustering$z_dirichlet, model_data = model_data
      )

      w_post_prob = NA
    }

  }else{
    w_post_prob = NA 
    z_post_prob = NA
  }

  
  out = list(
    theta = theta,
    alpha = alpha,
    psi = psi,
    z = z,
    w = w,
    mu = mu,
    beta = beta,
    sigma = sigma,
    pw = pw,
    pz = pz,
    alpha_precision = alpha_precision,
    omega = omega,
    H_active = H_active,
    H = H,
    w_post_prob = w_post_prob,
    z_post_prob = z_post_prob
  )

  return(out)

}


