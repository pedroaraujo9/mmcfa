update_chain = function(H,
                        epsilon,
                        alpha,
                        alpha_precision,
                        psi,
                        z,
                        mu,
                        sigma,
                        beta,
                        w,
                        pw,
                        model_data,
                        add_cluster = TRUE,
                        add_prob_spline = TRUE,
                        update_w_params = TRUE,
                        update_z = TRUE,
                        update_w = TRUE) {

  epsilon = update_epsilon(
    mu = mu,
    sigma = sigma,
    alpha = alpha,
    psi = psi,
    z = z,
    model_data,
    add_cluster = add_cluster
  )

  theta = update_theta(
    epsilon = epsilon,
    model_data = model_data
  )

  alpha = update_alpha(
    theta = cbind(theta[, 1:H]),
    psi = psi,
    prior_precision = cbind(alpha_precision[, 1:H]),
    model_data = model_data
  )

  psi = update_psi(
    theta = cbind(theta[, 1:H]),
    alpha = cbind(alpha[, 1:H]),
    model_data = model_data,
    a = 1,
    b = 1
  )

  mu = update_mu(
    epsilon = epsilon,
    sigma = sigma,
    z = z,
    model_data = model_data
  )



  if(update_z == TRUE) {

    z = update_z(
      epsilon = epsilon,
      mu = mu,
      sigma = sigma,
      w = w,
      beta = beta,
      model_data = model_data,
      add_prob_spline = add_prob_spline
    )

  }

  if(update_w_params == TRUE) {

    beta = update_beta(
      beta = beta,
      z = z,
      w = w,
      model_data = model_data
    )

    pw = update_pw(
      w = w,
      model_data = model_data
    )

  }

  if(update_w == TRUE) {

    w_out = update_w(
      beta = beta,
      z = z,
      pw = rbind(pw),
      model_data = model_data
    )

    w = w_out$w
    w_post_prob = w_out$w_post_prob

  }


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
    pw = pw
  )

  return(out)

}


