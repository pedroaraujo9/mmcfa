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
                        model_data) {

  epsilon = update_epsilon2(
    mu = mu,
    sigma = sigma,
    alpha = alpha,
    psi = psi,
    z = z,
    model_data
  )

  theta = update_theta(
    epsilon = epsilon,
    model_data = model_data
  )

  mu = update_mu(
    epsilon = epsilon,
    sigma = sigma,
    z = z,
    model_data = model_data
  )

  sigma = update_sigma(
    epsilon = epsilon,
    mu = mu,
    z = z,
    model_data = model_data,
    sigma_a = 1,
    sigma_b = 1
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

  z = update_z(
    epsilon = epsilon,
    mu = mu,
    sigma = sigma,
    w = w,
    beta = beta,
    model_data = model_data
  )

  beta = update_beta(
    beta = beta,
    z = z,
    w = w,
    model_data = model_data
  )

  w = update_w(
    beta = beta,
    z = z,
    pw = pw,
    model_data = model_data
  )

  pw = update_pw(
    w = w,
    model_data = model_data
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
    pw = pw
  )

  return(out)

}


