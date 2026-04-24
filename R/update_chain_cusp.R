update_chain_cusp = function(H,
                             H_max,
                             alpha,
                             alpha_precision,
                             psi,
                             sigma,
                             omega,
                             v,
                             model_data,
                             adapt_H = FALSE,
                             est_sigma = TRUE) {

  n = model_data$dims$n
  G = model_data$dims$G

  # No clustering in this chain update: use fixed placeholders.
  z_dummy = rep(1, n)
  mu_dummy = matrix(0, nrow = max(1, G), ncol = H_max)

  theta = update_theta(
    H = H,
    mu = mu_dummy[, 1:H],
    sigma = sigma,
    alpha = alpha[, 1:H],
    psi = psi,
    z = z_dummy,
    model_data = model_data,
    add_cluster = FALSE
  ) |> complete_dim(dimension = H_max)

  theta_sd = apply(theta, 2, sd)
  theta = scale(theta, center = FALSE, scale = theta_sd)

  if(est_sigma == TRUE) {

    sigma = update_sigma(
      H = H,
      theta = theta[, 1:H],
      mu = mu_dummy[, 1:H],
      z = z_dummy,
      model_data = model_data,
      add_cluster = FALSE
    )

  }

  alpha = update_alpha(
    H = H,
    theta = cbind(theta[, 1:H]),
    psi = psi,
    prior_precision = alpha_precision[1:H],
    model_data = model_data
  ) |> complete_dim(dimension = H_max)

  alpha = scale(alpha, center = FALSE, scale = 1/theta_sd)

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
  omega = out_cusp$omega[1:H] |> complete_dim(dimension = H_max)
  v = out_cusp$v[1:H] |> complete_dim(dimension = H_max)
  alpha_precision = out_cusp$prec |> complete_dim(dimension = H_max)

  if(adapt_H == TRUE) {

    out_adapt = adapt_cusp(
      H_active = H_active,
      H = H,
      H_max = H_max,
      ind = ind,
      alpha = alpha,
      mu = mu_dummy,
      z = z_dummy,
      theta = theta,
      sigma = sigma,
      omega = omega,
      alpha_precision = alpha_precision,
      v = v,
      nu = model_data$cusp$nu,
      min_var = model_data$cusp$min_var,
      model_data = model_data,
      add_cluster = FALSE
    )

    alpha = out_adapt$alpha
    theta = out_adapt$theta
    omega = out_adapt$omega
    alpha_precision = out_adapt$alpha_precision
    v = out_adapt$v
    H = out_adapt$H

  }

  psi = update_psi(
    theta = cbind(theta[, 1:H]),
    alpha = cbind(alpha[, 1:H]),
    model_data = model_data,
    a = 2,
    b = 2
  )

  out = list(
    theta = theta,
    alpha = alpha,
    psi = psi,
    sigma = sigma,
    alpha_precision = alpha_precision,
    omega = omega,
    v = v,
    H_active = H_active,
    H = H
  )

  return(out)

}
