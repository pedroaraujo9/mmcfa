update_epsilon = function(mu, sigma, alpha, psi, z, model_data) {

  # dimensions
  H = ncol(alpha)
  J = model_data$dims$J
  n_time = model_data$dims$n_time
  n = model_data$dims$n
  n_id = model_data$dims$n_id

  # data
  y = model_data$data$y
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  R = model_data$theta_spline$R
  Rty = model_data$theta_spline$Rty

  # RR = kronecker(diag(H), R)
  # one = matrix(1, H * n_time)

  # variability from data
  inv_psi_matrix = matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  alpha_scaled_psi = alpha * inv_psi_matrix
  alphat_psi_alpha = crossprod(alpha_scaled_psi, alpha)
  RtR = crossprod(R)
  prec_data = kronecker(alphat_psi_alpha, RtR)
  Rty_alpha_scaled = Rty %*% alpha_scaled_psi

  # sampling
  epsilon = matrix(nrow = n, ncol = H)

  for(i in 1:n_id) {

    id_i = id_unique[i]
    y_i = y[id == id_i, ]
    z_i = z[id == id_i]

    Sigma = diag(1/(sigma[z_i]^2))
    prec_prior = kronecker(diag(H), Sigma)

    post_cov = solve(prec_prior + prec_data)
    V = t(chol(post_cov))

    post_center = post_cov %*% as.vector(Sigma %*% mu[z_i, ] + Rty_alpha_scaled[id == id_i, ])
    epsilon_i_vec = post_center + V %*% rnorm(n_time * H)

    # epsilon_i_vec = epsilon_i_vec - post_cov %*% t(RR) %*% (one) %*% (1/(t(one) %*% RR %*% post_cov %*% t(RR) %*% one)) %*% t(one) %*% RR %*% epsilon_i_vec

    epsilon_i = matrix(epsilon_i_vec, nrow = n_time, ncol = H, byrow = FALSE)
    epsilon[id == id_i, ] = epsilon_i
  }

  return(epsilon)

}

update_epsilon2 = function(mu, sigma, alpha, psi, z, model_data) {

  # dimensions
  H = ncol(alpha)
  J = model_data$dims$J
  n_time = model_data$dims$n_time
  n = model_data$dims$n
  n_id = model_data$dims$n_id

  # data
  y = model_data$data$y
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  R = model_data$theta_spline$R
  Rty = model_data$theta_spline$Rty

  # variability from data
  inv_psi_matrix = matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  alpha_scaled_psi = alpha * inv_psi_matrix
  alphat_psi_alpha = crossprod(alpha_scaled_psi, alpha)
  RtR = crossprod(R)
  prec_data = kronecker(alphat_psi_alpha, RtR)
  Rty_alpha_scaled = Rty %*% alpha_scaled_psi

  SIG = matrix((1/sigma[z])^2, nrow = n, ncol = H, byrow = FALSE)
  MU_scaled = mu[z, ] * SIG

  # sampling
  epsilon = matrix(nrow = n, ncol = H)

  for(i in 1:n_id) {

    id_i = id_unique[i]
    f = (id == id_i)

    z_i = z[f]

    Sigma = diag(SIG[f, 1])
    prec_prior = kronecker(diag(H), Sigma)

    epsilon_i_vec = post_epsilon_cpp2(
      prec_prior = prec_prior,
      prec_data = prec_data,
      MU_scaled = MU_scaled[f, ],
      Rty_alpha_scaled =  Rty_alpha_scaled[f, ]
    )

    epsilon_i = matrix(epsilon_i_vec, nrow = n_time, ncol = H, byrow = FALSE)
    epsilon[f, ] = epsilon_i

  }

  return(epsilon)

}

