update_epsilon = function(mu, sigma, alpha, psi, z, model_data, add_cluster = TRUE) {

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
  RtR = model_data$theta_spline$RtR

  # variability from data
  inv_psi_matrix = matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  alpha_scaled_psi = alpha * inv_psi_matrix
  alphat_psi_alpha = crossprod(alpha_scaled_psi, alpha)
  prec_data = kronecker(alphat_psi_alpha, RtR)
  Rty_alpha_scaled = Rty %*% alpha_scaled_psi

  inv_sigma = (1/sigma[z])^2
  inv_SIGMA = matrix(inv_sigma, nrow = n, ncol = H, byrow = FALSE)

  # variability from clusters
  if(add_cluster == TRUE) {

    MU_scaled = mu[z, ] * inv_SIGMA

  }else{

    MU_scaled = matrix(0, nrow = n, ncol = H)

  }


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

    post_center = post_cov %*% as.vector(MU_scaled[z_i, ] + Rty_alpha_scaled[id == id_i, ])
    epsilon_i_vec = post_center + V %*% rnorm(n_time * H)

    # epsilon_i_vec = epsilon_i_vec - post_cov %*% t(RR) %*% (one) %*% (1/(t(one) %*% RR %*% post_cov %*% t(RR) %*% one)) %*% t(one) %*% RR %*% epsilon_i_vec

    epsilon_i = matrix(epsilon_i_vec, nrow = n_time, ncol = H, byrow = FALSE)
    epsilon[id == id_i, ] = epsilon_i
  }

  return(epsilon)

}

update_epsilon2 = function(mu, sigma, alpha, psi, z, model_data, add_cluster = TRUE) {

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
  RtR = model_data$theta_spline$RtR

  # variability from data
  alpha_scaled_psi = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)

  mats = utils_comp_epsilon(
    alpha = alpha,
    alpha_scaled_psi = alpha_scaled_psi,
    RtR = RtR,
    Rty = Rty
  )

  prec_data = mats$prec_data
  Rty_alpha_scaled = mats$Rty_alpha_scaled

  # variability from clusters
  inv_sigma = (1/sigma[z])^2
  inv_SIGMA = matrix(inv_sigma, nrow = n, ncol = H, byrow = FALSE)

  if(add_cluster == TRUE) {


    MU_scaled = mu[z, ] * inv_SIGMA

  }else{

    MU_scaled = matrix(0, nrow = n, ncol = H)

  }

  # sampling
  epsilon = matrix(nrow = n, ncol = H)

  for(i in 1:n_id) {

    id_i = id_unique[i]
    f = (id_i == id)
    z_i = z[f]

    prec_prior = diag(rep(inv_SIGMA[f, 1], times = H))

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

update_epsilon_exp1 = function(mu, sigma, alpha, psi, z, model_data, add_cluster = T) {

  n_time = model_data$dims$n_time
  R = model_data$theta_spline$R
  K = model_data$dims$K
  I = diag(K)
  y = model_data$data$y
  id = model_data$data$id
  id_unique = model_data$data$id_unique

  alpha_scaled_psi = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  A = crossprod(alpha_scaled_psi, alpha)

  I = diag(K)

  St = lapply(1:n_time, function(t){
    solve(I + R[t, t]^2 * A)
  })

  Lt = lapply(1:n_time, function(t){
    t(chol(St[[t]]))
  })

  for(i in 1:n_id) {

    id_i = id_unique[i]

    epsilon_i = epsilon[id == id_i, ]
    mu_i = mu[z[id == id_i], ]
    y_i = y[id == id_i, ]


    for(t in 1:n_time) {

      ct = rbind(R[t, -t]) %*% epsilon_i[-t, , drop = FALSE]
      rt = y_i[t, ] - ct %*% t(alpha)
      S = St[[t]]
      m = S %*% (mu_i[t,] + R[t, t] * t(alpha_scaled_psi) %*% t(rt))
      epsilon_i[t, ] = m + Lt[[t]] %*% rnorm(K)

    }

    epsilon[id == id_i, ] = epsilon_i

  }

  return(epsilon)

}

# microbenchmark::microbenchmark(
#
#   current = update_epsilon2(
#     mu, sigma, alpha, psi, z, model_data, add_cluster = T
#   ),
#
#   ind = update_epsilon_exp1(
#     mu, sigma, alpha, psi, z, model_data, add_cluster = T
#   ),
#
#   times = 100
# )


