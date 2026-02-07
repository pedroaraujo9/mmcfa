simulate_data = function(seed) {

  set.seed(seed)

  #### define dimensions and hyperparameters ####
  K = 2
  M = 3
  G = 4
  J = 20
  n_id = 30
  n_time = 50
  n = n_id * n_time
  id_unique = paste0("id", 1:n_id)
  id = rep(id_unique, each = n_time)
  time = rep(1:n_time, times = n_id)
  n_basis = 10
  theta_spline_penalty = 1
  theta_intercept_penalty = 1

  #### sample global group ####
  w = rep(1:M, each = n_id / M)
  names(w) = id_unique

  #### sample local group ####
  ZM = matrix(nrow = n_id, ncol = n_time)

  for(i in 1:n_id) {

    if(w[i] == 1) {

      thr = sample(1:10, size = 2)
      ZM[i,] = c(
        rep(1, 20 - thr[1]),
        rep(3, 20 - thr[2]),
        rep(4, 50 - 40 + thr[1] + thr[2])
      )

    }else if(w[i] == 2) {

      if(w[i] == 2) {
        thr = sample(1:10, size = 2)
        ZM[i,] = c(
          rep(1, 20 - thr[1]),
          rep(2, 20 - thr[2]),
          rep(4, 50 - 40 + thr[1] + thr[2])
        )
      }

    }else if(w[i] == 3) {
      thr = sample(1:10, size = 2)
      ZM[i,] = c(
        rep(2, 20 - thr[1]),
        rep(3, 20 - thr[2]),
        rep(3, 50 - 40 + thr[1] + thr[2])
      )
    }
  }

  z = as.vector(t(ZM))

  #### define parameters ####
  mu = rbind(
    "1" = c(-2, -2),
    "2" = c(0, -2),
    "3" = c(-1, 0),
    "4" = c(2, 2)
  )

  sigma = rbind(
    c(1, 1),
    c(0.5, 0.5),
    c(1, 1),
    c(0.5, 0.5)
  )

  epsilon = mu[z, ] + sigma[z, ] * gen_normal_mat(n, K)
  alpha = gen_normal_mat(J, K)
  psi = runif(J, 1, 5)
  psi_matrix  = matrix(psi, nrow = n, ncol = J, byrow = TRUE)

  #### generate smooth latent effets ####
  theta_basis_funcions = create_basis_matrix(
    n_basis = n_basis,
    time = 1:n_time
  )

  B_theta = theta_basis_funcions$model_matrix
  S_theta = theta_basis_funcions$nD
  S_theta = S_theta * theta_spline_penalty
  S_theta[1, 1] = theta_intercept_penalty

  R = B_theta %*% solve(t(B_theta) %*% B_theta + S_theta) %*% t(B_theta)
  theta = kronecker(diag(n_id), R) %*% epsilon

  theta_df = data.frame(theta, z = z, id = id, time = time, w = w[id]) |>
    tidyr::gather(dim, value, -z, -id, -time, -w) |>
    dplyr::mutate(dim = as.integer(gsub("X", "", dim)))

  #### generate data ####
  eta = theta %*% t(alpha)
  y = gen_normal_mat(n, J) * psi_matrix + eta

  #### return ####
  data_sim = list(
    y = y,
    id = id,
    time = time,
    true_params = list(
      z = z,
      w = w,
      mu = mu,
      sigma = sigma,
      theta = theta,
      theta_df = theta_df,
      epsilon = epsilon,
      alpha = alpha,
      psi = psi,
      B_theta = B_theta,
      S_theta = S_theta,
      n_basis = n_basis,
      theta_spline_penalty = theta_spline_penalty,
      theta_intercept_penalty = theta_intercept_penalty
    )
  )

  return(data_sim)

}



