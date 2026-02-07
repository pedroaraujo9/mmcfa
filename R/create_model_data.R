create_model_data = function(y,
                             id,
                             time,
                             K,
                             G,
                             M,
                             center = TRUE,
                             scale = FALSE,
                             n_basis = 10,
                             theta_intercept_penalty = 1,
                             theta_spline_penalty = 1,
                             z_spline_penalty = 1,
                             w_dirichlet = 1) {

  data = data.frame(
    id = id,
    time = time
  ) |>
    dplyr::bind_cols(as.data.frame(y)) |>
    dplyr::arrange(id, time)

  y = data |> dplyr::select(-id, -time) |> as.matrix()
  id = data$id
  time = data$time

  n = nrow(y)
  J = ncol(y)

  id_unique = unique(id)
  time_unique = unique(time)

  n_id = length(id_unique)
  n_time = length(time_unique)

  order = 1
  time_seq = time - min(time) + 1
  time_seq_unique = unique(time_seq)

  data_mean = colMeans(y)
  data_sd = apply(y, 2, sd)
  y = scale(y, center = center, scale = scale) |> as.matrix()

  id_time_df = tibble::tibble(
    id = id,
    time_seq = time_seq,
    time = time
  ) |>
    dplyr::mutate(id = factor(id, levels = id_unique))

  theta_basis_funcions = gen_basis_matrix(
    n_basis = n_basis,
    time = time_seq_unique
  )

  B_theta = theta_basis_funcions$model_matrix
  S_theta = theta_basis_funcions$nD
  S_theta = S_theta * theta_spline_penalty
  S_theta[1, 1] = theta_intercept_penalty

  R = B_theta %*% solve(t(B_theta) %*% B_theta + S_theta) %*% t(B_theta)

  out = list()

  out$data = list(
    y = y,
    id = id,
    time = time,
    time_seq = time_seq,
    time_unique = time_unique,
    id_unique = id_unique,
    data_mean = data_mean,
    data_sd = data_sd,
    center = center,
    scale = scale,
    var_names = colnames(y)
  )

  out$dims = list(
    n = n,
    J = J,
    n_id = n_id,
    n_time = n_time,
    G = G,
    K = K,
    M = M
  )

  out$theta_spline = list(
    n_basis = n_basis,
    theta_intercept_penalty = theta_intercept_penalty,
    theta_spline_penalty = theta_spline_penalty,
    B_theta = B_theta,
    S_theta = S_theta,
    R = R
  )

  return(out)

}
