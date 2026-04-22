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
                             z_intercept_penalty = 1,
                             z_spline_penalty = 1,
                             w_dirichlet = 1,
                             cusp_nu = 2,
                             cusp_a = 1,
                             cusp_b = 1,
                             cusp_min_var = 0.05) {

  data = data.frame(
    id = id,
    time = time
  ) |>
    mutate(id = factor(id, levels = unique(id))) |>
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

  idx = rep((1:n_id)- 1, each = n_time)
  idx_unique = 0:(n_id-1)

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

  theta_basis_funcions = create_basis_matrix(
    n_basis = n_basis,
    time = time_seq_unique
  )

  B_theta = theta_basis_funcions$model_matrix
  S_theta = theta_basis_funcions$nD
  S_theta = S_theta * theta_spline_penalty
  S_theta[1, 1] = theta_intercept_penalty
  S_expand = kronecker(diag(M), S_theta)
  noise = 1

  B_n = kronecker(diag(n_id), B_theta)
  S_n = kronecker(diag(n_id), S_theta)

  Ri = B_theta %*% solve(t(B_theta) %*% B_theta + S_theta) %*% t(B_theta)
  # Ri = 0.5 * Ri + (1-0.5) * diag(n_time)

  R = kronecker(diag(n_id), Ri)
  Ry = R %*% y
  RtR_i = crossprod(Ri)
  Ri_inv = solve(Ri + diag(noise, n_time))
  R_inv = kronecker(diag(n_id), Ri_inv)

  RRt_inv_i = solve(RtR_i + diag(noise, n_time))
  RRt_inv = kronecker(diag(n_id), RRt_inv_i)
  R = Matrix::bdiag(R)
  RRt_inv = Matrix::bdiag(RRt_inv)


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
    var_names = colnames(y),
    idx = idx,
    idx_unique = idx_unique
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

  out$cusp = list(
    nu = cusp_nu,
    a = cusp_a,
    b = cusp_b,
    min_var = cusp_min_var
  )

  out$theta_spline = list(
    n_basis = n_basis,
    theta_intercept_penalty = theta_intercept_penalty,
    theta_spline_penalty = theta_spline_penalty,
    B_theta = B_theta,
    S_theta = S_theta,
    Ri = Ri,
    R = R,
    Ry = Ry,
    RtR_i = RtR_i,
    Ri_inv = Ri_inv,
    R_inv = R_inv,
    RRt_inv_i = RRt_inv_i,
    RRt_inv = RRt_inv,
    B_n = B_n,
    S_n = S_n,
    S_expand = S_expand,
    w_dirichlet = w_dirichlet
  )

  return(out)

}
