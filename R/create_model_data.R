create_model_data = function(y,
                             id,
                             time,
                             K,
                             G,
                             M,
                             center = TRUE,
                             scale = FALSE,
                             n_basis = 10,
                             spline_intercept_penalty = 1,
                             spline_penalty = 1,
                             w_dirichlet = 1,
                             z_dirichlet = 1,
                             cusp_nu = 2,
                             cusp_a = 1,
                             cusp_b = 1,
                             cusp_min_var = 0.05) {

  data = data.frame(
    id = id,
    time = time
  ) |>
    dplyr::mutate(id = factor(id, levels = unique(id))) |>
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
  time_seq = as.integer(time - min(time) + 1)
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

  basis_funcions = create_basis_matrix(
    n_basis = n_basis,
    time = time_seq_unique
  )

  B = basis_funcions$model_matrix
  S = basis_funcions$nD
  S = S * spline_penalty
  S[1, 1] = spline_intercept_penalty
  S_expand = kronecker(diag(M), S)

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

  out$clustering = list(
    n_basis = n_basis,
    spline_intercept_penalty = spline_intercept_penalty,
    spline_penalty = spline_penalty,
    B = B,
    S = S,
    S_expand = S_expand,
    w_dirichlet = w_dirichlet,
    z_dirichlet = z_dirichlet
  )

  return(out)

}
