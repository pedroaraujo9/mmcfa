generate_R = function(time,
                      n_basis = 10,
                      add_penalty = TRUE,
                      spline_penalty = 1,
                      spline_intercept_penalty = 1) {

  theta_basis_funcions = create_basis_matrix(
    n_basis = n_basis,
    time = unique(time)
  )

  B_theta = theta_basis_funcions$model_matrix
  S_theta = theta_basis_funcions$nD
  S_theta = S_theta * spline_penalty
  S_theta[1, 1] = spline_intercept_penalty

  if(add_penalty == TRUE) {
    R = B_theta %*% solve(t(B_theta) %*% B_theta + S_theta) %*% t(B_theta)
  }else{
    R = B_theta %*% solve(t(B_theta) %*% B_theta) %*% t(B_theta)
  }

  return(R)

}

smooth_theta = function(theta, id, time, n_basis = 10, spline_penalty = 1) {

  id_unique = unique(id)
  n_id = length(id_unique)

  R = generate_R(time = time, n_basis = n_basis, spline_penalty = spline_penalty)

  theta_smooth = theta

  for(i in 1:n_id) {
    id_i = id_unique[i]
    theta_smooth[id == id_i, ] = R %*% theta[id == id_i, ]
  }

  return(theta_smooth)

}

find_init_FA = function(y, time, id, K_max, n_basis = 10) {

  y_centred = scale(y, center = TRUE, scale = TRUE)
  dec = eigen(cov(y_centred))
  V = cbind(dec$vectors[, 1:K_max])

  if(K_max == 1) {

    L = as.matrix(sqrt(dec$values[1:K_max]))

  }else{

    L = diag(sqrt(dec$values[1:K_max]))

  }

  alpha = V %*% L
  theta = y_centred %*% V %*% solve(L)
  theta = smooth_theta(theta, id, time, n_basis = n_basis)

  out = list(
    alpha = alpha,
    theta = theta
  )

  return(out)

}

find_init_z = function(theta, G_max) {

  d = dist(theta)

  asw = lapply(2:G_max, function(g){

    z = kmeans(theta, centers = g, iter.max = 100, nstart = 100)$cluster
    mean(cluster::silhouette(z, d)[, "sil_width"])

  }) %>% do.call(c, .)

  G = 1 + which.max(asw)
  z = kmeans(theta, centers = G, iter.max = 100, nstart = 100)$cluster

  out = list(z = z, G = G, asw = asw)
  return(out)

}

find_init_w = function(z_seq, M_max) {

  d_seq = TraMineR::seqdist(z_seq, method = "HAM") |> suppressMessages() |> as.dist()

  asw_seq = lapply(2:M_max, function(m){

    hc = hclust(d_seq, method = "ward.D")
    w = hc |> cutree(k = m)
    mean(cluster::silhouette(w, d_seq)[, "sil_width"])

  }) %>% do.call(c, .)

  M = 1 + which.max(asw_seq)
  w = hclust(d_seq, method = "ward.D") |> cutree(k = M)
  names(w) = unique(id)

  out = list(w = w, M = M, asw_seq = asw_seq)
  return(out)

}

create_init_values2 = function(y, id, time, K_max, G_max, M_max, n_basis) {

  n_id = length(unique(id))
  n_time = length(unique(time))

  fa_init = find_init_FA(
    y = y, K_max = K_max, time = time, id = id, n_basis = n_basis
  )

  alpha = fa_init$alpha
  theta = fa_init$theta

  # z_init = find_init_z(theta, G_max)
  # G = z_init$G
  # z = z_init$z
  # z_seq = z |>
  #   matrix(nrow = n_id, ncol = n_time, byrow = T) |>
  #   TraMineR::seqdef() |>
  #   suppressMessages()
  #
  # w_init = find_init_w(z_seq = z_seq, M_max = M_max)
  # M = w_init$M
  # w = w_init$w

  out = list(
    alpha = alpha,
    theta = theta #,
    # w = w,
    # z = z,
    # G = G,
    # M = M,
    # asw = z_init$asw,
    # asw_seq = w_init$asw_seq
  )

  return(out)

}
