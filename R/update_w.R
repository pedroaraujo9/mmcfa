update_w = function(beta, z, pw, model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$theta_spline$n_basis

  id = model_data$data$id
  id_unique = model_data$data$id_unique
  time_seq = model_data$data$time_seq
  B = model_data$theta_spline$B_theta

  Z = create_dummy(z, G)
  id = factor(id, levels = id_unique)

  ll = matrix(nrow = n_id, ncol = M)

  intercept = beta[1, , drop = FALSE]
  beta = beta[-1, , drop = FALSE]
  B = cbind(1, B)

  for(m in 1:M) {

    idx = ((m-1)*n_basis + 1):(m*n_basis)
    beta_group = beta[idx, , drop = FALSE]
    beta_group = rbind(intercept, beta_group)

    prob_group = compute_prob_group(B, beta_group, time_seq-1)
    log_pz = log(rowSums(Z * prob_group))

    ll[, m] = fast_aggregate_sum(log_pz, id)[, 1] + log(pw[m])

  }

  ll = ll - matrix(
    mclust::logsumexp(ll), nrow = n_id, ncol = M, byrow = F
  )

  w = as.integer(extraDistr::rcatlp(n = n_id, log_prob = ll) + 1)
  names(w) = model_data$data$id_unique

  return(w)

}


update_w2 = function(beta, z, pw, model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id

  id = model_data$data$id
  id_unique = model_data$data$id_unique
  B = model_data$theta_spline$B_theta

  Z = create_dummy(z, G)
  id = factor(id, levels = id_unique)

  ll = matrix(nrow = n_id, ncol = M)

  for(m in 1:M) {

    w = rep(m, n_id)
    W = create_dummy(w, M)
    X = kronecker(W, B)

    prob = mclust::softmax(X %*% beta)
    log_pz = log(rowSums(Z * prob))
    ll[, m] = aggregate(log_pz ~ id, FUN = sum)$log_pz + log(pw[m])

  }

  ll = ll - matrix(
    mclust::logsumexp(ll), nrow = n_id, ncol = M, byrow = F
  )

  w = as.integer(extraDistr::rcatlp(n = n_id, log_prob = ll) + 1)
  names(w) = id_unique

  return(w)
}

