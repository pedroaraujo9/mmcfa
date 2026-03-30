update_w = function(beta, z, pw, model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$theta_spline$n_basis

  id_unique = model_data$data$id_unique
  id = model_data$data$id
  id = factor(id, levels = id_unique)
  time_seq = model_data$data$time_seq
  B = model_data$theta_spline$B

  Z = create_dummy(z, G)

  w_post_prob = matrix(nrow = n_id, ncol = M)

  for(m in 1:M) {

    if(m == 1) {

      beta_group = beta[1:n_basis, ]
      prob_group = compute_prob_group(B, beta_group, time_seq-1)

    }else{

      idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
      beta_group = rbind(beta[1, ], beta[idx, ])
      prob_group = compute_prob_group(cbind(1, B), beta_group, time_seq-1)

    }

    log_pz = log(rowSums(Z * prob_group))

    w_post_prob[, m] = exp(fast_aggregate_sum(log_pz, id)[, 1]) * pw[, m]

  }

  w_post_prob = w_post_prob / rowSums(w_post_prob)

  w = as.integer(extraDistr::rcat(n = n_id, prob = w_post_prob))
  names(w) = id_unique
  out = list(w = w, w_post_prob = w_post_prob)

  return(out)

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

