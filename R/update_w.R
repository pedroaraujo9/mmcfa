update_w = function(beta,
                    z,
                    pw,
                    model_data,
                    w_prior) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$clustering$n_basis

  id_unique = model_data$data$id_unique
  id = model_data$data$id
  time_seq = model_data$data$time_seq
  B = model_data$clustering$B

  Z = create_dummy(z, G)

  if(!is.null(w_prior)) {

    pw = w_prior

  }else{

    pw = matrix(as.numeric(pw), nrow = 1, ncol = M, byrow = TRUE)

  }

  log_w_post_prob = matrix(nrow = n_id, ncol = M)

  for(m in 1:M) {

    if(m == 1) {

      beta_group = beta[1:n_basis, ]
      prob_group = compute_prob_group(B, beta_group, time_seq-1L)

    }else{

      idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
      beta_group = rbind(beta[1, ], beta[idx, ])
      prob_group = compute_prob_group(cbind(1, B), beta_group, time_seq-1L)

    }

    log_pz = log(rowSums(Z * prob_group))
    log_pz_id = as.numeric(rowsum(log_pz, id))
    log_w_post_prob[, m] = log_pz_id + log(pw[, m] + 1e-300)

  }

  log_w_post_prob = norm_mat(log_w_post_prob)
  w_post_prob = exp(log_w_post_prob)

  w = sample_cat(log_w_post_prob)
  names(w) = id_unique
  out = list(w = w, w_post_prob = w_post_prob)

  return(out)

}


