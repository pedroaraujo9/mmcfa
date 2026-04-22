update_zw = function(w,
                     theta,
                     mu,
                     sigma,
                     H,
                     pz,
                     beta,
                     model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  n = model_data$dims$n
  n_id = model_data$dims$n_id
  n_time = model_data$dims$n_time
  B = model_data$theta_spline$B_theta

  time = model_data$data$time
  time_seq = model_data$data$time_seq
  y = model_data$data$y
  id_unique = model_data$data$id_unique
  id = model_data$data$id
  R = model_data$theta_spline$R
  idx = as.integer(id)

  sigma_id = sigma[w]
  s = matrix(sigma_id[idx], nrow = n, ncol = H, byrow = F)

  prob = matrix(pz, nrow = n_time, ncol = G, byrow = T)

  z_w = matrix(nrow = M, ncol = n_time)

  for(m in 1:M) {

    nw = sum(w == m)

    if(nw > 0) {

      filter = id %in% id_unique[w == m]
      n_zw = sum(filter)

      ll = lapply(1:G, function(g){

        mu_g = matrix(mu[g, ], nrow = n_zw, ncol = H, byrow = TRUE)
        logp = dnorm(theta[filter, ], mean = mu_g, sd = s[filter], log = TRUE)
        rowsum(rowSums(logp), time[filter])

      }) %>% do.call(cbind, .)

      z_w[m, ] = sample_cat(ll + log(prob))
      accept = TRUE

    }else{

      z_w[m, ] = sample(1:G, size = n_time, replace = TRUE)

    }

  }

  out = list(
    z_w = z_w,
    z = as.numeric(t(z_w[w, ]))
  )

  return(out)


}

update_wz = function(beta, z_w, pw, theta, mu, sigma, model_data, smooth) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$theta_spline$n_basis

  id_unique = model_data$data$id_unique
  id = model_data$data$id
  id = factor(id, levels = id_unique)
  time_seq = model_data$data$time_seq
  B = model_data$theta_spline$B
  R = model_data$theta_spline$R

  w_post_prob = matrix(nrow = n_id, ncol = M)

  for(m in 1:M) {


    z_m = rep(z_w[m, ], times = n_id)
    if(smooth == TRUE) {
      MN = compute_kernel(theta, M = R %*% mu[z_m, ], S = sigma[m], id = id)
    }else{
      MN = compute_kernel(theta, M = mu[z_m, ], S = sigma[m], id = id)
    }

    log_pz_id = as.numeric(MN)
    w_post_prob[, m] = log_pz_id + log(pw[m])

  }

  log_w_post_prob = norm_mat(w_post_prob)

  w = sample_cat(log_w_post_prob)
  names(w) = id_unique
  out = list(w = w, w_post_prob = exp(log_w_post_prob))

  return(out)

}

