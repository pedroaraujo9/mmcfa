update_alpha_rep = function(theta, psi, model_data, prior_precision) {

  H = ncol(theta)
  J = model_data$dims$J
  y = model_data$data$y
  S0 = diag(prior_precision[1, ])
  Psi = diag(1/psi)

  OtO = crossprod(theta)
  V = solve(kronecker(OtO, Psi) + kronecker(S0, diag(J)))

  m = V %*% as.vector(Psi %*% t(y) %*% theta)
  alpha = m + t(chol(V)) %*% rnorm(J*H)
  alpha = matrix(alpha, nrow = J, ncol = H, byrow = FALSE)

  return(alpha)
}

update_theta_new = function(mu,
                            sigma,
                            alpha,
                            psi,
                            z,
                            H,
                            model_data,
                            sigma_state,
                            w,
                            clust_var,
                            smooth) {

  # dimensions
  J = model_data$dims$J
  n_time = model_data$dims$n_time
  n = model_data$dims$n
  n_id = model_data$dims$n_id
  G = model_data$dims$G

  # data
  y = model_data$data$y
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  y = model_data$data$y
  Ri = model_data$theta_spline$Ri

  I_time = diag(n_time)
  I_H = diag(H)

  # variability from data
  AS = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  tASA = crossprod(AS, alpha)
  YAS = y %*% AS

  theta = matrix(nrow = n, ncol = H)

  if(clust_var == TRUE) {
    tau = 1/(sigma_state^2)
  }else{
    tau = 1/(sigma^2)
  }

  for(i in 1:n_id) {

    filter = id == id_unique[i]
    zi = z[filter]

    if(smooth == TRUE) {
      Mi = Ri %*% mu[zi, ]
    }else{
      Mi = mu[zi, ]
    }

    if(clust_var == FALSE) {

      tau_i = tau[i]

      V = solve(tASA + I_H * tau_i)

      YASi = YAS[filter, , drop = FALSE]
      mi = (YASi + tau_i * Mi) %*% V

      theta[filter, ] = mi +  gen_normal_mat(n_time, H) %*% chol(V)
    }else{


      if(smooth == TRUE) {
        Mi = Ri %*% mu[zi, ]
      }else{
        Mi = mu[zi, ]
      }

      YASi = YAS[filter, , drop = FALSE]
      noise = gen_normal_mat(n, H)

      for(g in unique(zi)) {

        zig = zi == g
        n_zig = sum(zig)
        tau_i = tau[i, g]
        V = solve(tASA + I_H * tau_i)

        mi = (YASi[zig, ] + tau_i * Mi[zig, ]) %*% V
        theta[filter, ][zig, ] = mi +  noise[filter, ][zig, ] %*% chol(V)

      }

    }

  }

  return(theta)

}

update_epsilon_new = function(mu,
                              sigma,
                              alpha,
                              psi,
                              z,
                              H,
                              model_data,
                              add_cluster,
                              mixscat_prior,
                              w,
                              smooth) {

  # dimensions
  J = model_data$dims$J
  n_time = model_data$dims$n_time
  n = model_data$dims$n
  n_id = model_data$dims$n_id
  G = model_data$dims$G

  # data
  y = model_data$data$y
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  y = model_data$data$y
  Ri = model_data$theta_spline$Ri
  Ry = model_data$theta_spline$Ry
  RtR = model_data$theta_spline$RtR

  I_time = diag(n_time)
  I_H = diag(H)

  if(mixscat_prior == TRUE) {
    sigma = sigma[w]
  }

  # variability from data
  AS = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  tASA = crossprod(AS, alpha)
  RYAS = Ry %*% AS
  V1 = kronecker(tASA, RtR)

  epsilon = matrix(nrow = n, ncol = H)
  tau = 1/(sigma^2)

  for(i in 1:n_id) {

    filter = id == id_unique[i]
    zi = z[filter]
    tau_i = tau[i]

    if(smooth == TRUE) {
      Mi = mu[zi, ]
    }else{
      Mi = mu[zi, ]
    }

    V0 = kronecker(I_H,  I_time * tau_i)
    V = solve(V1 + V0)

    RYASi = RYAS[filter, , drop = FALSE]
    mi = V %*% as.vector(RYASi + tau_i * Mi)

    epsilon[filter, ] = matrix(mi +  t(chol(V)) %*% gen_normal_mat(n_time*H, 1),
                               nrow = n_time, ncol = H, byrow = FALSE)

  }

  return(epsilon)

}

update_mu_new = function(z,
                         sigma,
                         theta,
                         H,
                         model_data,
                         w,
                         mixscat_prior,
                         sigma_state,
                         clust_var,
                         sigma_mu,
                         smooth) {

  G = model_data$dims$G
  n = model_data$dims$n
  R = model_data$theta_spline$R
  idx = as.integer(model_data$data$id)


  if(clust_var == TRUE) {

    sigma = sigma_state[cbind(idx, z)]
    tau = 1/(sigma^2)

  }else{

    tau = 1/(sigma[idx]^2)

  }

  L_inv = diag(tau)


  Z = create_dummy(z, G)

  if(smooth == TRUE) {

    A = R %*% Z

  }else{

    A = Z

  }

  V0 = diag(1/(sigma_mu^2))
  tAL = crossprod(A, L_inv)
  V1 = tAL %*% A
  V = solve(V0 + V1)

  m = V %*% tAL %*% theta
  mu = as.matrix(m + t(chol(V)) %*% gen_normal_mat(G, H))

  return(mu)

}

update_sigma_mu = function(mu, model_data, H){

  G = nrow(mu)
  sse = rowSums(mu^2)
  tau = rgamma(G, shape = 2 + 0.5*H, rate = 2 + 0.5*sse)
  sigma_mu = 1/sqrt(tau)
  return(sigma_mu)

}

update_sigma_new = function(H, theta,
                            mu,
                            z,
                            model_data,
                            smooth,
                            w) {

  G = model_data$dims$G
  n = model_data$dims$n
  id = model_data$data$id
  n_id = model_data$dims$n_id
  n_time = model_data$dims$n_time
  R = model_data$theta_spline$R
  Ri = model_data$theta_spline$Ri
  M = model_data$dims$M

  if(smooth == TRUE) {
    se = as.matrix(theta - R %*% mu[z, ])^2
  }else{
    se = (theta - mu[z, ])^2
  }

  if(mixscat_prior == FALSE) {

    sse_id = rowSums(rowsum(as.matrix(se), id))
    tau = rgamma(n_id, shape = 2 + 0.5*n_time*H, rate = 2 + 0.5*sse_id)


  }else{

    w_vec = w[id]
    n_w = table(factor(w, levels = 1:M)) %>% as.numeric()
    see_cluster = rowsum(se, w_vec) %>% rowSums()
    see_cluster = see_cluster[as.character(1:M)]
    see_cluster[is.na(see_cluster)] = 0
    tau = rgamma(M, shape = 2 + 0.5*n_w*n_time*H, rate = 2 + 0.5*see_cluster)

  }

  sigma = 1/sqrt(tau)

  return(sigma)

}

update_sigma_clust = function(H,
                              theta,
                              mu,
                              z,
                              est_epsilon,
                              model_data,
                              smooth,
                              w) {

  G = model_data$dims$G
  n = model_data$dims$n
  id = model_data$data$id
  n_id = model_data$dims$n_id
  n_time = model_data$dims$n_time
  R = model_data$theta_spline$R
  Ri = model_data$theta_spline$Ri
  M = model_data$dims$M
  id_unique = model_data$data$id_unique

  if(smooth == TRUE) {

    se = as.matrix(theta - R %*% mu[z, ])^2

  }else{

    se = (theta - mu[z, ])^2

  }

  zf = factor(z, levels = 1:G)
  n_id_g = as.matrix(table(id, zf))

  sigma_state = lapply(1:G, function(g){

    g_filter = z == g

    if(sum(g_filter) > 0) {

      see_g = rowsum(x = rbind(se[g_filter, ]), group = id[g_filter]) %>% rowSums(na.rm = TRUE)
      see_g = see_g[as.character(id_unique)]
      see_g[is.na(see_g)] = 0
      tau = rgamma(n_id, shape = 2 + 0.5*n_id_g[, g]*H, rate = 2 + 0.5*see_g)

    }else{

      tau = rgamma(n_id, shape = 2, rate = 2)

    }

    1/sqrt(tau)


  }) %>% do.call(cbind, .)

  return(sigma_state)

}

update_z_new = function(z,
                        w,
                        theta,
                        mu,
                        sigma,
                        H,
                        pz,
                        beta,
                        model_data,
                        mixscat_prior,
                        add_sigma_w,
                        time_pz,
                        clust_var,
                        sigma_state,
                        logP_proposal = NULL,
                        global_update = TRUE) {

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

  if(mixscat_prior == TRUE) {

    prob = compute_probs(w = w, M = M, B = B, beta = beta)

  }else{

    prob = matrix(pz, nrow = n, ncol = G, byrow = T)

  }

  idx = as.integer(id)


  if(clust_var == FALSE) {
    s = matrix(sigma[idx], nrow = n, ncol = H, byrow = F)

  }else{
    s = matrix(sigma_state[cbind(idx, z)], nrow = n, ncol = H, byrow = F)
  }

  eta = 1

  if(!is.null(logP_proposal)) {

    if(global_update == TRUE) {

      # ll = lapply(1:G, function(g){
      #   mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
      #   logp = dnorm(theta, mean = mu_g, sd = eta * s, log = TRUE)
      #   rowSums(logp)
      # }) %>% do.call(cbind, .)

      # logP_proposal = norm_mat(ll + log(prob))

      logP_proposal = logP_proposal * eta
      z_new = sample_cat(logP_proposal)

      #R_new = compute_wR(z_new, logP = logP_proposal, model_data = model_data)
      #R_current = compute_wR(z, logP = logP_proposal, model_data = model_data)

      p_new = compute_kernel(theta, M = R %*% mu[z_new, ], S = s, id = id)
      p_old = compute_kernel(theta, M = R %*% mu[z, ], S = s, id = id)
      p_new = p_new + rowsum(log(prob[cbind(1:n, z_new)]), id)
      p_old = p_old + rowsum(log(prob[cbind(1:n, z)]), id)

      q_new = rowsum(logP_proposal[cbind(1:n, z_new)], id)
      q_old = rowsum(logP_proposal[cbind(1:n, z)], id)
      # q_new = compute_kernel(theta, M = mu[z_new, ], S = s, id = id)
      # q_old = compute_kernel(theta, M = mu[z, ], S = s, id = id)

      log_accept_prob = as.numeric((p_new - p_old) + (q_old - q_new))

      accept = log(runif(n_id)) < log_accept_prob
      z[accept[idx]] = z_new[accept[idx]]
      logP = NULL

    }else{

      z_new = sample_cat(logP_proposal)
      accept_vec = numeric(n)

      for(t in 1:n_time) {

        filter = t == time_seq

        if(any(z_new[filter] != z[filter])) {

          z_prop = z
          z_prop[filter] = z_new[filter]

          z_prop_t = z_prop[filter]
          z_curr_t = z[filter]
          logP_t = logP_proposal[filter, ]

          p_new = compute_kernel(theta, M = R %*% mu[z_prop, ], S = s, id = id)
          p_old = compute_kernel(theta, M = R %*% mu[z, ], S = s, id = id)
          p_new = p_new + rowsum(log(prob[cbind(1:n, z_prop)]), id)
          p_old = p_old + rowsum(log(prob[cbind(1:n, z)]), id)

          q_new = logP_t[cbind(1:n_id, z_prop_t)]
          q_old = logP_t[cbind(1:n_id, z_curr_t)]

          log_accept_prob = as.numeric((p_new - p_old) + (q_old - q_new))
          accept = log(runif(n_id)) < log_accept_prob
          z[filter][accept] = z_prop[filter][accept]

          accept_vec[filter] = accept
          logP = NULL

        }

      }

      accept = accept_vec


    }

  }else{

    logP = lapply(1:G, function(g){

      mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
      logp = dnorm(theta, mean = mu_g, sd = s, log = TRUE)
      rowSums(logp)

    }) %>% do.call(cbind, .)

    logP = norm_mat(logP + log(prob))

    if(any((is.na(logP) | is.infinite(logP)))) {
      print("problem is in Z")
    }
    z = sample_cat(logP)
    accept = TRUE
    log_accept_prob = rep(0, n)

  }

  out = list(
    z = z,
    accept = accept,
    log_accept_prob = log_accept_prob,
    logP = logP
  )

  return(out)


}

# update_z_new = function(theta,
#                         mu,
#                         z,
#                         sigma,
#                         pz,
#                         model_data,
#                         smooth) {
#
#   H = ncol(theta)
#   G = model_data$dims$G
#   M = model_data$dims$M
#   n = model_data$dims$n
#   y = model_data$data$y
#   id_unique = model_data$data$id_unique
#   id = model_data$data$id
#   n_id = model_data$dims$n_id
#   n_time = model_data$dims$n_time
#
#   R = model_data$theta_spline$R
#   idx = id %>% as.integer()
#
#   ll = lapply(1:G, function(g){
#
#     mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#     s_g = matrix(sigma[idx], nrow = n, ncol = H, byrow = TRUE)
#     logp = dnorm(theta, mean = mu_g, sd = s_g, log = TRUE)
#     rowSums(logp)
#
#   }) %>% do.call(cbind, .)
#
#   prob = matrix(pz, nrow = n, ncol = G, byrow = T)
#   ll = norm_mat(ll + log(prob))
#
#   z = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)
#
#   return(z)
#
# }

# ll = lapply(1:G, function(g){
#   mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#   logp = dnorm(theta, mean = mu_g, sd = eta * s, log = TRUE)
#   rowSums(logp)
# }) %>% do.call(cbind, .)



# z_current = z
# accept_vec = numeric(n)
#
# for(t in 1:n_time) {
#
#   filter = t == time_seq
#
#   # Compare against the rolling state, not the original state
#   if(any(z_new[filter] != z_current[filter])) {
#
#     z_prop = z_current
#     z_prop[filter] = z_new[filter]
#
#     z_prop_t = z_prop[filter]
#     z_curr_t = z_current[filter]
#     logP_t = logP[filter, ]
#
#     p_new = compute_kernel(theta, M = R %*% mu[z_prop, ], S = s, id = id) + rowsum(log(pz[z_prop]), id)
#     p_old = compute_kernel(theta, M = R %*% mu[z_current, ], S = s, id = id) + rowsum(log(pz[z_current]), id)
#
#     q_new = logP_t[cbind(1:n_id, z_prop_t)]
#     q_old = logP_t[cbind(1:n_id, z_curr_t)]
#
#     log_accept_prob = as.numeric((p_new - p_old) + (q_old - q_new))
#     accept = log(runif(n_id)) < log_accept_prob
#
#     z_current[filter][accept] = z_prop[filter][accept]
#     accept_vec[filter] = accept
#
#   } else {
#     accept_vec[filter] = 1
#   }
#
# }
#
# z = z_current



# p_new = compute_kernel(theta, M = R %*% mu[z_new, ], S = s, id = id)
# p = compute_kernel(theta, M = R %*% mu[z, ], S = s, id = id)
# p_new = p_new + rowsum(log(pz[z_new]), id)
# p = p + rowsum(log(pz[z]), id)
#
# q_new = rowsum(logP[cbind(1:n, z_new)], id)
# q = rowsum(logP[cbind(1:n, z)], id)
#
# log_accept_prob = as.numeric((p_new - p) + (q - q_new))
# # log_accept_prob = as.numeric((p_new - p))
# accept = log(runif(n_id)) < log_accept_prob
# z[accept[idx]] = z_new[accept[idx]]



# ll = lapply(1:G, function(g){
#   mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#   logp = dnorm(theta, mean = mu_g, sd = s, log = TRUE)
#   rowSums(logp)
# }) %>% do.call(cbind, .)
#
# prob = matrix(pz, nrow = n, ncol = G, byrow = TRUE)
#
# ll = norm_mat(ll + log(prob))
# P_indep = exp(ll)
#
# rho = 0.25
# stay_mask = runif(n_id) < rho
#
# z_drawn = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)
# z_new = z
# z_new[!stay_mask[idx]] = z_drawn[!stay_mask[idx]]
#
# P_indep_prop = rowSums(create_dummy(z_new, G) *  P_indep)
# P_indep_current = rowSums(create_dummy(z, G) *  P_indep)
#
# q_proposal = log(rho * (z_new == z) + (1 - rho) * P_indep_prop)
# q_current = log(rho * (z == z_new) + (1 - rho) * P_indep_current)
#
# E_proposal = rowSums(dnorm(as.matrix(theta - R %*% mu[z_new, ]), sd = s, log = TRUE))
# E_current  = rowSums(dnorm(as.matrix(theta - R %*% mu[z, ]), sd = s, log = TRUE))
#
# prior_proposal = log(prob)[cbind(1:n, z_new)]
# prior_current = log(prob)[cbind(1:n, z)]
#
# E_proposal = E_proposal + prior_proposal
# E_current  = E_current + prior_current
#
# log_alpha = rowsum(E_proposal - E_current, id) + rowsum(q_current - q_proposal, id)
#
# accept = (log(runif(n_id)) < log_alpha)[, 1]
# accept_vec = accept[idx]
#
# # Update state
# z[accept_vec] = z_new[accept_vec]




# if(smooth == TRUE) {
#
#   ll = lapply(1:G, function(g){
#
#     mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#     logp = dnorm(theta, mean = mu_g, sd = s, log = TRUE)
#     rowSums(logp)
#
#   }) %>% do.call(cbind, .)
#
#   prob = matrix(pz, nrow = n, ncol = G, byrow = T)
#   ll = norm_mat(ll + log(prob))
#
#   z_new = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)
#
#   E_proposal = rowSums(dnorm(as.matrix(theta - R %*% mu[z_new, ]), sd = s, log = T))
#   E_current = rowSums(dnorm(as.matrix(theta - R %*% mu[z, ]), sd = s, log = T))
#
#   D_proposal = rowSums(dnorm(theta - mu[z_new, ], sd = s, log = T))
#   D_current = rowSums(dnorm(theta - mu[z, ], sd = s, log = T))
#
#   log_alpha = rowsum(E_proposal - E_current, id) + rowsum(D_current - D_proposal, id)
#
#   accept = (log(runif(n_id)) < log_alpha)[, 1]
#   accept_vec = accept[idx]
#   z[accept_vec] = z_new[accept_vec]
#
# }





# Z = create_dummy(z, G)
#
# for(t in 1:n_time) {
#
#   filter = t == time
#   ll = matrix(nrow = n_id, ncol = G)
#
#   for(g in 1:G) {
#
#     Z[filter, ] = 0
#     Z[filter, g] = 1
#     E = (theta - crossprod(R, Z) %*% mu)
#     ll[, g] = rowSums(rowsum(dnorm(as.matrix(E), sd = s, log = T), group = id))
#
#   }
#
#   ll = norm_mat(ll + log(prob))
#   z[filter] = as.integer(extraDistr::rcatlp(n = n_id, log_prob = ll) + 1)
#   Z = create_dummy(z, G)
#
# }








#
#
#
#
# update_z_new2 = function(theta,
#                          mu,
#                          z,
#                          sigma,
#                          pz,
#                          model_data,
#                          smooth) {
#
#   H = ncol(theta)
#   G = model_data$dims$G
#   M = model_data$dims$M
#   n = model_data$dims$n
#   y = model_data$data$y
#   id_unique = model_data$data$id_unique
#   id = model_data$data$id
#   n_id = model_data$dims$n_id
#   n_time = model_data$dims$n_time
#
#   R = model_data$theta_spline$R
#   idx = id %>% as.integer()
#   time = model_data$data$time_seq
#
#   ll = lapply(1:G, function(g){
#
#     mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#     s_g = matrix(sigma[idx], nrow = n, ncol = H, byrow = TRUE)
#     logp = dnorm(theta, mean = mu_g, sd = s_g, log = TRUE)
#     rowSums(logp)
#
#   }) %>% do.call(cbind, .)
#
#   prob = matrix(pz, nrow = n, ncol = G, byrow = T)
#   ll = norm_mat(ll + log(prob))
#
#   z = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)
#
#   return(z)
#
# }
#
#
# update_mu_rep2 = function(z,
#                           sigma,
#                           epsilon,
#                           sigma_theta,
#                           alpha,
#                           psi,
#                           model_data,
#                           smooth) {
#
#   G = model_data$dims$G
#   H = model_data$dims$K
#
#   mu = matrix(0, nrow = G, ncol = H)
#
#   for(g in 1:G) {
#
#     ng = sum(z == g)
#
#     if(ng == 0) {
#
#       mu[g, ] = rnorm(H, mean = 0, sd = 1)
#
#     } else {
#
#       epsilon_bar = colMeans(epsilon[z == g, , drop = FALSE])
#       tau = 1/(sigma[g]^2)
#       nu = 1 / (1 + ng * tau)
#       m = nu * (ng * tau) * epsilon_bar
#       mu[g, ] = rnorm(H, m, sqrt(nu))
#
#     }
#   }
#
#   return(mu)
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# update_mu_sigma_ng = function(epsilon,
#                               z,
#                               model_data,
#                               m0 = 0,
#                               kappa0 = 1,
#                               a0 = 2,
#                               b0 = 2) {
#
#   G = model_data$dims$G
#   H = ncol(epsilon)
#
#   mu = matrix(0, G, H)
#   sigma = rep(0, G)
#
#   for(g in 1:G) {
#
#     idx = which(z == g)
#     ng = length(idx)
#
#     if(ng == 0) {
#
#       tau = rgamma(1, shape = a0, rate = b0)
#
#       mu[g, ] = rnorm(H, m0, sqrt(1/(kappa0 * tau)))
#       sigma[g] = 1 / sqrt(tau)
#
#     } else {
#
#       x = epsilon[idx, , drop = FALSE]
#
#       xbar = colMeans(x)
#
#       # cluster-level pooled SS (ALL dimensions share variance)
#       ss = sum((x - matrix(xbar, ng, H, byrow = TRUE))^2)
#
#       kappa_n = kappa0 + ng
#       m_n = (kappa0 * m0 + ng * xbar) / kappa_n
#
#       a_n = a0 + (ng * H) / 2
#       b_n = b0 + 0.5 * ss +
#         (kappa0 * ng * sum((xbar - m0)^2)) / (2 * kappa_n)
#
#       tau = rgamma(1, shape = a_n, rate = b_n)
#
#       mu[g, ] = rnorm(H, m_n, sqrt(1/(kappa_n * tau)))
#       sigma[g] = 1 / sqrt(tau)
#     }
#   }
#
#   list(mu = mu, sigma = sigma)
# }





# update_epsilon_rep2 = function(mu,
#                                sigma,
#                                sigma_theta,
#                                alpha,
#                                psi,
#                                z,
#                                model_data,
#                                add_cluster,
#                                smooth) {
#
#   # dimensions
#   H = ncol(alpha)
#   J = model_data$dims$J
#   n_time = model_data$dims$n_time
#   n = model_data$dims$n
#   n_id = model_data$dims$n_id
#   G = model_data$dims$G
#
#   # data
#   y = model_data$data$y
#   id = model_data$data$id
#   id_unique = model_data$data$id_unique
#   I_time = diag(n_time)
#   I_H = diag(H)
#   Ri = model_data$theta_spline$Ri
#   Li_inv = model_data$theta_spline$Li_inv
#
#   # variability from data
#   AS = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
#   tASA = crossprod(AS, alpha)
#   YAS = y %*% AS
#
#   # sampling
#   if(smooth == TRUE) {
#
#     V1 = kronecker(tASA,  I_time)
#
#     if(add_cluster == FALSE) {
#
#       V0 = kronecker(I_H, Li_inv)
#       V = solve(V0 + V1)
#       V_chol = t(chol(V))
#
#     }
#
#     epsilon = matrix(nrow = n, ncol = H)
#
#     for(i in 1:n_id) {
#
#       filter = id == id_unique[i]
#       y_i = y[filter, ]
#       z_i = z[filter]
#       YAS_i = YAS[filter, ]
#
#       if(add_cluster == TRUE) {
#
#         Si = diag(1/sigma[z_i]^2)
#         Li_inv = Si + Li_inv
#
#         V0 = kronecker(I_H, Li_inv)
#         V = solve(V0 + V1)
#         V_chol = t(chol(V))
#
#         M_i = Li_inv %*% Ri %*% mu[z_i, ]
#
#       }else{
#
#         M_i = matrix(0, nrow = n_time, ncol = H)
#
#       }
#
#       m = V %*% as.vector(YAS_i + M_i)
#       epsilon_i_vec = m + V_chol %*% rnorm(n_time * H)
#       epsilon[filter, ] = matrix(epsilon_i_vec, nrow = n_time, ncol = H, byrow = FALSE)
#
#     }
#
#   }else{
#
#     if(add_cluster == TRUE) {
#
#       epsilon = matrix(nrow = n, ncol = H)
#
#       for(g in 1:G) {
#
#         filter = z == g
#         n_g = sum(filter)
#
#         if(n_g > 0) {
#
#           YAS_g = YAS[filter, , drop = FALSE]
#           M_g = matrix(mu[g, ], n_g, H, byrow = TRUE)
#           V0 = diag(1/(sigma[g, ]^2))
#           V = solve(tASA + V0)
#           m = (YAS_g + M_g %*% V0) %*% V
#
#           epsilon[filter, ] = m + gen_normal_mat(n_g, H) %*% t(chol(V))
#
#         }
#
#       }
#
#
#     }else{
#
#       M = matrix(0, nrow = n, H)
#
#       V0 = diag(1/(sigma_theta^2))
#       V = solve(tASA + V0)
#       m = YAS %*% V
#
#       epsilon = m + gen_normal_mat(n, H) %*% t(chol(V))
#
#     }
#
#   }
#
#   return(epsilon)
#
# }
#
#
# update_mu_rep2 = function(z,
#                           sigma,
#                           epsilon,
#                           sigma_theta,
#                           alpha,
#                           psi,
#                           model_data,
#                           smooth) {
#
#   G = model_data$dims$G
#   H = model_data$dims$K
#
#   mu = matrix(0, nrow = G, ncol = H)
#
#   for(g in 1:G) {
#
#     ng = sum(z == g)
#
#     if(ng == 0) {
#
#       mu[g, ] = rnorm(H, mean = 0, sd = sigma[g, ])
#
#     } else {
#
#       epsilon_bar = colMeans(epsilon[z == g, , drop = FALSE])
#       tau = 1/(sigma[g, ]^2)
#       nu = 1 / ((1 + ng) * tau)
#       m = nu * (ng * tau) * epsilon_bar
#       mu[g, ] = rnorm(H, m, sqrt(nu))
#
#     }
#   }
#
#   return(mu)
#
# }
#
# # update_mu_rep = function(z, sigma, epsilon, sigma_theta, alpha, psi, model_data, smooth) {
# #
# #   # dimensions
# #   H = ncol(alpha)
# #   J = model_data$dims$J
# #   n_time = model_data$dims$n_time
# #   n = model_data$dims$n
# #   n_id = model_data$dims$n_id
# #   G = model_data$dims$G
# #
# #   L_inv = model_data$theta_spline$L_inv
# #   R = model_data$theta_spline$R
# #   Z = create_dummy(z, G)
# #
# #   # S = matrix(1/(sigma[z]^2), nrow = n, ncol = n, byrow = TRUE)
# #   if(smooth == TRUE) {
# #
# #     S = diag(1/(sigma[z]^2))
# #     L_inv = S + L_inv
# #
# #     RZ = R %*% Z
# #     tRZL = crossprod(RZ, L_inv)
# #     V1 = tRZL %*% RZ
# #
# #     V = solve(V1 + diag(G))
# #     m = V %*% tRZL %*% epsilon
# #     mu = m + t(chol(V)) %*% gen_normal_mat(G, H)
# #
# #   }else{
# #
# #     U = diag(1/(sigma[z]^2))
# #     V0 = diag(1/(sigma_theta^2))
# #
# #     tZU = crossprod(Z, U)
# #     V1 = tZU %*% Z
# #
# #     V = solve(kronecker(V0, V1) + diag(G*H))
# #     m = V %*% as.vector(tZU %*% epsilon %*% V0)
# #     mu = m + t(chol(V)) %*% rnorm(G*H)
# #     mu = matrix(mu, nrow = G, ncol = H, byrow = FALSE)
# #
# #   }
# #
# #   return(mu)
# #
# # }
#
# mf_update_z = function(theta, R, mu, sigma, n, G, max_iter = 100, tol = 1e-6) {
#   q_z = matrix(1/G, nrow = n, ncol = G)
#   R_diag = diag(R)
#   R_inv = solve(R)
#   diag(R) = 0
#
#   for (iter in 1:max_iter) {
#     q_z_old = q_z
#
#     E_mu = q_z %*% mu
#     neighbor_mu = R %*% E_mu
#
#     log_q = matrix(0, nrow = n, ncol = G)
#     for (g in 1:G) {
#       scaled_cross = (theta - neighbor_mu) %*% mu[g, ] / sigma[g]^2
#       scaled_quad  = 0.5 * R_diag * sum(mu[g, ]^2) / sigma[g]^2
#       var_penalty  = log(sigma[g])
#       resid_g      = theta - outer(rep(1, n), mu[g, ])
#       resid_prec   = 0.5 * rowSums((R_inv %*% resid_g) * resid_g) / sigma[g]^2
#
#       log_q[, g] = scaled_cross - scaled_quad - var_penalty - resid_prec
#     }
#
#     log_q = log_q - apply(log_q, 1, max)
#     q_z = exp(log_q) / rowSums(exp(log_q))
#
#     if (max(abs(q_z - q_z_old)) < tol) break
#   }
#
#   return(q_z)
# }
#
# update_pi = function(theta, R, L_inv, mu, sigma, q){
#
#   n = nrow(theta)
#   G = nrow(mu)
#   K = ncol(theta)
#
#   mu_bar = q %*% mu
#
#   R_mu_bar = R %*% mu_bar
#   R_diag = diag(R)
#
#   m_bar = R_mu_bar - R_diag * mu_bar
#
#   Rinv_theta = L_inv %*% theta
#
#   log_q = matrix(0, n, G)
#
#   for(g in 1:G){
#
#     mu_g = matrix(mu[g, ], n, K, byrow=TRUE)
#
#     diff1 = theta - R_diag * mu_g - m_bar
#     diff2 = Rinv_theta - mu_g
#
#     A = rowSums(diff1 * diff2)
#
#     log_q[, g] = -0.5 * (1 / sigma[g]^2) * A
#   }
#
#   log_q = log_q - apply(log_q, 1, max)
#   q_new = exp(log_q)
#   q_new = q_new / rowSums(q_new)
#
#   q_new
# }
#
# update_z_rep2 = function(epsilon,
#                          mu,
#                          z,
#                          sigma,
#                          sigma_theta,
#                          alpha,
#                          psi,
#                          pz,
#                          model_data,
#                          smooth) {
#
#   H = ncol(epsilon)
#   G = model_data$dims$G
#   M = model_data$dims$M
#   n = model_data$dims$n
#   id = model_data$data$id
#   y = model_data$data$y
#   id = model_data$data$id
#   id_unique = model_data$data$id_unique
#   n_id = model_data$dims$n_id
#   n_time = model_data$dims$n_time
#
#   ll = lapply(1:G, function(g){
#
#     mu_g = matrix(mu[g, ], nrow = n, ncol = H, byrow = TRUE)
#     s_g = matrix(sigma[g, ], nrow = n, ncol = H, byrow = TRUE)
#
#     dnorm(epsilon, mean = mu_g, sd = s_g, log = TRUE) %>% rowSums()
#
#   }) %>% do.call(cbind, .)
#
#   prob = matrix(pz, nrow = n, ncol = G, byrow = T)
#   ll = norm_mat(ll + log(prob))
#
#   z = as.integer(extraDistr::rcatlp(n = n, log_prob = ll) + 1)
#
#
#   return(z)
#
# }
#
# # update_sigma_rep = function(theta, mu, z, model_data, sigma_theta, smooth) {
# #
# #   G = model_data$dims$G
# #   K = model_data$dims$K
# #   R = model_data$theta_spline$R
# #   L_inv = model_data$theta_spline$L_inv
# #
# #   Z = create_dummy(z, G)
# #   n_g = colSums(Z)
# #
# #   if(smooth == TRUE) {
# #     E = theta - R %*% Z %*% mu
# #   }else{
# #     V = diag(1/(sigma_theta^2))
# #     E = (theta - Z %*% mu) %*% V
# #   }
# #
# #   sigma = numeric(G)
# #
# #   for(g in 1:G) {
# #
# #     n_g = sum(z == g)
# #
# #     if(n_g > 0) {
# #
# #       ss2 = sum(E[z  == g, ]^2)
# #
# #       shape_post = 2 + 0.5 * n_g * K
# #       rate_post = 2 + 0.5 * ss2
# #       sigma[g] = sqrt(1 / rgamma(1, shape = shape_post, rate = rate_post))
# #
# #     }else{
# #       sigma[g] = sqrt(1 / rgamma(1, shape = 2, rate = 2))
# #     }
# #
# #   }
# #
# #   return(sigma)
# #
# # }
#
# update_sigma_rep2 = function(epsilon, mu, z, model_data, sigma_theta, smooth) {
#
#   G = model_data$dims$G
#   H = ncol(epsilon)
#
#   res = (epsilon - mu[z, ])^2
#
#   sigma = matrix(nrow = G, ncol = H)
#
#   for(g in 1:G) {
#
#     if(sum(z == g) == 0) {
#
#       sigma[g, ] = 1/sqrt(rgamma(n = H, shape = 2, rate = 2))
#
#     }else{
#
#       ss = colSums(rbind(res[z == g, ]))
#       ng = sum(z == g)
#       sigma[g, ] = 1/sqrt(rgamma(n = H, shape = rep(2, H) + ng/2, rate = 2 + ss/2))
#
#     }
#
#   }
#
#   return(sigma)
#
# }
#
#
# # update_z_rep = function(epsilon,
# #                         z,
# #                         mu,
# #                         sigma,
# #                         alpha,
# #                         psi,
# #                         pz,
# #                         model_data) {
# #
# #   H = ncol(epsilon)
# #   G = model_data$dims$G
# #   M = model_data$dims$M
# #   n = model_data$dims$n
# #   id = model_data$data$id
# #   y = model_data$data$y
# #   RRt = model_data$theta_spline$RRt
# #   Rn = model_data$theta_spline$Rn
# #   RnRnt = model_data$theta_spline$RnRnt
# #   R = model_data$theta_spline$R
# #   id = model_data$data$id
# #   id_unique = model_data$data$id_unique
# #   n_id = model_data$dims$n_id
# #   n_time = model_data$dims$n_time
# #   Q = model_data$theta_spline$RRt
# #   Q_inv = model_data$theta_spline$RRt_inv
# #   time = model_data$data$time_seq
# #   Rn = model_data$theta_spline$Rn
# #
# #   prob = matrix(pz, nrow = n_id, ncol = G, byrow = T)
# #
# #   for(t in 1:n_time) {
# #
# #     Z = create_dummy(z, G = G)
# #
# #     filter = time == t
# #
# #     epsilon_t = epsilon[filter, ]
# #     center = Rn %*% Z %*% mu
# #     center_t = center[filter, ]
# #     ll = (epsilon_t - center_t) %*% t(mu)
# #
# #     ll = ll + log(prob)
# #
# #     ll = ll - matrix(
# #       mclust::logsumexp(ll), nrow = n_id, ncol = G, byrow = F
# #     )
# #
# #     z[filter] = as.integer(extraDistr::rcatlp(n = n_id, log_prob = ll) + 1)
# #
# #   }
# #
# #
# #   return(z)
# #
# # }
#
#
# update_mu_sigma_ng2 = function(epsilon,
#                                z,
#                                model_data,
#                                m0 = 0, kappa0 = 1,
#                                a0 = 2, b0 = 2) {
#
#   G = model_data$dims$G
#   H = ncol(epsilon)
#
#   mu = matrix(0, G, H)
#   sigma = matrix(0, G, H)
#
#   for(g in 1:G) {
#
#     idx = which(z == g)
#     ng = length(idx)
#
#     if(ng == 0) {
#
#       # prior draw
#       tau = rgamma(H, shape = a0, rate = b0)
#       mu[g, ] = rnorm(H, m0, sqrt(1/(kappa0 * tau)))
#       sigma[g, ] = 1 / sqrt(tau)
#
#     } else {
#
#       x = epsilon[idx, , drop = FALSE]
#
#       xbar = colMeans(x)
#       ss = colSums((x - matrix(xbar, ng, H, byrow = TRUE))^2)
#
#       kappa_n = kappa0 + ng
#       m_n = (kappa0 * m0 + ng * xbar) / kappa_n
#       a_n = a0 + ng/2
#       b_n = b0 + 0.5 * ss + (kappa0 * ng * (xbar - m0)^2) / (2 * kappa_n)
#
#       tau = rgamma(H, shape = a_n, rate = b_n)
#
#       mu[g, ] = rnorm(H, m_n, sqrt(1/(kappa_n * tau)))
#       sigma[g, ] = 1 / sqrt(tau)
#     }
#   }
#
#   list(mu = mu, sigma = sigma)
# }
#
#
#


# update_theta_rep = function(mu,
#                             sigma,
#                             sigma_theta,
#                             alpha,
#                             psi,
#                             z,
#                             model_data,
#                             add_cluster,
#                             smooth) {
#
#   # dimensions
#   H = ncol(alpha)
#   J = model_data$dims$J
#   n_time = model_data$dims$n_time
#   n = model_data$dims$n
#   n_id = model_data$dims$n_id
#   G = model_data$dims$G
#
#   # data
#   y = model_data$data$y
#   id = model_data$data$id
#   id_unique = model_data$data$id_unique
#   I_time = diag(n_time)
#   I_H = diag(H)
#   Ri = model_data$theta_spline$Ri
#   RtR_i = model_data$theta_spline$RtR_i
#   Ry = model_data$theta_spline$Ry
#
#   # variability from data
#   AS = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
#   tASA = crossprod(AS, alpha)
#   RYAS = Ry %*% AS
#   V1 = kronecker(tASA, RtR_i)
#
#   U =  I_H
#
#   epsilon = matrix(nrow = n, ncol = H)
#
#   for(i in 1:n_id) {
#
#     filter = id == id_unique[i]
#     z_i = z[filter]
#     m0_i = mu[z_i, ]
#     Li = diag(1/(sigma[z_i]^2))
#
#     V0 = kronecker(U, Li)
#     V = solve(V0 + V1)
#
#     RYAS_i = RYAS[filter, , drop = FALSE]
#     m_i = V %*% as.vector(RYAS_i +  Li %*% m0_i)
#     epsilon_i = m_i + t(chol(V)) %*% rnorm(n_time * H)
#
#     epsilon[filter, ] = matrix(
#       epsilon_i,
#       nrow = n_time, ncol = H, byrow = FALSE
#     )
#
#   }
#
#   return(epsilon)
#
# }
