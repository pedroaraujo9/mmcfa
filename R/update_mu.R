update_mu = function(H,
                     z,
                     sigma,
                     theta,
                     model_data) {

  G = model_data$dims$G
  idx = as.integer(model_data$data$id)

  sigma_long = sigma[cbind(idx, z)]
  tau = 1/(sigma_long^2)

  res_prec = tapply(tau, z, sum)
  prec_diag = rep(1, G)
  prec_diag[as.numeric(names(res_prec))] = 1 + res_prec

  m_num = matrix(0, nrow = G, ncol = H)
  present_states = as.numeric(names(res_prec))

  for(h in 1:H) {
    res_mean = tapply(theta[, h] * tau, z, sum)
    m_num[present_states, h] = res_mean
  }

  V_diag = 1 / prec_diag
  m = m_num * V_diag
  mu = m + (sqrt(V_diag) * gen_normal_mat(G, H))


  return(mu)

}
