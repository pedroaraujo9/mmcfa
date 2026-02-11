update_psi = function(alpha, theta, a, b, model_data) {

  n = model_data$dims$n
  J = model_data$dims$J
  y = model_data$data$y

  resid = y - (theta %*% t(alpha))

  sum_resid = colSums(resid^2)
  psi = 1/rgamma(J, shape = a + 0.5 * n, rate = b + 0.5*sum_resid)

  return(psi)

}
