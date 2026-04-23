update_alpha = function(H, theta, psi, model_data, prior_precision) {

  J = model_data$dims$J
  y = model_data$data$y

  #S0 = diag(prior_precision, nrow = H, ncol = H)
  #Psi = diag(1/psi, nrow = J, ncol = J)
  S0 = diag(prior_precision)
  Psi = diag(1/psi)

  OtO = crossprod(theta)
  V = solve(kronecker(OtO, Psi) + kronecker(S0, diag(J)))

  m = V %*% as.vector(Psi %*% t(y) %*% theta)
  alpha = m + t(chol(V)) %*% rnorm(J*H)
  alpha = matrix(alpha, nrow = J, ncol = H, byrow = FALSE)

  return(alpha)
}
