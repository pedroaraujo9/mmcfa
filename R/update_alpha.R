update_alpha = function(theta,
                        psi,
                        prior_precision,
                        model_data) {

  J = model_data$dims$J
  y = model_data$data$y
  n = model_data$dims$n
  H = ncol(theta)

  alpha = matrix(0, nrow = J, ncol = H)
  tau = 1/psi
  tau_matrix = matrix(1/psi, nrow = n, ncol = J, byrow = T)

  prior_mean = rep(0, H)

  for(j in 1:J) {

    if(H == 1) {

      prior_precision_j = matrix(prior_precision[j, ], nrow = 1, ncol = 1)

    }else{

      prior_precision_j = diag(prior_precision[j, ])

    }

    alpha[j, ] = propose_coef_rcpp(
      y = cbind(y[, j]),
      X = cbind(theta),
      X_prec = matrix(tau[j], nrow = n, ncol = H),
      y_prec = cbind(tau_matrix[, j]),
      prior_mean = prior_mean,
      prior_precision = prior_precision_j
    )

  }

  return(alpha)
}
