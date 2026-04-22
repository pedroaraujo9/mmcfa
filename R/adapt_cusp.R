adapt_cusp = function(H_active,
                      H,
                      H_max,
                      ind,
                      alpha,
                      mu,
                      z,
                      theta,
                      sigma,
                      omega,
                      alpha_precision,
                      v,
                      nu,
                      min_var,
                      model_data) {

  G = model_data$dims$G
  idx = as.integer(model_data$data$id)
  R = model_data$theta_spline$R

  if(H_active < H - 1) {

    # remove inactive
    active_col = c(ind[1:H] > 1:H, rep(FALSE, H_max - H))

    alpha = cbind(alpha[, active_col]) %>% complete_dim(dimension = H_max)
    mu = cbind(mu[, active_col]) %>% complete_dim(dimension = H_max) %>% cbind()
    theta = cbind(theta[, active_col]) %>% complete_dim(dimension = H_max) %>% cbind()
    omega = omega[active_col] %>% complete_dim(dimension = H_max) %>% cbind()
    alpha_precision = cbind(alpha_precision[, active_col]) %>% complete_dim(dimension = H_max)

    H = H_active + 1

    omega_H = 1 - sum(omega[1:(H-1)], na.rm = TRUE)
    prec_H = 1/min_var
    alpha_H = rnorm(nrow(alpha), 0, sqrt(min_var))
    mu_H = rnorm(G, 0, 1)
    theta_H = rnorm(nrow(theta), as.numeric(R %*% cbind(mu_H[z])), sigma[idx])
    alpha[, H] = alpha_H
    mu[, H] = mu_H
    theta[, H] = theta_H
    omega[H] = omega_H
    alpha_precision[, H] = prec_H

    v[H] = 1
    v[(H+1):H_max] = NA

  }else if(H < H_max){

    H = H + 1
    alpha_H = rnorm(nrow(alpha), 0, sqrt(min_var))

    v[H-1] = rbeta(1, 1, nu)
    v[H] = 1
    omega = v[1:H] %>% stick_breaking()

    mu_H = rnorm(G, 0, 1)
    theta_H = rnorm(nrow(theta), as.numeric(R %*% cbind(mu_H[z])), sigma[idx])

    # update
    alpha[, H] = alpha_H
    mu[, H] = mu_H
    theta[, H] = theta_H
    alpha_precision[, H] = 1/min_var
    omega = omega %>% complete_dim(dimension = H_max)

    if(H < H_max) {

      alpha[, (H+1):H_max] = NA
      theta[, (H+1):H_max] = NA
      mu[, (H+1):H_max] = NA
      v[(H+1):H_max] = NA

    }
  }

  out = list(
    H = H,
    alpha = alpha,
    mu = mu,
    theta = theta,
    omega = omega,
    alpha_precision = alpha_precision,
    v = v
  )

  return(out)

}
