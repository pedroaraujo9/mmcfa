update_theta = function(H,
                        mu,
                        sigma,
                        alpha,
                        psi,
                        z,
                        model_data,
                        add_cluster) {

  # dimensions
  J = model_data$dims$J
  n_time = model_data$dims$n_time
  n = model_data$dims$n
  n_id = model_data$dims$n_id
  G = model_data$dims$G

  # data
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  y = model_data$data$y

  I_H = diag(H)

  # variability from data
  AS = alpha * matrix(1/psi, nrow = J, ncol = H, byrow = FALSE)
  tASA = crossprod(AS, alpha)
  YAS = y %*% AS

  theta = matrix(nrow = n, ncol = H)
  noise = gen_normal_mat(n, H)

  if(add_cluster == TRUE) {

  tau = 1/(sigma^2)

  for(i in 1:n_id) {

      filter = id == id_unique[i]
      rows_i = which(filter)
      zi = z[filter]

      for(g in unique(zi)) {

        zig = (zi == g)
        idx_tg = rows_i[zig]

        tau_i = tau[i, g]
        V = solve(tASA + I_H * tau_i)

        mi = (YAS[idx_tg, , drop=FALSE] + tau_i * matrix(mu[g, ], nrow=sum(zig), ncol=H, byrow=TRUE)) %*% V

        theta[idx_tg, ] = mi + noise[idx_tg, , drop=FALSE] %*% chol(V)

      }
    }
    
  }else{

    if(is.matrix(sigma)) {
      sigma_id = sigma[, 1]
    }else{
      sigma_id = as.numeric(sigma)
    }

    tau_id = 1/(sigma_id^2)

    for(i in 1:n_id) {

      filter = id == id_unique[i]
      rows_i = which(filter)
      tau_i = tau_id[i]

      V = solve(tASA + I_H * tau_i)
      mi = YAS[rows_i, , drop = FALSE] %*% V

      theta[rows_i, ] = mi + noise[rows_i, , drop = FALSE] %*% chol(V)

    }

  }
  
  return(theta)

}
