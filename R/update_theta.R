update_theta = function(epsilon, model_data) {

  n_id = model_data$dims$n_id
  id = model_data$data$id
  id_unique = model_data$data$id_unique
  R = model_data$theta_spline$Ri

  theta = epsilon

  for(i in 1:n_id) {
    filter = id_unique[i] == id
    theta[filter, ] = R %*% epsilon[filter, ]
  }

  return(theta)
}








