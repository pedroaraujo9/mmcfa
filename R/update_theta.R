update_theta = function(epsilon, model_data) {

  n_id = model_data$dims$n_id
  id = model_data$data$id
  id_unique = model_data$data$id_unique

  R = model_data$theta_spline$R

  theta = epsilon

  for(i in 1:n_id) {
    id_i = id_unique[i]
    theta[id == id_i, ] = R %*% epsilon[id == id_i, ]
  }

  return(theta)
}








