create_init_values = function(model_data) {

  K = model_data$dims$K
  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_time = model_data$dims$n_time
  id_unique = model_data$data$id_unique
  time = model_data$data$time_seq

  #### FA inits ####
  y = model_data$data$y
  dec = eigen(cov(y))
  V = cbind(dec$vectors[, 1:K])

  if(K == 1) {

    L = as.matrix(sqrt(dec$values[1:K]))

  }else{

    L = diag(sqrt(dec$values[1:K]))

  }


  alpha = V %*% L
  theta = y %*% V %*% solve(L)
  theta_sp = update_theta(theta, model_data = model_data)

  #### local clustering mixture ####
  km = kmeans(theta_sp, centers = G, nstart = 1000, iter.max = 10000)
  z = km$cluster
  mu = km$centers

  #### global clustering mixture ####
  Zm = matrix(z, nrow = n_id, ncol = n_time, byrow = T) |>
    TraMineR::seqdef() |>
    suppressMessages()

  d = TraMineR::seqdist(Zm, method = "HAM") |> suppressMessages() |> as.dist()
  hc = hclust(d, method = "ward.D")
  w = hc |> cutree(k = M)
  names(w) = id_unique

  #### return ####
  out = list(
    alpha = alpha,
    theta = theta_sp,
    z = z,
    mu = mu,
    w = w
  )

  return(out)

}
