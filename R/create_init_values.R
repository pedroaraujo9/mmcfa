create_init_values = function(model_data) {

  #### FA inits ####
  y = model_data$y
  dec = eigen(cov(y))
  V = cbind(dec$vectors[, 1:K])

  if(K == 1) {

    L = as.matrix(sqrt(dec$values[1:K]))

  }else{

    L = diag(sqrt(dec$values[1:K]))

  }


  alpha = V %*% L
  theta = y %*% V %*% solve(L)

  #### local clustering mixture ####
  km = kmeans(theta_init, centers = G, nstart = 1000, iter.max = 10000)
  z = km$cluster
  mu = km$centers

  #### global clustering mixture ####
  Zm = matrix(z, nrow = model_data$n_id, ncol = model_data$n_time, byrow = T) |>
    TraMineR::seqdef()

  d = TraMineR::seqdist(Zm, method = "HAM") %>% as.dist()
  hc = hclust(d, method = "ward.D")
  w = hc |> cutree(k = M)
  names(w) = model_data$id_unique

  #### return ####
  out = list(
    alpha = alpha_init,
    theta = theta_init,
    z = z,
    mu = mu,
  )

  return(out)

}
