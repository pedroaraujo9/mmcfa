update_sigma = function(H,
                        theta,
                        mu,
                        z,
                        model_data,
                        add_cluster = TRUE) {

  G = model_data$dims$G
  n_id = model_data$dims$n_id

  id = model_data$data$id
  id_unique = model_data$data$id_unique

  if(add_cluster == TRUE) {

    se = (theta - mu[z, ])^2
    zf = factor(z, levels = 1:G)
    n_id_g = as.matrix(table(id, zf))

    sigma = lapply(1:G, function(g){

      g_filter = z == g

      if(sum(g_filter) > 0) {

        see_g = rowsum(x = rbind(se[g_filter, ]), group = id[g_filter]) %>% rowSums(na.rm = TRUE)
        see_g = see_g[as.character(id_unique)]
        see_g[is.na(see_g)] = 0
        tau = rgamma(n_id, shape = 2 + 0.5*n_id_g[, g]*H, rate = 2 + 0.5*see_g)

      }else{

        tau = rgamma(n_id, shape = 2, rate = 2)

      }

      1/sqrt(tau)


    }) %>% do.call(cbind, .)

  }else{

    se_id = rowsum(x = rowSums(theta^2), group = id)
    se_id = se_id[as.character(id_unique), , drop = FALSE]
    se_id[is.na(se_id)] = 0

    n_id_obs = table(id)
    n_id_obs = as.numeric(n_id_obs[as.character(id_unique)])

    tau = rgamma(
      n_id,
      shape = 2 + 0.5 * n_id_obs * H,
      rate = 2 + 0.5 * as.numeric(se_id)
    )

    sigma_id = 1/sqrt(tau)
    sigma = matrix(sigma_id, nrow = n_id, ncol = G)

  }

  return(sigma)

}
