update_pz = function(z, z_dir, model_data, time_pz, w, w_fixed) {

  G = model_data$dims$G
  id = model_data$data$id

  if(time_pz == FALSE) {
    nz = z |> factor(levels = 1:G) |> table()
    pz = extraDistr::rdirichlet(1, alpha = rbind(nz + z_dir))[1, ]

  }else{

    time = model_data$data$time_seq
    n_time = model_data$dims$n_time
    z_fac =  z |> factor(levels = 1:G)
    time_fac = time |> factor(1:n_time)

    pz = lapply(1:M, function(m){

      filter = w[id] == m
      count = as.matrix(table(time_fac, z_fac))
      count = count + z_dir

      extraDistr::rdirichlet(n_time, alpha = count)

    }) %>% do.call(rbind, .)

  }

  pz[pz < 1e-300] = 1e-300

  return(pz)
}
