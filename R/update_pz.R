update_pz = function(z, z_dir, model_data) {

  G = model_data$dims$G
  id = model_data$data$id

  nz = z |> factor(levels = 1:G) |> table()
  pz = extraDistr::rdirichlet(1, alpha = rbind(nz + z_dir))[1, ]
  pz[pz < 1e-300] = 1e-300

  return(pz)

}
