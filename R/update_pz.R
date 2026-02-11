update_pz = function(z, z_dir, G, model_data) {
  nz = z %>% factor(levels = 1:G) %>% table()
  pz = extraDistr::rdirichlet(1, alpha = rbind(nz + z_dir))[1, ]
  return(pz)
}
