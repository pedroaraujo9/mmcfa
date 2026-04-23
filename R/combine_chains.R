combine_chains = function(runs) {

  sample_list = runs[[1]]$sample_list

  if(length(runs) > 1) {

    for(i in 2:length(runs)) {

      params = c(
        "theta", "alpha", "mu", "sigma", "z", "pz", "z_post_prob",
        "H", "H_active",
        "beta", "pw", "w", "w_post_prob"
      )

      for(param in params) {

        sample_list[[param]] = abind::abind(
          sample_list[[param]], runs[[i]]$sample_list[[param]], along = 1
        )

      }

    }

  }

  return(sample_list)

}


build_block_perm_mat = function(perm, n_basis) {

  M = length(perm)
  P_sigma = diag(M)[perm, ]
  T_mat = kronecker(P_sigma, diag(n_basis))

  return(T_mat)
}

apply_relabel = function(sample_list, pivot_z, pivot_w, n_basis) {

  iters = dim(sample_list$w_post_prob)[1]
  G = dim(sample_list$z_post_prob)[3]
  M = dim(sample_list$w_post_prob)[3]

  ls_z = label.switching::label.switching(
    method = "ECR",
    z = sample_list$z,
    zpivot = pivot_z,
    K = G
  )

  ls_w = label.switching::label.switching(
    method = "ECR",
    z = sample_list$w,
    zpivot = pivot_w,
    K = M
  )

  for(i in 1:iters) {
    perm_z = labels = ls_z$permutations$`ECR`[i, ]
    perm_w = labels = ls_w$permutations$`ECR`[i, ]

    sample_list$z[i, ] = order(perm_z)[sample_list$z[i, ]]
    sample_list$w[i, ] = order(perm_w)[sample_list$w[i, ]]

    if(!is.null(sample_list$pz)) sample_list$pz[i, ] = sample_list$pz[i, perm_z]
    if(!is.null(sample_list$pw)) sample_list$pw[i, ] = sample_list$pw[i, perm_w]

    sample_list$w_post_prob[i,,] = sample_list$w_post_prob[i,, perm_w]
    sample_list$beta[i,,] = build_block_perm_mat(perm_w, n_basis) %*% sample_list$beta[i,,]

    sample_list$z_post_prob[i,,] = sample_list$z_post_prob[i,, perm_z]
    sample_list$mu[i,,] = sample_list$mu[i, perm_z, ]
    sample_list$sigma[i,,] = sample_list$sigma[i, , perm_z]

  }

  return(sample_list)

}


relabel = function(runs, n_basis, pivot_z = NULL, pivot_w = NULL) {

  last_iter = nrow(runs[[1]]$sample_list$w)

  if(is.null(pivot_z)) {
    pivot_z = runs[[1]]$sample_list$z[last_iter, ]
  }

  if(is.null(pivot_w)) {
    pivot_w = runs[[1]]$sample_list$w[last_iter, ]
  }

  for(chain in 1:length(runs)) {
    runs[[chain]]$sample_list = apply_relabel(
      sample_list = runs[[chain]]$sample_list,
      pivot_z = pivot_z,
      pivot_w = pivot_w,
      n_basis = n_basis
    )
  }

  return(runs)

}


