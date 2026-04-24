multiple_run_cusp = function(y,
                              id,
                              time,
                              K,
                              iters,
                              burn_in,
                              thin,
                              chains,
                              n_cores,
                              center_data = TRUE,
                              scale_data = TRUE,
                              init_list = NULL,
                              seed = NULL,
                              cusp_nu = 4,
                              cusp_a = 2,
                              cusp_b = 2,
                              cusp_min_var = 0.01,
                              cusp_adapt_H = FALSE,
                              cusp_alpha0 = 4,
                              cusp_alpha1 = 5*1e-3,
                              est_sigma = TRUE,
                              verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  G_dummy = 2
  M_dummy = 2
  n_basis_dummy = 2
  w_dirichlet_dummy = 1
  z_dirichlet_dummy = 1
  spline_intercept_penalty_dummy = 1
  spline_penalty_dummy = 1

  model_data = create_model_data(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G_dummy,
    M = M_dummy,
    center = center_data,
    scale = scale_data,
    n_basis = n_basis_dummy,
    spline_intercept_penalty = spline_intercept_penalty_dummy,
    spline_penalty = spline_penalty_dummy,
    w_dirichlet = w_dirichlet_dummy,
    z_dirichlet = z_dirichlet_dummy,
    cusp_nu = cusp_nu,
    cusp_a = cusp_a,
    cusp_b = cusp_b,
    cusp_min_var = cusp_min_var
  )

  n_cores = max(1L, as.integer(n_cores))
  old_plan = future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (n_cores > 1L) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

  single_run_copy = single_run_cusp

  if(n_cores > 1) {

    runs = future.apply::future_lapply(1:chains, FUN = function(i) {

      single_run_copy(
        model_data = model_data,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        cusp_adapt_H = cusp_adapt_H,
        cusp_alpha0 = cusp_alpha0,
        cusp_alpha1 = cusp_alpha1,
        est_sigma = est_sigma,
        init_list = init_list,
        seed = NULL,
        verbose = verbose
      )

    }, future.seed = TRUE)

  }else{

    runs = lapply(1:chains, function(i){

      single_run_copy(
        model_data = model_data,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        cusp_adapt_H = cusp_adapt_H,
        cusp_alpha0 = cusp_alpha0,
        cusp_alpha1 = cusp_alpha1,
        est_sigma = est_sigma,
        init_list = init_list,
        seed = NULL,
        verbose = verbose
      )

    })

  }

  names(runs) = paste0("chain=", 1:chains)

  out = list(
    chains = runs,
    model_data = model_data
  )

  return(out)

}
