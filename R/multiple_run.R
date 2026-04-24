multiple_run = function(y,
                        id,
                        time,
                        K,
                        G,
                        M,
                        n_basis,
                        iters,
                        burn_in,
                        thin,
                        chains,
                        n_cores,
                        mixscat_prior = c(TRUE, FALSE),
                        alpha_prior = c("normal", "cusp"),
                        center_data = TRUE,
                        scale_data = TRUE,
                        init_list = NULL,
                        z = NULL,
                        w = NULL,
                        w_prior = NULL,
                        seed = NULL,
                        w_dirichlet = 0.01,
                        z_dirichlet = 0.01,
                        spline_intercept_penalty = 1,
                        spline_penalty = 1,
                        cusp_nu = 4,
                        cusp_a = 2,
                        cusp_b = 2,
                        cusp_min_var = 0.01,
                        cusp_adapt_H = FALSE,
                        cusp_alpha0 = 4,
                        cusp_alpha1 = 5*1e-3,
                        add_cluster = TRUE,
                        verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  model_data = create_model_data(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M,
    center = center_data,
    scale = scale_data,
    n_basis = n_basis,
    spline_intercept_penalty = spline_intercept_penalty,
    spline_penalty = spline_penalty,
    w_dirichlet = w_dirichlet,
    z_dirichlet = z_dirichlet,
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


  single_run_copy = single_run


  if(n_cores > 1) {

    runs = future.apply::future_lapply(1:chains, FUN = function(i) {

      single_run_copy(
        model_data = model_data,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        alpha_prior = alpha_prior,
        mixscat_prior = mixscat_prior,
        cusp_adapt_H = cusp_adapt_H,
        cusp_alpha0 = cusp_alpha0,
        cusp_alpha1 = cusp_alpha1,
        add_cluster = add_cluster,
        init_list = init_list,
        z = z,
        w = w,
        w_prior = w_prior,
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
        alpha_prior = alpha_prior,
        mixscat_prior = mixscat_prior,
        cusp_adapt_H = cusp_adapt_H,
        cusp_alpha0 = cusp_alpha0,
        cusp_alpha1 = cusp_alpha1,
        add_cluster = add_cluster,
        init_list = init_list,
        z = z,
        w = w,
        w_prior = w_prior,
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
