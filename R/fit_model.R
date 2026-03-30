fit_model = function(y,
                     id,
                     time,
                     K,
                     G,
                     M,
                     chains,
                     iters,
                     burn_in,
                     thin,
                     n_basis = 10,
                     dirichlet_param = 0.01,
                     n_cores = 1,
                     add_prob_spline = TRUE,
                     seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  model_data = create_model_data(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M,
    center = TRUE,
    scale = FALSE,
    n_basis = n_basis,
    theta_spline_penalty = 1,
    theta_intercept_penalty = 1,
    z_spline_penalty = 1,
    z_intercept_penalty = 1,
    w_dirichlet = dirichlet_param
  )

  n_cores = max(1L, as.integer(n_cores))
  old_plan = future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (n_cores > 1L) {
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

  single_run_fn = single_run


  if(n_cores > 1) {

    runs = future.apply::future_lapply(1:chains, FUN = function(i) {

      single_run_fn(
        model_data = model_data,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        alpha_prior = "normal",
        seed = NULL,
        add_prob_spline = add_prob_spline
      )

    }, future.seed = TRUE)

  }else{

    runs = lapply(1:chains, FUN = function(i) {

      single_run_fn(
        model_data = model_data,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        alpha_prior = "normal",
        seed = NULL,
        add_prob_spline = add_prob_spline
      )

    })

  }

  out = list(
    model_data = model_data,
    runs = runs
  )

  return(out)

}
