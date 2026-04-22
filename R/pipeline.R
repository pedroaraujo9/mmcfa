pipeline = function(y,
                    id,
                    time,
                    G,
                    K,
                    M,
                    iters,
                    burn_in,
                    thin,
                    chains,
                    w_prior,
                    init_list = NULL,
                    z = NULL,
                    w = NULL,
                    alpha_prior = c("normal", "cusp"),
                    adapt_H = c(FALSE, TRUE),
                    seed = NULL) {

  model_data = create_model_data(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M,
    center = TRUE,
    scale = TRUE,
    n_basis = n_basis,
    theta_spline_penalty = 1,
    theta_intercept_penalty = 1,
    z_spline_penalty = 1,
    z_intercept_penalty = 1,
    w_dirichlet = 0.01,
    cusp_nu = 4,
    cusp_a = 2,
    cusp_b = 2,
    cusp_min_var = 0.01
  )

  runs = lapply(1:chains, function(i){

    single_run(
      model_data = model_data,
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      alpha_prior = alpha_prior,
      mixscat_prior = TRUE,
      tune = 0,
      init_list = init_list,
      z = NULL,
      w = NULL,
      seed = seed,
      w_prior = w_prior,
      logP_proposal = NULL,
      adapt_H = adapt_H,
      smooth = FALSE,
      alpha0 = 0,
      alpha1 = 0,
      clust_var = TRUE
    )

  })

  names(runs) = paste0("chain=", 1:chains)

  out = list(
    chains = runs,
    model_data = model_data
  )

  return(out)

}
