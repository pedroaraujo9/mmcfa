#' Fit the mixff model
#'
#' Fits the full model and returns chain-level samples, combined posterior
#' samples, model metadata, and runtime information.
#'
#' @param y A numeric matrix of observed outcomes with one row per observation
#'   and one column per variable.
#' @param id Subject identifier vector with one entry per row of y.
#' @param time Time index vector with one entry per row of y.
#' @param K Integer. Maximum number of latent factors.
#' @param G Integer. Number of local clusters.
#' @param M Integer. Number of global clusters.
#' @param n_basis Integer. Number of spline basis functions.
#' @param iters Integer. Total MCMC iterations per chain.
#' @param burn_in Integer. Number of warmup iterations.
#' @param thin Integer. Thinning interval.
#' @param chains Integer. Number of MCMC chains.
#' @param n_cores Integer. Number of CPU cores used for parallel chains.
#' @param mixscat_prior Logical. If TRUE, use the mixscat prior.
#' @param center_data Logical. If TRUE, center y before fitting.
#' @param scale_data Logical. If TRUE, scale y before fitting.
#' @param init_list Optional list of initial values.
#' @param z Optional fixed local cluster assignments.
#' @param w Optional fixed global cluster assignments.
#' @param w_prior Optional prior values for w.
#' @param seed Optional random seed.
#' @param w_dirichlet Numeric. Dirichlet concentration for w.
#' @param z_dirichlet Numeric. Dirichlet concentration for z.
#' @param spline_intercept_penalty Numeric. Penalty on spline intercept.
#' @param spline_penalty Numeric. Spline smoothness penalty.
#' @param add_cluster Logical. If FALSE, run without the clustering component.
#' @param verbose Logical. If TRUE, print iteration progress.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{chains}: Chain-specific outputs.
#'   \item \code{post_sample}: Combined posterior samples after relabeling.
#'   \item \code{model_data}: Internal data structures used for fitting.
#'   \item \code{args}: Input arguments used to call the function.
#'   \item \code{run_time}: Total runtime.
#' }
#'
#' @export
fit_model = function(y,
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
                     center_data = TRUE,
                     scale_data = TRUE,
                     init_list = NULL,
                     z = NULL,
                     w = NULL,
                     w_prior = NULL,
                     relabel = TRUE, 
                     seed = NULL,
                     w_dirichlet = 0.01,
                     z_dirichlet = 0.01,
                     spline_intercept_penalty = 1,
                     spline_penalty = 1,
                     add_cluster = TRUE,
                     verbose = TRUE) {

  init_time = Sys.time()

  alpha_prior = "normal"
  cusp_nu = 4
  cusp_a = 2
  cusp_b = 2
  cusp_min_var = 0.01
  cusp_adapt_H = FALSE
  cusp_alpha0 = 4
  cusp_alpha1 = 5 * 1e-3

  args = list(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M,
    n_basis = n_basis,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    chains = chains,
    n_cores = n_cores,
    mixscat_prior = mixscat_prior,
    center_data = center_data,
    scale_data = scale_data,
    init_list = init_list,
    z = z,
    w = w,
    w_prior = w_prior,
    seed = seed,
    w_dirichlet = w_dirichlet,
    z_dirichlet = z_dirichlet,
    spline_intercept_penalty = spline_intercept_penalty,
    spline_penalty = spline_penalty,
    alpha_prior = alpha_prior,
    add_cluster = add_cluster,
    verbose = verbose
  )

  runs = multiple_run(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M,
    n_basis = n_basis,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    chains = chains,
    n_cores = n_cores,
    mixscat_prior = mixscat_prior,
    alpha_prior = alpha_prior,
    center_data = center_data,
    scale_data = scale_data,
    init_list = init_list,
    z = z,
    w = w,
    w_prior = w_prior,
    seed = seed,
    w_dirichlet = w_dirichlet,
    z_dirichlet = z_dirichlet,
    spline_intercept_penalty = spline_intercept_penalty,
    spline_penalty = spline_penalty,
    cusp_nu = cusp_nu,
    cusp_a = cusp_a,
    cusp_b = cusp_b,
    cusp_min_var = cusp_min_var ,
    cusp_adapt_H = cusp_adapt_H,
    cusp_alpha0 = cusp_alpha0,
    cusp_alpha1 = cusp_alpha1,
    add_cluster = add_cluster,
    verbose = verbose
  )

  model_data = runs$model_data

  if(relabel == TRUE) {

    runs$chains = relabel(runs$chains, n_basis = model_data$clustering$n_basis)
    
  }

  post_sample = combine_chains(runs$chains)
  post_sample$spline_probs = compute_spline_probs(model_data, post_sample)

  if (isFALSE(mixscat_prior)) {
    drop_params = c("w", "w_post_prob", "beta", "pw")

    post_sample[drop_params] = NULL

    runs$chains = lapply(runs$chains, function(chain) {
      if (!is.null(chain$sample_list)) {
        chain$sample_list[drop_params] = NULL
      }
      if (!is.null(chain$init)) {
        chain$init[drop_params] = NULL
      }
      chain
    })
  }

  end_time = Sys.time()
  run_time = end_time - init_time

  out = list(
    chains = runs$chains,
    post_sample = post_sample,
    model_data = model_data,
    args = args,
    run_time = run_time
  )

  return(out)

}
