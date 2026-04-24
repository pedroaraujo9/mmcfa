#' Fit the cusp-only model
#'
#' Fits the non-clustering cusp model and returns chain-level samples,
#' combined posterior samples, model metadata, and runtime information.
#'
#' @param y A numeric matrix of observed outcomes with one row per observation
#'   and one column per variable.
#' @param id Subject identifier vector with one entry per row of y.
#' @param time Time index vector with one entry per row of y.
#' @param K Integer. Maximum number of latent factors.
#' @param iters Integer. Total MCMC iterations per chain.
#' @param burn_in Integer. Number of warmup iterations.
#' @param thin Integer. Thinning interval.
#' @param chains Integer. Number of MCMC chains.
#' @param n_cores Integer. Number of CPU cores used for parallel chains.
#' @param center_data Logical. If TRUE, center y before fitting.
#' @param scale_data Logical. If TRUE, scale y before fitting.
#' @param init_list Optional list of initial values.
#' @param seed Optional random seed.
#' @param cusp_nu Numeric. CUSP prior parameter nu.
#' @param cusp_a Numeric. CUSP prior parameter a.
#' @param cusp_b Numeric. CUSP prior parameter b.
#' @param cusp_min_var Numeric. Minimum variance used by the CUSP prior.
#' @param cusp_adapt_H Logical. If TRUE, adapt the active factor dimension.
#' @param cusp_alpha0 Numeric. Intercept in the adaptation probability schedule.
#' @param cusp_alpha1 Numeric. Slope in the adaptation probability schedule.
#' @param est_sigma Logical. If TRUE, sigma is updated during MCMC. If FALSE,
#'   sigma is kept at its initial value (1).
#' @param verbose Logical. If TRUE, print iteration progress.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{chains}: Chain-specific outputs.
#'   \item \code{post_sample}: Combined posterior samples.
#'   \item \code{model_data}: Internal data structures used for fitting.
#'   \item \code{args}: Input arguments used to call the function.
#'   \item \code{run_time}: Total runtime.
#' }
#'
#' @export
fit_cusp = function(y,
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

  init_time = Sys.time()

  args = list(
    y = y,
    id = id,
    time = time,
    K = K,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    chains = chains,
    n_cores = n_cores,
    center_data = center_data,
    scale_data = scale_data,
    init_list = init_list,
    seed = seed,
    cusp_nu = cusp_nu,
    cusp_a = cusp_a,
    cusp_b = cusp_b,
    cusp_min_var = cusp_min_var,
    cusp_adapt_H = cusp_adapt_H,
    cusp_alpha0 = cusp_alpha0,
    cusp_alpha1 = cusp_alpha1,
    est_sigma = est_sigma,
    verbose = verbose
  )

  runs = multiple_run_cusp(
    y = y,
    id = id,
    time = time,
    K = K,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    chains = chains,
    n_cores = n_cores,
    center_data = center_data,
    scale_data = scale_data,
    init_list = init_list,
    seed = seed,
    cusp_nu = cusp_nu,
    cusp_a = cusp_a,
    cusp_b = cusp_b,
    cusp_min_var = cusp_min_var,
    cusp_adapt_H = cusp_adapt_H,
    cusp_alpha0 = cusp_alpha0,
    cusp_alpha1 = cusp_alpha1,
    est_sigma = est_sigma,
    verbose = verbose
  )

  model_data = runs$model_data
  post_sample = combine_chains(runs$chains)

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
