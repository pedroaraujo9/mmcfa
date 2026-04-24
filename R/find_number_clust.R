#' Find the number of active clusters
#'
#' Runs multiple short initial chains and summarizes the posterior distribution
#' of the number of active global clusters. It also returns representative
#' cluster allocations and posterior probabilities for each candidate size.
#'
#' @param y A numeric matrix of observed outcomes with one row per observation
#'   and one column per variable.
#' @param id Subject identifier vector with one entry per row of y.
#' @param time Time index vector with one entry per row of y.
#' @param K Integer. Number of latent factors.
#' @param G Integer. Number of local clusters used in initialization runs.
#' @param M Integer. Maximum number of global clusters.
#' @param n_init Integer. Number of independent initialization runs.
#' @param init_iters Integer. Number of iterations for each initialization run.
#' @param n_basis Integer. Number of spline basis functions.
#' @param spline_intercept_penalty Numeric. Penalty on spline intercept.
#' @param spline_penalty Numeric. Spline smoothness penalty.
#' @param center_data Logical. If TRUE, center y before fitting.
#' @param scale_data Logical. If TRUE, scale y before fitting.
#' @param z Optional fixed local cluster assignments. If NULL, z is
#'   randomly initialized in each run.
#' @param init_list Optional list of initial values passed to initialization
#'   runs.
#' @param seed Optional random seed.
#' @param w_dirichlet Numeric. Dirichlet concentration for w.
#' @param verbose Logical. If TRUE, print progress.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{model_data}: Internal data structures used in runs.
#'   \item \code{runs}: List of initialization run outputs.
#'   \item \code{entropy}: Long-format entropy summaries across runs/iterations.
#'   \item \code{post_modes}: Representative posterior summaries by cluster size.
#'   \item \code{clust_size}: Number of active clusters per run and iteration.
#'   \item \code{clust_size_p}: Posterior proportions for final active sizes.
#'   \item \code{run_time}: Total runtime.
#' }
#'
#' @export
find_number_clust = function(y,
                             id,
                             time,
                             K,
                             G,
                             M,
                             n_init,
                             init_iters,
                             n_basis,
                             spline_intercept_penalty = 1,
                             spline_penalty = 1,
                             center_data = TRUE,
                             scale_data = TRUE,
                             z = NULL,
                             init_list = NULL,
                             seed = NULL,
                             w_dirichlet = 0.01,
                             verbose = TRUE
) {

  init_time = Sys.time()
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
    z_dirichlet = 1,
    cusp_nu = 0,
    cusp_a = 0,
    cusp_b = 0,
    cusp_min_var = 0
  )

  n_id = model_data$dims$n_id
  id_unique = model_data$data$id_unique

  runs = lapply(1:n_init, function(i){

    if(verbose == TRUE) cat("iter ", i, "\r")

    single_run(
      model_data = model_data,
      iters = init_iters,
      burn_in = 0,
      thin = 1,
      alpha_prior = "normal",
      mixscat_prior = TRUE,
      cusp_adapt_H = FALSE,
      cusp_alpha0 = 0,
      cusp_alpha1 = 0,
      init_list = init_list,
      z = z,
      w = NULL,
      w_prior = NULL,
      seed = NULL,
      verbose = FALSE
    )

  })

  entropy = lapply(runs, function(run){

    lapply(1:init_iters, function(i){

      run$sample_list$w_post_prob[i, ,] %>% compute_entropy() %>% sum()

    }) %>% do.call(c, .)

  }) %>% do.call(rbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(runs = 1:nrow(.)) %>%
    tidyr::gather(iter, entropy, -runs) %>%
    dplyr::mutate(iter = as.integer(gsub("V", "", iter))) %>%
    dplyr::arrange(runs)

  clust_size = lapply(runs, function(run){
    lapply(1:init_iters, function(j){
      run$sample_list$w[j, ] %>% unique() %>% length()
    }) %>% do.call(c, .)
  }) %>%
    do.call(rbind, .)

  w_last = lapply(runs, function(run){
    run$sample_list$w[init_iters, ]
  }) %>%
    do.call(rbind, .)

  z_last = lapply(runs, function(run){
    run$sample_list$z[init_iters, ]
  }) %>%
    do.call(rbind, .)

  #### number of active cluster distribution ####
  p = clust_size[, init_iters] %>% table() %>% prop.table()

  clust_size_p = data.frame(
    n_clust = as.numeric(names(p)),
    prop = as.numeric(p)
  )

  sizes = clust_size_p$n_clust

  w_size = apply(w_last, MARGIN = 1, FUN = function(x) {
    x %>% unique() %>% length()
  })

  z_size = apply(z_last, MARGIN = 1, FUN = function(x) {
    x %>% unique() %>% length()
  })

  post_modes = lapply(sizes, function(size){

    #### condionate on the chosen model #####
    w_m = rbind(w_last[w_size == size, ])
    z_m = rbind(z_last[w_size == size, ])

    zm_size = apply(z_m, MARGIN = 1, FUN = function(x) {
      x %>% unique() %>% length()
    })

    n_best = nrow(w_m)

    for(i in 1:n_best){
      w_m[i, ]  = w_m[i, ] %>% factor(labels = 1:size) %>% as.integer()
    }

    psm = mcclust::comp.psm(cls = w_m)
    colnames(psm) = rownames(psm) = id_unique

    if (n_best == 1) {
      w_pear = w_m[1, ]
    } else {
      cl = mcclust::maxpear(psm, cls.draw = w_m, max.k = size, method = "draws")
      w_pear = cl$cl
      if (length(w_pear) == 0) {
        w_pear = w_m[1, ]
      }
    }

    if (size >= 2) {
      ls = label.switching::label.switching(
        method = "ECR",
        z = w_m,
        zpivot = w_pear,
        K = size
      )

      for(i in 1:n_best) {
        perm = ls$permutations$`ECR`[i, ]
        w_m[i, ] = order(perm)[w_m[i, ]]
      }
    } else {
      # If size is 1, no label switching is possible or necessary.
      # w_m remains as all 1s.
      message("Note: size=1 detected for a model group. Skipping label switching.")
    }

    w_post_prob = lapply(1:n_id, function(i){
      w_m[, i] %>%
        factor(levels = 1:size) %>%
        table() %>%
        prop.table() %>%
        as.numeric()
    }) %>%
      do.call(rbind, .)

    rownames(w_post_prob) = id_unique
    w = w_post_prob %>% apply(MARGIN = 1, FUN = which.max)

    ### init z ###
    z_psm = mcclust::comp.psm(cls = z_m)
    z_cl = mcclust::maxpear(z_psm, cls.draw = z_m, max.k = G, method = "draws")
    z_pear = z_cl$cl

    out = list(
      w = w,
      z_size = zm_size,
      z = z_pear,
      w_pear = w_pear,
      w_sample = w_m,
      w_post_prob = w_post_prob
    )

    return(out)

  })

  names(post_modes) = sizes

  end_time = Sys.time()
  run_time = end_time - init_time

  out = list(
    model_data = model_data,
    entropy = entropy,
    post_modes = post_modes,
    clust_size = clust_size,
    clust_size_p = clust_size_p,
    run_time = run_time
  )

  return(out)

}
