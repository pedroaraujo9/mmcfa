find_number_clust_rep = function(M_max,
                                 G,
                                 K,
                                 id,
                                 time,
                                 n_basis,
                                 n_init = 50,
                                 init_iters = 10,
                                 lambda,
                                 dirichlet_param = 0.01,
                                 a_lambda = 1,
                                 b_lambda = 1,
                                 intercept_penalty = 1,
                                 seed,
                                 verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  model_data = create_model_data(
    y = y,
    id = id,
    time = time,
    K = K,
    G = G,
    M = M_max,
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

  n_id = model_data$dims$n_id
  id_unique = model_data$data$id_unique

  runs = lapply(1:n_init, function(i){

    if(verbose) cat("Init run:", i, "\r")

    single_run(
      model_data = model_data,
      iters = init_iters,
      burn_in = 0,
      thin = 1,
      alpha_prior = "normal",
      mixscat_prior = TRUE,
      tune = 0,
      init_list = NULL,
      z = NULL,
      w = NULL,
      seed = NULL,
      logP_proposal = NULL,
      adapt_H = FALSE,
      smooth = FALSE,
      alpha0 = 0,
      alpha1 = 0,
      clust_var = TRUE
    )
  })


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

  post_modes = lapply(sizes, function(size){

    #### condionate on the chosen model #####
    w_m = rbind(w_last[w_size == size, ])
    z_m = rbind(z_last[w_size == size, ])

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
      z = z_pear,
      w_pear = w_pear,
      w_sample = w_m,
      w_post_prob = w_post_prob
    )

    return(out)

  })

  names(post_modes) = sizes

  out = list(
    model_data = model_data,
    runs = runs,
    post_modes = post_modes,
    clust_size = clust_size,
    clust_size_p = clust_size_p
    # entropy = entropy,
    # prob_last = prob_last
  )

  return(out)

}
