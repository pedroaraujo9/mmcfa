create_sample_list = function(iters = 1000,
                              burn_in = iters/2,
                              thin = 1,
                              model_data,
                              alpha_prior = c("normal", "cusp"),
                              init_list = NULL,
                              seed = NULL) {

  n = model_data$dims$n
  J = model_data$dims$J
  K = model_data$dims$K
  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$theta_spline$n_basis

  id_unique = model_data$data$id_unique
  alpha_prior = match.arg(alpha_prior, choices = c("normal", "cusp"))

  if(!is.null(seed)) set.seed(seed)

  iters_vec = seq(from = burn_in + 1, to = iters, by = thin)
  iters = length(iters_vec)

  #### basic FA parameters ####
  alpha = gen_sample_array(
    iters = iters,
    dimension = c(J, K),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$alpha
  )

  alpha_precision = gen_sample_array(
    iters = iters,
    dimension = c(J, K),
    sampler = function(x) rgamma(x, shape = 1, rate = 1),
    init = init_list$alpha_precision
  )

  theta = epsilon = gen_sample_array(
    iters = iters,
    dimension = c(n, K),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$theta
  )

  psi = gen_sample_array(
    iters = iters,
    dimension = c(J),
    sampler = function(x){1},
    init = init_list$psi
  )

  sample_list = list(
    alpha = alpha,
    alpha_precision = alpha_precision,
    theta = theta,
    epsilon = epsilon,
    psi = psi
  )

  #### CUSP parameters ####

  if(alpha_prior == "cusp") {

    sample_list$omega = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) rep(1/K, x),
      init = init_list$omega
    )

    sample_list$ind = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) sample(1:K, size = x, replace = T),
      init = init_list$ind
    )

    sample_list$v = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) rbeta(x, 1, 1),
      init = init_list$v
    )

    sample_list$H_active = gen_sample_array(
      iters = iters,
      dimension = c(1),
      sampler = function(x) K,
      init = init_list$H_active
    )

    sample_list$H_effective = gen_sample_array(
      iters = iters,
      dimension = c(1),
      sampler = function(x) K,
      init = init_list$H_active
    )

  }

  #### local clustering parameters ####
  sample_list$mu = gen_sample_array(
    iters = iters,
    dimension = c(G, K),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$mu
  )

  sample_list$sigma = gen_sample_array(
    iters = iters,
    dimension = c(G),
    sampler = function(x) 1,
    init = init_list$sigma
  )

  sample_list$z = gen_sample_array(
    iters = iters,
    dimension = c(n),
    sampler = function(x) sample(1:G, size = n, replace = T),
    init = init_list$z
  )

  sample_list$pz = gen_sample_array(
    iters = iters,
    dimension = c(G),
    sampler = function(x) rep(1/G, G),
    init = init_list$pz
  )

  sample_list$z_post_prob = gen_sample_array(
    iters = iters,
    dimension = c(n, G),
    sampler = function(x) 1/G
  )

  #### global clustering parameters ####
  sample_list$beta = gen_sample_array(
    iters = iters,
    dimension = c(M*(n_basis), G),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$beta
  )

  sample_list$beta[1,,G] = 0

  sample_list$w = gen_sample_array(
    iters = iters,
    dimension = c(n_id),
    sampler = function(x) {
      w = sample(1:M, size = n_id, replace = TRUE)
      names(w) = id_unique
      return(w)
    },
    init = init_list$w
  )

  colnames(sample_list$w) = id_unique

  sample_list$w_post_prob = gen_sample_array(
    iters = iters,
    dimension = c(n_id, M),
    sampler = function(x) 1/M,
    init = init_list$w_post_prob
  )

  sample_list$pw = gen_sample_array(
    iters = iters,
    dimension = c(M),
    sampler = function(x) rep(1/M, M),
    init = init_list$pw
  )

  #### log posterior ####
  sample_list$logpost = gen_sample_array(
    iters = iters,
    dimension = c(9),
    sampler = function(x) 0,
    init = NULL
  )

  colnames(sample_list$logpost) = c(
    "y", "alpha", "theta", "mu",
    "psi", "pz", "beta", "pw", "logpost"
  )

  sample_list$iters_vec = iters_vec

  return(sample_list)

}

