create_sample_list = function(iters = 1000,
                              model_data,
                              alpha_prior = c("normal", "cusp"),
                              init_list = NULL,
                              seed = NULL) {

  n = model_data$n
  J = model_data$J
  K = model_data$K
  G = model_data$G
  M = model_data$M

  if(!is.null(seed)) set.seed(seed)

  #### basic FA parameters ####
  alpha = gen_sample_array(
    iters = iters,
    dimension = c(J, K),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$alpha
  )

  theta = gen_sample_array(
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

  sample = list(
    alpha = alpha,
    theta = theta,
    psi = psi
  )

  #### CUSP parameters ####

  if(alpha_prior == "cusp") {

    sample$omega = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) rep(1/K, x),
      init = init_list$omega
    )

    sample$ind = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) sample(1:K, size = x, replace = T),
      init = init_list$ind
    )

    sample$H_active = gen_sample_array(
      iters = iters,
      dimension = c(1),
      sampler = function(x) K,
      init = init_list$H_active
    )

    sample$H_effective = gen_sample_array(
      iters = iters,
      dimension = c(1),
      sampler = function(x) K,
      init = init_list$H_active
    )

    sample$omega = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) rep(1/K, K),
      init = init_list$omega
    )

    sample$v = gen_sample_array(
      iters = iters,
      dimension = c(K),
      sampler = function(x) rbeta(x, 1, 1),
      init = init_list$v
    )

  }

  #### local clustering parameters ####
  sample$mu = gen_sample_array(
    iters = iters,
    dimension = c(G, K),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$mu
  )

  sample$z = gen_sample_array(
    iters = iters,
    dimension = c(n),
    sampler = function(x) sample(1:G, size = n, replace = T),
    init = init_list$z
  )

  sample$pz = gen_sample_array(
    iters = iters,
    dimension = c(G),
    sampler = function(x) rep(1/G, G),
    init = init_list$pz
  )

  sample$z_prob = gen_sample_array(
    iters = iters,
    dimension = c(n, G),
    sample = function(x) 1/G
  )

  #### global clustering parameters ####
  sample_list$beta = gen_sample_array(
    iters = iters,
    dimension = c(M*n_basis, G),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$beta
  )

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

  sample_list$w_prob = gen_sample_array(
    iters = iters,
    dimension = c(n_id, M),
    sampler = function(x) 1/M,
    init = init_list$w_prob
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

  return(sample)

}

