compute_spline_probs = function(model_data, post_sample) {

  M = model_data$dims$M
  G = model_data$dims$G
  n_time = model_data$dim$n_time
  iters = dim(post_sample$beta)[1]

  spline_probs = array(NA, dim = c(iters, M * n_time, G))

  for(i in 1:iters) {
    spline_probs[i,,] = compute_probs(
      w = 1:M,
      M = M,
      B = model_data$clustering$B,
      beta = post_sample$beta[i,,]
    )
  }

  spline_probs
}
