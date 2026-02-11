test_that("function works across different dimensions", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  for(k in 1:2) {
    for(g in 1:2) {
      for(m in 1:2) {

        model_data = create_model_data(
          y = sim_data$y,
          id = sim_data$id,
          time = sim_data$time,
          G = g,
          K = k,
          M = m
        )

        sample_list = create_sample_list(
          iters = 10, model_data = model_data, alpha_prior = "cusp",
          init_list = NULL,
          seed = 1
        )

        epsilon = sample_list$epsilon[1,,]
        alpha = sample_list$alpha[1,,] |> cbind()
        alpha_precision = sample_list$alpha_precision[1,,] |> cbind()
        psi = sample_list$psi[1,]
        z = sample_list$z[1,]
        mu = sample_list$mu[1,,] |> matrix(nrow = g, ncol = k)
        sigma = sample_list$sigma[1,]

        out = update_chain(
          H = k,
          G = g,
          M = m,
          epsilon = epsilon,
          alpha = alpha,
          alpha_precision = alpha_precision,
          psi = psi,
          z = z,
          mu = mu,
          sigma = sigma,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output is a list with expected elements
        expect_type(out, "list")
        expect_named(out, c("theta", "epsilon", "alpha", "psi", "z", "mu", "sigma"))

      }
    }
  }

})
