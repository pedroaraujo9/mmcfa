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

        mu = sample_list$mu[1,,] |> matrix(nrow = g, ncol = k)
        alpha = sample_list$alpha[1,,] |> cbind()
        psi = sample_list$psi[1,]
        z = sample_list$z[1,]
        sigma = rep(1, g)  # one sigma per cluster

        epsilon = update_epsilon(
          mu = mu,
          sigma = sigma,
          alpha = alpha,
          psi = psi,
          z = z,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output properties
        expect_equal(nrow(epsilon), model_data$dims$n)
        expect_equal(ncol(epsilon), ncol(alpha))
        expect_type(epsilon, "double")

      }
    }
  }

})
