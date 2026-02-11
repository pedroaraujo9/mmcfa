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

        epsilon = sample_list$epsilon[1,,] |> cbind()
        mu = sample_list$mu[1,,] |> matrix(nrow = g, ncol = k)
        sigma = sample_list$sigma[1,]

        z = update_z(
          mu = mu,
          epsilon = epsilon,
          sigma = sigma,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output properties
        expect_length(z, model_data$dims$n)
        expect_true(all(z >= 1 & z <= g))
        expect_type(z, "integer")

      }
    }
  }

})
