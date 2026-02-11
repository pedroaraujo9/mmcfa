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
          G = k,
          K = g,
          M = m
        )

        sample_list = create_sample_list(
          iters = 10, model_data = model_data, alpha_prior = "cusp",
          init_list = NULL,
          seed = 1
        )

        theta = sample_list$theta[1,,] |> cbind()
        psi = sample_list$psi[1,]
        prior_precision = sample_list$alpha_precision[1,,] |> cbind()

        alpha = update_alpha(
          theta = theta,
          psi = psi,
          prior_precision = prior_precision,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

      }
    }
  }

})
