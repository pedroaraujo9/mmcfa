test_that("function works across different dimensions", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)
  iters = 100

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
          iters = iters,
          model_data = model_data,
          alpha_prior = "normal",
          init_list = NULL,
          seed = 1
        ) |>
          expect_no_error() |>
          expect_no_message() |>
          expect_no_warning()

        sample_list = create_sample_list(
          iters = iters,
          model_data = model_data,
          alpha_prior = "cusp",
          init_list = NULL,
          seed = 1
        ) |>
          expect_no_error() |>
          expect_no_message() |>
          expect_no_warning()

      }
    }
  }

})
