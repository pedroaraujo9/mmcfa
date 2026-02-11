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

        alpha = sample_list$alpha[1,,] |> cbind()
        theta = sample_list$theta[1,,] |> cbind()

        psi = update_psi(
          alpha = alpha,
          theta = theta,
          a = 1,
          b = 1,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output properties
        expect_length(psi, model_data$dims$J)
        expect_true(all(psi > 0))  # psi should be positive (variance)
        expect_type(psi, "double")

      }
    }
  }

})

test_that("function returns correct dimensions", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 2,
    K = 3,
    M = 2
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  alpha = sample_list$alpha[1,,] |> cbind()
  theta = sample_list$theta[1,,] |> cbind()

  psi = update_psi(
    alpha = alpha,
    theta = theta,
    a = 0.01,
    b = 0.01,
    model_data = model_data
  )

  expect_length(psi, model_data$dims$J)

})
