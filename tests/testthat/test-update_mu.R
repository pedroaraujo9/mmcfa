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
        z = sample_list$z[1,]

        mu = update_mu(
          epsilon = epsilon,
          z = z,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output properties
        expect_equal(nrow(mu), model_data$dims$G)
        expect_equal(ncol(mu), ncol(epsilon))
        expect_type(mu, "double")

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
    G = 4,
    K = 3,
    M = 2
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  epsilon = sample_list$epsilon[1,,] |> cbind()
  z = sample_list$z[1,]

  mu = update_mu(
    epsilon = epsilon,
    z = z,
    model_data = model_data
  )

  expect_equal(dim(mu), c(4, 3))  # G x K

})

test_that("function handles empty clusters", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 3,
    K = 2,
    M = 2
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  epsilon = sample_list$epsilon[1,,] |> cbind()

  # Create z where all observations are in cluster 1 (clusters 2 and 3 empty)
  n = model_data$dims$n
  z = rep(1, n)

  mu = update_mu(
    epsilon = epsilon,
    z = z,
    model_data = model_data
  ) |>
    expect_no_error()

  # Should still return mu for all G clusters
  expect_equal(nrow(mu), 3)
  expect_equal(ncol(mu), 2)

})

test_that("function handles single factor case (K = 1)", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 2,
    K = 1,
    M = 2
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  epsilon = sample_list$epsilon[1,,] |> cbind()
  z = sample_list$z[1,]

  mu = update_mu(
    epsilon = epsilon,
    z = z,
    model_data = model_data
  ) |>
    expect_no_error()

  expect_equal(dim(mu), c(2, 1))

})
