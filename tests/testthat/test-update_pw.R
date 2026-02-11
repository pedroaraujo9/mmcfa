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

        w = sample_list$w[1,]

        pw = update_pw(
          w = w,
          epsilon = 1,
          model_data = model_data
        ) |>
          expect_no_error() |>
          expect_no_message()

        # Check output properties
        expect_length(pw, model_data$dims$M)
        expect_true(all(pw > 0))
        expect_true(all(pw <= 1))
        expect_equal(sum(pw), 1, tolerance = 1e-10)  # probabilities sum to 1
        expect_type(pw, "double")

      }
    }
  }

})

test_that("function handles M = 1 case", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 2,
    K = 2,
    M = 1
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  w = sample_list$w[1,]

  pw = update_pw(
    w = w,
    epsilon = 1,
    model_data = model_data
  )

  expect_equal(pw, 1)

})

test_that("function handles empty clusters", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 2,
    K = 2,
    M = 3
  )

  # Create w where all observations are in cluster 1 (clusters 2 and 3 empty)
  n_id = model_data$dims$n_id
  w = rep(1, n_id)
  names(w) = model_data$data$id_unique

  pw = update_pw(
    w = w,
    epsilon = 1,
    model_data = model_data
  ) |>
    expect_no_error()

  # Should still return probabilities for all M clusters

  expect_length(pw, 3)
  expect_true(all(pw > 0))
  expect_equal(sum(pw), 1, tolerance = 1e-10)

})

test_that("function respects epsilon parameter", {

  set.seed(1)
  sim_data = simulate_data(seed = 1)

  model_data = create_model_data(
    y = sim_data$y,
    id = sim_data$id,
    time = sim_data$time,
    G = 2,
    K = 2,
    M = 3
  )

  sample_list = create_sample_list(
    iters = 10, model_data = model_data, alpha_prior = "normal",
    init_list = NULL,
    seed = 1
  )

  w = sample_list$w[1,]

  # Both should work without error
  pw_small = update_pw(w = w, epsilon = 0.1, model_data = model_data) |>
    expect_no_error()

  pw_large = update_pw(w = w, epsilon = 10, model_data = model_data) |>
    expect_no_error()

  expect_length(pw_small, model_data$dims$M)
  expect_length(pw_large, model_data$dims$M)

})
