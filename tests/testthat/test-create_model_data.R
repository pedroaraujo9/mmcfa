test_that("function works", {

  K = 2
  G = 4
  M = 3
  center = TRUE
  scale = FALSE
  n_basis = 10
  theta_spline_penalty = 1
  theta_intercept_penalty = 1
  z_spline_penalty = 1
  z_intercept_penalty = 1
  w_dirichlet = 1

  for(k in 1:3) {
    for(g in 1:2) {
      for(m in 1:2) {

        model_data = create_model_data(
          y = sim_data$y,
          id = sim_data$id,
          time = sim_data$time,
          K = k,
          G = g,
          M = m,
          center = center,
          scale = scale,
          n_basis = n_basis,
          theta_spline_penalty = theta_spline_penalty,
          theta_intercept_penalty = theta_intercept_penalty,
          z_spline_penalty = z_spline_penalty,
          z_intercept_penalty = z_intercept_penalty,
          w_dirichlet = w_dirichlet
        ) |>
          expect_no_error() |>
          expect_no_message() |>
          expect_type("list")

      }
    }
  }



})
