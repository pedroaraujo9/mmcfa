gen_normal_mat = function(n, k) {
  matrix(rnorm(n * k), nrow = n, ncol = k)
}
