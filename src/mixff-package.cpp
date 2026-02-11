// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec propose_coef_rcpp(const arma::vec& y,
                            const arma::mat& X,
                            const arma::vec& y_prec,
                            const arma::mat& X_prec,
                            const arma::vec& prior_mean,
                            const arma::mat& prior_precision) {

  arma::mat weighted_X = X % X_prec;
  arma::mat aux_cross = X.t() * weighted_X;
  arma::vec rhs = X.t() * (y % y_prec);

  arma::mat post_cov = arma::inv_sympd(prior_precision + aux_cross);
  arma::vec post_mean = post_cov * (rhs + prior_precision * prior_mean);

  arma::vec z = arma::randn<arma::vec>(X.n_cols);
  arma::mat R = arma::chol(post_cov);

  return post_mean + R.t() * z;

}

// [[Rcpp::export]]
arma::mat sample_theta_cpp(
    const arma::mat& y,
    const arma::mat& alpha,
    const arma::mat& prior_mean,
    const arma::vec& psi,
    const bool center
) {
  int n = y.n_rows;
  int K = alpha.n_cols;


  // 2. weighted alpha
  arma::mat alpha_w = alpha;
  alpha_w.each_col() /= psi;      // divide each row j by tau_j

  // 3. posterior covariance
  arma::mat V = inv_sympd(alpha_w.t() * alpha + arma::eye(K, K));
  arma::mat L = chol(V, "upper");

  // 4. posterior mean
  arma::mat m = (y * alpha_w + prior_mean) * V;

  // 5. random noise
  arma::mat e(n, K, arma::fill::randn);

  // 6. posterior draw
  arma::mat theta(n, K);

  if(center == false) {

    theta = m + e * L;

  }else{

    arma::vec one_vec = arma::ones<arma::vec>(n);
    arma::mat A = V * one_vec * inv_sympd(one_vec.t() * V * one_vec) * one_vec.t() * theta;
    theta = theta - A;

  }

  return theta;

}

// [[Rcpp::export]]
arma::mat sample_theta_i_cpp(
    const arma::mat& y,
    const arma::mat& alpha,
    const arma::mat& prior_mean,
    const arma::mat& prior_precision,
    const arma::vec& psi
) {
  int n = y.n_rows;
  int K = alpha.n_cols;


  // 2. weighted alpha
  arma::mat alpha_w = alpha;
  alpha_w.each_col() /= psi;      // divide each row j by tau_j

  // 3. posterior covariance
  arma::mat V = inv_sympd(alpha_w.t() * alpha + prior_precision);
  arma::mat L = chol(V, "upper");

  // 4. posterior mean
  arma::mat m = V * (y * alpha_w + prior_precision * prior_mean);

  // 5. random noise
  arma::mat e(n, K, arma::fill::randn);

  // 6. posterior draw
  arma::mat theta(n, K);
  theta = m + e * L;

  return theta;

}

// [[Rcpp::export]]
arma::mat cpp_compute_V(arma::mat X,
                        arma::vec omega,
                        arma::mat precision_matrix) {

  arma::mat XtOmegaX = X.each_col() % omega;
  XtOmegaX = XtOmegaX.t() * X;
  XtOmegaX += 1e-7 * arma::eye(X.n_cols, X.n_cols);
  arma::mat V = arma::inv(XtOmegaX + precision_matrix);

  return V;
}

// [[Rcpp::export]]
arma::mat cpp_compute_m(arma::mat V,
                        arma::mat X,
                        arma::vec z,
                        arma::vec omega,
                        arma::vec C,
                        arma::vec center,
                        arma::mat inv_cov) {
  arma::mat m = V * (X.t() * ((z - 0.5) + (omega % C)) + inv_cov * center);
  return m;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  // Ensure symmetry
  arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma_sym);
}

// [[Rcpp::export]]
arma::mat sample_beta(arma::mat X,
                      arma::vec omega,
                      arma::mat inv_cov,
                      arma::vec z,
                      arma::vec C,
                      arma::vec center) {

  arma::mat V = cpp_compute_V(X, omega, inv_cov);
  arma::mat mu = cpp_compute_m(V, X, z, omega, C, center, inv_cov);

  return mvrnormArma(1, mu, V).t();

}

// [[Rcpp::export]]
double logsumexp_cpp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}

// [[Rcpp::export]]
arma::mat update_theta2_cpp(arma::mat epsilon, Rcpp::List model_data) {

  // Extract nested list elements
  Rcpp::List theta_spline = model_data["theta_spline"];
  arma::mat R = Rcpp::as<arma::mat>(theta_spline["R"]);

  Rcpp::List dims = model_data["dims"];
  int n_id = dims["n_id"];

  Rcpp::List data_list = model_data["data"];
  arma::vec id = Rcpp::as<arma::vec>(data_list["id"]);
  arma::vec id_unique = Rcpp::as<arma::vec>(data_list["id_unique"]);

  // Initialize theta as a copy of epsilon
  arma::mat theta = epsilon;

  // Loop through each unique ID
  for(int i = 0; i < n_id; i++) {
    double id_i = id_unique(i);

    // Explicitly call arma::find and arma::uvec
    arma::uvec idx = arma::find(id == id_i);

    // Matrix multiplication and sub-matrix assignment
    theta.rows(idx) = R * epsilon.rows(idx);
  }

  return theta;
}

// [[Rcpp::export]]
arma::mat post_epsilon_cpp(arma::mat prec_prior,
                           arma::mat prec_data,
                           arma::mat MU_scaled,
                           arma::mat Rty_alpha_scaled) {

  int n = prec_prior.n_rows;
  arma::vec z = arma::randn<arma::vec>(n);

  arma::mat post_cov = arma::inv_sympd(prec_prior + prec_data);
  arma::mat V = arma::chol(post_cov, "lower");

  arma::mat post_center = post_cov * arma::vectorise(MU_scaled + Rty_alpha_scaled);
  arma::mat epsilon_vec = post_center + V * z;

  return epsilon_vec;

}

// [[Rcpp::export]]
arma::vec post_epsilon_cpp2(const arma::mat& prec_prior,
                            const arma::mat& prec_data,
                            const arma::mat& MU_scaled,
                            const arma::mat& Rty_alpha_scaled) {

  int n = prec_prior.n_rows;

  // Cholesky of precision (avoid explicit inverse)
  arma::mat L = arma::chol(prec_prior + prec_data, "lower");

  // RHS for mean
  arma::vec rhs = arma::vectorise(MU_scaled + Rty_alpha_scaled);

  // Solve for posterior mean via triangular solves: L * L^T * mu = rhs
  arma::vec post_mean = arma::solve(arma::trimatl(L), rhs);
  post_mean = arma::solve(arma::trimatu(L.t()), post_mean);

  // Sample noise: if prec = L*L^T, then Sigma = L^{-T}*L^{-1}, so L^{-T}*z ~ N(0, Sigma)
  arma::vec z = arma::randn<arma::vec>(n);
  arma::vec noise = arma::solve(arma::trimatu(L.t()), z);

  return post_mean + noise;

}


