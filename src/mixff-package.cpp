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
                        arma::mat inv_cov) {

  arma::mat XtOmegaX = X.each_col() % omega;
  XtOmegaX = XtOmegaX.t() * X;
  XtOmegaX += 1e-7 * arma::eye(X.n_cols, X.n_cols);
  arma::mat V = arma::inv(XtOmegaX + inv_cov);

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
