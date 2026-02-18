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
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  // Ensure symmetry
  arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma_sym);
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
                        arma::vec C) {
  arma::mat m = V * (X.t() * ((z - 0.5) + (omega % C)));
  return m;
}

// [[Rcpp::export]]
arma::mat sample_beta(arma::mat X,
                      arma::vec omega,
                      arma::mat inv_cov,
                      arma::vec z,
                      arma::vec C) {

  arma::mat V = cpp_compute_V(X, omega, inv_cov);
  arma::mat mu = cpp_compute_m(V, X, z, omega, C);

  return mvrnormArma(1, mu, V).t();

}

// [[Rcpp::export]]
arma::mat sample_beta2(arma::mat X,
                       arma::vec omega,
                       arma::mat precision_matrix,
                       arma::vec z,
                       arma::vec C) {

  arma::mat X_omega = X.each_col() % omega;
  arma::mat post_prec = (X_omega.t() * X) + precision_matrix;

  arma::mat L_prec = arma::chol(post_prec, "lower");

  arma::vec rhs = X.t() * ((z - 0.5) + (omega % C));
  arma::vec std_noise = arma::randn<arma::vec>(rhs.n_elem);

  return arma::solve(arma::trimatl(L_prec), rhs + std_noise);

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
Rcpp::List utils_comp_epsilon(arma::mat& alpha,
                              arma::mat& alpha_scaled_psi,
                              arma::mat& RtR,
                              arma::mat& Rty) {

  arma::mat alphat_psi_alpha = alpha_scaled_psi.t() * alpha;
  arma::mat prec_data = arma::kron(alphat_psi_alpha, RtR);
  arma::mat Rty_alpha_scaled = Rty * alpha_scaled_psi;

  return Rcpp::List::create(
    Rcpp::Named("prec_data") = prec_data,
    Rcpp::Named("Rty_alpha_scaled") = Rty_alpha_scaled
  );
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

  arma::vec z = arma::randn<arma::vec>(n);
  arma::vec noise = arma::solve(arma::trimatu(L.t()), z);

  return post_mean + noise;

}

// [[Rcpp::export]]
arma::mat update_epsilon_cpp(arma::mat& alpha,
                             arma::mat& alpha_scaled_psi,
                             arma::mat& RtR,
                             arma::mat& Rty,
                             arma::mat& MU_scaled,
                             arma::vec& inv_sigma,
                             arma::uvec& z,
                             arma::umat& idx) {

  int n = MU_scaled.n_rows;
  int K = alpha.n_cols;
  int n_id = idx.n_rows;
  int T = idx.n_cols;

  arma::mat alphat_psi_alpha = alpha_scaled_psi.t() * alpha;
  arma::mat prec_data = arma::kron(alphat_psi_alpha, RtR);
  arma::mat Rty_alpha_scaled = Rty * alpha_scaled_psi;

  arma::mat epsilon(n, K);

  for(int i=0; i < n_id; i++) {
    arma::uvec idx_i = idx.row(i).t();
    arma::uvec z_i = z(idx_i);
    arma::mat prec_prior_i = arma::diagmat(arma::repmat(inv_sigma(z_i), K, 1));

    arma::mat L = arma::chol(prec_prior_i + prec_data, "lower");
    arma::vec rhs = arma::vectorise(MU_scaled.rows(idx_i) + Rty_alpha_scaled.rows(idx_i));
    arma::vec post_mean = arma::solve(arma::trimatl(L), rhs);

    arma::vec std_noise = arma::randn<arma::vec>(T * K);
    arma::vec noise = arma::solve(arma::trimatu(L.t()), std_noise);
    arma::vec epsilon_i = post_mean + noise;

    epsilon.rows(idx_i) = arma::reshape(epsilon_i, T, K);

  }

  return epsilon;

}

// [[Rcpp::export]]
arma::mat update_epsilon_cpp_fast(arma::mat& alpha,
                                  arma::mat& alpha_scaled_psi,
                                  arma::mat& RtR,
                                  arma::mat& Rty,
                                  arma::mat& MU_scaled,
                                  arma::vec& inv_sigma,
                                  arma::uvec& z,
                                  arma::umat& idx) {

  int K = alpha.n_cols;
  int n_id = idx.n_rows;
  int T = idx.n_cols;
  int n = MU_scaled.n_rows;

  // Precompute constant data components
  // Instead of kron, we keep the components separate
  arma::mat alphat_psi_alpha = alpha_scaled_psi.t() * alpha;
  arma::mat Rty_alpha_scaled = Rty * alpha_scaled_psi;

  // Precompute Cholesky of data precision components if possible
  // But since prec_prior changes per 'i', we optimize the solve
  arma::mat epsilon(n, K);

  for(int i=0; i < n_id; i++) {
    arma::uvec idx_i = idx.row(i).t();

    arma::vec inv_sigma_i = inv_sigma(z(idx_i));
    arma::mat B = MU_scaled.rows(idx_i) + Rty_alpha_scaled.rows(idx_i);

    arma::mat prec_total = arma::kron(alphat_psi_alpha, RtR);
    prec_total.diag() += arma::repmat(inv_sigma_i, K, 1);

    arma::mat L = arma::chol(prec_total, "lower");

    arma::vec rhs = arma::vectorise(B);
    // Solve L*y = rhs -> L'*epsilon = y + noise
    arma::vec sol = arma::solve(arma::trimatl(L), rhs + arma::randn<arma::vec>(T * K));
    arma::vec epsilon_i = arma::solve(arma::trimatu(L.t()), sol);

    epsilon.rows(idx_i) = arma::reshape(epsilon_i, T, K);
  }

  return epsilon;
}

// [[Rcpp::export]]
arma::mat fast_dummy_dense(arma::ivec x, int G) {
  // x: vector of cluster/category assignments (0-indexed)
  // G: number of unique categories
  int N = x.n_elem;
  arma::mat out(N, G, arma::fill::zeros);

  for(int i = 0; i < N; ++i) {
    out(i, x(i)-1) = 1.0;
  }

  return out;
}

arma::mat sofmax_cpp(const arma::mat& x) {
  arma::mat ex = arma::exp(x);
  arma::vec row_sum = arma::sum(ex, 1);
  ex.each_col() /= row_sum;
  arma::mat P = ex;
  return P;
}

// [[Rcpp::export]]
arma::mat compute_prob_group(arma::mat& B,
                             arma::mat& beta_group,
                             arma::uvec& idx) {
  arma::mat prob = sofmax_cpp(B * beta_group);
  return prob.rows(idx);
}

// [[Rcpp::export]]
arma::mat predict_prob_cpp(int& M,
                           arma::ivec& w,
                           arma::mat& B,
                           arma::mat& beta) {

  arma::mat W = fast_dummy_dense(w, M);
  arma::mat X = arma::kron(W, B);
  arma::mat prob = sofmax_cpp(X * beta);
  return prob;
}

// [[Rcpp::export]]
arma::vec fast_aggregate_sum(arma::vec& log_pz, arma::ivec& id) {
  // 1. Find the range of IDs
  int min_id = id.min();
  int max_id = id.max();
  int range = max_id - min_id + 1;

  // 2. Initialize a result vector with zeros
  arma::vec sums(range, fill::zeros);

  // 3. Single-pass accumulation (The "O(n)" magic)
  for (unsigned int i = 0; i < id.n_elem; ++i) {
    // Offset by min_id so it starts at index 0
    sums(id(i) - min_id) += log_pz(i);
  }

  return sums;
}



// [[Rcpp::export]]
arma::mat update_theta_cpp(const arma::mat& epsilon,
                           const arma::mat& R,
                           const arma::vec& id,
                           const arma::vec& id_unique) {

  int n_id = id_unique.n_elem;
  arma::mat theta = epsilon;

  for(int i = 0; i < n_id; i++) {

    arma::uvec idx = arma::find(id == id_unique(i));
    theta.rows(idx) = R * epsilon.rows(idx);
  }

  return theta;
}

// arma::mat update_epsilon_t(const arma::mat& R,
//                            const arma::mat& epsilon,
//                            const arma::mat& alpha,
//                            const arma::mat& mu_t,
//                            const arma::mat& St,
//                            const arma::mat& y_t,
//                            const arma::mat& Rtt2,
//                            const arma::mat& A,
//                            const arma::mat& alpha_scaled_psi) {
//
//   arma::mat ct = R * epsilon;
//   arma::mat rt = y_t - ct * alpha.t();
//   arma::mat m = St * (mu_t + sqrt(Rtt2) * alpha_scaled_psi.t() * rt.t());
//   arma::mat epsilon_t = m + t(chol(S)) %*% rnorm(K)
//
// }
