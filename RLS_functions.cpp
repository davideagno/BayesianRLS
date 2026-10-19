#include <RcppArmadillo.h>
#include <thread>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// RLS function
// [[Rcpp::export]]
Rcpp::List rls(const vec& y, const mat& X, vec beta, mat P) {
  const vec beta_0 = beta;
  const mat P_0 = P;
  const int n = y.n_elem;
  for (int i = 0; i < n; i++) {
    vec x_i = X.row(i).t();
    vec Px = P * x_i;
    double a = 1.0 / (1.0 + as_scalar(x_i.t() * Px));
    beta += Px * (a * (y(i) - as_scalar(x_i.t() * beta)));
    P -= a * (Px * Px.t());
  }
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("P") = P,
    Rcpp::Named("y") = y,
    Rcpp::Named("X") = X,
    Rcpp::Named("beta_0") = beta_0,
    Rcpp::Named("P_0") = P_0
  );
}

// [[Rcpp::export]]
Rcpp::List rls_lambda(const arma::vec& y,
                      const arma::mat& X,
                      const arma::vec& beta_0,
                      const arma::vec& lambda,
                      const Rcpp::List& Omega) {
  
  // the P_0 on block j is lambda_j * Omega_j
  
  const uword p_blocks = Omega.size();
  if (lambda.n_elem != p_blocks) {
    Rcpp::stop("length(lambda) must equal length(Omega)");
  }
  
  // Build P_0 = blockdiag( 1e10, (lambda_1 Omega_1 + eps I)^{-1}, ... )
  mat P_0(X.n_cols, X.n_cols, fill::zeros);
  P_0(0, 0) = 1e10;
  
  const double eps = 1e-10;
  uword offset = 1;
  for (uword j = 0; j < p_blocks; ++j) {
    // Pull the j-th penalty matrix out of the R list and convert to arma::mat.
    mat Oj = Rcpp::as<mat>(Omega[j]);
    const uword k = Oj.n_rows;
    
    // Penalty matrices must be square.
    if (Oj.n_cols != k) Rcpp::stop("Omega[[%d]] is not square", (int)j + 1);
    
    // Form the regularized prior precision for this block: lambda_j * Omega_j + eps*I.
    mat A = lambda(j) * Oj;
    A.diag() += eps;
    
    // Invert to get the prior covariance for this block.
    // Try Cholesky-based sympd inversion first (faster & numerically nicer
    // for symmetric positive-definite matrices); fall back to a general
    // inverse if the factorization fails.
    mat Ainv;
    if (!inv_sympd(Ainv, A)) {
      Ainv = inv(A);          // fallback if sympd factorization fails
    }
    
    // Place this block's covariance on the diagonal of P_0.
    P_0.submat(offset, offset, offset + k - 1, offset + k - 1) = Ainv;
    offset += k;
  }
  
  // Sanity check: the blocks must exactly fill columns 1..(p-1) of X.
  if (offset != X.n_cols) {
    Rcpp::stop("block sizes in Omega don't match ncol(X) - 1");
  }
  
  // Use rls function to compute \beta_n and P_n
  arma::vec out = Rcpp::as<arma::vec>(rls(y, X, beta_0, P_0)["beta"]);
  arma::mat P = Rcpp::as<arma::mat>(rls(y, X, beta_0, P_0)["P"]);
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = out,
    Rcpp::Named("P") = P,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("Omega") = Omega,
    Rcpp::Named("y") = y,
    Rcpp::Named("X") = X,
    Rcpp::Named("beta_0") = beta_0,
    Rcpp::Named("P_0") = P_0
  );
}


// [[Rcpp::export]]
double pll_rls_t(const vec& y, const mat& X, vec beta, mat P, double sigma, double df) {
  double out = 0.0;
  const vec beta_0 = beta;
  const mat P_0 = P;
  const int n = y.n_elem;
  for (int i = 0; i < n; i++) {
    vec x_i = X.row(i).t();
    double mu = sum(dot(x_i, beta));
    out += R::dt((y(i) - mu) / sigma, df, 1) - std::log(sigma);
    
    vec Px = P * x_i;
    double a = 1.0 / (1.0 + as_scalar(x_i.t() * Px));
    beta += Px * (a * (y(i) - as_scalar(x_i.t() * beta)));
    P -= a * (Px * Px.t());
  }
  return out;
}


static inline double dt_log(double z, double df) {
  return R::lgammafn(0.5 * (df + 1.0))
  - R::lgammafn(0.5 * df)
  - 0.5 * std::log(df * M_PI)
  - 0.5 * (df + 1.0) * std::log1p(z * z / df);
}

// [[Rcpp::export]]
double pll_rls_t_lambda(double df,
                        double sigma,
                        const arma::vec& lambda,
                        const arma::vec& y,
                        const arma::mat& X,
                        const arma::vec& beta_0,
                        const Rcpp::List& Omega) {
  const uword p_blocks = Omega.size();
  if (lambda.n_elem != p_blocks) {
    Rcpp::stop("length(lambda) must equal length(Omega)");
  }
  if (sigma <= 0.0) Rcpp::stop("sigma must be positive");
  if (df    <= 0.0) Rcpp::stop("df must be positive");
  
  // Build P = blockdiag( 1e10,  (lambda_1 Omega_1 + eps I)^{-1},  ... )
  mat P(X.n_cols, X.n_cols, fill::zeros);
  P(0, 0) = 1e10;
  
  const double eps = 1e-10;
  uword offset = 1;
  for (uword j = 0; j < p_blocks; ++j) {
    mat Oj = Rcpp::as<mat>(Omega[j]);
    const uword k = Oj.n_rows;
    if (Oj.n_cols != k) Rcpp::stop("Omega[[%d]] is not square", (int)j + 1);
    
    mat A = lambda(j) * Oj;
    A.diag() += eps;
    
    // Symmetric positive-definite solve: Cholesky, then invert.
    mat Ainv;
    if (!inv_sympd(Ainv, A)) {
      // Fall back to general inverse if sympd factorization fails.
      Ainv = inv(A);
    }
    
    P.submat(offset, offset, offset + k - 1, offset + k - 1) = Ainv;
    offset += k;
  }
  if (offset != X.n_cols) {
    Rcpp::stop("block sizes in Omega don't match ncol(X) - 1");
  }
  
  vec beta = beta_0;
  const uword n = y.n_elem;
  const double log_sigma = std::log(sigma);
  
  double out = 0.0;
  for (uword i = 0; i < n; ++i) {
    vec x_i = X.row(i).t();
    
    const double mu  = dot(x_i, beta);
    const double res = y(i) - mu;
    
    out += dt_log(res / sigma, df) - log_sigma;
    
    vec Px = P * x_i;
    const double a = 1.0 / (1.0 + dot(x_i, Px));
    beta += Px * (a * res);
    P    -= a * (Px * Px.t());
  }
  return out;
}



// // ---------------------------------------------------------------------------
// // dt_log(): log density of the standard (location 0, scale 1) Student-t.
// // If dt_log is already defined in the same translation unit, delete this block
// // or compile with -DDT_LOG_DEFINED.
// // ---------------------------------------------------------------------------
// #ifndef DT_LOG_DEFINED
// inline double dt_log(double x, double df) {
//   return R::lgammafn(0.5 * (df + 1.0)) - R::lgammafn(0.5 * df)
//   - 0.5 * std::log(df * M_PI)
//   - 0.5 * (df + 1.0) * std::log1p(x * x / df);
// }
// #endif

namespace {

// Build P = blockdiag( 1e10, (lambda_1 Omega_1 + eps I)^{-1}, ... ) once.
// This depends only on (lambda, Omega), not on the permutation, so it is
// computed a single time and copied for each permutation.
mat build_P(const vec& lambda, const std::vector<mat>& Om, const uword p) {
  
  mat P(p, p, fill::zeros);
  P(0, 0) = 1e10;
  
  const double eps = 1e-10;
  uword offset = 1;
  
  for (uword j = 0; j < Om.size(); ++j) {
    const mat& Oj = Om[j];
    const uword k = Oj.n_rows;
    
    if (offset + k > p) Rcpp::stop("block sizes in Omega don't match ncol(X) - 1");
    
    mat A = lambda(j) * Oj;
    A.diag() += eps;
    
    mat Ainv;
    if (!inv_sympd(Ainv, A)) {
      Ainv = inv(A);              // fall back to a general inverse
    }
    
    P.submat(offset, offset, offset + k - 1, offset + k - 1) = Ainv;
    offset += k;
  }
  
  if (offset != p) Rcpp::stop("block sizes in Omega don't match ncol(X) - 1");
  return P;
}
  
  // Recursive-least-squares predictive log-likelihood over the rows listed in
  // idx (0-based). Xt is the *transpose* of X, so Xt.col(r) is a contiguous
  // column: this is much cache-friendlier than X.row(r) on a column-major matrix.
  // beta and P are taken by value so each permutation starts from a fresh copy.
  double pll_core(const double df,
                  const double sigma,
                  const vec&   y,
                  const mat&   Xt,
                  const uvec&  idx,
                  vec          beta,
                  mat          P) {
    
    const double log_sigma = std::log(sigma);
    const uword  n = idx.n_elem;
    const uword  p = Xt.n_rows;
    
    vec Px(p);
    double out = 0.0;
    
    for (uword i = 0; i < n; ++i) {
      const uword r = idx(i);
      const subview_col<double> x_i = Xt.col(r);
      
      const double res = y(r) - dot(x_i, beta);
      
      out += dt_log(res / sigma, df) - log_sigma;
      
      Px = P * x_i;
      const double a = 1.0 / (1.0 + dot(x_i, Px));
      
      beta += (a * res) * Px;
      P    -= (a * Px) * Px.t();
    }
    
    return out;
  }
  
} // anonymous namespace



// [[Rcpp::export]]
vec ps_rs_single_iteration_t_1(const mat& X, vec beta, mat P,
                               int N, double sigma, double df) {
  const int n = X.n_rows;
  for (int i = 0; i < N; i++) {
    int idx = randi<int>(distr_param(0, n - 1));
    vec x_new = X.row(idx).t();
    double t_new = R::rt(df);
    double y_new = as_scalar(x_new.t() * beta) + t_new * sigma;
    vec Px = P * x_new;
    double a = 1.0 / (1.0 + as_scalar(x_new.t() * Px));
    beta += Px * (a * (y_new - as_scalar(x_new.t() * beta)));
    P -= a * (Px * Px.t());
  }
  return beta;
}


//' P-splines Recursive Sampling
 //'
 //' @param X List of matrices (cubic splines regressors for each variable)
 //' @param beta List of vectors (one for each dimension)
 //' @param P List of matrices (one for each dimension)
 //' @param intercept Real number
 //' @param N Positive integer (length of forward samples)
 //' @param B Positive integer (sample size)
 //' @param sigma Positive number (standard deviation)
 //' @param seed Integer seed for random number generation
 //' @return List containing input parameters and beta_ps samples
 //' @export
 // [[Rcpp::export]]
 List ps_rs_cpp(List X, List beta, List P, double intercept, 
                int N, int B, double sigma, int seed = 1) {
   
   // Set random seed
   Rcpp::Environment base_env("package:base");
   Rcpp::Function set_seed = base_env["set.seed"];
   set_seed(seed);
   
   // Get dimensions
   mat X0 = as<mat>(X[0]);
   int n = X0.n_rows;
   int p = beta.size();
   
   // Initialize beta_ps - list of matrices to store samples
   List beta_ps(p);
   for (int j = 0; j < p; j++) {
     vec beta_j = as<vec>(beta[j]);
     beta_ps[j] = mat(B, beta_j.n_elem);
   }
   
   // Main sampling loop
   for (int b = 0; b < B; b++) {
     Rcout << b + 1 << std::endl;
     
     // Initialize beta_b and P_b for this iteration
     List beta_b = clone(beta);
     List P_b = clone(P);
     
     // Forward sampling loop
     for (int i = 0; i < N; i++) {
       
       // Sample new observation index
       int idx_new = as<int>(sample(n, 1)) - 1; // Convert to 0-based indexing
       
       // Extract x_new for all dimensions
       List x_new(p);
       for (int j = 0; j < p; j++) {
         mat X_j = as<mat>(X[j]);
         x_new[j] = trans(X_j.row(idx_new));
       }
       
       // Calculate mean for new observation
       double mean_new = 0.0;
       for (int j = 0; j < p; j++) {
         vec x_j = as<vec>(x_new[j]);
         vec beta_j = as<vec>(beta_b[j]);
         mean_new += dot(x_j, beta_j);
       }
       
       // Sample new y
       double y_new = mean_new + R::rnorm(0.0, sigma);
       
       // Create shuffled indices for updating
       IntegerVector idx_shuffle = sample(p, p, false);
       
       // Update coefficients in random order
       for (int idx_j = 0; idx_j < p; idx_j++) {
         int j = idx_shuffle[idx_j] - 1; // Convert to 0-based indexing
         
         // Calculate residual
         double r_ij = y_new;
         for (int idx_k = 0; idx_k < p; idx_k++) {
           int k = idx_shuffle[idx_k] - 1;
           if (k != j) {
             vec x_k = as<vec>(x_new[k]);
             vec beta_k = as<vec>(beta_b[k]);
             r_ij -= dot(x_k, beta_k);
           }
         }
         
         // Get current values
         vec x_ij = as<vec>(x_new[j]);
         mat P_ij = as<mat>(P_b[j]);
         vec beta_j = as<vec>(beta_b[j]);
         
         // Calculate update factor
         double a_ij = 1.0 / (1.0 + as_scalar(trans(x_ij) * P_ij * x_ij));
         
         // Calculate old prediction
         double old_f = dot(x_ij, beta_j);
         
         // Update beta
         beta_j = beta_j + a_ij * P_ij * x_ij * (r_ij - old_f);
         beta_b[j] = beta_j;
         
         // Update P
         P_ij = P_ij - a_ij * P_ij * x_ij * trans(x_ij) * P_ij;
         P_b[j] = P_ij;
       }
     }
     
     // Store results for this iteration
     for (int j = 0; j < p; j++) {
       vec beta_j = as<vec>(beta_b[j]);
       mat beta_ps_j = as<mat>(beta_ps[j]);
       beta_ps_j.row(b) = trans(beta_j);
       beta_ps[j] = beta_ps_j;
     }
   }
   
   // Return results
   return List::create(
     Named("X") = X,
     Named("beta") = beta,
     Named("P") = P,
     Named("intercept") = intercept,
     Named("B") = B,
     Named("N") = N,
     Named("sigma") = sigma,
     Named("seed") = seed,
     Named("beta_ps") = beta_ps
   );
 }


// [[Rcpp::export]]
List ps_rs_single_iteration(List X, List beta, List P, 
                            int N, double sigma, int seed) {
  
  // Get dimensions
  mat X0 = as<mat>(X[0]);
  int n = X0.n_rows;
  int p = beta.size();
  
  // Convert Lists to vectors of matrices for easier access
  std::vector<mat> X_vec(p);
  std::vector<vec> beta_b(p);
  std::vector<mat> P_b(p);
  
  for (int j = 0; j < p; j++) {
    X_vec[j] = as<mat>(X[j]);
    beta_b[j] = as<vec>(beta[j]);
    P_b[j] = as<mat>(P[j]);
  }
  
  // Set seed for this iteration
  arma::arma_rng::set_seed(seed);
  
  // Forward sampling loop
  for (int i = 0; i < N; i++) {
    
    // Sample new observation index
    int idx_new = arma::randi<int>(arma::distr_param(0, n - 1));
    
    // Extract x_new for all dimensions
    std::vector<vec> x_new(p);
    for (int j = 0; j < p; j++) {
      x_new[j] = trans(X_vec[j].row(idx_new));
    }
    
    // Calculate mean for new observation
    double mean_new = 0.0;
    for (int j = 0; j < p; j++) {
      mean_new += dot(x_new[j], beta_b[j]);
    }
    
    // Sample new y
    double y_new = mean_new + arma::randn() * sigma;
    
    // Create shuffled indices for updating
    uvec idx_shuffle = arma::randperm(p);
    
    // Update coefficients in random order
    for (int idx_j = 0; idx_j < p; idx_j++) {
      int j = idx_shuffle(idx_j);
      
      // Calculate residual
      double r_ij = y_new;
      for (int idx_k = 0; idx_k < p; idx_k++) {
        int k = idx_shuffle(idx_k);
        if (k != j) {
          r_ij -= dot(x_new[k], beta_b[k]);
        }
      }
      
      // Get current values
      vec x_ij = x_new[j];
      mat P_ij = P_b[j];
      vec beta_j = beta_b[j];
      
      // Calculate update factor
      double a_ij = 1.0 / (1.0 + as_scalar(trans(x_ij) * P_ij * x_ij));
      
      // Calculate old prediction
      double old_f = dot(x_ij, beta_j);
      
      // Update beta
      beta_b[j] = beta_j + a_ij * P_ij * x_ij * (r_ij - old_f);
      
      // Update P
      P_b[j] = P_ij - a_ij * P_ij * x_ij * trans(x_ij) * P_ij;
    }
  }
  
  // Convert back to List format
  List beta_result(p);
  for (int j = 0; j < p; j++) {
    beta_result[j] = beta_b[j];
  }
  
  return beta_result;
}




// [[Rcpp::export]]
List ps_rs_single_iteration_t(List X, List beta, List P, 
                              int N, double sigma, double df, int seed) {
  
  // Get dimensions
  mat X0 = as<mat>(X[0]);
  int n = X0.n_rows;
  int p = beta.size();
  
  // Convert Lists to vectors of matrices for easier access
  std::vector<mat> X_vec(p);
  std::vector<vec> beta_b(p);
  std::vector<mat> P_b(p);
  
  for (int j = 0; j < p; j++) {
    X_vec[j] = as<mat>(X[j]);
    beta_b[j] = as<vec>(beta[j]);
    P_b[j] = as<mat>(P[j]);
  }
  
  // Set seed for this iteration
  arma::arma_rng::set_seed(seed);
  
  // Forward sampling loop
  for (int i = 0; i < N; i++) {
    
    // Sample new observation index
    int idx_new = arma::randi<int>(arma::distr_param(0, n - 1));
    
    // Extract x_new for all dimensions
    std::vector<vec> x_new(p);
    for (int j = 0; j < p; j++) {
      x_new[j] = trans(X_vec[j].row(idx_new));
    }
    
    // Calculate mean for new observation
    double mean_new = 0.0;
    for (int j = 0; j < p; j++) {
      mean_new += dot(x_new[j], beta_b[j]);
    }
    
    // Sample new y from Student t distribution
    // Generate standard t and scale by sigma
    double t_sample = R::rt(df);
    double y_new = mean_new + t_sample * sigma;
    
    // Create shuffled indices for updating
    uvec idx_shuffle = arma::randperm(p);
    
    // Update coefficients in random order
    for (int idx_j = 0; idx_j < p; idx_j++) {
      int j = idx_shuffle(idx_j);
      
      // Calculate residual
      double r_ij = y_new;
      for (int idx_k = 0; idx_k < p; idx_k++) {
        int k = idx_shuffle(idx_k);
        if (k != j) {
          r_ij -= dot(x_new[k], beta_b[k]);
        }
      }
      
      // Get current values
      vec x_ij = x_new[j];
      mat P_ij = P_b[j];
      vec beta_j = beta_b[j];
      
      // Calculate update factor
      double a_ij = 1.0 / (1.0 + as_scalar(trans(x_ij) * P_ij * x_ij));
      
      // Calculate old prediction
      double old_f = dot(x_ij, beta_j);
      
      // Update beta
      beta_b[j] = beta_j + a_ij * P_ij * x_ij * (r_ij - old_f);
      
      // Update P
      P_b[j] = P_ij - a_ij * P_ij * x_ij * trans(x_ij) * P_ij;
    }
  }
  
  // Convert back to List format
  List beta_result(p);
  for (int j = 0; j < p; j++) {
    beta_result[j] = beta_b[j];
  }
  
  return beta_result;
}


// [[Rcpp::export]]
Rcpp::List invert_lambda_Omega(const arma::vec& lambda,
                               const Rcpp::List& Omega) {
  const R_xlen_t p = Omega.size();
  if (lambda.n_elem != static_cast<uword>(p)) {
    Rcpp::stop("length(lambda) must equal length(Omega)");
  }
  
  const double eps = 1e-10;
  Rcpp::List P(p);
  
  for (R_xlen_t j = 0; j < p; ++j) {
    mat Oj = Rcpp::as<mat>(Omega[j]);
    if (Oj.n_rows != Oj.n_cols) {
      Rcpp::stop("Omega[[%d]] is not square", (int)j + 1);
    }
    
    mat A = lambda(j) * Oj;
    A.diag() += eps;
    
    mat Ainv;
    if (!inv_sympd(Ainv, A)) {
      Ainv = inv(A);        // fallback if Cholesky-based inversion fails
    }
    
    P[j] = Ainv;
  }
  
  return P;
}


// [[Rcpp::export]]
List est_rs_cpp(NumericVector y, List X, List beta_0, List P_0) {
  const uword n = y.size();
  const uword p = beta_0.size();
  
  // Center y
  vec y_arma    = as<vec>(y);
  double intercept = mean(y_arma);
  vec y_centered   = y_arma - intercept;
  
  // Unpack lists into std::vector for fast indexing
  std::vector<mat> X_vec(p), P(p);
  std::vector<vec> beta(p);
  for (uword j = 0; j < p; ++j) {
    X_vec[j] = as<mat>(X[j]);
    beta[j]  = as<vec>(beta_0[j]);
    P[j]     = as<mat>(P_0[j]);
  }
  
  // Main loop: single pass over observations
  for (uword i = 0; i < n; ++i) {
    // Full prediction at i using current betas
    double f_hat = 0.0;
    for (uword j = 0; j < p; ++j) {
      f_hat += dot(X_vec[j].row(i), beta[j]);
    }
    
    // Random permutation of component indices (via R's RNG for reproducibility)
    IntegerVector perm = sample(p, p) - 1;
    
    for (uword idx_j = 0; idx_j < p; ++idx_j) {
      const uword j = perm[idx_j];
      const vec x_ij = X_vec[j].row(i).t();
      
      vec Px = P[j] * x_ij;
      const double a     = 1.0 / (1.0 + dot(x_ij, Px));
      const double resid = y_centered(i) - f_hat;
      
      const double f_j_old = dot(x_ij, beta[j]);
      beta[j] += a * resid * Px;
      P[j]    -= a * (Px * Px.t());
      const double f_j_new = dot(x_ij, beta[j]);
      
      f_hat += (f_j_new - f_j_old);
    }
  }
  
  // Pack back into R lists
  List beta_out(p), P_out(p);
  for (uword j = 0; j < p; ++j) {
    beta_out[j] = beta[j];
    P_out[j]    = P[j];
  }
  
  return List::create(
    _["y"]         = NumericVector(y_centered.begin(), y_centered.end()),
    _["X"]         = X,
    _["beta_0"]    = beta_0,
    _["P_0"]       = P_0,
    _["beta"]      = beta_out,
    _["P"]         = P_out,
    _["intercept"] = intercept
  );
}


// [[Rcpp::export]]
double pll_rs_t_cpp(double df,
                    double sigma,
                    List X,
                    NumericVector y,
                    List beta_0,
                    List P_0) {
  
  const int n = y.size();
  const int p = beta_0.size();
  
  if (sigma <= 0.0) stop("sigma must be positive");
  if (df    <= 0.0) stop("df must be positive");
  if ((int)X.size() != p || (int)P_0.size() != p) {
    stop("X, beta_0, and P_0 must all have length p");
  }
  
  // Center y once
  double intercept = 0.0;
  for (int i = 0; i < n; ++i) intercept += y[i];
  intercept /= n;
  
  vec y_centered(n);
  for (int i = 0; i < n; ++i) y_centered(i) = y[i] - intercept;
  
  // Unpack lists
  std::vector<mat> X_vec(p), P(p);
  std::vector<vec> beta(p);
  for (int j = 0; j < p; ++j) {
    X_vec[j] = as<mat>(X[j]);
    beta[j]  = as<vec>(beta_0[j]);
    P[j]     = as<mat>(P_0[j]);
  }
  
  const double log_sigma = std::log(sigma);
  double pll = 0.0;
  
  for (int i = 0; i < n; ++i) {
    
    // Full prediction at i (no intercept: y is already centered)
    double f_hat = 0.0;
    for (int j = 0; j < p; ++j) {
      f_hat += dot(X_vec[j].row(i), beta[j]);
    }
    
    // Predictive Student-t log-density BEFORE the update
    const double resid = y_centered(i) - f_hat;
    pll += R::dt(resid / sigma, df, 1) - log_sigma;
    
    // Random permutation of components via R's RNG
    // IntegerVector perm = sample(p, p) - 1;
    uvec perm = arma::randperm(p);
    
    for (int idx_j = 0; idx_j < p; ++idx_j) {
      const int j = perm[idx_j];
      const vec x_ij = X_vec[j].row(i).t();
      
      vec Px = P[j] * x_ij;
      const double a         = 1.0 / (1.0 + dot(x_ij, Px));
      const double resid_now = y_centered(i) - f_hat;
      
      const double f_j_old = dot(x_ij, beta[j]);
      beta[j] += a * resid_now * Px;
      P[j]    -= a * (Px * Px.t());
      const double f_j_new = dot(x_ij, beta[j]);
      
      // Keep running f_hat in sync with the just-updated beta[j]
      f_hat += (f_j_new - f_j_old);
    }
  }
  
  return pll;
}






















// [[Rcpp::export]]
List Sigma_n_lm(NumericMatrix V_n_, NumericMatrix X_future_) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = V_n_.nrow();
  
  if ((arma::uword) V_n_.ncol() != p)
    stop("V_n must be square.");
  if (N > 0 && (arma::uword) X_future_.ncol() != p)
    stop("ncol(X_future) must equal nrow(V_n).");
  
  // Read-only view: X_future is never modified.
  const arma::mat X(X_future_.begin(), N, p, false);
  
  // Updated in the loop, so these own their memory.
  arma::mat Sigma(p, p, arma::fill::zeros);
  arma::mat V_N = as<arma::mat>(V_n_);
  
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    
    Vx = V_N * x_i.t();
    const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
    
    const arma::mat VxVx = Vx * Vx.t();   // tcrossprod(Vx), formed once
    
    Sigma += (a_i * a_i) * VxVx;
    V_N   -= a_i * VxVx;
  }
  
  NumericMatrix Sigma_out = wrap(Sigma);
  NumericMatrix V_out     = wrap(V_N);
  
  SEXP dn = V_n_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    Sigma_out.attr("dimnames") = dn;
    V_out.attr("dimnames")     = dn;
  }
  
  return List::create(_["Sigma"] = Sigma_out,
                      _["V"]     = V_out,
                      _["V_n"]   = V_n_);
}




// [[Rcpp::export]]
NumericMatrix ps_lm_fixed(NumericVector beta_n_, NumericMatrix V_n_,
                          NumericMatrix X_future_, int B, double sigma,
                          double seed = 1.0) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = X_future_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X_future).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X_future).");
  if (B < 1)
    stop("B must be >= 1.");
  
  // Read-only view: X_future is never modified.
  const arma::mat X(X_future_.begin(), N, p, false);
  
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  arma::mat       V_N    = as<arma::mat>(V_n_);
  
  // --- gain sequence, computed once ---------------------------------
  arma::mat aVx(N, p);
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    
    Vx = V_N * x_i.t();
    const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
    
    aVx.row(i) = a_i * Vx.t();
    V_N -= a_i * (Vx * Vx.t());
    
    const arma::mat V_t = V_N.t();          // explicit temp: avoids aliasing
    V_N = 0.5 * (V_N + V_t);
  }
  
  // --- replicates ---------------------------------------------------
  arma::mat beta_out(B, p);
  arma::vec beta(p);
  
  Function set_seed = Environment::base_env()["set.seed"];
  
  for (int b = 1; b <= B; ++b) {
    beta = beta_n;
    set_seed(seed * b);
    
    for (arma::uword i = 0; i < N; ++i) {
      const double mu  = arma::dot(X.row(i), beta);
      const double y_i = R::rnorm(mu, sigma);
      const double e_i = y_i - mu;
      
      beta += e_i * aVx.row(i).t();
    }
    beta_out.row(b - 1) = beta.t();
  }
  
  NumericMatrix out = wrap(beta_out);
  
  SEXP dn = X_future_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  return out;
}





// [[Rcpp::export]]
NumericMatrix ps_lm(NumericVector beta_n_, NumericMatrix V_n_,
                    NumericMatrix X_, int N, int B, double sigma,
                    double seed = 1.0) {
  
  const arma::uword N_obs = X_.nrow();
  const arma::uword p     = X_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X).");
  if (N_obs < 1) stop("X must have at least one row.");
  if (N < 0)     stop("N must be >= 0.");
  if (B < 1)     stop("B must be >= 1.");
  
  // Read-only view: X is never modified.
  const arma::mat X(X_.begin(), N_obs, p, false);
  
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  const arma::mat V_start = as<arma::mat>(V_n_);
  
  arma::mat beta_out(B, p);
  
  arma::vec beta(p);
  arma::mat V_N(p, p);
  arma::vec Vx(p);
  std::vector<arma::uword> idx(N);
  
  Function set_seed = Environment::base_env()["set.seed"];
  set_seed(seed);
  
  const double dn = (double) N_obs;
  
  for (int b = 0; b < B; ++b) {
    beta = beta_n;
    V_N  = V_start;
    
    // Whole index path drawn first, as sample.int() does, so that the
    // uniform and normal draws are consumed in the same order as in R.
    for (int i = 0; i < N; ++i)
      idx[i] = (arma::uword) R_unif_index(dn);
    
    for (int i = 0; i < N; ++i) {
      const arma::rowvec x_i = X.row(idx[i]);
      
      Vx = V_N * x_i.t();
      const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
      
      const double mu  = arma::dot(x_i, beta);
      const double y_i = R::rnorm(mu, sigma);
      const double e_i = y_i - mu;
      
      beta += (a_i * e_i) * Vx;
      V_N  -= a_i * (Vx * Vx.t());
    }
    beta_out.row(b) = beta.t();
  }
  
  NumericMatrix out = wrap(beta_out);
  
  SEXP dnm = X_.attr("dimnames");
  if (!Rf_isNull(dnm)) {
    List dnl(dnm);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  return out;
}






// ---- thread-local RNG: xoshiro256++ seeded via splitmix64 -----------
static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

struct RNG {
  uint64_t s[4];
  double spare;
  bool has_spare;
  
  explicit RNG(uint64_t seed) : spare(0.0), has_spare(false) {
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) s[i] = splitmix64(x);
  }
  static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }
  inline uint64_t next() {
    const uint64_t r = rotl(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;    s[3]  = rotl(s[3], 45);
    return r;
  }
  inline double unif() { return ((next() >> 11) + 0.5) * 0x1.0p-53; }
  inline double norm() {
    if (has_spare) { has_spare = false; return spare; }
    const double u1 = unif(), u2 = unif();
    const double r  = std::sqrt(-2.0 * std::log(u1));
    const double th = 6.28318530717958647692 * u2;
    spare = r * std::sin(th);
    has_spare = true;
    return r * std::cos(th);
  }
  
  // Unbiased draw in [0, n). Mask-and-reject: no modulo bias.
  inline uint64_t bounded(uint64_t n) {
    uint64_t mask = n - 1;
    mask |= mask >> 1;  mask |= mask >> 2;  mask |= mask >> 4;
    mask |= mask >> 8;  mask |= mask >> 16; mask |= mask >> 32;
    uint64_t r;
    do { r = next() & mask; } while (r >= n);
    return r;
  }
};

// [[Rcpp::export]]
NumericMatrix ps_lm_fixed_par(NumericVector beta_n_, NumericMatrix V_n_,
                              NumericMatrix X_future_, int B, double sigma,
                              double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = X_future_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X_future).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X_future).");
  if (B < 1) stop("B must be >= 1.");
  
  const arma::mat X(X_future_.begin(), N, p, false);   // read-only view
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  arma::mat       V_N    = as<arma::mat>(V_n_);
  
  // --- gain sequence: serial ---------------------------------------
  arma::mat aVx(N, p);
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    Vx = V_N * x_i.t();
    const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
    aVx.row(i) = a_i * Vx.t();
    V_N -= a_i * (Vx * Vx.t());
    const arma::mat V_t = V_N.t();
    V_N = 0.5 * (V_N + V_t);
  }
  
  // aVx transposed once: the inner loop then walks contiguous columns.
  const arma::mat aVxT = aVx.t();                      // p x N
  
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta_n;
      for (arma::uword i = 0; i < N; ++i)
        beta += (sigma * rng.norm()) * aVxT.col(i);
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_future_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix ps_lm_fixed_par_alt(NumericVector beta_n_, NumericMatrix V_n_,
                                  NumericMatrix X_future_, int B, double sigma,
                                  double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = X_future_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X_future).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X_future).");
  if (B < 1) stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  
  const arma::mat X(X_future_.begin(), N, p, false);   // read-only view
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  arma::mat       V_N    = as<arma::mat>(V_n_);
  
  // --- gain sequence and predictive sds: serial ------------------------
  //   sigma_i^2 = sigma^2 / a_i = sigma^2 * (1 + x_i' V_{i-1} x_i)
  // Both depend only on X_future and V_n (not on b, not on beta), so they
  // are computed once. Replicate b is then
  //   beta^(b) = beta_n + sum_i z_i^(b) * sigma_i * a_i * V_{i-1} x_i,
  // with z_i^(b) ~ N(0,1) iid, i.e. e_i = y_i - x_i' beta_{i-1} ~ N(0, sigma_i^2).
  arma::mat gT(p, N);              // column i: sigma_i * a_i * V_{i-1} x_i
  arma::vec sigma_i(N);            // predictive sd of y_i
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    Vx = V_N * x_i.t();
    const double q_i = 1.0 + arma::dot(x_i, Vx);       // = 1 / a_i
    if (!(q_i > 0.0))
      stop("1 + x_i' V x_i <= 0 at row %d: V_n is not positive definite.",
           (int) i + 1);
    const double a_i = 1.0 / q_i;
    
    sigma_i[i]  = sigma * std::sqrt(q_i);
    gT.col(i)   = (sigma_i[i] * a_i) * Vx;
    
    V_N -= a_i * (Vx * Vx.t());
    const arma::mat V_t = V_N.t();
    V_N = 0.5 * (V_N + V_t);
  }
  
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta_n;
      for (arma::uword i = 0; i < N; ++i)
        beta += rng.norm() * gT.col(i);                 // e_i = sigma_i * z_i
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_future_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  // Predictive sds sigma_i (length N), common to all replicates.
  out.attr("sigma_i") = NumericVector(sigma_i.begin(), sigma_i.end());
  
  return out;
}



// [[Rcpp::export]]
NumericMatrix ps_lm_par(NumericVector beta_n_, NumericMatrix V_n_,
                        NumericMatrix X_, int N, int B, double sigma,
                        double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N_obs = X_.nrow();
  const arma::uword p     = X_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X).");
  if (N_obs < 1) stop("X must have at least one row.");
  if (N < 0)     stop("N must be >= 0.");
  if (B < 1)     stop("B must be >= 1.");
  
  // Transposed once: Xt.col(k) is contiguous, X.row(k) is strided.
  const arma::mat Xt = arma::mat(X_.begin(), N_obs, p, false).t();   // p x N_obs
  
  const arma::vec beta_n  = as<arma::vec>(beta_n_);
  const arma::mat V_start = as<arma::mat>(V_n_);
  
  // Stored p x B (one column per replicate): column writes are
  // contiguous and disjoint, so threads never share a cache line.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const uint64_t nobs = (uint64_t) N_obs;
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p), Vx(p);
    arma::mat V_N(p, p);
    
    for (int b = lo; b < hi; ++b) {
      // Seed from (seed, b) only: independent of thread count and of
      // scheduling, so the output is fully reproducible.
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      
      beta = beta_n;
      V_N  = V_start;
      
      for (int i = 0; i < N; ++i) {
        const arma::uword k = (arma::uword) rng.bounded(nobs);
        const arma::subview_col<double> x_i = Xt.col(k);
        
        Vx = V_N * x_i;
        const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
        const double e_i = sigma * rng.norm();   // y_i - x_i'beta
        
        beta += (a_i * e_i) * Vx;
        V_N  -= a_i * (Vx * Vx.t());
      }
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dnm = X_.attr("dimnames");
  if (!Rf_isNull(dnm)) {
    List dnl(dnm);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix ps_lm_par_alt(NumericVector beta_n_, NumericMatrix V_n_,
                            NumericMatrix X_, int N, int B, double sigma,
                            double seed = 1.0, int nthreads = 0,
                            bool keep_sigma = true) {
  
  const arma::uword N_obs = X_.nrow();
  const arma::uword p     = X_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X).");
  if (N_obs < 1) stop("X must have at least one row.");
  if (N < 0)     stop("N must be >= 0.");
  if (B < 1)     stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  
  // Transposed once: Xt.col(k) is contiguous, X.row(k) is strided.
  const arma::mat Xt = arma::mat(X_.begin(), N_obs, p, false).t();   // p x N_obs
  
  const arma::vec beta_n  = as<arma::vec>(beta_n_);
  const arma::mat V_start = as<arma::mat>(V_n_);
  
  // Stored p x B (one column per replicate): column writes are
  // contiguous and disjoint, so threads never share a cache line.
  arma::mat beta_out(p, B);
  
  // Predictive sds, N x B (column per replicate, same reason as above);
  // transposed to B x N at the end. Empty if keep_sigma = false.
  arma::mat sig_out(keep_sigma ? (arma::uword) N : 0,
                    keep_sigma ? (arma::uword) B : 0);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const uint64_t nobs = (uint64_t) N_obs;
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  // R's API (stop() included) must not be called from worker threads:
  // each thread records a failure in its own slot, reported after join.
  std::vector<int> bad_b(nthr, 0), bad_i(nthr, 0);
  
  auto worker = [&](int t, int lo, int hi) {
    arma::vec beta(p), Vx(p);
    arma::mat V_N(p, p);
    
    for (int b = lo; b < hi; ++b) {
      // Seed from (seed, b) only: independent of thread count and of
      // scheduling, so the output is fully reproducible.
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      
      beta = beta_n;
      V_N  = V_start;
      
      for (int i = 0; i < N; ++i) {
        const arma::uword k = (arma::uword) rng.bounded(nobs);
        const arma::subview_col<double> x_i = Xt.col(k);
        
        Vx = V_N * x_i;
        const double q_i = 1.0 + arma::dot(x_i, Vx);   // = 1 / a_i
        if (!(q_i > 0.0)) { bad_b[t] = b + 1; bad_i[t] = i + 1; return; }
        
        const double a_i = 1.0 / q_i;
        const double s_i = sigma * std::sqrt(q_i);      // sigma_i = sigma / sqrt(a_i)
        const double e_i = s_i * rng.norm();            // y_i - x_i'beta_{i-1} ~ N(0, sigma_i^2)
        
        if (keep_sigma) sig_out((arma::uword) i, (arma::uword) b) = s_i;
        
        beta += (a_i * e_i) * Vx;
        V_N  -= a_i * (Vx * Vx.t());
      }
      beta_out.col(b) = beta;
    }
  };
  
  if (nthr == 1) {
    worker(0, 0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, t, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  for (int t = 0; t < nthr; ++t)
    if (bad_b[t] > 0)
      stop("1 + x_i' V x_i <= 0 (replicate %d, step %d): V_n is not positive definite.",
           bad_b[t], bad_i[t]);
    
    NumericMatrix out = wrap(arma::mat(beta_out.t()));
    
    SEXP dnm = X_.attr("dimnames");
    if (!Rf_isNull(dnm)) {
      List dnl(dnm);
      if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
        out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
    }
    
    // sigma_i: B x N, row b = path of predictive sds of replicate b.
    if (keep_sigma)
      out.attr("sigma_i") = wrap(arma::mat(sig_out.t()));
    
    return out;
}





// ---------------------------------------------------------------------------
// Student-t analogue of ps_lm_fixed_par().
//
// The gain sequence a_i * V_i * x_i depends only on X_future, so it is still
// precomputed once and the b-th replicate is the linear form
//
//     beta^(b) = beta_n + sum_i (aVx_i) * e_i^(b),
//
// where e_i^(b) is the i-th predictive residual.  Only the law of e_i changes:
//
//     normal version : e_i ~ N(0, sigma^2)
//     this version   : e_i ~ sigma * t_df
//
// i.e. y_i | . ~ location-scale t with location x_i' beta_{i-1}, scale sigma,
// and df degrees of freedom.  Residuals are independent across i.
//
// NOTE on the parameterisation: sigma is the SCALE, not the sd.  For df > 2,
// Var(e_i) = sigma^2 * df/(df-2); for df <= 2 the variance is undefined.  If
// you want the innovations to have sd sigma, pass sigma * sqrt((df-2)/df).
// df = R_PosInf reproduces ps_lm_fixed_par() in distribution.
//
// Sampling: t = Z / sqrt(W/df) with W ~ chi2_df, W = 2G, G ~ Gamma(df/2, 1),
// so t = Z * sqrt(df / (2G)).  G is drawn with Marsaglia-Tsang (squeeze first,
// ~1.0-1.1 iterations per draw), extended to shape < 1 (i.e. df < 2, Cauchy
// included) by the boost identity G_a = G_{a+1} * U^(1/a).
// ---------------------------------------------------------------------------


// The only place where the RNG's uniform generator is named.  Adjust this one
// line if your RNG class calls it something other than unif() (unif01(), u01(),
// runif(), ...).  It must return U(0,1).
static inline double u01(RNG &rng) { return rng.unif(); }

// Sampler for a standardised t_df, constructed once per call: the Marsaglia-
// Tsang constants depend on df only.  Immutable, hence shareable across
// threads; all state lives in the RNG passed in.
struct StudentT {
  double df, half_df, shape, d, c, inv_shape;
  bool boost, gaussian;
  
  explicit StudentT(double df_)
    : df(df_), half_df(0.5 * df_), shape(0.5 * df_),
      boost(shape < 1.0), gaussian(!std::isfinite(df_)) {
    const double a = boost ? shape + 1.0 : shape;   // Marsaglia-Tsang needs a>=1
    d = a - 1.0 / 3.0;
    c = 1.0 / std::sqrt(9.0 * d);
    inv_shape = 1.0 / shape;
  }
  
  // G ~ Gamma(shape = df/2, scale = 1)
  double gamma1(RNG &rng) const {
    double g;
    for (;;) {
      double x, v;
      do {
        x = rng.norm();
        v = 1.0 + c * x;
      } while (v <= 0.0);
      v = v * v * v;
      const double u  = u01(rng);
      const double x2 = x * x;
      if (u < 1.0 - 0.0331 * x2 * x2) { g = d * v; break; }              // squeeze
      if (std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v))) { g = d * v; break; }
    }
    if (boost) {
      double u;
      do { u = u01(rng); } while (u <= 0.0);
      g *= std::pow(u, inv_shape);
    }
    return g > 1e-300 ? g : 1e-300;   // guard against underflow for tiny df
  }
  
  // standardised t_df draw
  double draw(RNG &rng) const {
    if (gaussian) return rng.norm();
    return rng.norm() * std::sqrt(half_df / gamma1(rng));
  }
};


// [[Rcpp::export]]
NumericMatrix ps_lm_fixed_par_t(NumericVector beta_n_, NumericMatrix V_n_,
                                NumericMatrix X_future_, int B, double sigma,
                                double df, double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = X_future_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X_future).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X_future).");
  if (B < 1) stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0)) stop("df must be > 0.");
  
  const arma::mat X(X_future_.begin(), N, p, false);   // read-only view
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  arma::mat       V_N    = as<arma::mat>(V_n_);
  
  // --- gain sequence: serial, unchanged (does not depend on the error law) --
  arma::mat aVx(N, p);
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    Vx = V_N * x_i.t();
    const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
    aVx.row(i) = a_i * Vx.t();
    V_N -= a_i * (Vx * Vx.t());
    const arma::mat V_t = V_N.t();
    V_N = 0.5 * (V_N + V_t);
  }
  
  // aVx transposed once: the inner loop then walks contiguous columns.
  const arma::mat aVxT = aVx.t();                      // p x N
  
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const StudentT tsamp(df);                            // shared, read-only
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta_n;
      for (arma::uword i = 0; i < N; ++i)
        beta += (sigma * tsamp.draw(rng)) * aVxT.col(i);
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_future_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix ps_lm_fixed_par_t_alt(NumericVector beta_n_, NumericMatrix V_n_,
                                    NumericMatrix X_future_, int B, double sigma,
                                    double df, double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N = X_future_.nrow();
  const arma::uword p = X_future_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X_future).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X_future).");
  if (B < 1) stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0)) stop("df must be > 0.");
  
  const arma::mat X(X_future_.begin(), N, p, false);   // read-only view
  const arma::vec beta_n = as<arma::vec>(beta_n_);
  arma::mat       V_N    = as<arma::mat>(V_n_);
  
  // --- gain sequence and predictive scales: serial ----------------------
  //   sigma_i^2 = sigma^2 / a_i = sigma^2 * (1 + x_i' V_{i-1} x_i)
  // Both depend only on X_future and V_n (not on b, not on beta, not on df),
  // so they are computed once. Replicate b is then
  //   beta^(b) = beta_n + sum_i t_i^(b) * sigma_i * a_i * V_{i-1} x_i,
  // with t_i^(b) ~ t_df iid, i.e. e_i = y_i - x_i' beta_{i-1} ~ sigma_i * t_df.
  arma::mat gT(p, N);              // column i: sigma_i * a_i * V_{i-1} x_i
  arma::vec sigma_i(N);            // predictive scale of y_i
  arma::vec Vx(p);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(i);
    Vx = V_N * x_i.t();
    const double q_i = 1.0 + arma::dot(x_i, Vx);       // = 1 / a_i
    if (!(q_i > 0.0))
      stop("1 + x_i' V x_i <= 0 at row %d: V_n is not positive definite.",
           (int) i + 1);
    const double a_i = 1.0 / q_i;
    
    sigma_i[i] = sigma * std::sqrt(q_i);
    gT.col(i)  = (sigma_i[i] * a_i) * Vx;
    
    V_N -= a_i * (Vx * Vx.t());
    const arma::mat V_t = V_N.t();
    V_N = 0.5 * (V_N + V_t);
  }
  
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const StudentT tsamp(df);                            // shared, read-only
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta_n;
      for (arma::uword i = 0; i < N; ++i)
        beta += tsamp.draw(rng) * gT.col(i);            // e_i = sigma_i * t_i
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_future_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  // Predictive scales sigma_i (length N), common to all replicates.
  out.attr("sigma_i") = NumericVector(sigma_i.begin(), sigma_i.end());
  
  return out;
}


// Valore restituito quando i parametri non sono ammissibili o la
// fattorizzazione fallisce: grande ma finito, cosi' nlminb accorcia il passo
// invece di fermarsi.
static const double BIG = 1e100;

// eps della funzione originale: rls_lambda() usa (lambda_j*Omega_j + eps*I)
// come precisione del blocco j e 1e10 come varianza dell'intercetta.
static const double EPS_ORIG = 1e-10;

// ---------------------------------------------------------------------------
// log-densita' t di Student standard.
// df >= DF_GAUSS viene trattato come limite normale: evita lgamma di argomenti
// enormi (che perde precisione) e rende la funzione continua fino a df = Inf,
// utile perche' senza penalizzazione l'ottimizzatore puo' spingersi molto in
// alto sul valore di df.
// ---------------------------------------------------------------------------
static const double DF_GAUSS = 1e8;


static const double SN_LOG2           = 0.69314718055994530942;  // log(2)
static const double SN_SQRT_2_OVER_PI = 0.79788456080286535588;  // sqrt(2/pi)

// log-densita' della skew normal standard SN(0, 1, alpha) in z.
// R::pnorm con log.p = TRUE resta accurato anche nella coda sinistra.
static inline double dsn_log_std(const double z, const double alpha) {
  return SN_LOG2 + R::dnorm(z, 0.0, 1.0, 1) + R::pnorm(alpha * z, 0.0, 1.0, 1, 1);
}

// Nucleo per il modello lineare: una Cholesky di V_n^{-1}, poi beta_n,
// leverage e residui LOO. Nulla dipende da (sigma, alpha): un errore qui e'
// un problema dei dati o di V_0, quindi stop() e non BIG.
static void loo_core_lm(const vec& y, const mat& X, const vec& beta_0,
                        const mat& V_0, vec& beta, vec& h, vec& r_loo) {
  
  mat S;                                        // V_0^{-1}
  const mat V0s = 0.5 * (V_0 + V_0.t());
  if (!inv_sympd(S, V0s))
    Rcpp::stop("V_0 deve essere simmetrica definita positiva");
  
  mat A = X.t() * X + S;                        // V_n^{-1}
  A = 0.5 * (A + A.t());
  
  mat R;                                        // A = R'R
  if (!chol(R, A))
    Rcpp::stop("V_0^{-1} + X'X non e' definita positiva");
  
  // beta_n = A^{-1} (X'y + V_0^{-1} beta_0)
  const vec b = X.t() * y + S * beta_0;
  beta = solve(trimatu(R), solve(trimatl(R.t()), b, solve_opts::fast),
               solve_opts::fast);
  
  // W = R^{-T} X'  =>  h_i = ||W.col(i)||^2 = x_i' V_n x_i
  const mat W = solve(trimatl(R.t()), X.t(), solve_opts::fast);
  h = trans(sum(square(W), 0));
  
  if (!h.is_finite() || h.max() >= 1.0 - 1e-12)
    Rcpp::stop("qualche leverage h_i >= 1: residui LOO non definiti");
  
  r_loo = (y - X * beta) / (1.0 - h);           // y_i - x_i' beta_{-i}
}








// ---------------------------------------------------------------------------
// Student-t analogue of ps_lm_par().
//
// Same predictive recursion, same row resampling; only the law of the one-step
// predictive residual changes:
//
//     normal version : e_i ~ N(0, sigma^2)
//     this version   : e_i ~ sigma * t_df
//
// i.e. y_i | . ~ location-scale t with location x_i' beta_{i-1}, scale sigma,
// df degrees of freedom; residuals independent across i.
//
// sigma is the SCALE, not the sd: Var(e_i) = sigma^2 * df/(df-2) for df > 2,
// undefined for df <= 2.  df = R_PosInf reproduces ps_lm_par() in distribution.
//
// The StudentT helper below is identical to the one used by
// ps_lm_fixed_par_t() -- if both functions live in the same translation unit,
// keep only one copy of the block between the two ==== markers.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix ps_lm_par_t(NumericVector beta_n_, NumericMatrix V_n_,
                          NumericMatrix X_, int N, int B, double sigma,
                          double df, double seed = 1.0, int nthreads = 0) {
  
  const arma::uword N_obs = X_.nrow();
  const arma::uword p     = X_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X).");
  if (N_obs < 1) stop("X must have at least one row.");
  if (N < 0)     stop("N must be >= 0.");
  if (B < 1)     stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0))     stop("df must be > 0.");
  
  // Transposed once: Xt.col(k) is contiguous, X.row(k) is strided.
  const arma::mat Xt = arma::mat(X_.begin(), N_obs, p, false).t();   // p x N_obs
  
  const arma::vec beta_n  = as<arma::vec>(beta_n_);
  const arma::mat V_start = as<arma::mat>(V_n_);
  
  // Stored p x B (one column per replicate): column writes are
  // contiguous and disjoint, so threads never share a cache line.
  arma::mat beta_out(p, B);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const uint64_t nobs = (uint64_t) N_obs;
  
  const StudentT tsamp(df);            // df-dependent constants, built once
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p), Vx(p);
    arma::mat V_N(p, p);
    
    for (int b = lo; b < hi; ++b) {
      // Seed from (seed, b) only: independent of thread count and of
      // scheduling, so the output is fully reproducible.
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      
      beta = beta_n;
      V_N  = V_start;
      
      for (int i = 0; i < N; ++i) {
        const arma::uword k = (arma::uword) rng.bounded(nobs);
        const arma::subview_col<double> x_i = Xt.col(k);
        
        Vx = V_N * x_i;
        const double a_i = 1.0 / (1.0 + arma::dot(x_i, Vx));
        const double e_i = sigma * tsamp.draw(rng);   // y_i - x_i'beta, ~ sigma * t_df
        
        beta += (a_i * e_i) * Vx;
        V_N  -= a_i * (Vx * Vx.t());
      }
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dnm = X_.attr("dimnames");
  if (!Rf_isNull(dnm)) {
    List dnl(dnm);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  return out;
}


// [[Rcpp::export]]
NumericMatrix ps_lm_par_t_alt(NumericVector beta_n_, NumericMatrix V_n_,
                              NumericMatrix X_, int N, int B, double sigma,
                              double df, double seed = 1.0, int nthreads = 0,
                              bool keep_sigma = true) {
  
  const arma::uword N_obs = X_.nrow();
  const arma::uword p     = X_.ncol();
  
  if ((arma::uword) beta_n_.size() != p)
    stop("length(beta_n) must equal ncol(X).");
  if ((arma::uword) V_n_.nrow() != p || (arma::uword) V_n_.ncol() != p)
    stop("V_n must be p x p with p = ncol(X).");
  if (N_obs < 1) stop("X must have at least one row.");
  if (N < 0)     stop("N must be >= 0.");
  if (B < 1)     stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0))     stop("df must be > 0.");
  
  // Transposed once: Xt.col(k) is contiguous, X.row(k) is strided.
  const arma::mat Xt = arma::mat(X_.begin(), N_obs, p, false).t();   // p x N_obs
  
  const arma::vec beta_n  = as<arma::vec>(beta_n_);
  const arma::mat V_start = as<arma::mat>(V_n_);
  
  // Stored p x B (one column per replicate): column writes are
  // contiguous and disjoint, so threads never share a cache line.
  arma::mat beta_out(p, B);
  
  // Predictive scales, N x B (column per replicate, same reason as above);
  // transposed to B x N at the end. Empty if keep_sigma = false.
  arma::mat sig_out(keep_sigma ? (arma::uword) N : 0,
                    keep_sigma ? (arma::uword) B : 0);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  const uint64_t nobs = (uint64_t) N_obs;
  
  const StudentT tsamp(df);            // df-dependent constants, built once
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  // R's API (stop() included) must not be called from worker threads:
  // each thread records a failure in its own slot, reported after join.
  std::vector<int> bad_b(nthr, 0), bad_i(nthr, 0);
  
  auto worker = [&](int t, int lo, int hi) {
    arma::vec beta(p), Vx(p);
    arma::mat V_N(p, p);
    
    for (int b = lo; b < hi; ++b) {
      // Seed from (seed, b) only: independent of thread count and of
      // scheduling, so the output is fully reproducible.
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      
      beta = beta_n;
      V_N  = V_start;
      
      for (int i = 0; i < N; ++i) {
        const arma::uword k = (arma::uword) rng.bounded(nobs);
        const arma::subview_col<double> x_i = Xt.col(k);
        
        Vx = V_N * x_i;
        const double q_i = 1.0 + arma::dot(x_i, Vx);   // = 1 / a_i
        if (!(q_i > 0.0)) { bad_b[t] = b + 1; bad_i[t] = i + 1; return; }
        
        const double a_i = 1.0 / q_i;
        const double s_i = sigma * std::sqrt(q_i);      // sigma_i = sigma / sqrt(a_i)
        const double e_i = s_i * tsamp.draw(rng);       // y_i - x_i'beta_{i-1} ~ sigma_i * t_df
        
        if (keep_sigma) sig_out((arma::uword) i, (arma::uword) b) = s_i;
        
        beta += (a_i * e_i) * Vx;
        V_N  -= a_i * (Vx * Vx.t());
      }
      beta_out.col(b) = beta;
    }
  };
  
  if (nthr == 1) {
    worker(0, 0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, t, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  for (int t = 0; t < nthr; ++t)
    if (bad_b[t] > 0)
      stop("1 + x_i' V x_i <= 0 (replicate %d, step %d): V_n is not positive definite.",
           bad_b[t], bad_i[t]);
    
    NumericMatrix out = wrap(arma::mat(beta_out.t()));
    
    SEXP dnm = X_.attr("dimnames");
    if (!Rf_isNull(dnm)) {
      List dnl(dnm);
      if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
        out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
    }
    
    // sigma_i: B x N, row b = path of predictive scales of replicate b.
    if (keep_sigma)
      out.attr("sigma_i") = wrap(arma::mat(sig_out.t()));
    
    return out;
}





// ---------------------------------------------------------------------------
// ps_rs_parallel_t_fixed
//
// Fully C++ version of ps_rs_parallel_t_1 / ps_rs_single_iteration_t_1,
// written in the style of ps_lm_fixed_par_t.
//
// In the recursive update y_new - x_new' beta = sigma * t_new, so every step
// reduces to   beta += sigma * t_new * (a * P x_new):   the gain a * P x
// depends only on the sampled rows and on the P recursion, not on b nor on
// beta.  Hence:
//   1. the N rows of X are sampled ONCE (serial, own RNG stream);
//   2. the gain sequence a_i * P_{i-1} x_i is computed ONCE (serial);
//   3. the loop over b runs in parallel (std::thread) and only draws the t's.
//
// Put it in the same .cpp file as ps_lm_fixed_par_t, after RNG and StudentT
// (same headers: RcppArmadillo, <thread>, <vector>, <cstdint>, <cmath>,
// <algorithm>).
//
// Returns a B x p matrix (one row per replicate, colnames taken from X) with
// attribute "idx" = 1-based indices of the N sampled rows of X.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix ps_rs_parallel_t_fixed(NumericMatrix X_, NumericVector beta_,
                                     NumericMatrix P_, int N, int B,
                                     double df, double sigma,
                                     double seed = 1.0, int nthreads = 0) {
  
  const arma::uword n = X_.nrow();
  const arma::uword p = X_.ncol();
  
  if (n < 1) stop("X must have at least one row.");
  if ((arma::uword) beta_.size() != p)
    stop("length(beta) must equal ncol(X).");
  if ((arma::uword) P_.nrow() != p || (arma::uword) P_.ncol() != p)
    stop("P must be p x p with p = ncol(X).");
  if (N < 0) stop("N must be >= 0.");
  if (B < 1) stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0)) stop("df must be > 0.");
  
  const arma::mat X(X_.begin(), n, p, false);          // read-only view
  const arma::vec beta0 = as<arma::vec>(beta_);
  arma::mat       P     = as<arma::mat>(P_);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  
  // --- 1. row sampling: once, serial ----------------------------------------
  // Own stream: seed 'base' corresponds to "b = -1" in the replicate scheme
  // below (multiplier 0), so it never coincides with the streams of the t's.
  arma::uvec idx((arma::uword) N);
  {
    RNG rng(base);
    for (int i = 0; i < N; ++i)
      idx[i] = (arma::uword) rng.bounded((uint64_t) n);
  }
  
  // --- 2. gain sequence a_i * P_{i-1} x_i: once, serial ---------------------
  // Written directly as p x N (one column per step): the inner loop walks
  // contiguous columns, and no N x p temporary + transpose is needed.
  arma::mat aPxT(p, (arma::uword) N);
  arma::vec Px(p);
  
  for (int i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(idx[i]);
    Px = P * x_i.t();
    const double a_i = 1.0 / (1.0 + arma::dot(x_i, Px));
    aPxT.col(i) = a_i * Px;
    P -= a_i * (Px * Px.t());
    const arma::mat P_t = P.t();
    P = 0.5 * (P + P_t);
  }
  
  // --- 3. replicates: parallel over b, only the t draws depend on b ---------
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  const StudentT tsamp(df);                            // shared, read-only
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta0;
      for (int i = 0; i < N; ++i)
        beta += (sigma * tsamp.draw(rng)) * aPxT.col(i);
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  IntegerVector idx_r(N);
  for (int i = 0; i < N; ++i) idx_r[i] = (int) idx[i] + 1;
  out.attr("idx") = idx_r;
  
  return out;
}










static inline double dt_log_2(const double z, const double df) {
  if (df >= DF_GAUSS)
    return -0.5 * z * z - 0.5 * std::log(2.0 * M_PI);
  return R::lgammafn(0.5 * (df + 1.0))
    - R::lgammafn(0.5 * df)
    - 0.5 * std::log(df * M_PI)
    - 0.5 * (df + 1.0) * std::log1p(z * z / df);
}

// ---------------------------------------------------------------------------
// S = P_0^{-1} = blockdiag( eps, lambda_1*Omega_1 + eps*I, ... )
//
// Nota: e' esattamente l'inversa della P_0 costruita da rls_lambda(), ma
// ottenuta senza invertire nulla. L'intercetta (colonna 0) ha varianza a
// priori 1e10, cioe' precisione 1e-10 = eps.
// ---------------------------------------------------------------------------
static mat build_S(const vec& lambda, const std::vector<mat>& Om, const uword p) {
  mat S(p, p, fill::zeros);
  S(0, 0) = EPS_ORIG;
  uword off = 1;
  for (uword j = 0; j < Om.size(); ++j) {
    const uword k = Om[j].n_rows;
    if (off + k > p)
      Rcpp::stop("le dimensioni dei blocchi di Omega non tornano con ncol(X)");
    mat B = lambda(j) * Om[j];
    B.diag() += EPS_ORIG;
    S.submat(off, off, off + k - 1, off + k - 1) = B;
    off += k;
  }
  if (off != p)
    Rcpp::stop("la somma dei blocchi di Omega deve dare ncol(X) - 1");
  return S;
}

static std::vector<mat> as_blocks(const Rcpp::List& Omega) {
  std::vector<mat> Om;
  Om.reserve(Omega.size());
  for (R_xlen_t j = 0; j < Omega.size(); ++j) {
    mat Oj = Rcpp::as<mat>(Omega[j]);
    if (Oj.n_rows != Oj.n_cols)
      Rcpp::stop("Omega[[%d]] non e' quadrata", (int) j + 1);
    Om.push_back(std::move(Oj));
  }
  return Om;
}

// ---------------------------------------------------------------------------
// Nucleo: una Cholesky, poi beta, leverage e residui LOO.
// Restituisce false se A non e' fattorizzabile o se qualche h_i sfiora 1.
// ---------------------------------------------------------------------------
static bool loo_core(const vec& y, const mat& X, const vec& beta_0,
                     const vec& lambda, const std::vector<mat>& Om,
                     vec& beta, vec& h, vec& r_loo) {
  
  const uword p = X.n_cols;
  const mat   S = build_S(lambda, Om, p);
  
  mat A = X.t() * X + S;
  A = 0.5 * (A + A.t());              // simmetria esatta prima di chol()
  
  mat R;                              // A = R'R, R triangolare superiore
  if (!chol(R, A)) return false;
  
  // beta = A^{-1} (X'y + S beta_0), risolto con due sostituzioni triangolari
  const vec b = X.t() * y + S * beta_0;
  beta = solve(trimatu(R), solve(trimatl(R.t()), b), solve_opts::fast);
  
  // W = R^{-T} X'   (p x n)   =>   h_i = ||W.col(i)||^2 = x_i' A^{-1} x_i
  const mat W = solve(trimatl(R.t()), X.t(), solve_opts::fast);
  h = trans(sum(square(W), 0));
  
  if (!h.is_finite() || h.max() >= 1.0 - 1e-12) return false;
  
  const vec e = y - X * beta;         // residui in-sample
  r_loo = e / (1.0 - h);              // residui leave-one-out, esatti
  return true;
}

// ---------------------------------------------------------------------------
//  Funzione obiettivo. Firma identica alla versione R.
//
//    log_par : c(log(df), log(sigma), log(lambda_1), ..., log(lambda_q))
//    y       : vettore risposta, lunghezza n
//    X       : matrice n x p, prima colonna = intercetta
//    beta_0  : media a priori, lunghezza p
//    Omega   : lista di q matrici di penalita' (i blocchi dopo l'intercetta)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double npll_loo(const arma::vec& log_par,
                const arma::vec& y,
                const arma::mat& X,
                const arma::vec& beta_0,
                const Rcpp::List& Omega) {
  
  const uword n = y.n_elem;
  const uword p = X.n_cols;
  const uword q = (uword) Omega.size();
  
  if (X.n_rows != n)       Rcpp::stop("nrow(X) deve essere uguale a length(y)");
  if (beta_0.n_elem != p)  Rcpp::stop("length(beta_0) deve essere uguale a ncol(X)");
  if (log_par.n_elem != q + 2)
    Rcpp::stop("length(log_par) deve essere length(Omega) + 2");
  if (!log_par.is_finite()) return BIG;
  
  const double df    = std::exp(log_par(0));
  const double sigma = std::exp(log_par(1));
  const vec    lambda = exp(log_par.subvec(2, log_par.n_elem - 1));
  
  if (!(df > 0.0) || !(sigma > 0.0) || !std::isfinite(sigma)
        || !lambda.is_finite() || lambda.min() <= 0.0) return BIG;
  
  const std::vector<mat> Om = as_blocks(Omega);
  
  vec beta, h, r;
  if (!loo_core(y, X, beta_0, lambda, Om, beta, h, r)) return BIG;
  
  // -sum_i [ log dt(r_i/sigma ; df) - log sigma ]
  // identico a  -sum(extraDistr::dlst(y, mu = mu_loo, df, sigma, log = TRUE))
  const double log_sigma = std::log(sigma);
  double out = 0.0;
  for (uword i = 0; i < n; ++i)
    out += dt_log_2(r(i) / sigma, df) - log_sigma;
  
  if (!std::isfinite(out)) return BIG;
  return -out;
}

// ---------------------------------------------------------------------------
//  Diagnostica: gli stessi ingredienti, restituiti invece che sommati.
//  Serve per controllare mu_loo contro il ciclo RLS e per guardare i residui.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List loo_fit(const arma::vec& lambda,
                   const arma::vec& y,
                   const arma::mat& X,
                   const arma::vec& beta_0,
                   const Rcpp::List& Omega) {
  
  const std::vector<mat> Om = as_blocks(Omega);
  if (lambda.n_elem != Om.size())
    Rcpp::stop("length(lambda) deve essere uguale a length(Omega)");
  
  vec beta, h, r;
  if (!loo_core(y, X, beta_0, lambda, Om, beta, h, r))
    Rcpp::stop("fattorizzazione fallita, oppure qualche leverage >= 1");
  
  return Rcpp::List::create(
    Rcpp::Named("beta")   = beta,        // coefficienti sul campione completo
    Rcpp::Named("h")      = h,           // leverage
    Rcpp::Named("edf")    = accu(h),     // gradi di liberta' effettivi
    Rcpp::Named("mu_loo") = y - r,       // previsioni leave-one-out
    Rcpp::Named("r_loo")  = r);          // residui leave-one-out
}




// ---------------------------------------------------------------------------
//  Variante di npll_loo con scala predittiva specifica per osservazione.
//
//  Nel ciclo RLS la previsione di y_i ha scala
//      sigma_i^2 = sigma^2 * (1 + x_i' V_{i-1} x_i).
//  In leave-one-out V_{i-1} diventa V_{-i} = (A - x_i x_i')^{-1}, con
//  A = X'X + S, e per Sherman-Morrison
//      x_i' V_{-i} x_i = h_i / (1 - h_i),   h_i = x_i' A^{-1} x_i,
//  quindi
//      sigma_i^2 = sigma^2 / (1 - h_i).
//  Il leverage h_i e' gia' calcolato da loo_core(); loo_core() garantisce
//  anche h_i < 1, per cui sigma_i e' sempre finito.
//
//    log_par : c(log(df), log(sigma), log(lambda_1), ..., log(lambda_q))
//    y       : vettore risposta, lunghezza n
//    X       : matrice n x p, prima colonna = intercetta
//    beta_0  : media a priori, lunghezza p
//    Omega   : lista di q matrici di penalita' (i blocchi dopo l'intercetta)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double npll_loo_alt(const arma::vec& log_par,
                    const arma::vec& y,
                    const arma::mat& X,
                    const arma::vec& beta_0,
                    const Rcpp::List& Omega) {
  
  const uword n = y.n_elem;
  const uword p = X.n_cols;
  const uword q = (uword) Omega.size();
  
  if (X.n_rows != n)       Rcpp::stop("nrow(X) deve essere uguale a length(y)");
  if (beta_0.n_elem != p)  Rcpp::stop("length(beta_0) deve essere uguale a ncol(X)");
  if (log_par.n_elem != q + 2)
    Rcpp::stop("length(log_par) deve essere length(Omega) + 2");
  if (!log_par.is_finite()) return BIG;
  
  const double df    = std::exp(log_par(0));
  const double sigma = std::exp(log_par(1));
  const vec    lambda = exp(log_par.subvec(2, log_par.n_elem - 1));
  
  if (!(df > 0.0) || !(sigma > 0.0) || !std::isfinite(sigma)
        || !lambda.is_finite() || lambda.min() <= 0.0) return BIG;
  
  const std::vector<mat> Om = as_blocks(Omega);
  
  vec beta, h, r;
  if (!loo_core(y, X, beta_0, lambda, Om, beta, h, r)) return BIG;
  
  // -sum_i [ log dt(r_i/sigma_i ; df) - log sigma_i ],
  //   sigma_i     = sigma / sqrt(1 - h_i)
  //   log sigma_i = log sigma - 0.5 * log(1 - h_i)
  const double log_sigma = std::log(sigma);
  double out = 0.0;
  for (uword i = 0; i < n; ++i) {
    const double log1mh      = std::log1p(-h(i));            // log(1 - h_i)
    const double sigma_i     = sigma / std::sqrt(1.0 - h(i));
    const double log_sigma_i = log_sigma - 0.5 * log1mh;
    out += dt_log_2(r(i) / sigma_i, df) - log_sigma_i;
  }
  
  if (!std::isfinite(out)) return BIG;
  return -out;
}


// [[Rcpp::export]]
NumericMatrix ps_rs_parallel_t_fixed_alt(NumericMatrix X_, NumericVector beta_,
                                         NumericMatrix P_, int N, int B,
                                         double df, double sigma,
                                         double seed = 1.0, int nthreads = 0) {
  
  const arma::uword n = X_.nrow();
  const arma::uword p = X_.ncol();
  
  if (n < 1) stop("X must have at least one row.");
  if ((arma::uword) beta_.size() != p)
    stop("length(beta) must equal ncol(X).");
  if ((arma::uword) P_.nrow() != p || (arma::uword) P_.ncol() != p)
    stop("P must be p x p with p = ncol(X).");
  if (N < 0) stop("N must be >= 0.");
  if (B < 1) stop("B must be >= 1.");
  if (!(sigma >= 0.0)) stop("sigma must be >= 0.");
  if (!(df > 0.0)) stop("df must be > 0.");
  
  const arma::mat X(X_.begin(), n, p, false);          // read-only view
  const arma::vec beta0 = as<arma::vec>(beta_);
  arma::mat       P     = as<arma::mat>(P_);
  
  const uint64_t base = (uint64_t) (int64_t) std::llround(seed);
  
  // --- 1. row sampling: once, serial ----------------------------------------
  // Own stream: seed 'base' corresponds to "b = -1" in the replicate scheme
  // below (multiplier 0), so it never coincides with the streams of the t's.
  arma::uvec idx((arma::uword) N);
  {
    RNG rng(base);
    for (int i = 0; i < N; ++i)
      idx[i] = (arma::uword) rng.bounded((uint64_t) n);
  }
  
  // --- 2. gain sequence and predictive scales: once, serial -----------------
  //   sigma_i^2 = sigma^2 / a_i = sigma^2 * (1 + x_i' P_{i-1} x_i)
  // Both depend only on the sampled rows and on the P recursion (not on b,
  // not on beta, not on df). Replicate b is then
  //   beta^(b) = beta_0 + sum_i t_i^(b) * sigma_i * a_i * P_{i-1} x_i,
  // with t_i^(b) ~ t_df iid, i.e. y_i - x_i' beta_{i-1} ~ sigma_i * t_df.
  // Written directly as p x N (one column per step): the inner loop walks
  // contiguous columns.
  arma::mat gT(p, (arma::uword) N);    // column i: sigma_i * a_i * P_{i-1} x_i
  arma::vec sigma_i((arma::uword) N);  // predictive scale of y_i
  arma::vec Px(p);
  
  for (int i = 0; i < N; ++i) {
    const arma::rowvec x_i = X.row(idx[i]);
    Px = P * x_i.t();
    const double q_i = 1.0 + arma::dot(x_i, Px);       // = 1 / a_i
    if (!(q_i > 0.0))
      stop("1 + x_i' P x_i <= 0 at step %d: P is not positive definite.", i + 1);
    const double a_i = 1.0 / q_i;
    
    sigma_i[i] = sigma * std::sqrt(q_i);
    gT.col(i)  = (sigma_i[i] * a_i) * Px;
    
    P -= a_i * (Px * Px.t());
    const arma::mat P_t = P.t();
    P = 0.5 * (P + P_t);
  }
  
  // --- 3. replicates: parallel over b, only the t draws depend on b ---------
  // Stored p x B (column per replicate) to avoid false sharing between
  // threads; transposed back at the end.
  arma::mat beta_out(p, B);
  const StudentT tsamp(df);                            // shared, read-only
  
  auto worker = [&](int lo, int hi) {
    arma::vec beta(p);
    for (int b = lo; b < hi; ++b) {
      RNG rng(base + 0x9E3779B97F4A7C15ULL * (uint64_t)(b + 1));
      beta = beta0;
      for (int i = 0; i < N; ++i)
        beta += tsamp.draw(rng) * gT.col(i);            // e_i = sigma_i * t_i
      beta_out.col(b) = beta;
    }
  };
  
  int nthr = nthreads > 0 ? nthreads : (int) std::thread::hardware_concurrency();
  if (nthr < 1) nthr = 1;
  nthr = std::min(nthr, B);
  
  if (nthr == 1) {
    worker(0, B);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthr);
    const int chunk = (B + nthr - 1) / nthr;
    for (int t = 0; t < nthr; ++t) {
      const int lo = t * chunk, hi = std::min(B, lo + chunk);
      if (lo < hi) pool.emplace_back(worker, lo, hi);
    }
    for (std::thread &th : pool) th.join();
  }
  
  NumericMatrix out = wrap(arma::mat(beta_out.t()));
  
  SEXP dn = X_.attr("dimnames");
  if (!Rf_isNull(dn)) {
    List dnl(dn);
    if (dnl.size() >= 2 && !Rf_isNull(dnl[1]))
      out.attr("dimnames") = List::create(R_NilValue, dnl[1]);
  }
  
  IntegerVector idx_r(N);
  for (int i = 0; i < N; ++i) idx_r[i] = (int) idx[i] + 1;
  out.attr("idx") = idx_r;
  
  // Predictive scales sigma_i (length N, aligned with "idx"), common to all
  // replicates.
  out.attr("sigma_i") = NumericVector(sigma_i.begin(), sigma_i.end());
  
  return out;
}




// ---------------------------------------------------------------------------
//  loo-pll (Section 4.3) for the LINEAR MODEL with P_n a location-scale
//  Student-t.
//
//  Same as npll_loo_lm_skew_alt, only the predictive density changes:
//    - no Omega / lambda: V_0 is fixed and passed in;
//    - P_n is a location-scale t centred at x_i' beta_{-i};
//    - the only parameters are (df, sigma).
//
//  RLS / generalized ridge estimate, eq. (5) of the paper:
//      V_n^{-1} = V_0^{-1} + X'X,     beta_n = V_n (V_0^{-1} beta_0 + X'y).
//
//  Leave-one-out predictions (Sec. 4.3):
//      x_i' beta_{-i} = y_i - (y_i - x_i' beta_n) / (1 - h_i),
//      h_i = x_i' V_n x_i.
//
//  "alt" predictive scale (sigma_n^2 = sigma^2 / a_n in LOO form):
//      1 + x_i' V_{-i} x_i = 1 / (1 - h_i)   =>   sigma_i = sigma / sqrt(1 - h_i).
//
//  P_n:  y_i ~ x_i' beta_{-i} + sigma_i * t_df,
//  i.e. exactly the law simulated by ps_lm_fixed_par_t_alt(). sigma is the
//  SCALE phi, as in the sampler and in npll_loo_alt; the sd of P_n is
//  sigma_i * sqrt(df / (df - 2)) for df > 2.
//
//  i-th contribution:
//      log p_i = log dt(r_i / sigma_i ; df) - log sigma_i,
//      r_i     = y_i - x_i' beta_{-i}.
//
//    log_par : c(log(df), log(sigma))      (same order as npll_loo_alt)
//    y       : response vector, length n
//    X       : n x p matrix (with or without an intercept column)
//    beta_0  : initial mean, length p
//    V_0     : p x p symmetric positive-definite matrix
//
//  Reuses, from the same file:
//    - BIG and dt_log_2()   (defined before npll_loo)
//    - loo_core_lm()        (defined with npll_loo_lm_skew_alt)
//  so it must be placed AFTER npll_loo_lm_skew_alt.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double npll_loo_lm_t_alt(const arma::vec& log_par,
                         const arma::vec& y,
                         const arma::mat& X,
                         const arma::vec& beta_0,
                         const arma::mat& V_0) {
  
  const uword n = y.n_elem;
  const uword p = X.n_cols;
  
  if (X.n_rows != n)       Rcpp::stop("nrow(X) must equal length(y)");
  if (beta_0.n_elem != p)  Rcpp::stop("length(beta_0) must equal ncol(X)");
  if (V_0.n_rows != p || V_0.n_cols != p)
    Rcpp::stop("V_0 must be p x p with p = ncol(X)");
  if (log_par.n_elem != 2) Rcpp::stop("log_par must be c(log(df), log(sigma))");
  if (!log_par.is_finite()) return BIG;
  
  const double df    = std::exp(log_par(0));
  const double sigma = std::exp(log_par(1));
  
  if (!(df > 0.0) || !(sigma > 0.0) || !std::isfinite(sigma)) return BIG;
  
  vec beta, h, r;
  loo_core_lm(y, X, beta_0, V_0, beta, h, r);
  
  // -sum_i [ log dt(r_i/sigma_i ; df) - log sigma_i ],
  //   sigma_i     = sigma / sqrt(1 - h_i)
  //   log sigma_i = log sigma - 0.5 * log(1 - h_i)
  const double log_sigma = std::log(sigma);
  double out = 0.0;
  for (uword i = 0; i < n; ++i) {
    const double log1mh      = std::log1p(-h(i));            // log(1 - h_i)
    const double sigma_i     = sigma / std::sqrt(1.0 - h(i));
    const double log_sigma_i = log_sigma - 0.5 * log1mh;
    out += dt_log_2(r(i) / sigma_i, df) - log_sigma_i;
  }
  
  if (!std::isfinite(out)) return BIG;
  return -out;
}
