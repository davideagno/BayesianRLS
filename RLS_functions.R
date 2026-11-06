# Shared R functions (the C++ functions are in RLS_functions.cpp).

# Simulated data for Section 5.3: additive model with Student-t errors
generate_am_t <- function(n = 500, p = 4, df = 10, sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  X <- matrix(runif(n * p, -2, 2), nrow = n, ncol = p)
  colnames(X) <- paste0("x", 1:p)
  
  # Define component functions
  f <- list(
    function(x) 0.5 * x^3 - x,              # cubic
    function(x) 3 * x * exp(-x^2),          # bump
    function(x) x^2 - 1,                    # quadratic
    function(x) 2 * sin(pi * x)             # smooth sinusoidal
  )
  
  # Evaluate each component
  F_mat <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
    fj <- f[[(j - 1) %% length(f) + 1]]  # recycle if p > length(f)
    F_mat[, j] <- fj(X[, j])
  }
  colnames(F_mat) <- paste0("f", 1:p)
  
  # Center each component (identifiability)
  F_mat <- scale(F_mat, center = TRUE, scale = FALSE)
  
  # Intercept + noise
  mu <- 1.0
  # eps <- rnorm(n, 0, sigma)
  eps <- rt(n, df)
  y <- mu + rowSums(F_mat) + sigma*eps
  
  list(
    y = y,
    X = X,
    F_true = F_mat,
    mu = mu,
    df = df,
    n = n,
    p = p,
    true_functions = f
  )
}

# True additive components at new covariates, centered as in the training data
true_f_test <- function(obj, Xnew) {
  out <- matrix(nrow = NROW(Xnew), ncol = NCOL(Xnew))
  
  loc <- attr(obj$F_true, "scaled:center")
  f_list <- obj$true_functions
  
  for (j in 1:NCOL(Xnew)) {
    fj <- f_list[[(j - 1) %% length(f_list) + 1]]
    out[,j] <- fj(Xnew[,j])
    out[,j] <- out[,j] - loc[j]
  }
  colnames(out) <- paste0("f", 1:NCOL(Xnew))
  out
}


# B forward paths of length N of predictive resampling (Algorithm 1) with a Student-t
# P_n and sigma_n^2 = sigma^2 / a_n, recording the squared Euclidean distance between
# beta_m and beta_n along each path (Figure S6 of the supplement)
ps_rs_conv_t_1_alt <- function(X, beta, P, N, B, df, sigma, seed = 1) {
  out <- matrix(nrow = N, ncol = B)
  n <- NROW(X)
  for (b in 1:B) {
    print(b)
    beta_b <- beta
    P_b <- P
    for (i in 1:N) {
      idx_new <- sample(1:n, size = 1)
      x_new <- X[idx_new,]
      Px <- P_b %*% x_new
      a <- c(1 / (1 + crossprod(x_new, Px)))
      sigma_new <- sigma / sqrt(a)
      y_new <- sum(x_new * beta_b) + rt(1, df)*sigma_new
      beta_b <- beta_b + Px * (a * (y_new - sum(x_new * beta_b)))
      P_b <- P_b - a * tcrossprod(Px)
      out[i,b] <- sum((unlist(beta) - unlist(beta_b))^2)
    }
  }
  out
}


# Posterior summaries of each row of a matrix of draws (rows: quantities, columns: draws):
# mean, standard deviation, 95% HPD interval (as robustBLME::hpd, vectorized) and the
# 2.5% and 97.5% quantiles (type 7, as quantile()).
summarize_draws <- function(draws, prob = 0.95) {
  draws <- as.matrix(draws)
  m <- nrow(draws)
  n <- ncol(draws)
  nn <- round(n * min(prob, 1 - prob))
  s <- t(apply(draws, 1, sort))
  if (m == 1) s <- matrix(s, nrow = 1)
  len <- s[, (n - nn + 1):n, drop = FALSE] - s[, 1:nn, drop = FALSE]
  k <- max.col(-len, ties.method = "first")
  q_type7 <- function(a) {
    h <- (n - 1) * a + 1
    lo <- floor(h)
    s[, lo] + (h - lo) * (s[, min(lo + 1, n)] - s[, lo])
  }
  cbind(mean = rowMeans(draws),
        sd = sqrt(rowSums((draws - rowMeans(draws))^2) / (n - 1)),
        hpd_lo = s[cbind(seq_len(m), k)],
        hpd_hi = s[cbind(seq_len(m), n - nn + k)],
        q025 = q_type7((1 - prob) / 2),
        q975 = q_type7(1 - (1 - prob) / 2))
}


# Maximizes the LOO-PLL with nlminb(), using its value and analytic gradient from
# npll_loo_alt_grad(). X'X and X'y depend only on the data, so they are computed once
# here, for the dataset at hand. If fixed_df is given, df is not optimized and start is
# c(log(sigma), log(lambda)). The first n_par columns of X are unpenalized.
optimize_loo_pll <- function(start, y, X, beta_0, Omega, fixed_df = 0, n_par = 1) {
  XtX <- crossprod(X)
  Xty <- c(crossprod(X, y))
  last_par <- NULL
  last <- NULL
  value_grad <- function(par) {
    if (!identical(par, last_par)) {
      last <<- npll_loo_alt_grad(par, y, X, XtX, Xty, beta_0, Omega, fixed_df, n_par)
      last_par <<- par
    }
    last
  }
  nlminb(start = start,
         objective = function(par) value_grad(par)$value,
         gradient = function(par) c(value_grad(par)$gradient))
}


# ---------------------------------------------------------------------------
# Simple interface for the linear model (see the tutorial in README.md)
# ---------------------------------------------------------------------------

# Posterior mean by recursive least squares, started from beta_0 and V_0.
# Returns the posterior mean beta (beta_n in the paper) and the matrix V (V_n).
rls_fit <- function(y, X, beta_0 = rep(0, ncol(X)), V_0 = diag(1000, ncol(X))) {
  out <- rls(y, X, beta_0, V_0)
  beta <- c(out$beta)
  names(beta) <- colnames(X)
  list(beta = beta, V = out$P, y = y, X = X, beta_0 = beta_0, V_0 = V_0)
}

# Predictive distribution selected by maximizing the leave-one-out prequential
# log-likelihood (LOO-PLL): a location-scale Student-t with degrees of freedom nu
# and scale phi, or a Gaussian (nu = Inf) with standard deviation phi.
rls_predictive <- function(fit, family = c("student_t", "gaussian")) {
  family <- match.arg(family)
  if (family == "student_t") {
    opt <- nlminb(start = c(log(5), log(sd(fit$y))), objective = npll_loo_lm_t_alt,
                  y = fit$y, X = fit$X, beta_0 = fit$beta_0, V_0 = fit$V_0)
    nu <- exp(opt$par[1])
    phi <- exp(opt$par[2])
  } else {
    # dt_log_2() uses the Gaussian density for df >= 1e8
    opt <- nlminb(start = log(sd(fit$y)), objective = function(x)
      npll_loo_lm_t_alt(c(log(1e8), x), fit$y, fit$X, fit$beta_0, fit$V_0))
    nu <- Inf
    phi <- exp(opt$par)
  }
  list(family = family, nu = nu, phi = phi, loo_pll = -opt$objective)
}

# Independent draws from the martingale posterior by predictive resampling
# (Algorithm 2 of the paper), with n_future future covariates resampled from the
# observed ones. The omitted part of the series, about a fraction n / (n + n_future)
# of the posterior variance, is added back by the Gaussian correction of
# Algorithm 2. The random numbers follow set.seed().
rls_posterior <- function(fit, pred, B = 10000, n_future = 9 * nrow(fit$X)) {
  X_future <- fit$X[sample(nrow(fit$X), n_future, replace = TRUE), , drop = FALSE]
  draws <- ps_lm_fixed_par_t_alt(fit$beta, fit$V, X_future, B = B, sigma = pred$phi,
                                 df = pred$nu, seed = sample.int(.Machine$integer.max, 1),
                                 correct = TRUE)
  colnames(draws) <- colnames(fit$X)
  draws
}
