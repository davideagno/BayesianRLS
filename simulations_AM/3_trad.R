
rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)
library(extraDistr)

Rcpp::sourceCpp("RLS_functions.cpp")
source("RLS_functions.R")

# Preliminary settings
p <- 15
n_obs <- 30*p
nSim <- 100
n_test <- 1e4

sigma_true <- 1
df_true <- 5
name_load <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, ".RData")
load(name_load)

time_trad <- mse_trad <- rep(NA, nSim)
ess_trad <- in_ci_trad <- width_ci_trad <- matrix(nrow = n_test, ncol = nSim)

mse_trad_q <- rep(NA, nSim)
in_ci_trad_q <- width_trad_q <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90

ncores = 8

for (s in 1:nSim) {
  
  print(s)
  
  # Train set & settings
  dat <- dat_sim[[s]]
  df <- data.frame(y = dat$y, dat$X)
  smooths <- lapply(1:p, function(j) {
    spec <- eval(call("s", as.name(paste0("x", j))))   # builds s(xj)
    mgcv::smoothCon(spec, data = df, knots = NULL, absorb.cons = TRUE)[[1]]
  })
  Omega <- Z <- list(); X_splines <- rep(1, length(dat$y))
  for (j in 1:p) {
    Z[[j]] <- smooths[[j]]$X
    Omega[[j]] <- smooths[[j]]$S[[1]]
    X_splines <- cbind(X_splines, Z[[j]])
  }
  
  # Test set
  dat_test <- dat_test_sim
  X_test <- dat_test$X
  df_test <- data.frame(y = dat_test$y, dat_test$X)
  X_splines_test <- cbind(1, do.call(cbind,
                                     lapply(smooths, function(sm) mgcv::PredictMat(sm, df_test))))
  mean_true <- dat_test$mu + rowSums(dat_test$F_true)
  
  d <- vapply(Z, ncol, integer(1))
  D <- d[1]; stopifnot(all(d == D))
  ranks <- vapply(Omega, function(x) sum(eigen(x, symmetric = TRUE, only.values = TRUE)$values > 1e-10),
                  integer(1))
  dat_splines <- data.frame(y = dat$y, X_splines)
  colnames(dat_splines) <- c("y", paste0("X", 1:((p*D)+1)))
  # Pack the penalty matrices into a 3-D array Omega_arr[j, , ]
  Omega_arr <- array(0, dim = c(p, D, D))
  for (j in 1:p) Omega_arr[j, , ] <- Omega[[j]]
  
  # custom Stan code
  stanvars <-
    stanvar(as.integer(p),     name = "P",     scode = "  int<lower=1> P;") +
    stanvar(as.integer(D),     name = "D",     scode = "  int<lower=1> D;") +
    stanvar(as.integer(ranks), name = "rk",    scode = "  array[P] int<lower=1> rk;") +
    stanvar(Omega_arr,         name = "Omega", scode = "  array[P] matrix[D, D] Omega;") +
    stanvar(scode = "  vector<lower=0>[P] lambda;  // smoothing parameters",
            block = "parameters") +
    stanvar(scode = "
    // weakly-informative prior on the smoothing parameters
    target += gamma_lpdf(lambda | 2, 0.5);
    
    // beta ~ N(0, sigma^2 * S^{-1}),  S = blockdiag(lambda_j * Omega_j)
    // => precision per block is (lambda_j / sigma^2) * Omega_j
    for (j in 1:P) {
      vector[D] beta_j = b[((j - 1) * D + 1):(j * D)];
      target += -0.5 * (lambda[j] / square(sigma)) * quad_form_sym(Omega[j], beta_j)
                + 0.5 * D * log(lambda[j]) - D * log(sigma);
    }",
            block = "model")
  
  # Fit
  t_trad <- system.time(fit_trad <- brm(formula = y ~ -1+., data = dat_splines, family = student(),
                                        stanvars = stanvars, prior = prior("", class = "b"), 
                                        chains = 5, cores = 5, iter = 20000, warmup = 19000, seed = p*s))
  post_samp <- coda::as.mcmc(as.matrix(fit_trad))
  beta_mcmc <- post_samp[,1:(9*p+1)]
  sigma_mcmc <- post_samp[,9*p+2]
  df_mcmc <- post_samp[,9*p+3]
  
  
  # Prediction
  mu_trad_test <- c(X_splines_test %*% colMeans(beta_mcmc))
  mse_trad[s] <- mean((mean_true - mu_trad_test)^2)
  
  # UQ
  in_ci_test <- rep(0, n_test)
  width_ci_test <- rep(NA, n_test)
  for (i in 1:n_test) {
    rp_i <- c(beta_mcmc %*% X_splines_test[i,])
    ci_i <- robustBLME::hpd(rp_i)
    width_ci_test[i] <- abs(ci_i[2] - ci_i[1])
    if (ci_i[1] <= mean_true[i] & mean_true[i] <= ci_i[2]) in_ci_test[i] <- 1
  }
  in_ci_trad[,s] <- in_ci_test
  width_ci_trad[,s] <- width_ci_test

  # Quantile
  sq_err <- rep(NA, n_test)
  in_ci_test_q <- rep(0, n_test)
  width_test_q <- rep(NA, n_test)
  for (i in 1:n_test) {
    pq_i <- c(beta_mcmc %*% X_splines_test[i,]) + sigma_mcmc*qt(lev, df = df_mcmc)
    ci_i <- robustBLME::hpd(pq_i)
    width_test_q[i] <- abs(ci_i[2] - ci_i[1])
    pq_true <- mean_true[i] + sigma_true*qt(lev, df_true)
    if (ci_i[1] <= pq_true & pq_true <= ci_i[2]) in_ci_test_q[i] <- 1
    sq_err[i] <- (mean(pq_i) - pq_true)^2
  }
  mse_trad_q[s] <- mean(sq_err)
  in_ci_trad_q[,s] <- in_ci_test_q
  width_trad_q[,s] <- width_test_q
  
  # ESS
  mu_mcmc_matrix <- X_splines_test %*% t(beta_mcmc)
  mu_mcmc_matrix_t <- t(mu_mcmc_matrix)
  ess_trad[,s] <- coda::effectiveSize(mu_mcmc_matrix_t)
  
  # Time
  time_trad[s] <- t_trad[3]
  
  
  name_out <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_trad_1.RData")
  save(mse_trad, in_ci_trad, width_ci_trad, ess_trad, time_trad, 
       mse_trad_q, in_ci_trad_q, width_trad_q, file = name_out)
}

print("COMPLETED")

