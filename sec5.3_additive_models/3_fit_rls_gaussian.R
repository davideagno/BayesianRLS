# Section 5.3, step 2: proposed method with Gaussian P_n (RLS-n).

rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)
library(extraDistr)

Rcpp::sourceCpp("../RLS_functions.cpp")
source("../RLS_functions.R")

# Preliminary settings
p <- 15
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 30  # n = n_mult * p
n_obs <- n_mult*p
dir.create("output", showWarnings = FALSE)
nSim <- 30
n_test <- 1e4

sigma_true <- 1
df_true <- 5
name_load <- file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, ".RData"))
load(name_load)

time_rls_norm <- time_est_norm <- time_opt_norm <- time_ps_norm <- rep(NA, nSim)
mse_rls_norm <- rep(NA, nSim)
ess_rls_norm <- in_ci_rls_norm <- width_ci_rls_norm <- matrix(nrow = n_test, ncol = nSim)
sigma_est_norm <- rep(NA, nSim)

mse_rls_q_norm <- rep(NA, nSim)
in_ci_rls_q_norm <- width_rls_q_norm <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90
# Posterior summaries (mean, sd, HPD, 2.5% and 97.5% quantiles) at each test point,
# with the corresponding true values, for any further analysis
summary_f <- summary_q <- array(NA, dim = c(n_test, 6, nSim))
truth_f <- truth_q <- matrix(NA, nrow = n_test, ncol = nSim)
lambda_est <- matrix(NA, nrow = p, ncol = nSim)

ncores <- 8

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
  b_0 <- rep(0, NCOL(X_splines))
  
  # Test set
  dat_test <- dat_test_sim
  X_test <- dat_test$X
  df_test <- data.frame(y = dat_test$y, dat_test$X)
  X_splines_test <- cbind(1, do.call(cbind,
                                     lapply(smooths, function(sm) mgcv::PredictMat(sm, df_test))))
  # Truth centered as the training data (see true_f_test)
  mean_true <- dat$mu + rowSums(true_f_test(dat, X_test))
  
  # Optimization - GAUSSIAN
  df_norm <- 10000
  t_opt <- system.time(opt <- optimize_loo_pll(c(log(1), rep(0.5, p)), dat$y, X_splines, b_0, Omega,
                                               fixed_df = df_norm))
  sigma_opt <- exp(opt$par[1])
  lambda_opt <- exp(opt$par[-1])
  
  lambda_est[,s] <- lambda_opt
  
  # Point Estimate
  t_est <- system.time(bb <- rls_lambda(dat$y, X_splines, b_0, lambda_opt, Omega))
  
  # Predictive resampling
  t_ps <- system.time(beta_ps <- ps_rs_parallel_t_fixed_alt(X_splines, bb$beta, bb$P, N = 9 * n_obs, B = 5000, 
                                                            df = df_norm, sigma = sigma_opt, seed = p*s, nthreads = ncores,
                                                            correct = TRUE))
  beta_ps <- as.mcmc(beta_ps)
  
  
  # Prediction
  mu_rls_test <- c(X_splines_test %*% bb$beta)
  mse_rls_norm[s] <- mean((mean_true - mu_rls_test)^2)
  
  # UQ for f(x) and Q_0.90(x): posterior summaries at each test point
  f_draws <- X_splines_test %*% t(beta_ps)   # n_test x number of draws
  q_draws <- f_draws + rep(sigma_opt*qnorm(lev), each = n_test)
  q_true <- mean_true + sigma_true*qt(lev, df_true)
  sum_f <- summarize_draws(f_draws)
  sum_q <- summarize_draws(q_draws)
  summary_f[,,s] <- sum_f
  summary_q[,,s] <- sum_q
  truth_f[,s] <- mean_true
  truth_q[,s] <- q_true
  in_ci_rls_norm[,s] <- 1*(sum_f[,"hpd_lo"] <= mean_true & mean_true <= sum_f[,"hpd_hi"])
  width_ci_rls_norm[,s] <- sum_f[,"hpd_hi"] - sum_f[,"hpd_lo"]
  mse_rls_q_norm[s] <- mean((sum_q[,"mean"] - q_true)^2)
  in_ci_rls_q_norm[,s] <- 1*(sum_q[,"hpd_lo"] <= q_true & q_true <= sum_q[,"hpd_hi"])
  width_rls_q_norm[,s] <- sum_q[,"hpd_hi"] - sum_q[,"hpd_lo"]
  
  # ESS
  ess_rls_norm[,s] <- coda::effectiveSize(t(f_draws))
  
  # Time
  time_opt_norm[s] <- t_opt[3]
  time_est_norm[s] <- t_est[3]
  time_ps_norm[s] <- t_ps[3]
  time_rls_norm[s] <- t_opt[3] + t_est[3] + t_ps[3]
  
  # Hyperparameters
  sigma_est_norm[s] <- sigma_opt
  
  
  name_out <- file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_norm_alt.RData"))
  save(mse_rls_norm, ess_rls_norm, in_ci_rls_norm, width_ci_rls_norm, 
       time_rls_norm, time_opt_norm, time_est_norm, time_ps_norm, sigma_est_norm, 
       mse_rls_q_norm, in_ci_rls_q_norm, width_rls_q_norm,
       summary_f, summary_q, truth_f, truth_q, lambda_est, file = name_out)
}

print("COMPLETED")

