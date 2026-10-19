
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

time_rls_norm <- time_est_norm <- time_opt_norm <- time_ps_norm <- rep(NA, nSim)
mse_rls_norm <- rep(NA, nSim)
ess_rls_norm <- in_ci_rls_norm <- width_ci_rls_norm <- matrix(nrow = n_test, ncol = nSim)
sigma_est_norm <- rep(NA, nSim)

mse_rls_q_norm <- rep(NA, nSim)
in_ci_rls_q_norm <- width_rls_q_norm <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90

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
  mean_true <- dat_test$mu + rowSums(dat_test$F_true)
  
  # Optimization - GAUSSIAN
  df_norm <- 10000
  t_opt <- system.time(opt <- nlminb(start = c(log(1), rep(0.5, p)), 
                                     objective = function(x) npll_loo_alt(c(log(df_norm), x),
                                                                      y = dat$y, 
                                                                      X = X_splines, 
                                                                      beta_0 = b_0, 
                                                                      Omega = Omega)))
  sigma_opt <- exp(opt$par[1])
  lambda_opt <- exp(opt$par[-1])
  
  # Point Estimate
  t_est <- system.time(bb <- rls_lambda(dat$y, X_splines, b_0, lambda_opt, Omega))
  
  # Predictive resampling
  t_ps <- system.time(beta_ps <- ps_rs_parallel_t_fixed_alt(X_splines, bb$beta, bb$P, N = 1e5, B = 5000, 
                                                            df = df_norm, sigma = sigma_opt, seed = p*s, nthreads = ncores))
  beta_ps <- as.mcmc(beta_ps)
  
  
  # Prediction
  mu_rls_test <- c(X_splines_test %*% bb$beta)
  mse_rls_norm[s] <- mean((mean_true - mu_rls_test)^2)
  
  # UQ
  in_ci_test <- rep(0, n_test)
  width_ci_test <- rep(NA, n_test)
  for (i in 1:n_test) {
    rp_i <- c(beta_ps %*% X_splines_test[i,])
    ci_i <- robustBLME::hpd(rp_i)
    width_ci_test[i] <- abs(ci_i[2] - ci_i[1])
    if (ci_i[1] <= mean_true[i] & mean_true[i] <= ci_i[2]) in_ci_test[i] <- 1
  }
  in_ci_rls_norm[,s] <- in_ci_test
  width_ci_rls_norm[,s] <- width_ci_test

  # Quantile
  sq_err <- rep(NA, n_test)
  in_ci_test_q <- rep(0, n_test)
  width_test_q <- rep(NA, n_test)
  for (i in 1:n_test) {
    pq_i <- c(beta_ps %*% X_splines_test[i,]) + sigma_opt*qnorm(lev)
    ci_i <- robustBLME::hpd(pq_i)
    width_test_q[i] <- abs(ci_i[2] - ci_i[1])
    pq_true <- mean_true[i] + sigma_true*qt(lev, df_true)
    if (ci_i[1] <= pq_true & pq_true <= ci_i[2]) in_ci_test_q[i] <- 1
    sq_err[i] <- (mean(pq_i) - pq_true)^2
  }
  mse_rls_q_norm[s] <- mean(sq_err)
  in_ci_rls_q_norm[,s] <- in_ci_test_q
  width_rls_q_norm[,s] <- width_test_q
  
  # ESS
  mu_ps_matrix <- X_splines_test %*% t(beta_ps)
  mu_ps_matrix_t <- t(mu_ps_matrix)
  ess_rls_norm[,s] <- coda::effectiveSize(mu_ps_matrix_t)
  
  # Time
  time_opt_norm[s] <- t_opt[3]
  time_est_norm[s] <- t_est[3]
  time_ps_norm[s] <- t_ps[3]
  time_rls_norm[s] <- t_opt[3] + t_est[3] + t_ps[3]
  
  # Hyperparameters
  sigma_est_norm[s] <- sigma_opt
  
  
  name_out <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_norm_alt.RData")
  save(mse_rls_norm, ess_rls_norm, in_ci_rls_norm, width_ci_rls_norm, 
       time_rls_norm, time_opt_norm, time_est_norm, time_ps_norm, sigma_est_norm, 
       mse_rls_q_norm, in_ci_rls_q_norm, width_rls_q_norm, file = name_out)
}

print("COMPLETED")

