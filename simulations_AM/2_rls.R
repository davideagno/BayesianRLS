
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

time_rls <- time_est <- time_opt <- time_ps <- rep(NA, nSim)
mse_rls <- rep(NA, nSim)
ess_rls <- in_ci_rls <- width_ci_rls <- matrix(nrow = n_test, ncol = nSim)
df_est <- sigma_est <- rep(NA, nSim)

mse_rls_q <- rep(NA, nSim)
in_ci_rls_q <- width_rls_q <- matrix(nrow = n_test, ncol = nSim)
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
  
  # Optimization
  t_opt <- system.time(opt <- nlminb(start = c(log(2), log(1), rep(0.5, p)), objective = npll_loo_alt,
                                     y = dat$y, X = X_splines, beta_0 = b_0, Omega = Omega))
  t_opt
  opt
  df_opt <- exp(opt$par[1])
  sigma_opt <- exp(opt$par[2])
  lambda_opt <- exp(opt$par[-c(1:2)])
  
  # Point Estimate
  t_est <- system.time(bb <- rls_lambda(dat$y, X_splines, b_0, lambda_opt, Omega))
  
  # Predictive resampling
  t_ps <- system.time(beta_ps <- ps_rs_parallel_t_fixed_alt(X_splines, bb$beta, bb$P, N = 1e5, B = 5000, 
                                                            df = df_opt, sigma = sigma_opt, seed = p*s, nthreads = ncores))
  beta_ps <- as.mcmc(beta_ps)
  
  
  # Prediction
  mu_rls_test <- c(X_splines_test %*% bb$beta)
  mse_rls[s] <- mean((mean_true - mu_rls_test)^2)
  
  # UQ
  in_ci_test <- rep(0, n_test)
  width_ci_test <- rep(NA, n_test)
  for (i in 1:n_test) {
    rp_i <- c(beta_ps %*% X_splines_test[i,])
    ci_i <- robustBLME::hpd(rp_i)
    width_ci_test[i] <- abs(ci_i[2] - ci_i[1])
    if (ci_i[1] <= mean_true[i] & mean_true[i] <= ci_i[2]) in_ci_test[i] <- 1
  }
  in_ci_rls[,s] <- in_ci_test
  width_ci_rls[,s] <- width_ci_test

  # Quantile
  sq_err <- rep(NA, n_test)
  in_ci_test_q <- rep(0, n_test)
  width_test_q <- rep(NA, n_test)
  for (i in 1:n_test) {
    pq_i <- c(beta_ps %*% X_splines_test[i,]) + sigma_opt*qt(lev, df_opt)
    ci_i <- robustBLME::hpd(pq_i)
    width_test_q[i] <- abs(ci_i[2] - ci_i[1])
    pq_true <- mean_true[i] + sigma_true*qt(lev, df_true)
    if (ci_i[1] <= pq_true & pq_true <= ci_i[2]) in_ci_test_q[i] <- 1
    sq_err[i] <- (mean(pq_i) - pq_true)^2
  }
  mse_rls_q[s] <- mean(sq_err)
  in_ci_rls_q[,s] <- in_ci_test_q
  width_rls_q[,s] <- width_test_q
  
  # ESS
  mu_ps_matrix <- X_splines_test %*% t(beta_ps)
  mu_ps_matrix_t <- t(mu_ps_matrix)
  ess_rls[,s] <- coda::effectiveSize(mu_ps_matrix_t)
  
  # Time
  time_opt[s] <- t_opt[3]
  time_est[s] <- t_est[3]
  time_ps[s] <- t_ps[3]
  time_rls[s] <- t_opt[3] + t_est[3] + t_ps[3]
  
  # Hyperparameters
  df_est[s] <- df_opt
  sigma_est[s] <- sigma_opt
                                            
  
  name_out <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_rls_alt.RData")
  save(mse_rls, ess_rls, in_ci_rls, width_ci_rls,
       time_rls, time_opt, time_est, time_ps, df_est, sigma_est, 
       mse_rls_q, in_ci_rls_q, width_rls_q, file = name_out)
}

print("COMPLETED")
