
rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)
library(extraDistr)

Rcpp::sourceCpp("AM_sim.cpp")
source("functions_AM_sim.R")

# Preliminary settings
p <- 15
n_obs <- 30*p
nSim <- 100
n_test <- 1e4

sigma_true <- 1
df_true <- 5
name_load <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, ".RData")
load(name_load)

time_bart <- mse_bart <- rep(NA, nSim)
ess_bart <- in_ci_bart <- width_ci_bart <- matrix(nrow = n_test, ncol = nSim)

mse_bart_q <- rep(NA, nSim)
in_ci_bart_q <- width_bart_q <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90

ncores = 8

for (s in 1:nSim) {
  
  print(s)
  
  # Train set & settings
  dat <- dat_sim[[s]]
  df <- data.frame(y = dat$y, dat$X)
  
  # Test set
  dat_test <- dat_test_sim
  df_test <- data.frame(y = dat_test$y, dat_test$X)
  mean_true <- dat_test$mu + rowSums(dat_test$F_true)
  
  # BART
  t_bart <- system.time(fit_bart <- BART::gbart(x.train = dat$X, y.train = dat$y,
                                                x.test = dat_test$X, ndpost = 1000, nskip = 2000, mc.cores = ncores))
  xx <- predict(fit_bart, dat_test$X)
  
  # Prediction
  mu_bart_test <- apply(fit_bart$yhat.test, 2, mean)
  mse_bart[s] <- mean((mean_true - mu_bart_test)^2)
  
  # UQ
  in_ci_test <- rep(0, n_test)
  width_ci_test <- rep(NA, n_test)
  for (i in 1:n_test) {
    rp_i <- fit_bart$yhat.test[,i]
    ci_i <- robustBLME::hpd(rp_i)
    width_ci_test[i] <- abs(ci_i[2] - ci_i[1])
    if (ci_i[1] <= mean_true[i] & mean_true[i] <= ci_i[2]) in_ci_test[i] <- 1
  }
  in_ci_bart[,s] <- in_ci_test
  width_ci_bart[,s] <- width_ci_test

  # Quantile
  sq_err <- rep(NA, n_test)
  in_ci_test_q <- rep(0, n_test)
  width_test_q <- rep(NA, n_test)
  for (i in 1:n_test) {
    pq_i <- fit_bart$yhat.test[,i] + fit_bart$sigma*qnorm(lev)
    ci_i <- robustBLME::hpd(pq_i)
    width_test_q[i] <- abs(ci_i[2] - ci_i[1])
    pq_true <- mean_true[i] + sigma_true*qt(lev, df_true)
    if (ci_i[1] <= pq_true & pq_true <= ci_i[2]) in_ci_test_q[i] <- 1
    sq_err[i] <- (mean(pq_i) - pq_true)^2
  }
  mse_bart_q[s] <- mean(sq_err)
  in_ci_bart_q[,s] <- in_ci_test_q
  width_bart_q[,s] <- width_test_q
  
  # ESS
  ess_bart[,s] <- coda::effectiveSize(fit_bart$yhat.test)
  
  # Time
  time_bart[s] <- t_bart[3]
  
  
  name_out <- paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_bart.RData")
  save(mse_bart, in_ci_bart, width_ci_bart, ess_bart, time_bart, 
       mse_bart_q, in_ci_bart_q, width_bart_q, file = name_out)
}

print("COMPLETED")

