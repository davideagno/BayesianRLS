# Section 5.3, step 4: BART.

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

time_bart <- mse_bart <- rep(NA, nSim)
ess_bart <- in_ci_bart <- width_ci_bart <- matrix(nrow = n_test, ncol = nSim)

mse_bart_q <- rep(NA, nSim)
in_ci_bart_q <- width_bart_q <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90
ndpost <- 5000   # BART draws kept
nskip <- 50000   # BART burn-in
# Posterior summaries (mean, sd, HPD, 2.5% and 97.5% quantiles) at each test point,
# with the corresponding true values, for any further analysis
summary_f <- summary_q <- array(NA, dim = c(n_test, 6, nSim))
truth_f <- truth_q <- matrix(NA, nrow = n_test, ncol = nSim)

ncores = 8

for (s in 1:nSim) {
  
  print(s)
  
  # Train set & settings
  dat <- dat_sim[[s]]
  df <- data.frame(y = dat$y, dat$X)
  
  # Test set
  dat_test <- dat_test_sim
  df_test <- data.frame(y = dat_test$y, dat_test$X)
  # Truth centered as the training data (see true_f_test)
  mean_true <- dat$mu + rowSums(true_f_test(dat, dat_test$X))
  
  # BART
  t_bart <- system.time(fit_bart <- BART::gbart(x.train = dat$X, y.train = dat$y,
                                                x.test = dat_test$X, ndpost = ndpost, nskip = nskip, mc.cores = ncores))
  
  # Prediction
  mu_bart_test <- apply(fit_bart$yhat.test, 2, mean)
  mse_bart[s] <- mean((mean_true - mu_bart_test)^2)
  
  # UQ for f(x) and Q_0.90(x): posterior summaries at each test point
  f_draws <- t(fit_bart$yhat.test)   # n_test x number of draws
  q_draws <- f_draws + rep(fit_bart$sigma[-(1:nskip)]*qnorm(lev), each = n_test)
  q_true <- mean_true + sigma_true*qt(lev, df_true)
  sum_f <- summarize_draws(f_draws)
  sum_q <- summarize_draws(q_draws)
  summary_f[,,s] <- sum_f
  summary_q[,,s] <- sum_q
  truth_f[,s] <- mean_true
  truth_q[,s] <- q_true
  in_ci_bart[,s] <- 1*(sum_f[,"hpd_lo"] <= mean_true & mean_true <= sum_f[,"hpd_hi"])
  width_ci_bart[,s] <- sum_f[,"hpd_hi"] - sum_f[,"hpd_lo"]
  mse_bart_q[s] <- mean((sum_q[,"mean"] - q_true)^2)
  in_ci_bart_q[,s] <- 1*(sum_q[,"hpd_lo"] <= q_true & q_true <= sum_q[,"hpd_hi"])
  width_bart_q[,s] <- sum_q[,"hpd_hi"] - sum_q[,"hpd_lo"]
  
  # ESS
  ess_bart[,s] <- coda::effectiveSize(t(f_draws))
  
  # Time
  time_bart[s] <- t_bart[3]
  
  
  name_out <- file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_bart.RData"))
  save(mse_bart, in_ci_bart, width_ci_bart, ess_bart, time_bart, 
       mse_bart_q, in_ci_bart_q, width_bart_q,
       summary_f, summary_q, truth_f, truth_q, file = name_out)
}

print("COMPLETED")

