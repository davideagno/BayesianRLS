# Section 5.2, step 3: traditional Bayesian additive model (BLM), fitted with brms.

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

time_trad <- mse_trad <- rep(NA, nSim)
ess_trad <- in_ci_trad <- width_ci_trad <- matrix(nrow = n_test, ncol = nSim)

mse_trad_q <- rep(NA, nSim)
in_ci_trad_q <- width_trad_q <- matrix(nrow = n_test, ncol = nSim)
lev <- 0.90
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
  # Truth centered as the training data (see true_f_test)
  mean_true <- dat$mu + rowSums(true_f_test(dat, X_test))
  
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
                + 0.5 * rk[j] * log(lambda[j]) - rk[j] * log(sigma);
    }",
            block = "model")
  
  # Compile the Stan model once (not timed); later fits reuse it through update()
  if (s == 1) {
    cat("Compiling the Stan model\n")
    fit_template <- brm(formula = y ~ -1+., data = dat_splines, family = student(),
                        stanvars = stanvars, prior = prior("", class = "b"), chains = 0)
  }

  # Fit (sampling only)
  t_trad <- system.time(fit_trad <- update(fit_template, newdata = dat_splines, stanvars = stanvars,
                                           recompile = FALSE, chains = 5, cores = 5,
                                           iter = 20000, warmup = 19000, seed = p*s))
  post_samp <- coda::as.mcmc(as.matrix(fit_trad))
  beta_mcmc <- post_samp[,1:(9*p+1)]
  sigma_mcmc <- post_samp[,9*p+2]
  df_mcmc <- post_samp[,9*p+3]
  
  
  # Prediction
  mu_trad_test <- c(X_splines_test %*% colMeans(beta_mcmc))
  mse_trad[s] <- mean((mean_true - mu_trad_test)^2)
  
  # UQ for f(x) and Q_0.90(x): posterior summaries at each test point
  f_draws <- X_splines_test %*% t(beta_mcmc)   # n_test x number of draws
  q_draws <- f_draws + rep(sigma_mcmc*qt(lev, df = df_mcmc), each = n_test)
  q_true <- mean_true + sigma_true*qt(lev, df_true)
  sum_f <- summarize_draws(f_draws)
  sum_q <- summarize_draws(q_draws)
  summary_f[,,s] <- sum_f
  summary_q[,,s] <- sum_q
  truth_f[,s] <- mean_true
  truth_q[,s] <- q_true
  in_ci_trad[,s] <- 1*(sum_f[,"hpd_lo"] <= mean_true & mean_true <= sum_f[,"hpd_hi"])
  width_ci_trad[,s] <- sum_f[,"hpd_hi"] - sum_f[,"hpd_lo"]
  mse_trad_q[s] <- mean((sum_q[,"mean"] - q_true)^2)
  in_ci_trad_q[,s] <- 1*(sum_q[,"hpd_lo"] <= q_true & q_true <= sum_q[,"hpd_hi"])
  width_trad_q[,s] <- sum_q[,"hpd_hi"] - sum_q[,"hpd_lo"]
  
  # ESS
  ess_trad[,s] <- coda::effectiveSize(t(f_draws))
  
  # Time
  time_trad[s] <- t_trad[3]
  
  
  name_out <- file.path("output", paste0("sim_p", p, "_n", n_obs, "_df", df_true, "_sigma", sigma_true, "_trad.RData"))
  save(mse_trad, in_ci_trad, width_ci_trad, ess_trad, time_trad, 
       mse_trad_q, in_ci_trad_q, width_trad_q,
       summary_f, summary_q, truth_f, truth_q, file = name_out)
}

print("COMPLETED")

