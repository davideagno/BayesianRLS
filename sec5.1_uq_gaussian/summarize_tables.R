# Summaries for Table 1 (output of simulate.R). For Table S1, load instead the
# output of ../secS2.2_uq_student_t/simulate.R (second argument "t").
#
# For each choice of sigma_n, the three entries of each line refer to predictive
# resampling with Algorithm 1 (PS-1), Algorithm 2 (PS-2), and the exact posterior
# (Monte Carlo reference for the errors; Gaussian approximation for Student-t responses).
#   log ||mean - beta_n||^2 / p   : Monte Carlo error of the posterior mean
#   log ||cov - Sigma_n||_F       : error of the posterior covariance (lower triangle)
#   width, cover                  : average width and coverage of the 95% intervals
#   time                          : seconds (PS-1 and PS-2 only)

rm(list = ls())

p <- 50
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 3  # n = n_mult * p
n_obs <- n_mult*p
model <- if (length(commandArgs(TRUE)) > 1) commandArgs(TRUE)[2] else "norm"  # "norm" or "t"
dir_load <- if (model == "t") "../secS2.2_uq_student_t/output" else "output"
name_load <- file.path(dir_load, paste0("sim_th_p", p, "_n", n_obs, "_", model, ".RData"))
sink(file.path(dir_load, paste0("summary_", model, "_n", n_mult, "p.txt")), split = TRUE)
load(name_load)


# sigma_n^2 = sigma^2 / a_n
round(c(mean(log(mse_A)), mean(log(mse_B)), mean(log(mse_ex_Vn))), 2)
round(c(mean(log(fr_Vn_A)), mean(log(fr_Vn_B)), mean(log(fr_Vn_ex))), 2)
round(c(mean(width_A), mean(width_B), mean(width_Vn)), 3)
round(c(mean(rowMeans(in_ci_A)), mean(rowMeans(in_ci_B)), mean(rowMeans(in_ci_Vn))), 3)
round(c(mean(time_A), mean(time_B)), 2)

round(c(sd(log(mse_A)), sd(log(mse_B)), sd(log(mse_ex_Vn))), 2)
round(c(sd(log(fr_Vn_A)), sd(log(fr_Vn_B)), sd(log(fr_Vn_ex))), 2)
round(c(sd(width_A), sd(width_B), sd(width_Vn)), 3)
round(c(sd(colMeans(in_ci_A)), sd(colMeans(in_ci_B)), sd(colMeans(in_ci_Vn))), 3)
round(c(sd(time_A), sd(time_B)), 2)


# sigma_n = sigma
round(c(mean(log(mse_1)), mean(log(mse_2)), mean(log(mse_ex_Sn))), 2)
round(c(mean(log(fr_Sn_1)), mean(log(fr_Sn_2)), mean(log(fr_Sn_ex))), 2)
round(c(mean(width_1), mean(width_2), mean(width_Sn)), 3)
round(c(mean(rowMeans(in_ci_1)), mean(rowMeans(in_ci_2)), mean(rowMeans(in_ci_Sn))), 3)
round(c(mean(time_1), mean(time_2)), 2)

round(c(sd(log(mse_1)), sd(log(mse_2)), sd(log(mse_ex_Sn))), 2)
round(c(sd(log(fr_Sn_1)), sd(log(fr_Sn_2)), sd(log(fr_Sn_ex))), 2)
round(c(sd(width_1), sd(width_2), sd(width_Sn)), 3)
round(c(sd(colMeans(in_ci_1)), sd(colMeans(in_ci_2)), sd(colMeans(in_ci_Sn))), 3)
round(c(sd(time_1), sd(time_2)), 2)


# Fraction of posterior variance omitted by stopping at N = 10 n, before the correction
# (mean and range over datasets), close to n / N = 0.1
round(c(mean(frac_omit), range(frac_omit)), 3)

sink()
