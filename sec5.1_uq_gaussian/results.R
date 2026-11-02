# Summaries for Table 1 (output of sim_normal.R). For Table S1, load instead the
# output of ../secS2.2_uq_student_t/sim_t.R (file suffix _t).

rm(list = ls())

p <- 50
n_mult <- if (length(commandArgs(TRUE)) > 0) as.numeric(commandArgs(TRUE)[1]) else 3  # n = n_mult * p
n_obs <- n_mult*p
model <- if (length(commandArgs(TRUE)) > 1) commandArgs(TRUE)[2] else "norm"  # "norm" or "t"
dir_load <- if (model == "t") "../secS2.2_uq_student_t/output" else "output"
name_load <- file.path(dir_load, paste0("sim_th_p", p, "_n", n_obs, "_", model, ".RData"))
sink(file.path(dir_load, paste0("summary_", model, "_n", n_mult, "p.txt")), split = TRUE)
load(name_load)


round(c(mean(log(mse_A)), mean(log(mse_B))), 3)
round(c(mean(log(fr_Vn_A)), mean(log(fr_Vn_B))), 2)
round(c(mean(rowMeans(in_ci_A)), mean(rowMeans(in_ci_B)), mean(rowMeans(in_ci_Vn))), 3) 
round(c(mean(time_A), mean(time_B)), 1)

round(c(sd(log(mse_A)), sd(log(mse_B))), 3)
round(c(sd(log(fr_Vn_A)), sd(log(fr_Vn_B))), 2)
round(c(sd(rowMeans(in_ci_A)), sd(rowMeans(in_ci_B)), sd(rowMeans(in_ci_Vn))), 3) 
round(c(sd(time_A), sd(time_B)), 1)



round(c(mean(log(mse_1)), mean(log(mse_2))), 3)
round(c(mean(log(fr_Sn_1)), mean(log(fr_Sn_2))), 2)
round(c(mean(rowMeans(in_ci_1)), mean(rowMeans(in_ci_2)), mean(rowMeans(in_ci_Sn))), 3) 
round(c(mean(time_1), mean(time_2)), 1)

round(c(sd(log(mse_1)), sd(log(mse_2))), 3)
round(c(sd(log(fr_Sn_1)), sd(log(fr_Sn_2))), 2)
round(c(sd(rowMeans(in_ci_1)), sd(rowMeans(in_ci_2)), sd(rowMeans(in_ci_Sn))), 3) 
round(c(sd(time_1), sd(time_2)), 1)

sink()
