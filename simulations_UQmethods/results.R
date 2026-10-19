
rm(list = ls())

p <- 50
n_obs <- 10*p
name_load <- paste0("sim_th_p", p, "_n", n_obs, "_norm.RData")
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
