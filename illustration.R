
rm(list = ls())

library(mgcv)
library(coda)
library(brms)
library(parallel)

Rcpp::sourceCpp("RLS_functions.cpp")
source("RLS_functions.R")

p <- 5
n_obs <- 15
n_all <- 1e5
sigma_true <- 0.2
df_true <- 3
var_true <- sigma_true^2 * df_true / (df_true-2)

set.seed(1997); beta_true <- rnorm(p, sd = 2)
set.seed(2)
X_all <- matrix(rnorm(n_all*p), nrow = n_all, ncol = p)
X <- X_all[1:n_obs,]
X_new <- X_all[(n_obs+1):n_all,]
mu_true <- c(X %*% beta_true)
eps <- rt(n_obs, df = df_true)
y <- mu_true + sigma_true*eps

# RLS
b_0 <- rep(0, p)
V_0 <- diag(p)*1000
bb <- rls(y, X, b_0, V_0)

# Fixed covariates for predictive resamplings
set.seed(5)
X_ps <- matrix(nrow = NROW(X_new), ncol = p)
for (i in 1:NROW(X_ps)) {
  idx <- sample(1:n_obs, size = 1)
  X_ps[i,] <- X[idx,]
}

# Student-t P_n
opt <- nlminb(start = c(2,2), objective = npll_loo_lm_t_alt,
              y = y, X = X, beta_0 = b_0, V_0 = V_0)
opt
(df_opt <- exp(opt$par[1]))
(sigma_opt <- exp(opt$par[2])) 
var_opt <- sigma_opt^2 * df_opt / (df_opt - 2)
B_ps <- 10000
beta_ps <-  ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_opt, df_opt, nthreads = 7)

Gaussian P_n
opt_n <- nlminb(start = 1, objective = function(x) npll_loo_lm_t_alt(c(log(1e5), x), y, X, b_0, V_0))
(sigma_opt_n <- exp(opt_n$par))
beta_ps_n <-  ps_lm_fixed_par_t_alt(bb$beta, bb$P, X_ps, B_ps, sigma_opt_n, 1e5, nthreads = 7)

library(ggcorrplot)
library(latex2exp)
library(ggthemes)

# Posterior
idx <- 4
df_d1 <- data.frame(rls = as.numeric(beta_ps[,idx]), 
                    blm = as.numeric(beta_ps_n[,idx]), 
                    title = "Posterior distribution")
pl <- ggplot(df_d1) + 
  facet_wrap(~ title) +
  stat_density(aes(x = rls, col = "rls"),
               geom = "line", position = "identity", adjust = 1.8) +
  stat_density(aes(x = blm, col = "blm"),
               geom = "line", position = "identity", adjust = 1.8, linetype = "dashed") +
  scale_color_tableau( name = "", palette = "Color Blind",
                       labels = c("rls" = "Student-t", 
                                  "blm" = "Gaussian")) +
  xlab("Beta") + ylab("Density") + xlim(c(-1, 0.5)) +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl

# Plug-in predictive
y_val <- seq(-1.8, 1.8, length = 2000)
df_plot <- data.frame(
  y     = rep(y_val, 3),
  dens  = c(dnorm(y_val, sd = sigma_opt_n),
            dlst(y_val, df = df_opt,  sigma = sigma_opt),
            dlst(y_val, df = df_true, sigma = sigma_true)),
  distr = factor(rep(c("rls_n", "rls_t", "true"), each = length(y_val)),
                 levels = c("rls_n", "rls_t", "true")),
  title = "Data generating distribution"
)
lab_distr <- c("rls_t" = "Student-t", "rls_n" = "Gaussian", "true" = "True")
pl_p <- ggplot(df_plot, aes(x = y, y = dens, colour = distr, linetype = distr)) +
  facet_wrap(~ title) +
  geom_line(linewidth = 0.5) +
  scale_color_tableau(name = "", palette = "Color Blind", labels = lab_distr) +
  scale_linetype_manual(name = "", labels = lab_distr,
                        values = c("rls_t" = "solid", "rls_n" = "dashed", "true" = "dotdash")) +
  labs(x = "Y", y = "Density") +
  theme_light() +
  theme(legend.text.align = 0,
        strip.text = element_text(size = 11, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
pl_p

pl_tot <- ggpubr::ggarrange(pl_p, pl)
pl_tot
# ggsave(filename = "pl_ill.pdf", plot = pl_tot, width = 9, height = 5)

