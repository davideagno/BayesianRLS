# Section 5.2: posterior distributions of a regression coefficient in the two examples
# (Student-t and skew Student-t errors), with the matching and the Gaussian P_n.
# Run after student_t_errors.R and skew_t_errors.R, which save the draws in output/.

rm(list = ls())

library(ggplot2)

examples <- c("draws_student_t.RData" = "Student-t errors", "draws_skew_t.RData" = "Skew Student-t errors")
# Kernel density estimates (adjust = 1.8, as in the supplement figures), plotted on the
# range between the 0.3% and 99.7% quantiles of the draws of each example
df_post <- do.call(rbind, lapply(names(examples), function(f) {
  load(file.path("output", f))
  rng <- range(sapply(draws, quantile, probs = c(0.003, 0.997)))
  do.call(rbind, lapply(c("fit", "gauss"), function(d) {
    dens <- density(draws[[d]], adjust = 1.8, n = 1024)
    keep <- dens$x >= rng[1] & dens$x <= rng[2]
    data.frame(beta = dens$x[keep], dens = dens$y[keep], distr = d, title = unname(examples[f]))
  }))
}))
df_post$title <- factor(df_post$title, levels = examples)

lab_distr <- c("fit" = expression("Matching" ~ P[n]), "gauss" = expression("Gaussian" ~ P[n]))
pl_post <- ggplot(df_post, aes(x = beta, y = dens, colour = distr, linetype = distr)) +
  facet_wrap(~title, scales = "free") +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(name = "", labels = lab_distr, values = c("fit" = "#fc7d0b", "gauss" = "#1170aa")) +
  scale_linetype_manual(name = "", labels = lab_distr, values = c("fit" = "solid", "gauss" = "dashed")) +
  labs(x = "Regression coefficient", y = "Density") +
  theme_light(base_size = 9) +
  theme(
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "gray82"),
    panel.grid.major = element_line(linewidth = 0.3, colour = "gray93"),
    panel.grid.minor = element_line(linewidth = 0.15, colour = "gray93"),
    axis.text.y = element_blank(),
    legend.position = "top"
  )

# Paper figures are written directly to PAPER/figures when the script is run within the
# paper project, and to output/ otherwise
fig_dir <- if (dir.exists("../../PAPER/figures")) "../../PAPER/figures" else "output"
ggsave(file.path(fig_dir, "pl_ill_post.pdf"), plot = pl_post, width = 6.3, height = 2.6, device = cairo_pdf)
