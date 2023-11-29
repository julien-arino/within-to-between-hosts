library(ggplot2)
library(latex2exp)
# For a little bit of fun
# library(xkcd)
# library(extrafont)
# download.file("http://simonsoftware.se/other/xkcd.ttf",
#              dest="xkcd.ttf", mode="wb")

source("functions_all.R")

# Plots of sensitivity values computed using the sensitivity package
vals = value_indicators(params_change = pars.sobol, params_fixed = params_all)
# Compute the correlation between the indicators and the parameters
writeLines("Computing partial rank correlations - V_max")
tmp = as.numeric(vals$V_max)
x_V_max= pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
             rank = TRUE, semi = FALSE)
ggplot(x_V_max)
writeLines("Computing partial rank correlations - F_U_max")
tmp = as.numeric(vals$F_U_max)
x_F_U_max = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                rank = TRUE, semi = FALSE)
ggplot(x_F_U_max)
writeLines("Computing partial rank correlations - F_B_max")
tmp = as.numeric(vals$F_B_max)
x_F_B_max = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                rank = TRUE, semi = FALSE)
ggplot(x_F_B_max)
writeLines("Computing partial rank correlations - T_max_V")
tmp = as.numeric(vals$T_max_V)
x_T_max_V = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                rank = TRUE, semi = FALSE)
ggplot(x_T_max_V)
writeLines("Computing partial rank correlations - T_max_F_U")
tmp = as.numeric(vals$T_max_F_U)
x_T_max_F_U = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                  rank = TRUE, semi = FALSE)
ggplot(x_T_max_F_U)
writeLines("Computing partial rank correlations - T_max_F_B")
tmp = as.numeric(vals$T_max_F_B)
x_T_max_F_B = pcc(pars.sobol[,2:dim(pars.sobol)[2]], tmp,
                  rank = TRUE, semi = FALSE)
ggplot(x_T_max_F_B)
# Combine values for aggregate plots. All values, then by type
PRCC = data.frame(
  names = names(pars.list),
  V_max = x_V_max$PRCC,
  F_U_max = x_F_U_max$PRCC,
  F_B_max = x_F_B_max$PRCC,
  T_max_V = x_T_max_V$PRCC,
  T_max_F_U = x_T_max_F_U$PRCC,
  T_max_F_B = x_T_max_F_B$PRCC
)
colnames(PRCC) = c("names", "V_max", "F_U_max", "F_B_max", "T_max_V", "T_max_F_U", "T_max_F_B")
PRCC_values = PRCC[,c("names", "V_max", "F_U_max", "F_B_max")]
PRCC_times = PRCC[,c("names", "T_max_V", "T_max_F_U", "T_max_F_B")]
# Order the parameters by decreasing absolute value of PRCC
order_PRCC_values = order(apply(abs(PRCC_values[2:4]), 1, max), 
                          decreasing = TRUE)
order_PRCC_times = order(apply(abs(PRCC_times[2:4]), 1, max), 
                         decreasing = TRUE)
PRCC_values = PRCC_values[order_PRCC_values,]
PRCC_times = PRCC_times[order_PRCC_times,]
# Reshape for plotting 
PRCC_values = reshape2::melt(PRCC_values, id.vars = "names")
PRCC_times = reshape2::melt(PRCC_times, id.vars = "names")
# Needed for x axis labels
x_labels = c(TeX("$\\beta$"), 
             TeX("$d_D$"), 
             TeX("$d_I$"), 
             TeX("$d_V$"), 
             TeX("$\\delta$"), 
             TeX("$\\epsilon_{FI}$"), 
             TeX("$\\eta_{FI}$"), 
             TeX("$k_{B_F}$"), 
             TeX("$k_{int_f}$"), 
             TeX("$k_{lin_f}$"), 
             TeX("$k_{U_f}$"), 
             TeX("$\\lambda_{S}$"),
             TeX("$p$"),
             TeX("$p_{FI}$"),
             TeX("$\\psi_{F_{prod}}$"),
             TeX("$S_{max}$"),
             TeX("$\\tau_{I}$"),
             TeX("$T^*$"))
# Do some scaling for the size of points
scaled_size_values = PRCC_values$value
max_value_values = max(abs(PRCC_values$value))
half_saturation_values = 0.01
for (i in 1:dim(PRCC_values)[1]) {
  if (PRCC_values$value[i] < 0) {
    scaled_size_values[i] = 
      PRCC_values$value[i]/(-half_saturation_values-max_value_values)
  } else {
    scaled_size_values[i] = 
      PRCC_values$value[i]/(half_saturation_values+max_value_values)
  }
}
scaled_size_times = PRCC_times$value
max_value_times = max(abs(PRCC_times$value))
half_saturation_times = 0.01
for (i in 1:dim(PRCC_times)[1]) {
  if (PRCC_times$value[i] < 0) {
    scaled_size_times[i] = 
      PRCC_times$value[i]/(-half_saturation_times-max_value_times)
  } else {
    scaled_size_times[i] = 
      PRCC_times$value[i]/(half_saturation_times+max_value_times)
  }
}
# Plot the values
ggplot(data = PRCC_values, aes(names, value, col = variable)) + 
  geom_point(size = scaled_size_values*12) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  ylim(-1,1) +
  scale_color_discrete(labels = c(TeX("$V_{max}$"), 
                                  TeX("$F_{U_{max}}$"), 
                                  TeX("$F_{B_{max}}$"))) +
  labs(color="Indicator") +
  scale_x_discrete(labels = x_labels[order_PRCC_values],
                   limits = names(pars.list)[order_PRCC_values]) +
  theme(legend.position = c(0.9,0.9)) +
  #theme_xkcd()
  theme_minimal()
ggsave(filename = "FIGS/sensitivities-PRCC-values.png",
       width = 20, height = 15, units = "cm")
# Plot the times
ggplot(data = PRCC_times, aes(names, value, col = variable)) + 
  geom_point(size = scaled_size_times*12) +
  xlab("Parameter") +
  ylab("Partial rank correlation coefficients") +
  ylim(-1,1) +
  scale_color_discrete(labels = c(TeX("$t_{V_{max}}$"), 
                                  TeX("$t_{F_{U_{max}}}$"), 
                                  TeX("$t_{F_{B_{max}}}$"))) +
  labs(color="Indicator") +
  scale_x_discrete(labels = x_labels[order_PRCC_times],
                   limits = names(pars.list)[order_PRCC_times]) +
  theme(legend.position = c(0.9,0.9)) +
  theme_minimal()
ggsave(filename = "FIGS/sensitivities-PRCC-times.png",
       width = 20, height = 15, units = "cm")

