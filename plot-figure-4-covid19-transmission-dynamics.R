# library(reshape)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(scales)

# Set equal time range for all figures
time_range = c(0,30)
# Legend position for small figures
pos_legend = c(0.7,0.6)

# Plot beta_hat versus time. 
# For the regions, in order to make them appear with
# the same colour as the curve, we make the lower limit for
# one class be the upper limit for the next, so the regions
# don't overlap in practice. (In reality, they do overlap, 
# of course, but that means the colours mix.)
df_low <- format_df(time = STATE_VARS$time, data = SUMMARIES$beta_hat_outcome[[1]])
df_medium <- format_df(time = STATE_VARS$time, data = SUMMARIES$beta_hat_outcome[[2]])
df_high <- format_df(time = STATE_VARS$time, data = SUMMARIES$beta_hat_outcome[[3]])
df_medium$lower = df_low$upper
df_high$lower = df_medium$upper
# Now plot
p_beta_vs_time = ggplot(df_low) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="dodgerblue4", alpha = 0.1) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  geom_ribbon(data = df_medium, aes(x = time, ymin = lower, ymax = upper), fill="darkolivegreen4", alpha = 0.1) +
  geom_line(data = df_medium, aes(x = time, y = line), colour = "darkolivegreen4", size = 1.5) +
  geom_line(data = df_high, aes(x= time, y = line), colour = "firebrick4", size = 1.5) +
  geom_ribbon(data = df_high, aes(x = time, ymin = lower, ymax = upper), fill="firebrick4", alpha = 0.1) +
  xlim(time_range) +
  xlab("Time (days)") +
  ylab(TeX("$\\hat{\\beta}$ (day$^{-1}$)")) +
  ggtitle("Transmissibility") +
  theme(plot.title = element_text(hjust = 0.5))
  
# Boxplot of max transmissibility
max_beta_hat = data.frame(Status = c(), value = c())
for (s in 1:3) {
  max_beta_hat = 
    rbind(max_beta_hat,
          data.frame(Status = rep(s, length(SUMMARIES$beta_hat_max_value_outcome[[s]])[1]),
                     value = SUMMARIES$beta_hat_max_value_outcome[[s]]))
}
max_beta_hat$Status[which(max_beta_hat$Status == 1)] = "Mild"
max_beta_hat$Status[which(max_beta_hat$Status == 2)] = "ICU"
max_beta_hat$Status[which(max_beta_hat$Status == 3)] = "Dead"
max_beta_hat$Status <- 
  factor(max_beta_hat$Status , 
         levels=c("Mild", "ICU", "Dead"))
max_beta_hat$Status = as.factor(max_beta_hat$Status)
p_boxplot = ggplot(max_beta_hat, 
                   aes(x = Status, 
                       y = value, 
                       color = Status)) +
  geom_violin() + 
  scale_color_manual(values = c("dodgerblue4", 
                       "darkolivegreen4",
                       "firebrick4")) + 
  theme(legend.position = "none") +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  stat_summary(fun=median, geom="point", shape=22, size=4) +
  ylab(TeX("$\\hat{\\beta}$ (day$^{-1}$)")) +
  ggtitle("Peak transmissibility") +
  theme(plot.title = element_text(hjust = 0.5))

# Boxplot of time of max transmissibility
time_max_beta_hat = data.frame(Status = c(), value = c())
for (s in 1:3) {
  time_max_beta_hat = 
    rbind(time_max_beta_hat,
          data.frame(Status = rep(s, length(SUMMARIES$beta_hat_max_time_outcome[[s]])[1]),
                     value = SUMMARIES$beta_hat_max_time_outcome[[s]]))
}
time_max_beta_hat$Status[which(time_max_beta_hat$Status == 1)] = "Mild"
time_max_beta_hat$Status[which(time_max_beta_hat$Status == 2)] = "ICU"
time_max_beta_hat$Status[which(time_max_beta_hat$Status == 3)] = "Dead"
#time_max_beta_hat$Status = as.factor(time_max_beta_hat$Status)
time_max_beta_hat$Status <- 
  factor(time_max_beta_hat$Status , 
         levels=c("Mild", "ICU", "Dead"))
p_boxplot_time = ggplot(time_max_beta_hat, 
                        aes(x = Status, 
                            y = value, 
                            color = Status)) +
  geom_violin() + 
  scale_color_manual(values=c("dodgerblue4", 
                              "darkolivegreen4",
                              "firebrick4")) + 
  theme(legend.position = "none") +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  stat_summary(fun=median, geom="point", shape=22, size=4) +
  ylab("Time (days)") +
  ggtitle("Time peak transmissibility") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(0,30)


# Set up figure
pp = ggarrange(ggarrange(p_beta_vs_time, p_boxplot,
                         p_boxplot_time,
                         nrow = 1, ncol =3),
               # ggarrange(p_beta_vs_time, p_boxplot,
               #           nrow = 1, ncol =2),
               nrow = 1, ncol = 1,
               widths = c(1, 1.5))
pp
ggsave(filename = "FIGS/figure-4-covid-19-transmission.png",
       width = 20, height = 15, units = "cm")


