# library(reshape)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(scales)

# Function that sets up the data frame for plotting (changes names, formats, etc.)
# By default, prepares data with mean and 2.5 and 97.5 percentiles. Change as needed.
format_df = function(time,
                     data,
                     line_plotted = "mean", 
                     lower = "2.5%",
                     upper = "97.5%") {
  df = data.frame(time = time,
                  lower = data[,lower],
                  line = data[,line_plotted],
                  upper = data[,upper])
  return(df)
}

# Set equal time range for all figures
time_range = c(0,80)
# Legend position for small figures
pos_legend = c(0.7,0.6)

# Plot the damage
to_plot = "Phi"
df <- format_df(time = STATE_VARS$time, data = SUMMARIES[[to_plot]])
p_damage = ggplot(df) +
  annotate('rect', 
           xmin = 0, 
           xmax = time_range[2], 
           ymin = -85, ymax = -75, 
           alpha=.2, fill='orange') + 
  annotate('text',
           x = time_range[2]/2, y = - 80, label = "ICU patients") +
  annotate('rect', 
           xmin = 0, 
           xmax = time_range[2], 
           ymin = -100, ymax = -85, 
           alpha=.2, fill='red') + 
  annotate('text',
           x = time_range[2]/2, y = - 92.5, label = "Dead patients") +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="grey", alpha = 0.5) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  xlim(time_range) +
  xlab("Time (days)") +
  ylab(TeX("Percentage $\\Psi$ of tissue damage")) +
  ggtitle("Lung tissue damage") +
  theme(plot.title = element_text(hjust = 0.5))
  
# Plot the viral load
to_plot = "V"
df <- format_df(time = STATE_VARS$time, data = SUMMARIES[[to_plot]])
p_viral_load = ggplot(df) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="grey", alpha=0.5) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  xlim(time_range) +
  labs(x = "") +
  ylab(TeX("$\\log_{10}$(copies/ml)")) +
  ggtitle("Viral load") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot the number of infected calls
to_plot = "I"
df <- format_df(time = STATE_VARS$time, data = SUMMARIES[[to_plot]])
p_infected = ggplot(df) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="grey", alpha=0.5) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  xlim(time_range) +
  labs(x = "") +
  ylab(TeX("$10^8$ cell/ml")) +
  ggtitle("Infected cells") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot the bound IFN
# (Function from https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales)
scientific_10 = function(x) {
  ifelse(
    x==0, "0",
    parse(text = sub("e[+]?", " %*% 10^ ", scientific_format()(x)))
  )
} 
to_plot = "F_B"
df <- format_df(time = STATE_VARS$time, data = SUMMARIES[[to_plot]])
p_bound_IF = ggplot(df) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="grey", alpha=0.5) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  xlim(time_range) +
  xlab("Time (days)") +
  ylab("pg/ml") +
  scale_y_continuous(label = scientific_10) +
  ggtitle("Bound IFN") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot the unbound IFN
to_plot = "F_U"
df <- format_df(time = STATE_VARS$time, data = SUMMARIES[[to_plot]])
p_unbound_IF = ggplot(df) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill="grey", alpha=0.5) +
  geom_line(aes(x = time, y = line), colour = "dodgerblue4", size = 1.5) +
  xlim(time_range) +
  xlab("Time (days)") +
  ylab("pg/ml") +
  ggtitle("Unbound IFN") +
  theme(plot.title = element_text(hjust = 0.5))

# Set up figure
pp = ggarrange(p_damage,
               ggarrange(p_viral_load, p_infected, 
                         p_bound_IF, p_unbound_IF,
                         nrow = 2, ncol =2),
               nrow = 1, ncol = 2,
               widths = c(1, 1.5))
pp
ggsave(filename = "FIGS/figure-3-state-and-damage.png",
       width = 20, height = 15, units = "cm")


