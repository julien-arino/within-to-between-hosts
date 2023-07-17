library(reshape)
library(ggplot2)

# What state variable we are plotting
to_plot = "V"

# Plot all the values (V1)
df <- data.frame(x = STATE_VARS$time, STATE_VARS[[to_plot]])
# Long format
df <- melt(df, id.vars = "x")
ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "none")
ggsave(filename = "FIGS/viral-load-v1.png",
       width = 20, height = 15, units = "cm")

# Plot all the values (V2)
if (FALSE) {
  df <- data.frame(x = STATE_VARS$time, STATE_VARS[[to_plot]])
  # Long format
  df <- melt(df, id.vars = "x")
  ggplot(df, aes(x = x, y = value, color = variable)) +
    geom_point() +
    geom_smooth() +
    theme(legend.position = "none")
  ggsave(filename = "FIGS/viral-load-v2.png",
         width = 20, height = 15, units = "cm")
}

# Plot the percentiles (2.5, 25, 50, 75 and 97.5%)
if (FALSE) {
  df <- data.frame(x = STATE_VARS$time, SUMMARIES[[to_plot]][,1:5])
  df <- melt(df, id.vars = "x")
  ggplot(df, aes(x = x, y = value, color = variable)) +
    geom_line() +
    theme(legend.position = "none")
  ggsave(filename = "FIGS/summary-viral-load.png",
         width = 20, height = 15, units = "cm")
}

# Plot the mean with error bars
df <- data.frame(time = STATE_VARS$time, SUMMARIES[[to_plot]])
colnames(df) = c("time", "p2d5", "p25","p50","p75","p97d5","mean")
df = df[,c("time", "p25", "mean", "p75")]
df = melt(df, id.vars = "time")
ggplot(df, aes(x = time, y = value, color = variable)) +
  geom_line(aes(linetype = variable, size = variable, color = variable)) +
  scale_linetype_manual(values = c("dashed","solid","dashed")) +
  scale_size_manual(values = c(0.5, 1, 0.5)) +
  scale_color_manual(values = c("grey",
                                "dodgerblue4",
                                "grey20")) +
  theme(legend.position = c(0.8,0.5)) +
  scale_fill_manual(name = "Dose", labels = c("A", "B", "C")) +
  xlim(0,10) +
  xlab("Time (days)") +
  ylab("Viral load")
ggsave(filename = "FIGS/summary-viral-load.png",
       width = 20, height = 15, units = "cm")


