library(reshape)
library(ggplot2)
library(latex2exp)


# Variable we are plotting
to_plot = "V"

# Plot all the values of the viral load
if (FALSE) {
  df <- data.frame(x = STATE_VARS$time, STATE_VARS[[to_plot]])
  # Long format
  df <- melt(df, id.vars = "x")
  ggplot(df, aes(x = x, y = value, color = variable)) +
    geom_line() +
    xlab("Time (days)") +
    ylab("Viral load") +
    theme(legend.position = "none")
  ggsave(filename = "FIGS/viral-load.png",
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
  xlim(0,60) +
  xlab("Time (days)") +
  ylab("Viral load")
ggsave(filename = "FIGS/summary-viral-load.png",
       width = 20, height = 15, units = "cm")


# Variable we are plotting
to_plot = "Phi"

# Plot all the values (V1)
if (FALSE) {
  df <- data.frame(x = STATE_VARS$time, STATE_VARS[[to_plot]])
  # Long format
  df <- melt(df, id.vars = "x")
  # Add category for type of outcome
  tmp = as.numeric(gsub("X", "", df$variable))
  df$outcome = SUMMARIES$lung_loss$outcome[tmp]
  # Now plot..
  ggplot(df, aes(x = x, y = value, color = outcome)) +
    geom_line() +
    xlab("Time (days)") +
    ylab(TeX("$\\Phi$")) +
  theme(legend.position = "none")
  ggsave(filename = "FIGS/lung-damage.png",
         width = 20, height = 15, units = "cm")
}

# Plot the mean with error bars
df <- data.frame(time = STATE_VARS$time, SUMMARIES[[to_plot]])
colnames(df) = c("time", "p2d5", "p25","p50","p75","p97d5","mean")
df = df[,c("time", "p2d5", "mean", "p97d5")]
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
  xlim(0,80) +
  xlab("Time (days)") +
  ylab(TeX("$\\Phi$"))
ggsave(filename = "FIGS/summary-lung-damage.png",
       width = 20, height = 15, units = "cm")

# Plot time of death and hospitalisation info
df = data.frame(value = SUMMARIES$time_of_death)
ggplot(df, aes(x = value)) + 
  geom_histogram(binwidth=.25) +
  xlab("Time (days)") +
  ylab("# deaths")
ggsave(filename = "FIGS/deaths.png",
       width = 20, height = 15, units = "cm")
df = data.frame(value = SUMMARIES$time_of_hospitalisation)
ggplot(df, aes(x = value)) + 
  geom_histogram(binwidth=.25) +
  xlab("Time (days)") +
  ylab("# hospitalisations without death")
ggsave(filename = "FIGS/hospitalisations.png",
       width = 20, height = 15, units = "cm")

# Plot both distributions
df1 = data.frame(value = SUMMARIES$time_of_death,
                 Events = rep("Deaths", 
                            length(SUMMARIES$time_of_death)))
df2 = data.frame(value = SUMMARIES$time_of_recovery,
                 Events = rep("Recoveries", 
                            length(SUMMARIES$time_of_recovery)))
df = rbind(df1, df2)
df.mean <- plyr::ddply(df, "Events", 
                       plyr::summarise, 
                       rating.mean=mean(value))
ggplot(df, aes(x = value, fill = Events)) + 
  geom_density(alpha=.3)+
  xlim(0,30) +
  xlab("Time (days)") +
  ylab("Density") +
#  geom_vline(data=df, aes(xintercept=df.mean,  colour=Event),
#             linetype="dashed", size=1) +
  theme(legend.position = c(0.8,0.5))
ggsave(filename = "FIGS/deaths-recoveries.png",
       width = 20, height = 15, units = "cm")
