library(reshape)
library(ggplot2)

df <- data.frame(x = t, V)
# Long format
df <- melt(df, id.vars = "x")

ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "none")
ggsave(filename = "FIGS/viral-load.png",
       width = 20, height = 15, units = "cm")

df <- data.frame(x = STATE_VARS$time, SUMMARIES$V[,1:4])
df <- melt(df, id.vars = "x")
ggplot(df, aes(x = x, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "none")
ggsave(filename = "FIGS/summary-viral-load.png",
       width = 20, height = 15, units = "cm")

