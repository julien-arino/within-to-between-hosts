library(parallel)
library(reshape)
library(ggplot2)

source("functions_all.R")

# t = COHORT[[1]][,"time"]
# V = mat.or.vec(nr = length(t), nc = N)
# for (i in 1:N) {
#   V[,i] = COHORT[[i]][,"V"]
# }
# 
# df <- data.frame(x = t, V)
# # Long format
# df <- melt(df, id.vars = "x")
# 
# ggplot(df, aes(x = x, y = value, color = variable)) +
#   geom_line() + 
#   theme(legend.position = "none")
# ggsave(filename = "FIGS/viral-load.png",
#        width = 20, height = 15, units = "cm")
