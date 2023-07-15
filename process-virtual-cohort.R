library(parallel)

source("functions_all.R")

t = COHORT[[1]][,"time"]
V = mat.or.vec(nr = length(t), nc = N)
for (i in 1:N) {
  V[,i] = COHORT[[i]][,"V"]
}

