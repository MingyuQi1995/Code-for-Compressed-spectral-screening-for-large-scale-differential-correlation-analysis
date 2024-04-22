library(irlba)
library(Matrix)
library(DGCA)
library(plyr)
library(ROCR)
library(parallel)
library(iterators)
library(doParallel)
registerDoParallel(cores=26)

load("data/Diffgroups.RData")
load("data/functions.rda")


set.seed(163666) 
#set.seed(3438)

B = 60

n1 = 40
n2 = 40
d =  20
g =  200
r =  0.1
n_nv1 = 3
n_nv2 = 3

mu = 0.5
sig = 0.2


data <- foreach(b = 1:B) %dopar% {
  truegroup <- results[[b]][[1]]
  nv1 <- truegroup[1,]
  nv2 <- truegroup[2,]
  dat <- Dat_generation(n1,n2,d,g,r, nv1 = nv1, nv2 = nv2, mu = mu, sig = sig)
  list(dat)
}



save(data, file = "data/datalow.RData")
#save(data, file = "data/data4sd.RData")