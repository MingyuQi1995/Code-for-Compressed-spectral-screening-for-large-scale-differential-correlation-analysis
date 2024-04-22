library(irlba)
library(Matrix)
library(DGCA)
library(plyr)
library(ROCR)
library(parallel)
library(iterators)
library(doParallel)

# Number of cores to use
numCores <- 50

# Create a cluster object with the specified number of cores
cl <- makeCluster(numCores)
registerDoParallel(cl)


set.seed(1111)

B = 100

g = 200

n_nv1 = 3
n_nv2 = 3

results <- foreach(b = 1:B) %dopar% {
  nv1 <- sample(1:g, n_nv1)
  nv2 <- sample(1:g, n_nv2)
  A <- rbind(nv1, nv2)
  list(A) 
}


save(results, file = "data/Diffgroups.RData")

# After exiting the outer loop, stop the parallel cluster
stopCluster(cl)