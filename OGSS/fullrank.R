library(irlba)
library(Matrix)
library(DGCA)
library(plyr)
library(ROCR)
library(MASS)
library(Matrix)
library(foreach)
library(doParallel)
library(parallel)
no_cores <- 6
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#set.seed(991726)

load("data/datalow.RData")
load("data/functions.rda")
load("data/Diffgroups.RData")

n1 = 40
n2 = 40
d =  20
g =  200
r =  0.1
n_nv1 = 3
n_nv2 = 3

mu = 0.5
sig = 0.2

nlambda = 25


B = 5
#Roc
N = 200


evalres <- foreach(b = 1:B, .combine = rbind, .packages = c("irlba", "Matrix", "DGCA", "plyr", "ROCR")) %dopar% {
  
     truegroup       <- results[[b]][[1]]
     nv1             <- truegroup[1,]
     nv2             <- truegroup[2,]
     dat <- data[[b]][[1]]
     G               <- dat$G
     p               <- dim(G)[2]
     g               <- dim(G)[1]
     Sigma1          <- dat$Sigma1
     Sigma2          <- dat$Sigma2
     Real_D          <- dat$Real_D
     v_true          <- unique(unlist(dat$ind))
     v_false         <- setdiff(c(1:p), v_true)
     data1           <- dat$data1
     data2           <- dat$data2
     colnames(data1) <- 1:ncol(data1)
     colnames(data2) <- 1:ncol(data2)
     penalty         <- as.vector(sqrt(rowSums(G)))
     D_hat           <- cov(data1) - cov(data2)
    
     k <- 2
     kreal <- k
     khat <- k
     
     sv              <- F
     SVD             <- irlba(A = D_hat, nv = k, maxit = 10000)
     U_hat           <- if (sv)  SVD$u[, 1:k] %*% diag(sqrt( SVD$d[1:k])) else  SVD$u[, 1:k]   
    
     ## generate lambda sequence
     lambda_seq      <- GenLambdaSeq(U_hat, G, nlambda = nlambda, penalty =  penalty, type = "1")
     
     tuning <- T
     
     SR_tuning       <- data.frame(b = 0, sr = 0, method = "")
     GR_tuning       <- data.frame(b = 0, gr = 0, method = "")
     
     if(tuning == T){
     ## evaluate lambda
     U_aug           <- as.matrix(Ufaug(U_hat ,G))
     temp<- calcEvalMetric(data1, data2, D_hat, U_aug, SVD = SVD, G, Params = list(lambda_seq = lambda_seq, k = k, sv = F, cv =  T),  pf = as.vector(sqrt(rowSums(G))))
  
     lambda_sel      <- sel_lambda(temp)
   
     
     for(l in 1:length(lambda_sel)){
       res_aug <- myog_sel(U_aug, G, pf = as.vector(sqrt(rowSums(G))), lambda = lambda_sel[l], eps = 1e-5, type = "1")
       res <- augtoorifast(res_aug, G)
       IoUres <- CalcIoU(res,G, nv1 = nv1, nv2 = nv2, v_true = v_true)
       SR_tuning[l,]  <- c(b = b, IoUres$SR, method = names(lambda_sel)[l])
       GR_tuning[l,]  <- c(b = b, IoUres$GR, method = names(lambda_sel)[l])
     }
     
     }
     
       
     ####report res
     klist <- list(kreal,khat)
  
     list(SR = SR_tuning, GR = GR_tuning)
   
}


save(evalres, file = "res/evalrestest.RData")



stopCluster(cl)

