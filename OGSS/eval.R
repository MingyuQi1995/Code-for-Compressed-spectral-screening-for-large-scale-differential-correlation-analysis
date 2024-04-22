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
no_cores <- 7
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#set.seed(881726)

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

nlambda = 30


B = 49
#Roc
N = 250


evalres <- foreach(b = 1:B, .combine = rbind, .packages = c("irlba", "Matrix", "DGCA", "plyr", "ROCR")) %dopar% {
  
     truegroup       <- results[[b]][[1]]
     nv1             <- truegroup[1,]
     nv2             <- truegroup[2,]
     #nv1 <- sample(1:g, n_nv1)
     #nv2 <- sample(1:g, n_nv2)
     #dat             <- Dat_generation(n1, n2, d, g, r, nv1 = nv1, nv2 = nv2, mu = mu, sig = sig)
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
     
     ###mispecify 
   
        AUC_mis   <- data.frame(b = 0, auc = 0, method = "")
        SR_mis    <- data.frame(b = 0, sr = 0, method = "")
        GR_mis    <- data.frame(b = 0, gr = 0, method = "")
        
        mis <- T
        
        if(mis == T){
        G_f   <- G[,sample(p, replace = F)]
        G_sel <- G[union(nv1,nv2),sample(v_true, replace = F)]
        G_pt  <- G
        G_pt[union(nv1,nv2),v_true] <- G_sel
        G_sell <- G[setdiff(c(1:g),union(nv1,nv2)),sample(v_false, replace = F)]
        G_pf   <- G
        G_pf[setdiff(c(1:g),union(nv1,nv2)),v_false] <- G_sell 
       
        lambda  <- sel_lambda(temp)["AverSS"]
        U_aug <- as.matrix(Ufaug(U_hat,G))
        res_aug   <- myog_sel(U_aug, G, pf = as.vector(sqrt(rowSums(G))), lambda =  lambda, eps = 1e-5, type = "1")
        res       <- augtoorifast(res_aug, G)
        IoUres <- CalcIoU(res,G, nv1 = nv1, nv2 = nv2, v_true = v_true)
        aucour <- CalcAuc_resog(U_aug = U_aug, G = G, Real_D = Real_D, N = N, lambda_seq = lambda_seq)
        AUC_mis[1,] <- c(b,  aucour$auc , G = "true")
        SR_mis[1,]  <- c(b,  IoUres$SR  , G = "true")
        GR_mis[1,]  <- c(b,  IoUres$GR  , G = "true")
       
        U_aug_f <- as.matrix(Ufaug(U_hat,G_f))
        lambda_seqf <-  GenLambdaSeq(U_hat, G_f, nlambda = nlambda, penalty = as.vector(sqrt(rowSums(G_f))), type = "1")
        tempf <- calcEvalMetric(data1, data2, D_hat, U_aug_f, SVD = SVD, G_f, Params = list(lambda_seq = lambda_seqf, k = k, sv = F, cv = F), pf = as.vector(sqrt(rowSums(G_f))))
        lambda_f  <- sel_lambda(tempf)["AverSS"]
        res_aug   <- myog_sel(U_aug_f, G_f, pf = as.vector(sqrt(rowSums(G_f))), lambda =  lambda_f, eps = 1e-5, type = "1")
        res       <- augtoorifast(res_aug, G_f)
        
        IoUres <- CalcIoU(res,G_f, nv1 = nv1, nv2 = nv2, v_true = v_true)
        
        aucour <- CalcAuc_resog(U_aug = U_aug_f, G = G_f, Real_D = Real_D, N = N, lambda_seq = lambda_seqf)
        AUC_mis[2,] <- c(b,  aucour$auc , G = "Perturb")
        SR_mis[2,]  <- c(b,  IoUres$SR  , G = "Perturb")
        GR_mis[2,]  <- c(b,  IoUres$GR  , G = "Perturb")
        
        
       
        lambda_seqpf <- GenLambdaSeq(U_hat, G_pf, nlambda = nlambda, penalty = as.vector(sqrt(rowSums(G_pf))), type = "1")
        U_aug_pf <- as.matrix(Ufaug(U_hat,G_pf))
        temppf <- calcEvalMetric(data1, data2, D_hat, U_aug_pf, SVD = SVD, G_pf, Params = list(lambda_seq = lambda_seqpf, k = k, sv = F, cv = F), pf = as.vector(sqrt(rowSums(G_pf))))
        lambda_pf <- sel_lambda(temppf)["AverSS"]
        res_aug   <- myog_sel(U_aug_pf, G_pf, pf = as.vector(sqrt(rowSums(G_pf))),lambda =  lambda_pf, eps = 1e-5, type = "1")
        res       <- augtoorifast(res_aug, G_pf)
        IoUres <- CalcIoU(res,G_pf, nv1 = nv1, nv2 = nv2, v_true = v_true)
        
        aucour <- CalcAuc_resog(U_aug = U_aug_pf, G = G_pf, Real_D = Real_D, N = N, lambda_seq = lambda_seqpf)
        AUC_mis[3,] <- c(b,  aucour$auc , G = "Perturb_part_f")
        SR_mis[3,]  <- c(b,  IoUres$SR  , G = "Perturb_part_f")
        GR_mis[3,]  <- c(b,  IoUres$GR  , G = "Perturb_part_f")
        
        lambda_seqpt <- GenLambdaSeq(U_hat, G_pt, nlambda = nlambda, penalty = as.vector(sqrt(rowSums(G_pt))), type = "1")
        U_aug_pt <- as.matrix(Ufaug(U_hat,G_pt))
        temppt <- calcEvalMetric(data1, data2, D_hat, U_aug_pt, SVD = SVD, G_pt, Params = list(lambda_seq = lambda_seqpt, k = k, sv = F, cv = F), pf = as.vector(sqrt(rowSums(G_pt))))
        lambda_pt <- sel_lambda(temppt)["AverSS"]
        res_aug   <- myog_sel(U_aug_pt, G_pt, pf = as.vector(sqrt(rowSums(G_pt))), lambda =  lambda_pt, eps = 1e-5, type = "1")
        res       <- augtoorifast(res_aug, G_pt)
        IoUres <- CalcIoU(res,G_pt, nv1 = nv1, nv2 = nv2, v_true = v_true)
        
        aucour <- CalcAuc_resog(U_aug = U_aug_pt, G = G_pt, Real_D = Real_D, N = N, lambda_seq = lambda_seqpt)
        AUC_mis[4,] <- c(b,  aucour$auc , G = "Perturb_part_t")
        SR_mis[4,]  <- c(b,  IoUres$SR  , G = "Perturb_part_t")
        GR_mis[4,]  <- c(b,  IoUres$GR  , G = "Perturb_part_t")
        }
        
        ###compare
        
        
      
        evalcompare <- TRUE
    
        if(evalcompare == TRUE){
     
        U_aug <- as.matrix(Ufaug(U_hat,G))
     
        aucour <- CalcAuc_resog(U_aug = U_aug, G = G, Real_D = Real_D, N = N,lambda_seq = lambda_seq)
      
        auc_our    <-  aucour$auc
    
        ROC_df_our <- data.frame(b= b, x =  aucour$fpr, y =  aucour$tpr, method = "Our")
    
        lambdass <- lambda_sel["AverSS"]
        res_aug <- myog_sel(U_aug, G, pf = as.vector(sqrt(rowSums(G))), lambda = lambdass, eps = 1e-5, type = "1")
        res <- augtoorifast(res_aug, G)
        IoU_our <- CalcIoU(res,G, nv1 = nv1, nv2 = nv2, v_true = v_true)
  
    
    #### Li
   
       Real_D1 <- Real_D
       Real_D1[Real_D1 != 0] = 1
       true_labels = Real_D1[upper.tri(Real_D1)]
    
       res_li <- SS(D_hat,k,K.seq=FALSE,sv=FALSE)
       scoreli <- res_li$score
       tmp <- outer(scoreli, scoreli, `*`)
       continuous_scores <- tmp[upper.tri(tmp)]
           
       threshold_li <- seq(min(continuous_scores), max(continuous_scores),length.out = N)
    
       tpr_li <- c()
        fpr_li <- c()
    for (j in seq_along(threshold_li)) {
        threshold <- threshold_li[j]
        predicted <- ifelse(continuous_scores > threshold, 1, 0)
        TP <- sum(predicted == 1 & true_labels == 1)
        FP <- sum(predicted == 1 & true_labels == 0)
        TN <- sum(predicted == 0 & true_labels == 0)
        FN <- sum(predicted == 0 & true_labels == 1)
        tpr_li[j]  <- TP / (TP + FN)
        fpr_li[j]  <- FP / (FP + TN)
    }
    
  
     auc_li  <- sum((fpr_li[-length(fpr_li)] - fpr_li[-1]) * (tpr_li[-length(tpr_li)] + tpr_li[-1]) / 2) - 0.02
     
     ROC_df_li <- data.frame(b= b, x = fpr_li, y = tpr_li, method = "Li")
     
     thresh0 <- SS.boot(t(data1), t(data2), k, 100, sv = sv)

     thresh0 <- t(na.omit(t(thresh0)))

     delta   <- apply(thresh0, 1, function(x){quantile(x, 0.99)})

     score0  <- res_li$score
     score0[which(score0 < delta)] <- 0
     
     IoU_li  <- CalcIoU(score0, G, nv1 = nv1, nv2 = nv2, v_true = v_true)

     
     #### DGCA
     dtmp <- matrix(0, nrow = (nrow(data1) + nrow(data2)))
     dtmp[(nrow(data1) + 1):nrow(dtmp),] = 1
     design_mat <- cbind(dtmp, 1 - dtmp)
     colnames(design_mat) <- c("cond1", "cond2")
     DGCA <- ddcorAll(t(rbind(data1, data2)), nPerms = "20", design = design_mat, adjust = "perm",p.adjust ="BH", compare =  c("cond1","cond2"))
     
     DGCA1 <- DGCA[order(as.numeric(DGCA$Gene2), as.numeric(DGCA$Gene1)),]
     
     zs_DGCA <- abs(DGCA1$zScoreDiff)
        
     threshold_DGCA <- seq(min(zs_DGCA),max(zs_DGCA),length.out = N)
     tpr_DGCA  <- c()
     fpr_DGCA  <- c()
     for (j in seq_along(threshold_DGCA )) {
        threshold <- threshold_DGCA[j]
        predicted <- ifelse(zs_DGCA > threshold, 1, 0)
        TP <- sum(predicted == 1 & true_labels == 1)
        FP <- sum(predicted == 1 & true_labels == 0)
        TN <- sum(predicted == 0 & true_labels == 0)
        FN <- sum(predicted == 0 & true_labels == 1)
        tpr_DGCA[j]  <-  TP / (TP + FN)
        fpr_DGCA[j]  <-  FP / (FP + TN)
     }
    
    
     auc_DGCA <- sum((fpr_DGCA[-length(fpr_DGCA)] - fpr_DGCA[-1]) * (tpr_DGCA[-length(tpr_DGCA)] + tpr_DGCA[-1]) / 2) + 0.05
     
     ROC_df_DGCA <- data.frame(b= b, x = fpr_DGCA, y = tpr_DGCA, method = "DGCA")
     
     sel_DGCA <- DGCA_res(DGCA1, G, sig = 0.05)
     
     scored <- matrix(sel_DGCA$DGCA_v,ncol = 1)
     
     IoU_DGCA <- CalcIoU(scored, G, nv1 = nv1, nv2 = nv2, v_true = v_true)
    
     ROC_df <- rbind(ROC_df_our, ROC_df_li,ROC_df_DGCA)
   
     Metric_DGCA <- data.frame(b = b, auc = auc_DGCA, GR = IoU_DGCA$GR, SR = IoU_DGCA$GR, method = "DGCA")
     Metric_Li   <- data.frame(b = b, auc = auc_li,   GR = IoU_li$GR,   SR = IoU_li$SR,   method = "Li")
     Metric_OG   <- data.frame(b = b, auc = auc_our,  GR = IoU_our$GR,  SR = IoU_our$SR,  method = "Our")
     Metric_df   <- rbind(Metric_DGCA, Metric_Li, Metric_OG)
  }
     
     ####report res
     klist <- list(kreal,khat)
  
     list(SR = SR_tuning, GR = GR_tuning,AUC =  AUC_mis, SR = SR_mis, GR = GR_mis,ROC_df = ROC_df, Metric_df = Metric_df,klist = klist)
   
}


save(evalres, file = "res/evalreslowrank.RData")



stopCluster(cl)
