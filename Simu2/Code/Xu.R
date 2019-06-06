rm(list = ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu2/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"

library(mvtnorm)
library(msm)
library(expm)
library(pROC)
library(MBSGS)
set.seed(20160125)
ns<-100
AUC<-group_FDR<-group_FOR<-feature_FDR<-feature_FOR<-MSE<-rep(NA,ns)
t0<-proc.time()
for(s in 1:ns){
  # simulate data
  seed <- 2016+100000*s
  source(paste0(code_folder, "R/GeneSimu2_v2.R"))
  data <- geneSimu2(seed)
  X <- data$X
  Y <- data$Y
  U1 <- data$U1
  types <- data$types
  true_beta <- data$true_beta

  # split training and testing
  train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
  X_train <- X[train_id, ]
  Y_train <- Y[train_id]
  X_test <- X[-train_id, ]
  Y_test <- Y[-train_id]

  # duplicate features
  U1_colsum <- apply(U1, 2, sum)
  feature_dup_index <- unlist(lapply(1:ncol(U1), 
    function(g) which(U1[, g]==1)))
  X_train_dup <- X_train[, feature_dup_index]
  X_test_dup <- X_test[, feature_dup_index]
  g_index <- unlist(sapply(1:ncol(U1), function(g) 
    rep(g, U1_colsum[g])))

  # MCMC
  start_nsim<-2001
  nsim<-3000

  res <- BSGSSS(Y_train, X_train_dup, niter=nsim, 
    burnin=nsim - start_nsim +1, group_size=U1_colsum, 
    num_update = 100, niter.update = 100)
 
  #####################
  # Check results
  #####################
  
  ##### feature level selection
  BETA_dup <- res$coef
  BETA <-  sapply(1:nrow(U1), function(x) 
    apply(BETA_dup[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
  BETA <- t(BETA)
  ALL_IND <- BETA != 0

  f_ind <- apply(ALL_IND,1,mean)  
  f_null <- 1-f_ind
  t <- seq(0,1,by=0.01)
  f_qvalue <- rep(NA,length(f_ind))
  for(j in 1:length(f_ind)){
    t_g<-t[t>=f_null[j]]
    f_qvalue[j]<-min(sapply(1:length(t_g), 
      function(x)sum(f_null[f_null<=t_g[x]])/sum(f_null<=t_g[x])))
  }
  feature_FDR[s]<-sum(true_beta==0 & f_qvalue<=0.1)/sum(f_qvalue<=0.1)
  feature_FOR[s]<-sum(true_beta!=0 & f_qvalue>0.1)/sum(f_qvalue>0.1)
  
  #### beta ####
  beta_median<-apply(BETA,1,median)
  MSE[s]<-mean((X_test %*% beta_median - Y_test)^2)

  ##### compute AUC for features
  AUC[s]<-pROC::auc(roc(factor(true_beta!=0), f_ind))

  cat(paste("repeat=",s,"\n",sep=""))
}  
t1<-proc.time()-t0  #9 hours
save.image("Results/Xu.RData")


