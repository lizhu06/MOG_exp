rm(list = ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu3/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"

library(MASS)
library(GIGrvg)
library(pROC)

C<-10
G<-100

U_seq<-c(0.2,0.5)
ns<-100
AUC<-feature_FDR<-feature_FOR<-MSE<-matrix(NA,length(U_seq), ns)

t0<-proc.time()
for(du in 1:length(U_seq)){
  for(s in 1:ns){
    U<-U_seq[du]

    seed <- 2016+100000*s+1000*du
   
    source(paste0(code_folder, "R/GeneSimu3_v2.R"))
    data <- geneSimu3(U=U_seq[du], seed)
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
    U1_rowsum <- apply(U1, 1, sum)
    feature_dup_index <- unlist(lapply(1:nrow(U1), 
      function(j) rep(j, U1_rowsum[j])))
    X_train_dup <- X_train[, feature_dup_index]
    X_test_dup <- X_test[, feature_dup_index]
    g_index <- unlist(sapply(1:nrow(U1), function(j) which(U1[j, ]==1)))

    # MCMC
    start_nsim<-2001
    nsim<-3000

    source("/mnt/glusterfs/liz86/MOG_Regression/R/HSVS/HSVS_modify.R")
    #t0<-proc.time() 
    res <- HSVS(burnin=(start_nsim-1),  N=(nsim - start_nsim + 1),  
      Y_train,  X_train_dup,  rep(seq(1,100), each=3))  

    #####################
    # Check results
    #####################
    BETA_dup <- res[[1]]
    BETA <-  sapply(1:nrow(U1), function(x) 
      apply(BETA_dup[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
    BETA <- t(BETA)

    get_95ci <- function(x)quantile(x,probs=c(0.025,0.975))
    beta_CI <- t(apply(BETA,1,get_95ci))
    beta_ind <- 1-(beta_CI[,1]<=0 & beta_CI[,2]>=0)*1
    feature_FDR[du, s] <- sum(true_beta==0 & beta_ind==1)/sum(beta_ind==1)
    feature_FOR[du, s] <- sum(true_beta!=0 & beta_ind==0)/sum(beta_ind==0)

    #### MSE
    beta_median <- apply(BETA, 1, median)
    MSE[du, s]<-mean((X_test %*% beta_median - Y_test)^2)
    
    ##### compute AUC for features
    PM <- (BETA!=0)*1
    pm<-apply(PM,1,mean)
    pos_frac<-apply(BETA>0,1,mean)
    neg_frac<-apply(BETA<0,1,mean)
    max_frac<-sapply(1:length(pos_frac),function(x) 
      max(pos_frac[x],neg_frac[x]))
    
    AUC[du, s]<-pROC::auc(roc(factor((true_beta!=0)*1), max_frac))
    
    cat(paste("repeat=",s,"\n",sep=""))

  }
}
t1<-proc.time()-t0  # 8 hours

save.image(file="Results/Veera.RData")
