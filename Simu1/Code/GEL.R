rm(list = ls())
library(pROC)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu1/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"
library(grpregOverlap)
res_name <- "GEL"

ns <- 100
AUC <- feature_FDR <- feature_FOR <- MSE <- AUC <- rep(NA,ns)
t0 <- proc.time()
for(s in 1:ns){

  seed <- 2016+100000*s
  source(paste0(code_folder, "R/GeneSimu1.R"))

  data <- geneSimu1(seed)
  X <- data$X
  Y <- data$Y
  U1 <- data$U1
  types <- data$types
  true_beta <- data$true_beta

  train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
  X_train <- X[train_id, ]
  Y_train <- Y[train_id]
  X_test <- X[-train_id, ]
  Y_test <- Y[-train_id]

  group_list <- lapply(1:10, function(x) paste0("f_", 
    seq(((x-1)*(20)+1), (x*20))))
  colnames(X_train) <- paste0("f_", seq(1, ncol(X_train)))

  cv_res <- cv.grpregOverlap(X_train, Y_train, 
    group_list, penalty="gel")
  #cv_res <- cv.grpregOverlap(X_train, Y_train, 
  #  group_list, penalty="cMCP")
  beta <- coef(cv_res) # include intercept
  #beta <- cv_res$fit$beta[, which.min(cv_res$cve)]
  beta_est <- beta[-1] # remove intercept

  #####################
  # Check results
  #####################
  feature_FDR[s] <- sum(true_beta==0 & beta_est!=0)/sum(beta_est!=0)
  feature_FOR[s] <- sum(true_beta!=0 & beta_est==0)/sum(beta_est==0)
  
  #### beta ####
  MSE[s] <- mean((X_test %*% beta_est + beta[1] - Y_test)^2)

  ##### compute AUC for features
  prob <- apply(cv_res$fit[[1]][-1,]!=0,1,mean)
  AUC[s] <- pROC::auc(roc(factor(true_beta!=0), prob))

  cat(paste("repeat=",s,"\n",sep=""))
}  
t1<-proc.time()-t0  

save(t1, feature_FDR, feature_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))



