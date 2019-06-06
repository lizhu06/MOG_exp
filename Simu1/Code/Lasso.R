rm(list = ls())
library(pROC)
library(glmnet)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu1/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"
set.seed(20160125)
ns<-100
AUC<-group_FDR<-group_FOR<-feature_FDR<-feature_FOR<-MSE<-rep(NA,ns)
t0<-proc.time()
for(s in 1:ns){
  res_name <- "Lasso_intercept"

  seed <- 2016+100000*s
  source(paste0(code_folder, "R/GeneSimu1.R"))

  data <- geneSimu1(seed)
  X <- data$X
  Y <- data$Y
  U1 <- data$U1
  types <- data$types
  true_beta <- data$true_beta

  ## split the data to training and testing
  burnInIter <- 2000
  totalIter <- 3000
  train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
  X_train <- X[train_id, ]
  Y_train <- Y[train_id]
  X_test <- X[-train_id, ]
  Y_test <- Y[-train_id]

  ##### Apply the function
  cvfit <- cv.glmnet(X_train, Y_train,family="gaussian", alpha=1)
  coef <- coef(cvfit,s="lambda.min")
  intercept <- as.vector(coef)[1]
  coef_vec <- as.vector(coef)[-1]

  ##### feature level
  feature_FDR[s]<-sum(true_beta==0 & coef_vec!=0)/sum(coef_vec!=0)
  feature_FOR[s]<-sum(true_beta!=0 & coef_vec==0)/sum(coef_vec==0)
  
  ##### beta 
  # deal with overlapping
  MSE[s]<-mean((X_test %*% coef_vec+coef[1] - Y_test)^2)
  #pred <- predict.cv.glmnet(cvfit, X_test, s="lambda.min") # same

  # feature AUC
  BETA<-sapply(1:length(cvfit$lambda),function(x)
    coef(cvfit,s=cvfit$lambda[x])[-1])
  prob<-apply(BETA!=0,1,mean)
  AUC[s]<-pROC::auc(roc(factor(true_beta!=0), prob))

  cat(paste("repeat=",s,"\n",sep=""))
}  
t1<-proc.time()-t0
save(t1, feature_FDR, feature_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))

