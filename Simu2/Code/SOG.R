rm(list = ls())
library(pROC)
library(Rcpp)
library(RcppEigen)
library(coda)
set.seed(20160125)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu2/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"

ns<-100
AUC<-group_FDR<-group_FOR<-feature_FDR<-feature_FOR<-MSE<-rep(NA,ns)
t0<-proc.time()
for(s in 1:ns){

  lassoInit <- FALSE
  MH_ind <- 1
  res_name <- "SOG_MH_intercept_autoConv"
  burnInIter <- 5000
  keepIter <- 5000
  maxIter <- 50000
  pi2_prop_n <- 10
  BernoulliWeighted <- TRUE

  seed <- 2016+100000*s

  source(paste0(code_folder, "R/GeneSimu2_v2.R"))

  data <- geneSimu2(seed)
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

  ## apply SOG
  sourceCpp(paste0(code_folder, "src/SOG_continuous_new.cpp"))
  source(paste0(code_folder, "R/SOG_continuous_new.R")) 
 
  res <- SOG_continuous(X_train, Y_train, U1, types=types,
    center=TRUE, lassoInit=lassoInit, 
    BernoulliWeighted=BernoulliWeighted, 
    pi2_prop_n=pi2_prop_n, MH_ind=MH_ind,
    burnInIter=burnInIter, keepIter=keepIter, maxIter=maxIter, 
    printInt=1000, seed=20170529,  
    PreEstBeta=FALSE, alpha2=1, beta2=1,s2Weighted=FALSE)

  ##### beta 
  # deal with overlapping
  BETA <- res$BETA
  beta_median <- apply(BETA,1,median)
 
  feature_FDR[s]<-sum(true_beta==0 & res$f_qvalue<=0.1)/sum(res$f_qvalue<=0.1)
  feature_FOR[s]<-sum(true_beta!=0 & res$f_qvalue>0.1)/sum(res$f_qvalue>0.1)

  
  beta_median <- apply(BETA,1,median)
  intercept_median <- median(res$intercept)
  Y_hat <- intercept_median + X_test %*% beta_median
  MSE[s]<-mean((Y_hat - Y_test)^2)
  ##### compute AUC for features
  AUC[s]<-pROC::auc(pROC::roc(response=factor(true_beta!=0), 
    predictor=res$f_prob))

  cat(paste("repeat=",s,"\n",sep=""))
}  
t1<-proc.time()-t0  #35.040 

save(t1, feature_FDR, feature_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))

