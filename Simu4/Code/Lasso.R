rm(list=ls())
library(Rcpp)
library(RcppEigen)
library(pROC)
library(coda)
library(glmnet)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu4/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"

U_seq<-c(0.2,0.5)
ns<-100
AUC<-F_FDR<-F_FOR<-MSE<-matrix(NA,length(U_seq),ns)

t0<-proc.time()
for(du in 1:length(U_seq)){
  for(s in 1:ns){
#du <- 1
#s <- 1
    res_name <- "Lasso"

    seed <- 2016+100000*s
    source(paste0(code_folder, "R/GeneSimu3_v2_nonOverlapping.R"))

    data <- geneSimu3_nonOverlapping(U=U_seq[du], seed)
    X <- data$X
    Y <- data$Y
    U1 <- data$U1
    U2 <- data$U2
    types <- data$types
    true_beta <- data$true_beta

    ## split the data to training and testing
    train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
    X_train <- X[train_id, ]
    Y_train <- Y[train_id]
    X_test <- X[-train_id, ]
    Y_test <- Y[-train_id]

    ##### Apply the function
    cvfit<-cv.glmnet(X_train,Y_train,family="gaussian", alpha=1, intercept=TRUE)
    coef<-coef(cvfit,s="lambda.min")
    intercept <- as.vector(coef)[1]
    coef_vec<-as.vector(coef)[-1]
    
    ##### feature level
    F_FDR[du,s]<-sum(true_beta==0 & coef_vec!=0)/sum(coef_vec!=0)
    F_FOR[du,s]<-sum(true_beta!=0 & coef_vec==0)/sum(coef_vec==0)
    
    ##### beta 
    MSE[du,s]<-mean((X_test %*% coef_vec+intercept - Y_test)^2)

    # feature AUC
    BETA<-sapply(1:length(cvfit$lambda),function(x)
      coef(cvfit,s=cvfit$lambda[x])[-1])
    prob<-apply(BETA!=0,1,mean)
    AUC[du,s]<-pROC::auc(pROC::roc(factor(true_beta!=0),prob))

    cat(paste("repeat=",s,"\n",sep=""))

  }
  cat(paste("du=",du,"\n",sep=""))
}
t1<-proc.time()-t0

save(t1, F_FDR, F_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))


