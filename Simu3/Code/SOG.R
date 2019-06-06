rm(list=ls())
library(Rcpp)
library(RcppEigen)
library(pROC)
library(coda)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu3/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"


U_seq<-c(0.2,0.5)
ns<-100
AUC<-feature_FDR<-feature_FOR<-MSE<-matrix(NA,length(U_seq),ns)

t0<-proc.time()
for(du in 1:length(U_seq)){
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
    source(paste0(code_folder, "R/GeneSimu3_v2.R"))

    data <- geneSimu3(U=U_seq[du], seed)
    X <- data$X
    Y <- data$Y
    U1 <- data$U1
    types <- data$types
    true_beta <- data$true_beta

    ## split the data to training and testing
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
      printInt=1000, seed=20160813,  
      PreEstBeta=FALSE, alpha2=1, beta2=1,s2Weighted=FALSE)

    ##### beta 
    # deal with overlapping
    BETA <- res$BETA
    beta_median <- apply(BETA,1,median)
    
    feature_FDR[du, s]<-sum(true_beta==0 & res$f_qvalue<=0.1) / 
        sum(res$f_qvalue<=0.1)
    feature_FOR[du, s]<-sum(true_beta!=0 & res$f_qvalue>0.1) / 
        sum(res$f_qvalue>0.1)

    beta_median <- apply(BETA,1,median)
    intercept_median <- median(res$intercept)
    Y_hat <- intercept_median + X_test %*% beta_median
    MSE[du, s]<-mean((Y_hat - Y_test)^2)
    ##### compute AUC for features
    AUC[du, s]<-pROC::auc(pROC::roc(response=factor(true_beta!=0), 
      predictor=res$f_prob))


    cat(paste("du=",du,", s=",s,"\n",sep=""))
  }
}

t1<-proc.time()-t0

save(t1, feature_FDR, feature_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))
