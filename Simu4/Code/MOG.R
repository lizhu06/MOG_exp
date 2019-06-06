rm(list=ls())
library(Rcpp)
library(RcppEigen)
library(pROC)
library(coda)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu4/")
code_folder <- "/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/"

U_seq<-c(0.2,0.5)
ns<-100
AUC<-feature_FDR<-feature_FOR<-MSE<-matrix(NA,length(U_seq),ns)

t0<-proc.time()
for(du in 1:length(U_seq)){
  for(s in 1:ns){
#du <- 1
#s <- 1

    lassoInit <- FALSE
    MH_ind <- 0
    res_name <- "MOG_MH_intercept_autoConv"
    burnInIter <- 5000
    keepIter <- 5000
    maxIter <- 10000
    pi2_prop_n <- 10
    pi3_prop_n <- 10
    weight_Tj <- FALSE
    weight_Dk <- FALSE
    weight_Rj <- FALSE


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

    ## apply MOG
    sourceCpp(paste0(code_folder, "src/MOG_continuous_new.cpp"))
    source(paste0(code_folder, "R/MOG_continuous_new.R")) 

    res <- MOG_continuous(X, Y, U1, U2, types=NULL, center=TRUE, 
      lassoInit=lassoInit,  weight_Tj=weight_Tj, weight_Dk=weight_Dk, 
      weight_Rj=weight_Rj,
      pi2_prop_n=pi2_prop_n,  pi3_prop_n=pi3_prop_n, MH_ind=MH_ind,
      burnInIter=burnInIter, keepIter=keepIter, maxIter=maxIter, 
      printInt=1000, seed=seed, 
      alpha1=1, beta1=1, alpha2=1, beta2=1, alpha3=1, beta3=1, 
      weight_s2=FALSE)

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

    cat(paste("repeat=",s,"\n",sep=""))

  }
  cat(paste("du=",du,"\n",sep=""))
}
t1<-proc.time()-t0

save(t1, feature_FDR, feature_FOR, MSE, AUC, 
  file=paste0("Results/", res_name, ".RData"))


