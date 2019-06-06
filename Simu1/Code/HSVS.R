#########################################################
# Simulation 1 in paper (obe-layer group structure)
# 09/29/2016
# Li Zhu
########################################################
rm(list = ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1")

library(MASS)
library(GIGrvg)
library(pROC)
set.seed(20160125)
ns<-100
AUC<-group_FDR<-group_FOR<-feature_FDR<-feature_FOR<-MSE<-rep(NA,ns)
t0<-proc.time()
for(s in 1:ns){
  set.seed(2016+100000*s)
    #####################
    # Generate Data 
    #####################
    N<-125    
    P<-200           
    G<-10     
    mg<-rep(20,G)
    X<-matrix(NA,N,P)
    Y<-rep(NA,N)
    true_group_ind<-c(rep(1,5),rep(0,5))
    true_beta_nonzero<-rnorm(50,0,sqrt(5))
    true_beta<-c(c(true_beta_nonzero[1:20]),
                 c(true_beta_nonzero[21:30],rep(0,10)),
                 c(true_beta_nonzero[31:40],rep(0,10)),
                 c(true_beta_nonzero[41:45],rep(0,15)),
                 c(true_beta_nonzero[46:50],rep(0,15)),rep(0,100))
    for(n in 1:N){
      z<-rnorm(G,0,1)
      e<-rnorm(P,0,1)
      for(g in 1:G){
        X[n,((g-1)*mg[g]+1):(g*mg[g])]<-(z[g]+e[((g-1)*mg[g]+1):(g*mg[g])])/sqrt(1+1)
      }
    }
    # add overlapping features
    #overlap_index1<-c(1,21)  # overlapping feature 1
    #overlap_index2<-c(41,141) # overlapping feature 2
    #X[,overlap_index1]<-cbind(apply(X[,overlap_index1],1,mean),
    #  apply(X[,overlap_index1],1,mean)) 
    #X[,overlap_index2]<-cbind(apply(X[,overlap_index2],1,mean),
    #  apply(X[,overlap_index2],1,mean))
    Y<-X%*%true_beta+rnorm(N)  # subject to change

    # sum-up overlapping features get true_beta2 vector
    true_beta2<-true_beta
    #true_beta2[overlap_index1[1]]<-sum(true_beta[overlap_index1])
    #true_beta2[overlap_index2[1]]<-sum(true_beta[overlap_index2])
    #true_beta2<-true_beta2[-c(overlap_index1[2],overlap_index2[2])]
    F<-length(true_beta)
    F2<-length(true_beta2)

    ##### Apply the function
    g_index<-rep(seq(1,G),times=mg) 
    start_nsim<-2001
    nsim<-3000
    train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
    X_train <- X[train_id, ]
    Y_train <- Y[train_id]
    X_test <- X[-train_id, ]
    Y_test <- Y[-train_id]

    source("/mnt/glusterfs/liz86/MOG_Regression/R/HSVS/HSVS_modify.R")
    #t0 <- proc.time()
    res <- HSVS(burnin=(start_nsim-1),  N=(nsim - start_nsim + 1),  Y_train,  X_train,  g_index)
    #t1 <- proc.time()-t0   #206.737 
    
    #####################
    # Check results
    #####################
    ##### group level selection
    ##### feature level
    BETA2<-res[[1]]
    #BETA2[overlap_index1[1],]<-apply(res$BETA[overlap_index1,],2,sum)
    #BETA2[overlap_index2[1],]<-apply(res$BETA[overlap_index2,],2,sum)
    #BETA2<-BETA2[-c(overlap_index1[2],overlap_index2[2]),]
    get_95ci<-function(x)quantile(x,probs=c(0.025,0.975))
    beta_CI<-t(apply(BETA2,1,get_95ci))
    beta_ind<-1-(beta_CI[,1]<=0 & beta_CI[,2]>=0)*1
    feature_FDR[s]<-sum(true_beta2==0 & beta_ind==1)/sum(beta_ind==1)
    feature_FOR[s]<-sum(true_beta2!=0 & beta_ind==0)/sum(beta_ind==0)

    #### MSE
    beta_median<-apply(BETA2,1,median)
    MSE[s]<-mean((X_test %*% beta_median - Y_test)^2)
    
    ##### compute AUC for features
    PM<-(BETA2!=0)*1
    pm<-apply(PM,1,mean)
    pos_frac<-apply(BETA2>0,1,mean)
    neg_frac<-apply(BETA2<0,1,mean)
    max_frac<-sapply(1:length(pos_frac),function(x) max(pos_frac[x],neg_frac[x]))
    
    AUC[s]<-pROC::auc(roc(factor((true_beta2!=0)*1), max_frac))
    

    cat(paste("repeat=",s,"\n",sep=""))

}  
t1<-proc.time()-t0
save.image("Results/Veera.RData")  #602.099

