#########################################################
# Simulation 1 for paper (one-layer group structure)
# 09/28/2016
# Li Zhu
########################################################
rm(list = ls())
library(pROC)
library(SGL)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1")

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
      start_nsim<-1001
      nsim<-2000
      train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
      X_train <- X[train_id, ]
      Y_train <- Y[train_id]
      X_test <- X[-train_id, ]
      Y_test <- Y[-train_id]

  g_index<-rep(seq(1,G),times=mg) 

  ##### Apply the function
  cvfit<-cvSGL(list(x=X_train,y=Y_train),index=g_index,type="linear",
    standardize =FALSE)
  coef_vec<-cvfit$fit$beta[,which.min(cvfit$lldiff)]
  
  ##### group level selection
  selected_group<-unique(g_index[coef_vec!=0])
  g_ind<-seq(1,G) %in% selected_group
  group_FDR[s]<-sum(true_group_ind==0 & g_ind==1)/
    sum(g_ind==1)
  group_FOR[s]<-sum(true_group_ind==1 & g_ind==0)/
    sum(g_ind==0)

  ##### feature level
  coef_vec2<-coef_vec
  #coef_vec2[overlap_index1[1]]<-sum(coef_vec2[overlap_index1])
  #coef_vec2[overlap_index2[1]]<-sum(coef_vec2[overlap_index2])
  #coef_vec2<-coef_vec2[-c(overlap_index1[2],overlap_index2[2])]
  feature_FDR[s]<-sum(true_beta2==0 & coef_vec2!=0)/sum(coef_vec2!=0)
  feature_FOR[s]<-sum(true_beta2!=0 & coef_vec2==0)/sum(coef_vec2==0)
  MSE[s]<-mean((X_test %*% coef_vec2 - Y_test)^2)

  #### feature AUC
  BETA<-cvfit$fit$beta
  IND2<-IND<-BETA!=0
  #IND2[overlap_index1[1],]<-apply(IND2[overlap_index1,],2,max)
  #IND2[overlap_index2[1],]<-apply(IND2[overlap_index2,],2,max)
  #IND2<-IND2[-c(overlap_index1[2],overlap_index2[2]),]
  prob<-apply(IND2,1,mean)
  AUC[s]<-pROC::auc(roc(factor(true_beta2!=0), prob))
  cat(paste("repeat=",s,"\n",sep=""))
}  
t1<-proc.time()-t0

save.image("Results/SGL.RData")

