#########################################
# Simulation 2 for paper
# Li Zhu
# 09/28/2016
##################################
rm(list = ls())
library(AUC)
library(glmnet)
library(grpreg)
setwd("/home06/liz86/BayesGL/exp/paper_simulation_3_v3_train_test")

C<-10
G<-100
Pf<-300

U_seq<-c(0.2,0.5)
#U_seq<-c(0.2,0.3,0.4,0.5)
ns<-100
AUC<-F_FDR<-F_FOR<-MSE<-matrix(NA,length(U_seq),ns)

t0<-proc.time()
for(du in 1:length(U_seq)){
  for(s in 1:ns){
    U<-U_seq[du]
    U<-U_seq[du]
    set.seed(2016+100000*s+1000*du)
    #####################
    # Generate Data 
    #####################
    N<-200    # 200 subjects
    P<-306     # 300 genes
    C<-10      # 10 clusters 
    G<-102     # 100 groups 
    m_gc<-c(10,11,10,11,rep(10,6))   # number of groups in each cluster
    m_fg<-rep(3,G)   # number of features in each group
    m_fc<-c(30,33,30,33,rep(30,6))   # number of features in each cluster--
    X<-matrix(NA,N,P)
    Y<-rep(NA,N)
    c_index<-rep(seq(1,C),times=m_fc)   # cluster index for each feature
    c_index_group<-rep(seq(1,C),times=m_gc)   # cluster index for each group
    g_index<-rep(seq(1,G),times=m_fg)      # group index for each feature
    
    true_c_ind<-c(1,1,1,1,1,0,0,0,0,0)
    true_g_ind<-c(c(rep(1,6),rep(0,4)),c(rep(1,7),rep(0,4)),
      c(rep(1,6),rep(0,4)),c(rep(1,7),rep(0,4)),c(rep(1,6),rep(0,4)),
                  rep(0,50))
    true_beta_ind<-c(c(rep(1,18),rep(0,12)),c(rep(1,21),rep(0,12)),
      c(rep(1,18),rep(0,12)),c(rep(1,21),rep(0,12)),c(rep(1,18),rep(0,12)),
                     rep(0,150))
    level_ind<-c(c(rep(3,12),rep(2,6),rep(1,12)),
      c(rep(3,15),rep(2,6),rep(1,12)),c(rep(3,12),rep(2,6),rep(1,12)),
      c(rep(3,15),rep(2,6),rep(1,12)),c(rep(3,12),rep(2,6),rep(1,12)),
                 rep(1,150))
    true_beta<-rep(0,P)
    true_beta[level_ind==1]<-0
    true_beta[level_ind==2]<-sample(c(-1,1),sum(level_ind==2),replace=TRUE,
      prob=c(0.5,0.5))*runif(sum(level_ind==2),U,2*U)
    true_beta[level_ind==3]<-sample(c(-1,1),sum(level_ind==3),replace=TRUE,
      prob=c(0.5,0.5))*runif(sum(level_ind==3),2*U,3*U)
    for(n in 1:N){
      zc<-rnorm(C,0,sqrt(0.2))
      zc_rep<-rep(zc,times=m_fc)
      zg<-rnorm(G,0,sqrt(0.3))
      zg_rep<-rep(zg,times=m_fg)
      e<-rnorm(P,0,sqrt(0.5))
      X[n,]<-zc_rep+zg_rep+e
    }
    ##### create overlapping genes
    over_group_index1<-c(1,11)
    over_group_index2<-c(22,32)
    X[,g_index==over_group_index1[1]]<-X[,g_index==over_group_index1[2]]<-
      (X[,g_index==over_group_index1[1]]
      +X[,g_index==over_group_index1[2]])/2
    X[,g_index==over_group_index2[1]]<-X[,g_index==over_group_index2[2]]<-
      (X[,g_index==over_group_index2[1]]
      +X[,g_index==over_group_index2[2]])/2
    Xf<-X[,-c(which(g_index %in% c(over_group_index1[2],
      over_group_index2[2])))]
    Y<-X%*%true_beta+rnorm(N)
    P<-ncol(X)
    Pf<-ncol(Xf)
    
    # sum-up overlapping features
    true_beta2<-true_beta
    true_beta2[which(g_index==over_group_index1[1])]<-(true_beta[which(
      g_index==over_group_index1[1])]+true_beta[which(
      g_index==over_group_index1[2])])
    true_beta2[which(g_index==over_group_index2[1])]<-(true_beta[which(
      g_index==over_group_index2[1])]+true_beta[which(
      g_index==over_group_index2[2])])
    true_beta2<-true_beta2[-c(which(g_index %in% c(over_group_index1[2],
      over_group_index2[2])))]

    ## apply function
    start_nsim<-1001
    nsim<-2000
    train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
    X_train <- X[train_id, ]
    Y_train <- Y[train_id]
    X_test <- X[-train_id, ]
    Y_test <- Y[-train_id]
    Xf_train <- Xf[train_id, ]
    Xf_test <- Xf[-train_id, ]

    ##### Apply the function
    cvfit<-cv.grpreg(Xf_train,Y_train,group=rep(seq(1,100), each=3), penalty="grLasso",family="gaussian")
    coef<-coef(cvfit,s="lambda.min")
    coef_vec<-as.vector(coef)[-1]

    ##### feature level
    coef_vec2<-coef_vec
    #coef_vec2[which(g_index %in% over_group_index1[1])]<-(coef_vec[which(g_index %in% over_group_index1[1])]+coef_vec[which(g_index %in% over_group_index1[2])])
    #coef_vec2[which(g_index %in% over_group_index2[2])]<-(coef_vec[which(g_index %in% over_group_index2[1])]+coef_vec[which(g_index %in% over_group_index2[2])])
    #coef_vec2<-coef_vec2[-c(which(g_index %in% c(over_group_index1[2],over_group_index2[2])))]
    F_FDR[du,s]<-sum(true_beta2==0 & coef_vec2!=0)/sum(coef_vec2!=0)
    F_FOR[du,s]<-sum(true_beta2!=0 & coef_vec2==0)/sum(coef_vec2==0)

    ##### beta 
    MSE[du,s]<-mean((Xf_test %*% coef_vec2 + coef[1] - Y_test)^2)

    #### feature AUC
    BETA<-matrix(NA,length(true_beta2),length(cvfit$lambda))
    for(t in 1:length(cvfit$lambda)){
      coef_vec<-as.vector(coef(cvfit,s=cvfit$lambda[t]))[-1]
      coef_vec2<-coef_vec
      #coef_vec2[which(g_index %in% over_group_index1[1])]<-(coef_vec[which(g_index %in% over_group_index1[1])]+coef_vec[which(g_index %in% over_group_index1[2])])
      #coef_vec2[which(g_index %in% over_group_index2[2])]<-(coef_vec[which(g_index %in% over_group_index2[1])]+coef_vec[which(g_index %in% over_group_index2[2])])
      #coef_vec2<-coef_vec2[-c(which(g_index %in% c(over_group_index1[2],over_group_index2[2])))]
      BETA[,t]<-coef_vec2
    }
    prob<-apply(BETA!=0,1,mean)
    AUC[du,s]<-AUC::auc(roc(prob, factor(true_beta2!=0)))

    cat(paste("du=",du,", s=",s,"\n",sep=""))
  }
}  
t1<-proc.time()-t0

save.image("Data/GroupLasso.RData")




