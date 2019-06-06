#########################################################
# Simulation 1 in paper (one-layer group structure)
# 
# 1/17/2018
# Li Zhu @ pitt
########################################################
rm(list = ls())
library(pROC)
library(BSGS)
packageVersion("BSGS") # 2.2

## code to generate data
source("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2/Code/GeneSimu1.R")

simu1_BSGS_single_cpu <- function(selected_tau2, mcse, start_ite, 
    stop_ite, burnInIter=20000, totalIter=30000){

  ns <- stop_ite - start_ite + 1
  MCSE <- AUC <- feature_FDR <- feature_FOR <- MSE <- rep(NA,ns)

  for(s in seq(start_ite, stop_ite)){

    # for debugging (removed)
    if(FALSE){
        selected_tau2 <- 5
        mcse <- 9999
        start_ite <- 1
        stop_ite <- 1
        s <- 6
        burnInIter <- 20000
        totalIter <- 30000
    }

    # keep the same seed as other methods
    seed <- 2016+100000*s
    data <- geneSimu1(seed=seed)
    
    # if groups overlap, features will be duplicated
    Y <- data$Y
    X_raw <- data$X
    U1 <- data$U1
    N <- length(Y)
    m1 <- ncol(U1)
    U1_rowsum <- rowSums(U1)
    pk <- colSums(U1)   # number of features belonging to each group
    feature_dup_index <- unlist(lapply(1:nrow(U1), function(j) rep(j, U1_rowsum[j])))
    k_index <- unlist(lapply(1:nrow(U1), function(j) which(U1[j,]==1)))
    X <- X_raw[, feature_dup_index]
    P <- ncol(X)

    Group.Info <- list(Group.Index=k_index, num.of.groups=m1, group.lables=unique(k_index), 
      num.of.var.in.groups=pk)
    true_beta <- data$true_beta

    ##### split training and testing 
    train_id <- sample(seq(1, nrow(X)), size = nrow(X)*0.8, replace=FALSE)
    X_train <- X[train_id, ]
    Y_train <- Y[train_id]
    X_test <- X[-train_id, ]
    Y_test <- Y[-train_id]

    ## hyperparametersF
    tau2.value <- rep(selected_tau2, P) 
    rho.value <- rep(0.5, P) 
    theta.value <- rep(0.5, m1)

    ## Initial values and stopping point 
    r.value <- rbinom(P, 1, 0.5) 
    eta.value <- rbinom(m1, 1, 0.5) 
    beta.value <- cbind(c(t(solve(t(X_train) %*% X_train + diag(1/5, P)) %*% t(X_train) %*% Y_train) )) # beta.true 
    sigma2.value <- 1 
    MCSE.Sigma2.Given <- mcse

    res <- BSGS::BSGSSampleC(Y_train, X_train, Group.Info, r.value, eta.value, 
                                beta.value, tau2.value, rho.value, theta.value, 
                                sigma2.value, nu=0, lambda=1, Num.of.Iter.Inside.CompWise=10, 
                                Num.Of.Iteration=totalIter, MCSE.Sigma2.Given)

    BETA <- res$beta.samples
    ## sum up duplicated features
    BETA_true = sapply(1:nrow(U1), function(x) 
        apply(BETA[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
    #all(BETA==t(BETA_true))  # TRUE
    BETA <- t(BETA_true)

    # feature selection
    pos_prob <- apply((BETA!=0)[, -seq(1,burnInIter)], 1, mean)
    feature_FDR[s-start_ite+1]<-sum(true_beta==0 & pos_prob >= 0.5)/sum(pos_prob >= 0.5)
    feature_FOR[s-start_ite+1]<-sum(true_beta!=0 & pos_prob < 0.5)/sum(pos_prob < 0.5)
    AUC[s-start_ite+1]<-pROC::auc(pROC::roc(response=factor(true_beta!=0), predictor=pos_prob))

    # estimate and prediction
    #all(((res$beta.samples!=0)[, -seq(1,burnInIter)]) == 
    #    res$r.samples[, -seq(1,burnInIter)])  # TRUE
    BETA[which(BETA == 0)] <- NA
    mean_na <- function(x) mean(x,na.rm=TRUE)
    beta_mean <- apply(BETA[, -seq(1,burnInIter)],1,mean_na)
    beta_mean[is.na(beta_mean)] <- 0 # all r are zero
    MSE[s-start_ite+1] <- mean((X_test %*% beta_mean - Y_test)^2)
    
    # check MCSE of sigma2
    MCSE[s-start_ite+1] <- bm(res$sigma2.samples[-seq(1,burnInIter),])$se
    print(s)
  }  
  return(list(feature_FDR=feature_FDR, feature_FOR=feature_FOR, MSE=MSE, AUC=AUC, MCSE=MCSE))
}


simu1_BSGS <- function(selected_tau2, mcse, totalSimu=100, num_cpu=1, burnInIter=20000, totalIter=30000){

    setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2/")
    library(snowfall)

    t0 <- proc.time()
    interv <- totalSimu/num_cpu
    start_ite_seq <- seq(1,totalSimu, by=interv)
    stop_ite_seq <- seq(interv, totalSimu, by=interv)
    snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=num_cpu)
    snowfall::sfLibrary(pROC) 
    snowfall::sfLibrary(BSGS) 
    snowfall::sfLibrary(batchmeans) 
    snowfall::sfExport("simu1_BSGS_single_cpu", "geneSimu1", "selected_tau2", "mcse", "start_ite_seq", 
        "stop_ite_seq", "burnInIter", "totalIter")
    res_list <- snowfall::sfClusterApply(1:num_cpu,function(x) 
        simu1_BSGS_single_cpu(selected_tau2=selected_tau2, mcse=mcse, start_ite=start_ite_seq[x], 
            stop_ite=stop_ite_seq[x], burnInIter, totalIter))
    snowfall::sfStop()
    delta_t <- proc.time() - t0  #1.73 hours

    MCSE <- AUC <- MSE <- feature_FDR <- feature_FOR <- NULL
    for(i in 1:num_cpu){
        MCSE <- c(MCSE, res_list[[i]]$MCSE)
        AUC <- c(AUC, res_list[[i]]$AUC)
        MSE <- c(MSE, res_list[[i]]$MSE)
        feature_FDR <- c(feature_FDR, res_list[[i]]$feature_FDR)
        feature_FOR <- c(feature_FOR, res_list[[i]]$feature_FOR)
    }
    res <- list(feature_FDR=feature_FDR, feature_FOR=feature_FOR, AUC=AUC, MSE=MSE, MCSE=MCSE, delta_t=delta_t)
    save(res, file=paste("Results/BSGS_tau2_", selected_tau2, "_mcse_", mcse, "_res.RData", sep=""))
}








