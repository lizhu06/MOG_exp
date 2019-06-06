###############################################################
# Simulation 3 in paper (two-layer overlapping group structure)
# 
# 1/17/2018 (corrected on 1/19/2018)
# Li Zhu @ pitt
################################################################
rm(list = ls())
library(pROC)
library(BSGS)
library(batchmeans)
packageVersion("BSGS") # 2.2

## code to generate data
source("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu3/Code/GeneSimu3.R")

## inner cross-validation to select tau2
inner_cv <- function(tau2_cv, inner_test_id, P, 
    m1, X_train, X_raw_train, Y_train, mcse, feature_dup_index, 
    Group.Info, totalIter, burnInIter){
    MSE_inner <- rep(NA, length(inner_test_id))
    for(i in 1:length(inner_test_id)){
        X_train_inner <- X_train[-inner_test_id[[i]], ]
        Y_train_inner <- Y_train[-inner_test_id[[i]]]
        X_raw_test_inner <- X_raw_train[inner_test_id[[i]], ]
        Y_test_inner <- Y_train[inner_test_id[[i]]]

        ## hyperparametersF
        tau2.value <- rep(tau2_cv, P) 
        rho.value <- rep(0.5, P) 
        theta.value <- rep(0.5, m1)

        ## Initial values and stopping point 
        r.value <- rbinom(P, 1, 0.5) 
        eta.value <- rbinom(m1, 1, 0.5) 
        beta.value <- cbind(c(t(solve(t(X_train_inner) %*% X_train_inner + 
            diag(1/5, P)) %*% t(X_train_inner) %*% Y_train_inner) )) # beta.true 
        sigma2.value <- 1 
        MCSE.Sigma2.Given <- mcse

        ## fit the model
        res <- BSGS::BSGSSampleC(Y_train_inner, X_train_inner, 
            Group.Info, r.value, eta.value, 
                  beta.value, tau2.value, rho.value, theta.value, 
                  sigma2.value, nu=0, lambda=1, 
                  Num.of.Iter.Inside.CompWise=10, 
                  Num.Of.Iteration=totalIter, MCSE.Sigma2.Given)

        BETA <- res$beta.samples
        ## sum up duplicated features
        BETA_true = sapply(1:ncol(X_raw_train), function(x) 
            apply(BETA[which(feature_dup_index==x), ,drop=FALSE], 2, sum))
        BETA <- t(BETA_true)

        BETA[which(BETA == 0)] <- NA
        mean_na <- function(x) mean(x,na.rm=TRUE)
        beta_mean <- apply(BETA[, -seq(1,burnInIter)],1,mean_na)
        beta_mean[is.na(beta_mean)] <- 0 # all r are zero
        MSE_inner[i] <- mean((X_raw_test_inner %*% beta_mean - Y_test_inner)^2)
    }
    return(mean(MSE_inner))
}

## single CPU 
simu3_BSGS_single_cpu <- function(U, du, tau2_choice, mcse, start_ite, 
    stop_ite, burnInIter=20000, totalIter=30000){

  ns <- stop_ite - start_ite + 1
  SELECTED_TAU2 <- MCSE <- AUC <- feature_FDR <- feature_FOR <- MSE <- rep(NA,ns)

  for(s in seq(start_ite, stop_ite)){

    # for debugging 
    if(FALSE){
        U <- 0.5
        mcse <- 9999
        start_ite <- 1
        stop_ite <- 1
        s <- 1
        du <- 2
        burnInIter <- 20000
        totalIter <- 30000
        tau2_choice <- c(1,2, 3,4,5)
    }

    # keep the same seed as other methods
    seed <- 2016+100000*s+1000*du
    data <- geneSimu3(U=U, seed=seed)
    
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
    X_raw_train <- X_raw[train_id, ]
    X_raw_test <- X_raw[-train_id,]
    Y_train <- Y[train_id]
    Y_test <- Y[-train_id]

    ## split inner train ID (3 fold)
    eachfold <- floor(nrow(X_train)/3)
    inner_test_id <- list()
    for(t in 1:2){
        inner_test_id[[t]] <- sample(setdiff(seq(1, nrow(X_train)), 
            unlist(inner_test_id)), eachfold, replace=FALSE)
    }
    inner_test_id[[3]] <- setdiff(seq(1, nrow(X_train)), 
        unlist(inner_test_id))

    ## Inner cross-validation to select tau2
    inner_cv_mse <- sapply(1:length(tau2_choice), function(inner_choice) 
        inner_cv(tau2_choice[inner_choice], inner_test_id, P, 
        m1, X_train, X_raw_train, Y_train, mcse, feature_dup_index, 
        Group.Info, totalIter, burnInIter))

    ## select the best tau2
    SELECTED_TAU2[s-start_ite+1] <- selected_tau2 <- tau2_choice[which.min(inner_cv_mse)][1]

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
    BETA <- t(BETA_true)

    # feature selection
    pos_prob <- apply((BETA!=0)[, -seq(1,burnInIter)], 1, mean)
    feature_FDR[s-start_ite+1]<-sum(true_beta==0 & pos_prob >= 0.5)/sum(pos_prob >= 0.5)
    feature_FOR[s-start_ite+1]<-sum(true_beta!=0 & pos_prob < 0.5)/sum(pos_prob < 0.5)
    AUC[s-start_ite+1]<-pROC::auc(pROC::roc(response=factor(true_beta!=0), predictor=pos_prob))

    # prediction
    BETA[which(BETA == 0)] <- NA
    mean_na <- function(x) mean(x,na.rm=TRUE)
    beta_mean <- apply(BETA[, -seq(1,burnInIter)],1,mean_na)
    beta_mean[is.na(beta_mean)] <- 0 # all r are zero
    #mean((X_train %*% beta_mean - Y_train)^2)
    MSE[s-start_ite+1] <- mean((X_raw_test %*% beta_mean - Y_test)^2)
    
    # check MCSE of sigma2
    MCSE[s-start_ite+1] <- bm(res$sigma2.samples[-seq(1,burnInIter),])$se
    print(s)
  }  
  return(list(feature_FDR=feature_FDR, feature_FOR=feature_FOR, 
    MSE=MSE, AUC=AUC, MCSE=MCSE, SELECTED_TAU2=SELECTED_TAU2))
}


simu3_BSGS <- function(U, du, tau2_choice, mcse, totalSimu=100, 
    num_cpu=1, burnInIter=20000, totalIter=30000){

    setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu3/")
    library(snowfall)

    t0 <- proc.time()
    interv <- totalSimu/num_cpu
    start_ite_seq <- seq(1,totalSimu, by=interv)
    stop_ite_seq <- seq(interv, totalSimu, by=interv)
    snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=num_cpu)
    snowfall::sfLibrary(pROC) 
    snowfall::sfLibrary(BSGS) 
    snowfall::sfLibrary(batchmeans) 
    snowfall::sfExport("simu3_BSGS_single_cpu", "geneSimu3", "inner_cv", 
        "U", "du", "tau2_choice", "mcse", "start_ite_seq", 
        "stop_ite_seq", "burnInIter", "totalIter")
    res_list <- snowfall::sfClusterApply(1:num_cpu,function(x) 
        simu3_BSGS_single_cpu(U, du, tau2_choice, mcse=mcse, 
            start_ite=start_ite_seq[x], 
            stop_ite=stop_ite_seq[x], burnInIter, totalIter))
    snowfall::sfStop()
    delta_t <- proc.time() - t0  #1.73 hours

    SELECTED_TAU2 <- MCSE <- AUC <- MSE <- feature_FDR <- feature_FOR <- NULL
    for(i in 1:num_cpu){
        SELECTED_TAU2 <- c(SELECTED_TAU2, res_list[[i]]$SELECTED_TAU2)
        MCSE <- c(MCSE, res_list[[i]]$MCSE)
        AUC <- c(AUC, res_list[[i]]$AUC)
        MSE <- c(MSE, res_list[[i]]$MSE)
        feature_FDR <- c(feature_FDR, res_list[[i]]$feature_FDR)
        feature_FOR <- c(feature_FOR, res_list[[i]]$feature_FOR)
    }
    res <- list(feature_FDR=feature_FDR, feature_FOR=feature_FOR, AUC=AUC, 
        MSE=MSE, MCSE=MCSE, SELECTED_TAU2=SELECTED_TAU2, delta_t=delta_t)
    save(res, file=paste("Results/BSGS_U_", U, "_mcse_", mcse, "_res.RData", sep=""))
}








