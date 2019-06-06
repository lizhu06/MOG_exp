rm(list=ls())
#library(sparcl)
#library(gplots)
library(pROC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App2/")
load("Data/pathway_data_filtered.RData")
load("Data/Data_clean.RData")  

res_cv <-list()
for(i in 1:5){
    load(paste("Results/MOG_ER_CV", i, "_balanced_LassoInit.RData", sep=""))
    res_cv[[i]] <- res
}

AUC1_cv <- AUC2_cv <- rep(NA, 5)
for(i in 1:5){
    X_s_test <- X_s[fold_id == i, ]
    Y_test <- Y[fold_id == i]
    #auc(roc(X_s_test[,feature_types=="ESR1_gene"],as.factor(Y_test)))  #0.948
    PROB <- matrix(NA, length(Y_test), ncol(res_cv[[i]]$BETA))
    for(s in 1:ncol(res_cv[[i]]$BETA)) {
        PROB[, s] <- c(pnorm(X_s_test %*% 
            res_cv[[i]]$BETA[-1, s]+
            res_cv[[i]]$BETA[1,s]))
    }
    prob<-apply(PROB,1,mean)
    AUC2_cv[i] <- pROC::auc(roc(factor(Y_test), prob)) 

    ### plug-in posterior median of beta
    beta_median <- apply(res_cv[[i]]$BETA, 1, median)
    prob<-pnorm(X_s_test%*%beta_median[-1]+beta_median[1])
    AUC1_cv[i] <- pROC::auc(roc(factor(Y_test), prob))  
}
mean(AUC1_cv) 
sd(AUC1_cv)/sqrt(5)

mean(AUC2_cv)  
sd(AUC2_cv)/sqrt(5) 

################################
##### pathway impact score #####
################################
C <- length(pathway_data)
pm_per_pathway_nodup <- pm_per_pathway <- matrix(NA, C, 5)

DFIS <- matrix(NA, nrow(res_cv[[1]]$BETA_dup)-1, 5)
FIS <- matrix(NA, length(unique(feature_types)), 5)
uni_features <- unique(feature_types)
rownames(FIS) <- uni_features

feature_dup_index <- res_cv[[1]]$feature_dup_index
feature_dup_names <- feature_types[feature_dup_index]
rownames(DFIS) <- feature_dup_names

for(i in 1:5){
    # selection probability (already no-duplicated)
    PM <- res_cv[[i]]$BETA[-1, ] != 0
    FIS[,i] <- pm <- apply(PM, 1, mean)
    names(pm) <- colnames(X_s)
    # (pm[order(pm, decreasing=TRUE)])[1:20]
    # beta_median <- apply(res_cv[[i]]$BETA, 1, median)[-1]
    # beta_median

    PM_dup <- res_cv[[i]]$BETA_dup[-1, ]!=0
    k_index_dup <- res_cv[[i]]$k_index[-1]
    l_index_dup <- res_cv[[i]]$l_index[-1]

    ## sort pathways  
    DFIS[,i] <- pm_dup <- apply(PM_dup, 1, mean)
    names(pm_dup) <- feature_dup_names
    #pm_dup[order(pm_dup, decreasing=TRUE)]
    pm_per_pathway[, i] <- sapply(1:C, function(x) mean(
        pm_dup[l_index_dup==x]))

    pm_per_pathway_nodup[, i] <- sapply(1:C, function(x) mean(
        pm[names(pm) %in% pathway_data[[x]]]))

    pathway_length <- sapply(1:C, function(x) sum(l_index_dup==x))
    pm_per_feature <- sapply(1:C, function(x)
        paste(round(pm_dup[l_index_dup==x][order(
            pm_dup[l_index_dup==x], decreasing=T)], 
            2), collapse="//"))
    feature_names <- sapply(1:C, function(x)
        paste(names(pm_dup[l_index_dup==x])[order(
            pm_dup[l_index_dup==x], decreasing=TRUE)],
             collapse="//"))
       
    pathway_impact <- cbind(names(pathway_data), 
        round(pm_per_pathway[,i], 3),round(pm_per_pathway_nodup[,i], 3),
      pathway_length, pm_per_feature, feature_names)[
        order(pm_per_pathway[,i], decreasing=T), ]
    colnames(pathway_impact) <- c("pathway_name", "impact_score", 
        "impact_score_nodup","pathway_size", "PM_features", "feature_names")
    write.csv(pathway_impact, 
        file=paste("Results/MOG_ER_pathway_CV", i ,
            "_balanced_lassoInit.csv", sep=""), row.names=FALSE)
}

## sort pathways
rownames(pm_per_pathway) <- names(pathway_data)
mean_PIS <- apply(pm_per_pathway, 1, mean)
names(mean_PIS) <- names(pathway_data)
mean_PIS_sort <- mean_PIS[order(mean_PIS, decreasing=TRUE)]
mean_PIS_sort[1:5]
#      Estrogen signaling pathway_Homo sapiens_hsa04915
#                                           0.043716825 

## sort dup features (not used)
mean_DFIS <- apply(DFIS, 1, mean)
names(mean_DFIS) <- rownames(DFIS)
mean_DFIS_sort <- mean_DFIS[order(mean_DFIS, 
    decreasing=TRUE)]
mean_DFIS_sort[1:5]

num_dup_times <- sapply(1:length(names(mean_DFIS_sort)), 
    function(x) 
    sum(names(mean_DFIS_sort) == names(mean_DFIS_sort)[x]))

## sort features
mean_FIS <- apply(FIS, 1, mean)
names(mean_FIS) <- rownames(FIS)
mean_FIS_sort <- mean_FIS[order(mean_FIS, decreasing=TRUE)]
mean_FIS_sort[1:5]

num_top_gene <- 20
FIS_table <- matrix(NA, num_top_gene, 3)
for(i in 1:num_top_gene){
    gene_name <- strsplit(names(mean_FIS_sort)[i], split="_")[[1]][1]
    pathway_index <- sapply(1:length(pathway_data), function(x) gene_name %in% 
        pathway_data[[x]])
    in_pathway_name <- names(pathway_data)[pathway_index]
    FIS_table[i, 1] <- names(mean_FIS_sort)[i]
    FIS_table[i, 2] <- mean_FIS_sort[i]
    FIS_table[i, 3] <- paste(in_pathway_name, collapse="//")
    print(i)
}
write.csv(FIS_table, file="Results/FIS_table_top20.csv")

index <- which(names(pathway_data)=="Estrogen signaling pathway_Homo sapiens_hsa04915")
mean_FIS_sort_inER <- sapply(1:length(mean_FIS_sort), function(i) 
    strsplit(names(mean_FIS_sort)[i], split="_")[[1]][1] %in% pathway_data[[index]])
sum(mean_FIS_sort_inER[1:10]) #10
sum(mean_FIS_sort_inER[1:100]) #88
sum(mean_FIS_sort_inER[1:200]) #125
num_nonzero_ER_MOG_nodup <- sapply(1:length(mean_FIS_sort_inER), function(x) sum(mean_FIS_sort_inER[1:x]))
num_nonzero_MOG_nodup <- seq(1, length(mean_FIS_sort_inER))
save(num_nonzero_MOG_nodup, file="Results/num_nonzero_MOG_balanced_nodup_lassoInit.RData")
save(num_nonzero_ER_MOG_nodup, file="Results/num_nonzero_ER_MOG_balanced_nodup_lassoInit.RData")


## sort dup features
index <- which(names(pathway_data)=="Estrogen signaling pathway_Homo sapiens_hsa04915")
mean_dfis <- apply(DFIS, 1, mean)
mean_dfis_sort <- mean_dfis[order(mean_dfis, decreasing=TRUE)]


uni_mean_dfis <- unique(mean_dfis)
uni_mean_dfis_sort <- uni_mean_dfis[order(uni_mean_dfis, decreasing=TRUE)]
num_nonzero_MOG <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(mean_dfis >= uni_mean_dfis_sort[x]))
num_nonzero_ER_MOG <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(c_index[mean_dfis >= uni_mean_dfis_sort[x]]==index))

save(num_nonzero_MOG, file="Results/num_nonzero_MOG_balanced_lassoInit.RData")
save(num_nonzero_ER_MOG, file="Results/num_nonzero_ER_MOG_balanced_lassoInit.RData")



