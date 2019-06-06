rm(list=ls())
library(glmnet)
library(pROC)
library(RcppEigen)
library(Rcpp)
library(coda)
library(snowfall)
library(grpregOverlap)

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App1/")
load("Data/Data_clean.RData")  
load("Data/pathway_data_no_singleton.RData")
load("Data/type_nodup.RData")
set.seed(20161120)

####### TSG ##########
AUC_TSG <- rep(NA, 5)                         
AUC_TSG_Var <- rep(NA, 5)
coef_TSG <- matrix(NA, ncol(X_s), 5)
nonzero_TSG <- matrix(NA, 5, 100)
nonzero_ER_TSG <- matrix(NA, 5, 100)
for(i in 1:5){
  fold <- i
  X_s_train <- X_s[which(fold_id != fold), ]
  Y_train <- Y[which(fold_id != fold)]

  level0_group <- lapply(1:nrow(U1), function(x) paste0("f_",x))
  level1_group <- lapply(1:ncol(U1), function(x) paste0("f_", 
    which(U1[,x]==1)))
  level2_group <- lapply(1:ncol(U2), function(x) 
    unlist(level1_group[which(U2[,x]==1)]))

  group_list <- c(level0_group, level1_group, level2_group)
  #group_list[[2]] <- c(group_list[[1]][1], group_list[[2]])
  #group_list[[4]] <- c(group_list[[3]][1], group_list[[4]])
  colnames(X_s_train) <- paste0("f_", seq(1, ncol(X_s_train)))

  cv_res <- cv.grpregOverlap(X_s_train, Y_train, 
    group_list, penalty="grLasso",family="binomial")
  beta <- coef(cv_res) # include intercept
  #beta <- cv_res$fit$beta[, which.min(cv_res$cve)]
  coef_TSG[,i] <- beta_est <- beta[-1] # remove intercept
  select_gene_names <- colnames(X_s)[beta_est != 0]
  select_types <- type_nodup[beta_est != 0]

  Nonzero_BETA <- cv_res$fit[[1]][-1,]!=0
  nonzero_TSG[i, ] <- apply(Nonzero_BETA, 2, sum)
  nonzero_ER_TSG[i, ] <- sapply(1:ncol(Nonzero_BETA), function(x)
      sum(colnames(X_s)[Nonzero_BETA[,x]] %in% 
          pathway_data$'Estrogen signaling pathway_Homo sapiens_hsa04915'))

  ## write pathway tables
  selected_prop <- sapply(1:length(pathway_data), function(x)
      sum(select_gene_names %in% pathway_data[[x]]) / sum(colnames(X_s) %in% 
          pathway_data[[x]]))
  gene_names_per_path <- sapply(1:length(pathway_data), 
      function(x) paste(paste(
          select_gene_names[select_gene_names %in% pathway_data[[x]]], 
          select_types[select_gene_names %in% pathway_data[[x]]],
          sep="_"
          ), collapse="//"))
  pathway_tab <- cbind(names(pathway_data), round(selected_prop, 3),gene_names_per_path)[order(selected_prop, decreasing=TRUE), ]
  colnames(pathway_tab) <- c("Pathway_name", "Selected_prop", "gene_names")
  #write.csv(pathway_tab, file=paste("Results/ER_KEGG_TSG_pathway_table_CV", i, "_balanced.csv", sep=""),
  #    row.names=FALSE)

  ## predict
  pred_prob <- as.vector(X_s[fold_id==i, ] %*% 
    matrix(beta_est, length(beta_est),1))
  AUC_TSG[i] <-  pROC::auc(roc(factor(Y[fold_id == i]), pred_prob))
}

mean(AUC_TSG)  #0.9459368
sd(AUC_TSG)/sqrt(5) #0.01010044
save(nonzero_TSG, file="Results/nonzero_TSG_balanced.RData")
save(nonzero_ER_TSG, file="Results/nonzero_ER_TSG_balanced.RData")

## calculate number of genes and number of genes in ER pathway
apply(coef_TSG!=0,2,sum)
pm <- apply(coef_TSG!=0,1,sum)
names(pm) <- feature_types
pm[order(pm, decreasing=TRUE)][1:10]
genes <- colnames(X_s)[pm>0]

## Pathway analysis
prepare <- function(one_database,significant,whole)
{
  len.gene=length(whole)
  len.DEgene=length(significant)
  len.NonDEgene=len.gene-len.DEgene
  
  DE_IN=length(intersect(significant,one_database))
  DE_OUT=len.DEgene-DE_IN
  NonDE_IN=length(one_database)-DE_IN
  NonDE_OUT=len.NonDEgene-NonDE_IN
  FisherTable=matrix(c(DE_IN,DE_OUT,NonDE_IN,NonDE_OUT),ncol=2)
  return(FisherTable)
}

pathway_fisher_detail_gene <- function(significant,whole,fdr=NULL,database=NULL){
  if(is.null(fdr)) fdr=0.05
  if(is.null(database)){
    return ("Please specify database=Mm.gmtl.c2 or database=Mm.gmtl.c5")
  }
  ### Update the gene sets by dropping genes that do not appear in whole
  gene.overlap2=lapply(database,function(x) intersect(x,whole))
  gene.overlap=levels(factor(unlist(gene.overlap2)))
  
  ### Filter out gene sets that contain less than 5 genes or more than 200 genes
  genesets <- lapply(database,function(x) intersect(x,gene.overlap))
  index=which(sapply(genesets,function(x) length(x)<200 & length(x)>5))
  genesets_pro <- genesets[index]
  names(genesets_pro) <-names(genesets)[index]
  
  path_pval <- lapply(genesets_pro,function(x) fisher.test(prepare(x,significant,whole),alternative="greater")$p.value)
  path_qval <- p.adjust(path_pval, method = "BH")
  sigIndex = which(path_qval<fdr)
  gene_set <- names(genesets_pro)
  print(length(gene_set))
  Pathway = gene_set
  Path_pval = path_pval
  Path_qval = path_qval
  TotalNumOfDEGenes = length(significant)
  NumOfGenesInPath = sapply(genesets_pro,length)
  NumOfDEGenesInPath = sapply(genesets_pro,function(x) length(intersect(x,significant)))
  DEgenes = (sapply(genesets_pro,function(x) paste(intersect(x,significant),collapse=", ")))
  res = cbind(Pathway,Path_pval,Path_qval,TotalNumOfDEGenes,NumOfGenesInPath,NumOfDEGenesInPath,DEgenes )
  return (res)
}

#### apply the function for pathway enrichment analysis ######
# load pathway database
#load("/home06/liz86/BayesGL/exp/BRCA_ER_KEGG/Data/pathway_data_no_singleton.RData")
# add expr, cnv, methyl
pathway_data2<-pathway_data

if(FALSE){
    for(i in 1:length(pathway_data)){
      init<-NA
      for(j in 1:length(pathway_data[[i]])){
        init<-c(init,paste(pathway_data[[i]][j],"gene",sep="_"),paste(pathway_data[[i]][j],"CNV",sep="_"),paste(pathway_data[[i]][j],"methyl",sep="_"))
      }
      pathway_data2[[i]]<-init[-1]
    }
}

# load genelist
pathwayRes1 <- pathway_fisher_detail_gene(significant=genes, 
  whole=unique(colnames(X_s)),fdr=0.05, 
  database=pathway_data2)
res <- pathwayRes1[order(unlist(pathwayRes1[,
  "Path_pval"])),]
#AMPK signaling pathway
#0.2036223
