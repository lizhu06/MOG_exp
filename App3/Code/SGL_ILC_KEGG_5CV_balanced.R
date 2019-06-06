#####################################
# Predict ER using lasso 5-fold CV
# Li Zhu 
# 11/18/2016
####################################
rm(list=ls())
library(glmnet)
library(grpreg)
library(SGL)
library(AUC)
#library(pROC)

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/app2/")
load("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_morePathways/Data/paper_app1_more_pathways_input.RData")
load("Data/MOG_ILC_CV_Pre_balanced.RData")

##### process outcome
patients<-rownames(d_nodup)
pat_names<-substr(patients,start=1,stop=12)
load("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/Data/tcgaClinical_10_14_15_Processed.Rda")
clin<-tcga_clinical_10_14_15
sum(pat_names %in% rownames(clin))  #770
Hist<-as.character(clin[match(pat_names,rownames(clin)),"Histology_textCode"])

d_nodup<-d_nodup[which(Hist %in% c("ILC","IDC")),]
Y <- Hist[which(Hist %in% c("ILC","IDC"))]
X<-d_nodup
# transfrom to log(odds) for methylation
for(n in 1:nrow(X)){
    X[n,which(type_nodup == "methyl")] <- log(X[n, which(type_nodup == 
        "methyl")] / (1 - X[n, which(type_nodup == "methyl")]))
}
F<-ncol(X)
feature_types<-sapply(1:ncol(X),function(i)paste(colnames(d_nodup)[i],type_nodup[i],sep="_"))

## standardize X
X_s <- scale(X)
Y <- (Y == "ILC")*1

g_index <- as.numeric(factor(featureNames_nodup))

####### Lasso ##########
set.seed(20161120)
AUC_SGL <- rep(NA, 5)                         
coef_SGL <- matrix(NA, ncol(X), 5)
nonzero_SGL <- matrix(NA, 5, 100)
nonzero_ER_SGL <- matrix(NA, 5, 100)
for(i in 1:5){
    cvfit <- cvSGL(list(x=X_s[fold_id != i, ] ,y=Y[fold_id != i]),
        g_index, type="logit", standardize=FALSE, nlam = 100)
    coef_vec <- cvfit$fit$beta[, which.min(cvfit$lldiff)]
    coef_SGL[,i] <- coef_vec
    select_gene_names <- colnames(X_s)[coef_vec != 0]
    select_types <- type_nodup[coef_vec != 0]

    Nonzero_BETA <- cvfit$fit$beta!=0
    nonzero_SGL[i, ] <- apply(Nonzero_BETA, 2, sum)
    nonzero_ER_SGL[i, ] <- sapply(1:ncol(Nonzero_BETA), function(x)
        sum(featureNames_nodup[Nonzero_BETA[,x]] %in% 
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
    pathwayTable <- cbind(names(pathway_data), round(selected_prop, 3),gene_names_per_path)[order(selected_prop, decreasing=TRUE), ]
    colnames(pathwayTable) <- c("Pathway_name", "Selected_prop", "gene_names")
    write.csv(pathwayTable, 
        file=paste("Results/ILC_KEGG_SGL_pathway_table_CV", i, "_balanced.csv", 
            sep=""), row.names=FALSE)

    ## predict
    eta <- X_s[fold_id == i, ] %*% coef_vec + 
        cvfit$fit$intercept[which.min(cvfit$lldiff)]
    pred_prob <- exp(eta) / (1 + exp(eta))
    AUC_SGL[i] <- AUC::auc(roc(pred_prob, factor(Y[fold_id == i])))  
}
mean(AUC_SGL)  #0.683 (balanced: 0.901)
sd(AUC_SGL)/sqrt(5) #0.113 (balanced: 0.011)

save(nonzero_SGL, file="Results/Lasso_nonzero_SGL_balanced.RData")
save(nonzero_ER_SGL, file="Results/Lasso_nonzero_ER_SGL_balanced.RData")

## calculate number of genes and number of genes in ER pathway
apply(coef_SGL!=0,2,sum)
pm <- apply(coef_SGL!=0,1,sum)
names(pm) <- feature_types
pm[order(pm, decreasing=TRUE)][1:10]
genes <- featureNames_nodup[pm>0]

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
pathwayRes1 <- pathway_fisher_detail_gene(significant=genes,whole=unique(featureNames_nodup),fdr=0.05,database=pathway_data2)
res <- pathwayRes1[order(unlist(pathwayRes1[,"Path_pval"])),]
#Endometrial cancer 0.017

save.image("Results/SGL_ILC_5CV_balanced.RData")

