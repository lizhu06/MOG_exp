rm(list=ls())
#gc()

#GCMdir <- '/home08/xiaoguang/IMSKM/result/20151202/exampleGCM'
#system(paste('mkdir -p',GCMdir))
data_dir <- "/home05/liz86/SOG_MOG_proj/RawData_from_Caleb/TCGA_Integration_Bayesian/"

## load gene expression data
load(paste0(data_dir, 'TCGA_Log2_RNA_Exp_4_3_15.Rda'))
Exprs <- tcgaLog2_rna_exp
dim(Exprs) ## 20531 * 1095 (0,21)

#TCGA methylation: downloaded by Caleb. 2015/10/11
load(paste0(data_dir, "BRCA_methylation450.Rdata"))

dim(BRCA_methylation450)
## 485577    894 (0,1)

#TCGA CNV: downloaded by Caleb. 2015/10/11
## from PGRR
load(paste0(data_dir, "BRCA_CNVc.Rdata"))

dim(BRCA_CNVc)
## 24776  1079 range(-1.29,3.65)


name_Expr <- colnames(Exprs)
name_methy <- substring(colnames(BRCA_methylation450),1,15)
name_methy[!substring(name_methy,1,1)=='T'] <- paste('T',substring(name_methy[!substring(name_methy,1,1)=='T'],1,14),sep='')
colnames(BRCA_methylation450) <- name_methy

name_CNV <- gsub(pattern='[.]',replacement='-',x=substring(colnames(BRCA_CNVc),1,15))
colnames(BRCA_CNVc) <- name_CNV

commonSamples <- intersect(intersect(name_Expr,name_methy),name_CNV)
length(commonSamples)  #770

ExprC <- Exprs[,commonSamples]
CNVC <- BRCA_CNVc[,commonSamples]
methylC <- BRCA_methylation450[,commonSamples]
methy_sd<-apply(methylC,1,sd)
ESR1<-Exprs["ESR1",commonSamples]



############################################################
############################################################
## methylation data cleaning
############################################################
############################################################

## get methylation annotation (from Caleb, can't find it any more)
#setwd("/home08/xiaoguang/IMSKM/data/BRCA/TCGA/Data20150912/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3")
#files <- list.files(pattern="*.txt")
#ProbeName = NULL
#afile <- files[1]

# downloaded from TCGA data portal myself
afile <- '/home05/liz86/SOG_MOG_proj/RawData_from_TCGA_DATA_portal/jhu-usc.edu_BRCA.HumanMethylation450.10.lvl-3.TCGA-A2-A1FV-01A-11D-A13K-05.txt'
ProbeName = NULL
annotation = read.table(afile,sep="\t",skip=1,as.is=TRUE,header=TRUE)
annotation.order = annotation[order(annotation$Composite.Element.REF),] 
Methy_ProbeName = annotation.order$Composite.Element.REF
Methy_genes = annotation.order$Gene_Symbol

all(Methy_ProbeName==rownames(methylC))


## methylation data proning
mehtyMissingProbes = apply(methylC,1,function(x) any(is.na(x)))
sum(mehtyMissingProbes) ## 120922
sum(!mehtyMissingProbes) ## 364655

## remove missing value
BRCA_methylation450_rmNA <- methylC[!mehtyMissingProbes,]
Methy_genes_rmNA <- Methy_genes[!mehtyMissingProbes]
annotation_order_rmNA<-annotation.order[!mehtyMissingProbes,]

## remove na gene name
Methy_genes_na <- is.na(Methy_genes_rmNA)
sum(Methy_genes_na) ## 63
sum(!Methy_genes_na) ## 364592
BRCA_methylation450_rmna <- BRCA_methylation450_rmNA[!Methy_genes_na,]
Methy_genes_rmna <- Methy_genes_rmNA[!Methy_genes_na]
annotation_order_rmna<-annotation_order_rmNA[!Methy_genes_na,]

## remove '' gene name
Methy_genes_NonName <- Methy_genes_rmna==''
sum(Methy_genes_NonName) ## 78210
sum(!Methy_genes_NonName) ## 286382
BRCA_methylation450_clean <- BRCA_methylation450_rmna[!Methy_genes_NonName,]
Methy_genes_clean <- Methy_genes_rmna[!Methy_genes_NonName]
annotation_order_clean<-annotation_order_rmna[!Methy_genes_NonName,]

## clean the gene name
Methy_genes_clean2 <- sapply(strsplit(Methy_genes_clean,';'),function(x) x[1])
dim(BRCA_methylation450_clean)  #286382    770
dim(annotation_order_clean)  #286382      5
length(Methy_genes_clean)  #286382
length(unique(Methy_genes_clean))  #24051
sd_methy<-apply(BRCA_methylation450_clean,1,sd)

## only keep 50kb up and down the start position of that gene
save(Methy_genes_clean2,file="/home06/liz86/BayesGL/exp/Application1/Data/Methy_genes_clean2.RData")
load("/home06/liz86/BayesGL/exp/Application1/Data/my.regions.RData")
load("/home06/liz86/BayesGL/exp/Application1/Data/missing_genes.RData")

sum(missing_genes %in% rownames(ExprC))
uni_genes<-unique(Methy_genes_clean2)
dup_length<-sapply(1:length(uni_genes),function(x)sum(Methy_genes_clean2==uni_genes[x]))
length(dup_length)#20147

num_in_window<-rep(NA,length(dup_length))
new_methy<-matrix(NA,length(uni_genes),ncol(BRCA_methylation450_clean))
for(i in 1:length(uni_genes)){
	if(dup_length[i]==1){
		new_methy[i,]<-BRCA_methylation450_clean[which(Methy_genes_clean2==uni_genes[i]),]
		num_in_window[i]<-"nodup"
	}else if(dup_length[i]>1 & !(uni_genes[i] %in% missing_genes)){
		sub_probe<-annotation_order_clean[which(Methy_genes_clean2==uni_genes[i]),]
		sub<-BRCA_methylation450_clean[which(Methy_genes_clean2==uni_genes[i]),]
		inwin_ind<-(abs(sub_probe[,"Genomic_Coordinate"]-
			mean(my.regions[which(my.regions[,1]==uni_genes[i]),"start_position"]))<=5000)*1
		num_in_window[i]<-sum(inwin_ind)
		if(sum(inwin_ind)>=1){
			new_methy[i,]<-apply(sub[inwin_ind==1,,drop=FALSE],2,mean)
		}else{
			new_methy[i,]<-apply(sub,2,mean)
		}
	}else{
		num_in_window[i]<-"missing"
		new_methy[i,]<-apply(BRCA_methylation450_clean[which(Methy_genes_clean2==uni_genes[i]),],2,mean)
	}
}
dim(new_methy) #20147   770
sd_methy<-apply(new_methy,1,sd)
sum(sd_methy==0)
colnames(new_methy)<-colnames(BRCA_methylation450_clean)
rownames(new_methy)<-uni_genes

## check one gene
s<-100
gene_a<-uni_genes[s]
probe<-annotation_order_clean[which(
	Methy_genes_clean2==gene_a),"Composite.Element.REF"]
new_methy[gene_a,1:10]
BRCA_methylation450_clean[probe,1:10]	

############################################################
############################################################
## CNV data cleaning
############################################################
############################################################


dim(CNVC) ## 24776   770
any(is.na(CNVC)) ## FALSE
index_7SK <- grep('7SK',rownames(CNVC))

CNVCrmq <- CNVC[-index_7SK,]
length(index_7SK) ## 298
dim(CNVCrmq) ## 24478   770

############################################################
############################################################
## gene expression data cleaning
############################################################
############################################################

dim(ExprC) ## 20531
any(is.na(ExprC)) ## FALSE
ExprCrmq <- ExprC[!rownames(ExprC) == '?',]
sum(rownames(ExprC) == '?') ## 29
sum(!rownames(ExprC) == '?') ## 20502

geneCutoffs <- c(0.9,0.8,0.7,0.6,0.5)[3]
for(ff in 1:length(geneCutoffs)){
	#geneCutoff <- geneCutoffs[ff]
	geneCutoff <- 0.5
	ageneFilteringCutoff <- geneCutoff
		

	## fileter by mean 
	gene_mean <- rowMeans(ExprCrmq)
	quantile(gene_mean,geneCutoff)
	ExprC_m <- ExprCrmq[gene_mean>quantile(gene_mean,geneCutoff),]

	dim(ExprC_m) ## 10251   770


	gene_sd <- apply(ExprC_m,1,sd)
	quantile(gene_sd,geneCutoff)

	ExprC_ms <- ExprC_m[gene_sd>quantile(gene_sd,geneCutoff),]
	dim(ExprC_ms) ## 5125  770


	############################################################
	############################################################
	## filter CNV and methylation by gene expression
	############################################################
	############################################################

	mRNAnames <- rownames(ExprC_ms)
	CNVindex <- rownames(CNVCrmq)%in%mRNAnames
	methylindex <- uni_genes%in%mRNAnames

	dim(ExprC_ms) # 5125  770
	CNVC_ms <- CNVCrmq[CNVindex,]
	methylC_ms <- new_methy[methylindex,]
	Methy_genes_ms <- uni_genes[methylindex]
	dim(CNVC_ms) # 4816  770
	dim(methylC_ms) # 5035  770


	omicsTyps <- c(rep('gene',nrow(ExprC_ms)),rep('CNV',nrow(CNVC_ms)),rep('methyl',nrow(methylC_ms)))
	featureNames <- c(rownames(ExprC_ms),rownames(CNVC_ms),Methy_genes_ms)

	group0 <- split(omicsTyps,featureNames)  # for each gene, list omics types
	table(sapply(group0,length))

	length(group0)

	dim_ExprC_ms <- dim(ExprC_ms)
	dim_CNVC_ms <- dim(CNVC_ms)
	dim_methylC_ms <- dim(methylC_ms)
	dimension <- c(dim_ExprC_ms=dim_ExprC_ms, dim_CNVC_ms=dim_CNVC_ms, dim_methylC_ms=dim_methylC_ms)
	write.csv(dimension,'dimension.csv')

	d <- t(rbind(ExprC_ms,CNVC_ms,methylC_ms))
	group1 <- matrix(0,ncol=ncol(d),nrow=length(group0))
	rownames(group1) <- names(group0)
	all(rownames(group1)==names(group0)) 
	groupNames <- names(group0)

	for(i in 1:nrow(group1)){
		print(i)
		group1[i,featureNames == groupNames[i]] = 1
	}


	dim(group1) ## 5125 26855

} ## end of filtering loop

colnames(d)<-featureNames
save(d,file="/home06/liz86/BayesGL/Data/TCGA/expr_meth_cnv.RData")
save(omicsTyps,file="/home06/liz86/BayesGL/Data/TCGA/omicsTyps.RData")
save(featureNames,file="/home06/liz86/BayesGL/Data/TCGA/featureNames.RData")
save(EGFR,file="/home06/liz86/BayesGL/Data/TCGA/EGFR.RData")
save(ESR1,file="/home06/liz86/BayesGL/Data/TCGA/ESR1.RData")





