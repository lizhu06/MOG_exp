rm(list=ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App2/") # transfered from wong06 (code in /Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/Bayesian_variable_Selection/Exp_revision1/Preprocess_RawData_for_App/Code/Preprocess_RawData_BRCA_ER.R)
load("Data/expr_meth_cnv.RData") 
load("Data/omicsTyps.RData")
load("Data/featureNames.RData")
load("Data/KEGG_2016.RData")

dim(d)  #770 14976
length(featureNames)  #14976
length(omicsTyps)   #14976

d2 <- d
# transfrom to log(odds) for methylation (not used)
for(n in 1:nrow(d)){
	d2[n,which(omicsTyps=="methyl")]<-log(d[n,which(omicsTyps=="methyl")]/
		(1-d[n,which(omicsTyps=="methyl")]))
}

### patients ####
patients <- rownames(d2)
pat_names <- substr(patients,start=1,stop=12)
rownames(d2) <- pat_names
load("Data/tcgaClinical_10_14_15_Processed.Rda")
clin<-tcga_clinical_10_14_15
sum(pat_names %in% rownames(clin))  #770
ER <- as.character(clin[match(pat_names,rownames(clin)),"er_status_by_ihc"])
names(ER) <- pat_names

# exclude outcome
d3 <- d2[which(ER %in%c("Positive","Negative")),]
ER <- ER[which(ER %in%c("Positive","Negative"))]
all(names(ER) == rownames(d3))

featureNames3<-featureNames
omicsTyps3<-omicsTyps

dim(d3)  #727 14976
length(featureNames3)  #14976
length(omicsTyps3)   #14976


# pathway data
pathwaydb <- KEGG_2016
num_gene_path<-sapply(1:length(pathwaydb),function(x)length(pathwaydb[[x]]))
#path_name<-sapply(1:length(names(pathwaydb)),function(x)strsplit(names(pathwaydb)[x],split="_")[[1]][1])
length(pathwaydb)  #293

uni_gene <- unique(featureNames) # unique genes in dataset
## exclude genes not in dataset
pathwaydb2<-lapply(1:length(pathwaydb),function(x)intersect(pathwaydb[[x]],uni_gene))
names(pathwaydb2)<-names(pathwaydb)
pathway_data <- pathwaydb2
## Calculate gene numbers in each pathway
pathway_length<-sapply(1:length(pathwaydb2),function(x)length(pathwaydb2[[x]]))

sum(pathway_length > 0)  #292
sum(pathway_length<=50)  #250
sum(pathway_length<=50 & pathway_length >=20)  #123
pathway_data<-pathwaydb2[which(pathway_length <= 50 & pathway_length>=20)]
uni_genes<-unique(unlist(pathway_data))
length(uni_genes)   #1316

## data without duplication
d_nodup <- d3[,which(featureNames %in% uni_genes)]
featureNames_nodup <- featureNames[which(featureNames %in% uni_genes)]
type_nodup <- omicsTyps[which(featureNames %in% uni_genes)]
dim(d_nodup) #727 3910
length(featureNames_nodup) #3910
length(type_nodup) #3910
save(d_nodup,file="Data/d_nodup.RData")
save(featureNames_nodup,file="Data/featureNames_nodup.RData")
save(type_nodup,file="Data/type_nodup.RData")
save(ER, file="Data/ER.RData")
save(pathway_data, file="Data/pathway_data_filtered.RData")

