

######################################
rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(maftools)
library(survminer)
library(Hmisc)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\mutation\\TCGA\\")
load("BRCAmaf.rda")
load("BRCA_predict_cluster.rda")

TCGA_cluster<-BRCA_predict_cluster$TCGA
dat=BRCA@data
sample1=unique(as.character(dat$Tumor_Sample_Barcode))
sample2=as.character(row.names(TCGA_cluster))
sample2<-substr(sample2,1,12)
rownames(TCGA_cluster)<-sample2
samp<-intersect(sample1,sample2)

TCGA_cluster1<-TCGA_cluster[samp,]
TCGA_maf1=subsetMaf(BRCA,tsb = samp)

clu1<-rownames(TCGA_cluster1[which(TCGA_cluster1[,3]==1),])
clu2<-rownames(TCGA_cluster1[which(TCGA_cluster1[,3]==2),])

TCGAclu1_maf=subsetMaf(BRCA,tsb = clu1)
TCGAclu2_maf=subsetMaf(BRCA,tsb = clu2)

#----------------------------------------------------------------------------
save(TCGAclu1_maf,TCGAclu2_maf,file="TCGA_cluster_mat.rda")

