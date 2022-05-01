

##########################################################
rm(list=ls())
library(limma)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\data\\")
load("BRCA_viper_HumanTF_new.rda")
dat<-BRCA_viper_HumanTF$`TCGA-new.rda`
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\ceRNA-TF\\")
load("BRCA_predict_cluster.rda")
TCGA_cluster<-BRCA_predict_cluster$TCGA


com_sam<-intersect(colnames(dat),rownames(TCGA_cluster))
dat1<-dat[,com_sam]
TCGA_cluster1<-TCGA_cluster[com_sam,]

clu1<-rownames(TCGA_cluster1[which(TCGA_cluster1$cluster==1),])
clu2<-rownames(TCGA_cluster1[which(TCGA_cluster1$cluster==2),])

c1<-dat1[,clu1]
c2<-dat1[,clu2]

data<-cbind(c1,c2)

sampleType<-c(rep("Cluster 1",ncol(c1)),rep("Cluster 2",ncol(c2)))
sampleType<-factor(sampleType,levels=c("Cluster 1","Cluster 2"))
design<-model.matrix(~sampleType)


fit<-lmFit(data,design)
fit<-eBayes(fit)
DEGs=topTable(fit,coef=2,n=Inf)

sig_up_TF<-rownames(DEGs[which(DEGs$logFC> 0.2 & DEGs$P.Value<0.05),])
sig_down_TF<-rownames(DEGs[which(DEGs$logFC< -0.2 & DEGs$P.Value<0.05),])

write.table(sig_up_TF,file = "sig_up_TF.txt",row.names = F)
write.table(sig_down_TF,file = "sig_down_TF.txt",row.names = F)

save(sig_down_TF,sig_up_TF,file="sig_TF.rda")


#---------------------------------------------------------------------------
library(readxl)
TF_GENE<-read_xlsx("3-TRANSFAC-version 12.4 (TF-Gene).xlsx")

up_TF<-TF_GENE[which(TF_GENE$`TF NAME` %in% sig_up_TF),c(2,4)]
down_TF<-TF_GENE[which(TF_GENE$`TF NAME` %in% sig_down_TF),c(2,4)]

up_TF<-data.frame(up_TF)
down_TF<-data.frame(down_TF)

save(up_TF,file = "up_TF.rda")
save(down_TF,file = "down_TF.rda")

write.table(up_TF,"up_TF.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(down_TF,"down_TF.txt",col.names = T,row.names = F,sep = "\t",quote = F)

target_up<-data.frame(c(up_TF[,1],up_TF[,2]),c(rep("TF",times = nrow(up_TF)),rep("gene",times = nrow(up_TF))))
target_down<-data.frame(c(down_TF[,1],down_TF[,2]),c(rep("TF",times = nrow(down_TF)),rep("gene",times = nrow(down_TF))))
target_up<-unique(target_up)
target_down<-unique(target_down)

colnames(target_up)<-c("TF","Gene")        
colnames(target_down)<-c("TF","Gene")

com1<-intersect(up_TF[,1],up_TF[,2])

com2<-intersect(down_TF[,1],down_TF[,2])

target_up1<-target_up[-(which(target_up$TF %in% com1 & target_up$Gene == "gene")),]
target_down1<-target_down[-(which(target_down$TF %in% com2 & target_down$Gene == "gene")),]

write.table(target_up1,"target_up.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(target_down1,"target_down.txt",col.names = T,row.names = F,sep = "\t",quote = F)

##########################################








####################################################
rm(list=ls())
library(limma)
library(multiMiR)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\ceRNA-TF\\")
load("lncRNA_miRNA.rda")
load("multimir_results.rda")
load("sig_TF.rda")
####################################

a<-multimir_results@data
b<-a[which(a$mature_mirna_id %in% lncRNA_miRNA$miRNA),]
TF<-c(sig_up_TF,sig_down_TF)
miRNA_TF<-b[which(b$target_symbol %in% TF),]
miRNA_TF1<-tidyr::unite(miRNA_TF, "mi_m", mature_mirna_id, target_symbol,remove = FALSE)


mirtarbase_TF<-miRNA_TF1[which(miRNA_TF1$database == "mirtarbase"),]
tarbase_TF<-miRNA_TF1[which(miRNA_TF1$database == "tarbase"),]


com<-Reduce(intersect, list(mirtarbase_TF$mi_m,
                            tarbase_TF$mi_m))

miRNA_TF<-data.frame(strsplit2(com,"_"))

colnames(miRNA_TF)<-c("miRNA","TF")

save(miRNA_TF,file = "miRNA_TF.rda")

######################################


lnc_mi_TF<-merge(miRNA_TF,lncRNA_miRNA,by = "miRNA")
lnc_mi_TF<-unique(lnc_mi_TF)

save(lnc_mi_TF,file = "lnc_mi_TF.rda")


##############################

lnc_mi_TF1<-unique(lnc_mi_TF[,c(1,3)])
lnc_mi_TF2<-unique(lnc_mi_TF[,c(1,2)])
colnames(lnc_mi_TF1)<-c("gene1","gene2")
colnames(lnc_mi_TF2)<-c("gene1","gene2")

lnc_mi_TF3<-rbind(lnc_mi_TF1,lnc_mi_TF2)

write.table(lnc_mi_TF3,file = "lnc_mi_TF.txt",row.names = F,quote = F,sep = "\t")


###########################################
##########################构建target信息文件


load("up-different-lncRNA.rda")
load("down-different-lncRNA.rda")
load("up-different-miRNA.rda")
load("down-different-miRNA.rda")

a<-multimir_results@data
b<-a[which(a$target_ensembl %in% rownames(up_lncRNA)),]
c<-a[which(a$target_ensembl %in% rownames(down_lncRNA)),]

target<-data.frame(c(lnc_mi_TF3[,1],lnc_mi_TF3[,2]))
colnames(target)<-"gene"
target[which(target$gene %in% sig_up_TF),2] = "up_TF"
target[which(target$gene %in% sig_down_TF),2] = "down_TF"
target[which(target$gene %in% b$target_symbol),2] = "up_lncRNA"
target[which(target$gene %in% c$target_symbol),2] = "down_lncRNA"
target[which(target$gene %in% rownames(up_miRNA)),2] = "up_miRNA"
target[which(target$gene %in% rownames(down_miRNA)),2] = "down_miRNA"


write.table(target,file = "target.txt",row.names = F,quote = F,sep = "\t")
