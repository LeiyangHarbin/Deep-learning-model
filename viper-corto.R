
###################################
rm(list=ls())
library(viper)
library(corto)
library(HGNChelper)
library(GEOquery)

setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BRCA-ZQ\\dataset\\")
load("BRCA_dataset.rda")

setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\viper\\")
load("TCGA_cluster_dat.rda")
load("Human_TF.rda")
############################
exprSet<-exprs(all_osdata$`TCGA-new.rda`)

tf<-Human_TF
##########################################
#regulon<-corto(a,centroids=tf,nthreads=1,nbootstraps=1000,verbose=TRUE,p=1e-3)
regulon<-corto(exprSet,centroids=tf,nthreads=1,nbootstraps=1000,verbose=TRUE,p=1e-3)
save(regulon,file="regulon_TCGA.rda")
#write.table(a[c(1:200),],file="TCGA_exp.txt")
#B<-round(a,8)
#write.table(a[c(1:1000),],file="TCGA_exp_3.txt")
#write(tf,file="TCGA_tf.txt")
###################################################
exp<-exprSet
c1<-TCGA_cluster_dat[which(TCGA_cluster_dat$cluster == 1),]
c2<-TCGA_cluster_dat[which(TCGA_cluster_dat$cluster == 2),]
exp1<-exp[,colnames(exp)%in%rownames(c1)]
exp2<-exp[,colnames(exp)%in%rownames(c2)]

signature <- rowTtest(exp1,exp2)
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *
                sign(signature$statistic))[, 1]
nullmodel <- ttestNull(exp1,exp2, per = 1000,
                       repos = TRUE, verbose = FALSE)

mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
save(mrs,file = "mrs.rda")
summary(mrs)
p1<-plot(mrs,cex=1.0)


NES<-mrs$es$nes
NES<-sort(NES,decreasing = T)
NES<-data.frame(NES)

pdf("viper-plot.pdf",width=8,height = 6)
plot(mrs, cex = .7,mrs = rownames(NES)[1:10])
dev.off()


################################
mrs_led<-ledge(mrs)
p1<-plot(mrsmrs_led,cex=1.0)
#########################################
##marina Inference of Master Regulators
library(diggit)
cores <- 3*(Sys.info()[1] != "Windows")+1
dgo <- marina(exp1,exp2, regulon=regulon, per=1000, cores=cores)
diggitMR(dgo)
save(dgo,file="TCGA_dgo.rda")





##############################################

############################################
rm(list=ls())
library(viper)
library(corto)
library(HGNChelper)
library(GEOquery)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\master regulator\\TCGA\\viper\\")
load("mrs.rda")
load("sig_MR_Cox.rda")

NES<-mrs$es$nes
NES<-sort(NES,decreasing = T)
NES<-data.frame(NES)

row.names(sig_MR_Cox)
pdf("viper_sig_MR.pdf",width=8.5,height = 6)
#plot(mrs, cex = .7,mrs = rownames(NES)[row.names(sig_MR_Cox),])
B<-c("KCMF1","KIN","NKRF","STAT6","XPA","ZHX1","ZNF226","ZNF396","ZNF491","ZNF497","ZNF589")
plot(mrs, cex = 1.0,mrs = B)
dev.off()




