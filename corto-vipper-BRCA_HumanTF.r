
##########################################
rm(list=ls())
library(corto)
library(viper)
library(HGNChelper)
library(GEOquery)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\BRCA-ZQ\\dataset\\")
load("BRCA_dataset.rda")
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\TF-TMK\\vipper-dataset\\vipper-BRCA\\HumanTF\\")
load("Human_TF.rda")
load("BRCA_regulon_HumanTF.rda")
BRCA_viper_HumanTF<-list()
#BRCA_regulon_HumanTF<-list()
for (i in 1:length(all_osdata))
{

exprSet<-exprs(all_osdata[[i]])
a<-t(exprSet)
an<-colnames(a)
tn<-Human_TF
tn<-checkGeneSymbols(tn)
an<-checkGeneSymbols(an)
an<-na.omit(an)
tn<-na.omit(tn)
a<-a[,an$x]
colnames(a)<-an$Suggested.Symbol
tf<-intersect(an$Suggested.Symbol,tn$Suggested.Symbol)
##########################################

a<-t(a)
a<-na.omit(a)
regulon<-BRCA_regulon_HumanTF[[i]]
#regulon<-corto(a,centroids=tf,nthreads=1,nbootstraps=1000,verbose=TRUE,p=1e-3)
#a是经过转换的表达矩阵，tf是转换后取交集得出的转录因子名

#vpres<-viper(a,regulon, dnull = NULL,pleiotropy = FALSE, nes = TRUE,
#      method = c("scale"),bootstraps =1000,#
#      minsize=10,adaptive.size = FALSE, eset.filter = TRUE,
#      mvws=1,pleiotropyArgs=list(regulators = 0.05, shadow = 0.05,
#                                      targets = 10, penalty=20, method = "adaptive"),cores=1,
#      verbose = TRUE)

vpres<-viper(a,regulon, dnull = NULL,pleiotropy = FALSE, nes =F,
             method = c("scale"),#
             minsize=10,adaptive.size = FALSE, eset.filter = TRUE,
             mvws=1,pleiotropyArgs=list(regulators = 0.05, shadow = 0.05,
                                        targets = 10, penalty=20, method = "adaptive"),cores=1,
             verbose = TRUE)
#这里是算Protein-activity的
BRCA_viper_HumanTF[i]=list(vpres)
#BRCA_regulon_HumanTF[i]<-list(regulon)
names(BRCA_viper_HumanTF)[i]=names(all_osdata)[i]
#names(BRCA_regulon_HumanTF)[i]=names(all_osdata)[i]
print(i)

}

save(BRCA_viper_HumanTF,file="BRCA_viper_HumanTF_new.rda")
#save(BRCA_regulon_HumanTF,file="BRCA_regulon_HumanTF.rda")

#B<-BRCA_viper_HumanTF[[22]]
#C<-B$nes
