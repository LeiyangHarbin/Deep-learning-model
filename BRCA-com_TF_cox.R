
###########################################
rm(list=ls())
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-com\\")
load("Sta_Time_BRCA.rda")
load("BRCA_viper_com_TF.rda")

###############################################################TCGA
dat_TCGA<-BRCA_viper_com_TF$TCGA
dat1<-cbind(Sta_Time_BRCA[colnames(dat_TCGA),c(2,3)],t(dat_TCGA))

library(survival)
Coxoutput_TCGA <- NULL
for (i in 3:ncol(dat1)) {
  cox<-coxph(Surv(times, status)~as.numeric(as.character(dat1[,i]))>median(as.numeric(as.character(dat1[,i]))),data=dat1)
  coxSummary = summary(cox)
  Coxoutput_TCGA=rbind.data.frame(Coxoutput_TCGA,
                                   data.frame(TF = colnames(dat1)[i],
                                              HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                              z=as.numeric(coxSummary$coefficients[,"z"]),
                                              pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                              lower=as.numeric(coxSummary$conf.int[,3]),
                                              upper=as.numeric(coxSummary$conf.int[,4]),
                                              stringsAsFactors = F))
}
head(Coxoutput_TCGA)

cox.pcutoff <- 0.05
sig_TCGA <- Coxoutput_TCGA[which(Coxoutput_TCGA$pvalue < cox.pcutoff),"TF"]
save(sig_TCGA,file = "sigTF_TCGA.rda")


###############################################################-----GSE42568
dat_GSE42568<-BRCA_viper_com_TF$GSE42568
dat1<-cbind(Sta_Time_BRCA[colnames(dat_GSE42568),c(2,3)],t(dat_GSE42568))
dat1<-na.omit(dat1)

library(survival)
Coxoutput_GSE42568 <- NULL
for (i in 3:ncol(dat1)) {
  cox<-coxph(Surv(times, status)~as.numeric(as.character(dat1[,i]))>median(as.numeric(as.character(dat1[,i]))),data=dat1)
  coxSummary = summary(cox)
  Coxoutput_GSE42568=rbind.data.frame(Coxoutput_GSE42568,
                                  data.frame(TF = colnames(dat1)[i],
                                             HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                             z=as.numeric(coxSummary$coefficients[,"z"]),
                                             pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                             lower=as.numeric(coxSummary$conf.int[,3]),
                                             upper=as.numeric(coxSummary$conf.int[,4]),
                                             stringsAsFactors = F))
}
head(Coxoutput_GSE42568)

cox.pcutoff <- 0.05
sig_GSE42568 <- Coxoutput_GSE42568[which(Coxoutput_GSE42568$pvalue < cox.pcutoff),"TF"]
save(sig_GSE42568,file = "sigTF_GSE42568.rda")


###############################################################---------GSE9893
dat_GSE9893<-BRCA_viper_com_TF$GSE9893
dat1<-cbind(Sta_Time_BRCA[colnames(dat_GSE9893),c(2,3)],t(dat_GSE9893))
dat1<-na.omit(dat1)

library(survival)
Coxoutput_GSE9893 <- NULL
for (i in 3:ncol(dat1)) {
  cox<-coxph(Surv(times, status)~as.numeric(as.character(dat1[,i]))>median(as.numeric(as.character(dat1[,i]))),data=dat1)
  coxSummary = summary(cox)
  Coxoutput_GSE9893=rbind.data.frame(Coxoutput_GSE9893,
                                      data.frame(TF = colnames(dat1)[i],
                                                 HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                 z=as.numeric(coxSummary$coefficients[,"z"]),
                                                 pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                 lower=as.numeric(coxSummary$conf.int[,3]),
                                                 upper=as.numeric(coxSummary$conf.int[,4]),
                                                 stringsAsFactors = F))
}
head(Coxoutput_GSE9893)

cox.pcutoff <- 0.05
sig_GSE9893 <- Coxoutput_GSE9893[which(Coxoutput_GSE9893$pvalue < cox.pcutoff),"TF"]
save(sig_GSE9893,file = "sigTF_GSE9893.rda")


###############################################################-----METABRIC
dat_METABRIC<-BRCA_viper_com_TF$METABRIC
dat1<-cbind(Sta_Time_BRCA[colnames(dat_METABRIC),c(2,3)],t(dat_METABRIC))
dat1<-na.omit(dat1)

library(survival)
Coxoutput_METABRIC <- NULL
for (i in 3:ncol(dat1)) {
  cox<-coxph(Surv(times, status)~as.numeric(as.character(dat1[,i]))>median(as.numeric(as.character(dat1[,i]))),data=dat1)
  coxSummary = summary(cox)
  Coxoutput_METABRIC=rbind.data.frame(Coxoutput_METABRIC,
                                     data.frame(TF = colnames(dat1)[i],
                                                HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                z=as.numeric(coxSummary$coefficients[,"z"]),
                                                pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                lower=as.numeric(coxSummary$conf.int[,3]),
                                                upper=as.numeric(coxSummary$conf.int[,4]),
                                                stringsAsFactors = F))
}
head(Coxoutput_METABRIC)

cox.pcutoff <- 0.01
sig_METABRIC <- Coxoutput_METABRIC[which(Coxoutput_METABRIC$pvalue < cox.pcutoff),"TF"]
save(sig_METABRIC,file = "sigTF_METABRIC.rda")



###############################################################---------GSE96058
dat_GSE96058<-BRCA_viper_com_TF$GSE96058
dat1<-cbind(Sta_Time_BRCA[colnames(dat_GSE96058),c(2,3)],t(dat_GSE96058))
dat1<-na.omit(dat1)

library(survival)
Coxoutput_GSE96058 <- NULL
for (i in 3:ncol(dat1)) {
  cox<-coxph(Surv(times, status)~as.numeric(as.character(dat1[,i]))>median(as.numeric(as.character(dat1[,i]))),data=dat1)
  coxSummary = summary(cox)
  Coxoutput_GSE96058=rbind.data.frame(Coxoutput_GSE96058,
                                      data.frame(TF = colnames(dat1)[i],
                                                 HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                 z=as.numeric(coxSummary$coefficients[,"z"]),
                                                 pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                 lower=as.numeric(coxSummary$conf.int[,3]),
                                                 upper=as.numeric(coxSummary$conf.int[,4]),
                                                 stringsAsFactors = F))
}
head(Coxoutput_GSE96058)

cox.pcutoff <-0.01
sig_GSE96058 <- Coxoutput_GSE96058[which(Coxoutput_GSE96058$pvalue < cox.pcutoff),"TF"]
save(sig_GSE96058,file = "sigTF_GSE96058.rda")


###################################################################
sig_cox_TF_com<-list()
sig_cox_TF_com[1]<-list(sig_TCGA)
sig_cox_TF_com[2]<-list(sig_METABRIC)
sig_cox_TF_com[3]<-list(sig_GSE9893)
sig_cox_TF_com[4]<-list(sig_GSE96058)
sig_cox_TF_com[5]<-list(sig_GSE42568)

names(sig_cox_TF_com)[1]<-c("sig_TCGA")
names(sig_cox_TF_com)[2]<-c("sig_METABRIC")
names(sig_cox_TF_com)[3]<-c("sig_GSE9893")
names(sig_cox_TF_com)[4]<-c("sig_GSE96058")
names(sig_cox_TF_com)[5]<-c("sig_GSE42568")

save(sig_cox_TF_com,file="sig_cox_TF_com.rda")


A1<-intersect(sig_TCGA,sig_METABRIC)
A2<-intersect(A1,sig_GSE9893)
A3<-intersect(A2,sig_GSE96058)
com_sig_TF<-intersect(A3,sig_GSE42568)
save(com_sig_TF,file ="com_sig_TF.rda")
