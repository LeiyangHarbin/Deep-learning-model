

#------------------------------------------------------------------3 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
library(rms)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")
load("BRCAtrain_TCGA_cluster_dat3.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 3 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=2.62e-02",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_3<-p6


ggsave("TCGA-predict-cluster-survival-3.pdf",width =6,height = 6)
save(TCGA_predict_cluster_survival_3,file="TCGA_predict_cluster_survival_3.rda")



#------------------------------------------------------------------4 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
library(rms)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")

load("BRCAtrain_TCGA_cluster_dat4.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 4 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=2.83e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_4<-p6

ggsave("TCGA-predict-cluster-survival-4.pdf",width =6,height = 6)
save(TCGA_predict_cluster_survival_4,file="TCGA_predict_cluster_survival_4.rda")


#------------------------------------------------------------------5 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
library(rms)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")

load("BRCAtrain_TCGA_cluster_dat5.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 5 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=8.35e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_5<-p6

ggsave("TCGA-predict-cluster-survival-5.pdf",width =6,height = 6)
save(TCGA_predict_cluster_survival_5,file="TCGA_predict_cluster_survival_5.rda")



#------------------------------------------------------------------6 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")

load("BRCAtrain_TCGA_cluster_dat6.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 6 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=1.28e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_6<-p6

ggsave("TCGA-predict-cluster-survival-6.pdf",width =7,height = 7)
save(TCGA_predict_cluster_survival_6,file="TCGA_predict_cluster_survival_6.rda")


#------------------------------------------------------------------7 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")

load("BRCAtrain_TCGA_cluster_dat7.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"
                            ,"Cluster 7"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 7 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=1.63e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,2.2),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_7<-p6

ggsave("TCGA-predict-cluster-survival-7.pdf",width =7,height = 7)
save(TCGA_predict_cluster_survival_7,file="TCGA_predict_cluster_survival_7.rda")


#------------------------------------------------------------------8 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")
load("BRCAtrain_TCGA_cluster_dat8.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"
                            ,"Cluster 7","Cluster 8"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 8 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=4.60e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,2.2),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_8<-p6

ggsave("TCGA-predict-cluster-survival-8.pdf",width =8,height = 8)
save(TCGA_predict_cluster_survival_8,file="TCGA_predict_cluster_survival_8.rda")


#------------------------------------------------------------------9 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")
load("BRCAtrain_TCGA_cluster_dat9.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"
                            ,"Cluster 7","Cluster 8","Cluster 9"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 9 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=1.44e-02",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,2.4),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_9<-p6

ggsave("TCGA-predict-cluster-survival-9.pdf",width =8,height = 8)
save(TCGA_predict_cluster_survival_9,file="TCGA_predict_cluster_survival_9.rda")



#------------------------------------------------------------------10 cluster
##########################################################################
rm(list=ls())
library(survminer)
library(survival)
library(survcomp)
library(ggplot2)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\cox-cluster-10\\")

load("BRCAtrain_TCGA_cluster_dat10.rda")
data1<-BRCAtrain_TCGA_cluster_dat1

fit2 <- survfit(Surv(times, status)~cluster, data = BRCAtrain_TCGA_cluster_dat1)
Bcox<-coxph(Surv(times, status)~cluster, data=BRCAtrain_TCGA_cluster_dat1)

summcph<-summary(Bcox)
summcph
summcph$sctest[3]
p1=ggsurvplot(fit2,legend.title="Subtype",
              main="Overall survival",
              legend.labs=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6"
                            ,"Cluster 7","Cluster 8","Cluster 9","Cluster 10"),
              pval=F, conf.int = T,#
              risk.table =T, # Add risk table
              #risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme=theme_survminer(base_size=12,base_family="serif",
                                      font.main=c(14, "bold","black"),
                                      font.submain= c(13, "bold", "black"),
                                      font.caption = c(14, "bold", "black"),
                                      font.x=c(14, "bold", "black"),font.legend = c(14,"bold"),
                                      font.y=c(14, "bold", "black"),
                                      font.tickslab=c(13, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Kaplan-Meier Curve for 10 Subtypes")+#设置plot图形
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.1))
p4<-p3+annotate("text", x=30, y=0.25, label="P=7.77e-01",color="black",cex=5)
p5<-p2$table+theme(legend.position="none")#设置table图形
p6=ggarrange(p4,p5,heights=c(5,2.4),ncol=1,nrow=2, align = "v")
p6
TCGA_predict_cluster_survival_10<-p6

ggsave("TCGA-predict-cluster-survival-10.pdf",width =8.5,height = 8.5)
save(TCGA_predict_cluster_survival_10,file="TCGA_predict_cluster_survival_10.rda")

