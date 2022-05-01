

################################
rm(list=ls())
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\GSEA_hallmark\\")
load('BRCA_TCGA_edgeR.rda')
hall<-read.gmt("h.all.v7.4.symbols.gmt")
hall<-separate(data = hall, col = term, into = c("Hallmark", "term"), sep = "HALLMARK_")
hall<-hall[,c(2:3)]

ge = nrDEG_edgeR$logFC
names(ge) = rownames(nrDEG_edgeR)
ge = sort(ge,decreasing = T)
head(ge)

GSEA_Hallmark_TCGA<-GSEA(ge, TERM2GENE = hall,
            nPerm=1000,
            verbose=FALSE,by="fgsea",
            pAdjustMethod="BH",
            pvalueCutoff=1)

save(GSEA_Hallmark_TCGA,file="GSEA_Hallmark_TCGA.rda")


################################
rm(list=ls())
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\GSEA_hallmark\\")
load('GSEA_Hallmark_TCGA.rda')


a<-order(GSEA_Hallmark_TCGA$NES,decreasing = T)[1:8]

gseaplot<-list()
for (i in 1:length(a)) {
  gseaplot[[i]]<-gseaplot2(GSEA_Hallmark_TCGA, geneSetID = a[i], title =GSEA_Hallmark_TCGA$Description[a[i]])
}


Hallmark_GSEA_BRCA1<-plot_grid(gseaplot[[1]],gseaplot[[2]],
          gseaplot[[3]],gseaplot[[4]],nrow=1)
ggsave("Hallmark_GSEA_BRCA1_TCGA.pdf",width =16,height =4.2)

save(Hallmark_GSEA_BRCA1,file="Hallmark_GSEA_BRCA1_TCGA.rda")

###########################
b<-order(GSEA_Hallmark_TCGA$NES,decreasing =F)[1:8]

gseaplot<-list()
for (i in 1:length(b)) {
  gseaplot[[i]]<-gseaplot2(GSEA_Hallmark_TCGA, geneSetID = b[i], title =GSEA_Hallmark_TCGA$Description[b[i]])
}


Hallmark_GSEA_BRCA2<-plot_grid(gseaplot[[1]],gseaplot[[2]],
                            gseaplot[[3]],gseaplot[[4]],nrow=1)
ggsave("Hallmark_GSEA_BRCA2_TCGA.pdf",width =16,height =4.2)

save(Hallmark_GSEA_BRCA2,file="Hallmark_GSEA_BRCA2_TCGA.rda")


Hallmark_GSEA_TCGA<-plot_grid(Hallmark_GSEA_BRCA1,
                         Hallmark_GSEA_BRCA2,nrow=2)
ggsave("Hallmark_GSEA_BRCA_TCGA.pdf",width =16,height=8.5)

save(Hallmark_GSEA_TCGA,file="Hallmark_GSEA_TCGA.rda")

##################################################





################################
rm(list=ls())
library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggrepel)
library(cowplot)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\GSEA_hallmark\\")
load('GSEA_Hallmark_TCGA.rda')

xx<-GSEA_Hallmark_TCGA
res<-xx@result
res<-res[order(res[,2],decreasing=F),]
res$Significant<-ifelse(res$pvalue<0.05,"S","NS")
res$Enriched<-ifelse(res$NES>0,"E","D")
res<-tidyr::unite(res, "Sig_Enriched",Significant,Enriched)

p<-ggplot(res,aes(x=NES, y=-log10(pvalue)))+scale_x_continuous(limits=c(-5,5))+
  geom_point(aes(color=Sig_Enriched),size=5.0)+scale_color_manual(values=c("green","darkmagenta","orange","red"))+
  theme_bw(base_size=12)+theme(legend.position="None")+
  geom_text_repel(data=subset(res,NES>-10&pvalue<1),family="serif",
                  aes(label=Description),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30),
  )#only p<0.05

p1<-p+geom_vline(xintercept=0)+geom_hline(yintercept=1.30103)
p2<-p1+theme(axis.title=element_text(size=15,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=14,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=14,family="serif",color="black",face="bold"),
             legend.title=element_text(size=14,family="serif",color="black",face="bold"),
             legend.text=element_text(size=14,family="serif",color="black",face="bold")
)
#p3<-p2+ggtitle("TCGA")+theme(plot.title=element_text(size=20,family="serif",color="black",face="bold",hjust=0.5))
p4A<-p2+labs(x="NES",y="-log10(P-value)")
p4A

GSEA_Hallmark_vol_TCGA<-p4A

ggsave("GSEA_Hallmark_vol_TCGA.pdf",width=8,height=8,units="in")
save(GSEA_Hallmark_vol_TCGA,file="GSEA_Hallmark_vol_TCGA.rda")

