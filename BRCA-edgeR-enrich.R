
################################
rm(list=ls())
library(edgeR)
library(ggpubr)
library(GEOquery)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\")
load("TCGA-BRCA-count.rda")
load("BRCA_predict_cluster.rda")
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\edgeR\\")
TCGA_data<-BRCA_predict_cluster$TCGA
TCGA<-exprs(eset3)
count<-TCGA[apply(TCGA, MARGIN=2,FUN=function(xxx)
{
  (sum(xxx==0)/length(xxx))<=0.5
}),]
TCGA_data<-TCGA_data[colnames(TCGA),]

data = count
cluster2<-rownames(TCGA_data[which(TCGA_data$cluster == 2),])
cluster1<-rownames(TCGA_data[which(TCGA_data$cluster == 1),])

data_2<-data[,cluster2]
data_1<-data[,cluster1]

data<-cbind(data_2,data_1)
group_list <- c(rep("Cluster2",time = length(cluster2)),rep("Cluster1",time = length(cluster1)))
group_list = factor(group_list,levels=c("Cluster2","Cluster1"))
design <- model.matrix(~0+group_list)
rownames(design) = colnames(data)
colnames(design)<-levels(group_list)


#差异表达矩阵
DGElist <- DGEList( counts = data, group = group_list)
## Counts per Million or Reads per Kilobase per Million
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 自定义
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)


fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1)) 
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
save(nrDEG_edgeR,file = 'BRCA_TCGA_edgeR.rda')
head(nrDEG_edgeR)


#提取差异显著的差异矩阵
padj = 0.05# 自定义
foldChange= 1.0 # 自定义
nrDEG_edgeR_signif=nrDEG_edgeR[(nrDEG_edgeR$FDR < padj & 
                                     (nrDEG_edgeR$logFC>foldChange | nrDEG_edgeR$logFC<(-foldChange))),]
nrDEG_edgeR_signif=nrDEG_edgeR_signif[order(nrDEG_edgeR_signif$logFC),]

up<-nrDEG_edgeR_signif[which(nrDEG_edgeR_signif $logFC>0),]
down<-nrDEG_edgeR_signif[which(nrDEG_edgeR_signif $logFC<0),]

save(nrDEG_edgeR_signif,file="BRCA_TCGA_edgeR_signif.rda")
save(up,file="up-edgeR.rda")
save(down,file="down-edgeR.rda")



########################################
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
x<-rownames(nrDEG_edgeR_signif)
eg <- bitr(x, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")

#?enrichGO
#go
go <- enrichGO(gene = eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='BP',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID',
               readable = T)
head(go)


p1<-barplot(go, drop = TRUE, showCategory =10)+
  #barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  #facet_grid(ONTOLOGY~., scale='free')+
  theme(text = element_text(size = 15,family="serif"),
        axis.title=element_text(size=15,family="serif",color="black",face="bold"),
         axis.text.x=element_text(size=15,family="serif",color="black"),
         axis.text.y=element_text(size=15,family="serif",color="black"),
         legend.text=element_text(size=15,family="serif",color="black"),
         legend.title=element_text(size=15,family="serif",color="black",face="bold"),
  )
p3<-p1+ggtitle("GO Biological Process Enrichment")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p5<-p3+labs(x="Gene number",y="Biological Process")
p5
GO_enrich<-p5
ggsave(GO_enrich,filename = "GO_BRCA_TCGA.pdf",width =8.5,height =6)
save(GO_enrich,file="GO_BRCA_TCGA_enrich.rda")

################################
#kegg
kegg<-enrichKEGG(gene = eg$ENTREZID,
                   organism = 'hsa', #KEGG可以用organism = 'hsa'
                   pvalueCutoff = 0.05)
head(kegg,2)

#######################
k1<-barplot(kegg,drop=TRUE,showCategory=10)
k2<-k1+theme(axis.title=element_text(size=15,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=15,family="serif",color="black"),
             axis.text.y=element_text(size=15,family="serif",color="black"),
             legend.text=element_text(size=15,family="serif",color="black"),
             legend.title=element_text(size=15,family="serif",color="black",face="bold")
)

k3<-k2+ggtitle("KEGG Enrichment")+theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
k5B<-k3+labs(x="Gene number",y="KEGG Pathway")
k5B
kegg_enrich<-k5B
ggsave(kegg_enrich,filename = "KEGG_BRCA_TCGA.pdf",width =8.5,height =6)
########################################
save(kegg_enrich,file="kegg_enrich_BRCA_TCGA.rda")

