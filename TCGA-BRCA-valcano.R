
########################################################
rm(list= ls()) 
library(ggpubr)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\Valcano\\")
load('BRCA_TCGA_edgeR.rda')

##########################
#火山图
deg_data<-nrDEG_edgeR
head(deg_data)
deg_data$logP<- -log10(deg_data$PValue)#对差异基因校正后p值进行log10转换
deg_data$group <- "not-significant"#新加一列group
deg_data$group[which((deg_data$PValue < 0.05) & (deg_data$logFC > 1))] <-"up-regulated"
deg_data$group[which((deg_data$PValue < 0.05) & (deg_data$logFC < -1))] <-"down-regulated"
table(deg_data$group)

deg_data$label<-""
deg_data<-deg_data[order(deg_data$PValue),]
deg_data$Symbol<-rownames(deg_data)
up_genes<-head(deg_data$Symbol[which(deg_data$group == "up-regulated")],10)
down_genes<-head(deg_data$Symbol[which(deg_data$group == "down-regulated")],10)
deg_top10_genes<-c(as.character(up_genes), as.character(down_genes))
deg_data$label[match(deg_top10_genes,deg_data$Symbol)] <- deg_top10_genes
colnames(deg_data)[7]<-c("Group")

p1<-ggscatter(deg_data,x = "logFC", y="logP",
          color = "Group",
          palette = c("red", "grey", "green"),
          size = 1,
          label = deg_data$label,
          font.label = 8,
          repel = T,
          xlab = "Log2FoldChange",
          ylab = "-Log10(P-value)",)  +
   geom_hline(yintercept = 1.30, linetype="dashed") +
   geom_vline(xintercept = c(-1,1), linetype="dashed")

p2<-p1+theme(axis.title=element_text(size=13,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=12,family="serif",color="black",face="bold"),
             legend.text=element_text(size=12,family="serif",color="black",face="bold")
)
p2
BRCA_edgeR_valcano_TCGA<-p2
save(BRCA_edgeR_valcano_TCGA,file="BRCA_edgeR_valcano_TCGA.rda")
ggsave("BRCA_DEG_valcano_TCGA.pdf",width =6,height = 6)



#########################################


########################################################
rm(list= ls()) 
library(ggpubr)
library(ggrepel)
library(ggplot2)
library(ggrepel)
setwd("D:\\OneDrive - hrbmu.edu.cn\\yanglei004\\Deeplearning\\BRCA-new\\DGE\\TCGA\\Valcano\\")
load('BRCA_TCGA_edgeR.rda')

##########################
#火山图
deg_data<-nrDEG_edgeR
head(deg_data)
deg_data$logP<- -log10(deg_data$PValue)#对差异基因校正后p值进行log10转换
deg_data$group <- "not-significant"#新加一列group
deg_data$group[which((deg_data$PValue < 0.05) & (deg_data$logFC > 1))] <-"up-regulated"
deg_data$group[which((deg_data$PValue < 0.05) & (deg_data$logFC < -1))] <-"down-regulated"
table(deg_data$group)

deg_data$label<-""
deg_data<-deg_data[order(deg_data$PValue),]
deg_data$Symbol<-rownames(deg_data)
up_genes<-head(deg_data$Symbol[which(deg_data$group == "up-regulated")],10)
down_genes<-head(deg_data$Symbol[which(deg_data$group == "down-regulated")],10)
deg_top10_genes<-c(as.character(up_genes), as.character(down_genes))
deg_data$label[match(deg_top10_genes,deg_data$Symbol)] <- deg_top10_genes
colnames(deg_data)[7]<-c("Group")


#################################
p1<-ggplot(deg_data,aes(logFC, logP))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=logP, color=logP))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(0.5,4))+
  # 主题调整：
  theme_bw()+
  # 调整主题和图例位置：
  #theme(panel.grid = element_blank(),
   #     legend.position = c(0.01,0.7),
   #     legend.justification = c(0,1)
  #)+
  theme(panel.grid = element_blank(),
       legend.position ="none")+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10_P-value"),
         size = "none")+
  # 添加标签：
 # geom_text(aes(label=label,color=logP),size=3,vjust=1.5,hjust=1)+
  geom_text_repel(data=deg_data,family="serif",
                  aes(label=label),size=2.5, box.padding=unit(0.25, "lines"),
                  point.padding=unit(0.30,"lines"),
                  max.overlaps = getOption("ggrepel.max.overlaps", default =30))+
   # 修改坐标轴：
  xlab("Log2(Fold change)")+
  ylab("-Log10(P-value)")

p1
p2<-p1+theme(axis.title=element_text(size=13,family="serif",color="black",face= "bold"),
             axis.text.y=element_text(size=12,family="serif",color="black",face="bold"),
             axis.text.x=element_text(size=12,family="serif",color="black",face="bold"),
             legend.title=element_text(size=12,family="serif",color="black",face="bold"),
             legend.text=element_text(size=12,family="serif",color="black",face="bold")
)
p2
BRCA_edgeR_valcano_TCGA_A<-p2
save(BRCA_edgeR_valcano_TCGA_A,file="BRCA_edgeR_valcano_TCGA_A.rda")
ggsave("BRCA_DEG_valcano_TCGA_A.pdf",width =6,height = 6)

