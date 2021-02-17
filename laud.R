####loading pacakges####
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(clusterProfiler)
library(AnnotationHub)  
library(org.Hs.eg.db)   
library(ggplot2)
library(DOSE) 
library(PerformanceAnalytics)
library(WGCNA)
library(pdftools)
set.seed(123456)
options(stringsAsFactors = F)
memory.limit(102400)
####loading data####
ladc.data <- read.table("./GSE69405_PROCESSED_GENE_TPM_ALL.txt",header = T,row.names = 1)

ladc_uni <- ladc.data[!duplicated(ladc.data$gene_name),]
rownames(ladc_uni) <- ladc_uni$gene_name
pt<- ladc_uni[,grep("^LC.PT",colnames(ladc_uni))]
mbt <- ladc_uni[,grep("^LC.MBT",colnames(ladc_uni))]

my_data<- cbind(pt,mbt)
dat <- my_data[,grep("_SC",colnames(my_data))]
colnames(dat)<- gsub("_",".",colnames(dat))

n1 <- sum(grepl(".PT",colnames(dat)))
dat_pt <- dat[,1:n1]
n2 <- ncol(dat)-n1
dat_mat <- dat[,n1+1:n2]
colnames(dat_pt) <- paste(colnames(dat_pt),"pt",sep = "_")
colnames(dat_mat) <- paste(colnames(dat_mat),"mbt",sep = "_")

sam.name <- "ladc_001"
dat_all <- cbind(dat_pt,dat_mat)
dir.create(sam.name)

ladc <- CreateSeuratObject(
  dat_all,
  project = "pt_mbt", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

save(ladc,file=paste0("./ladc_result/ladc_SeuratObject.RData"))
####quality control####
load("./ladc_result/ladc_SeuratObject.RData")
ladc[["percent.mt"]] <- PercentageFeatureSet(ladc, pattern = "^MT-")
p1_1 <- VlnPlot(ladc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+
  theme(legend.position = "none") +
  labs(tag = "A")
p1_1
plot1 <- FeatureScatter(ladc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ladc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

ladc <- subset(ladc, subset= nFeature_RNA > 6000 & percent.mt < 35)

p1_2 <- FeatureScatter(ladc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                       group.by = "orig.ident",pt.size = 1.3) +
  labs(tag = "B")
p1_2
#VariableFeatures
ladc[["RNA"]] <- CreateAssayObject(counts = log(ladc[["RNA"]]@counts+1))

ladc <- FindVariableFeatures(ladc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(ladc), 10)
plot1 <- VariableFeaturePlot(ladc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

p1_3 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) +
  theme(legend.position = "none") +
  labs(tag = "C")
p1_3
p1_1 | p1_2 | p1_3

####Dimension reduction####
all.genes <- rownames(ladc)
ladc <- ScaleData(ladc, features = all.genes)
ladc <- RunPCA(ladc, features = VariableFeatures(ladc),dims = 1:50) 
VizDimLoadings(ladc, dims = 1:2, reduction = "pca")
s1 <- DimHeatmap(ladc, dims = c(1:6), cells = 500, balanced = TRUE)

p2_1 <- DimPlot(ladc, reduction = "pca", group.by="orig.ident",label=T,label.size = 8,pt.size = 2)+
  labs(tag = "D")
p2_1

p2_2 <- ElbowPlot(ladc, ndims=30, reduction="pca") +
  theme(legend.position="bottom") +
  labs(tag = "E")
p2_2
p2_1| p2_2 

pc.num=1:20
ladc <- FindNeighbors(ladc, dims = pc.num) 

ladc <- FindClusters(ladc, resolution = 0.4)
table(ladc@meta.data$seurat_clusters)
ladc = RunTSNE(ladc, dims = pc.num)

new.cluster.ids <- c("0"="pt",
  "1"="mbt"
                     )
names(new.cluster.ids) <- levels(ladc)
ladc <- RenameIdents(ladc,new.cluster.ids)  

p3_1 <- DimPlot(ladc, reduction = "tsne",label=T,label.size = 8,pt.size = 2)+
  labs(tag = "F")
p3_1

diff.wilcox = FindAllMarkers(ladc,
                             min.pct = 0.1, 
                             logfc.threshold = 0.25
                             )

write.table(diff.wilcox,
            file=paste0("./ladc_result/total_marker_genes_tsne_20_PC.txt"),
            sep="\t",quote = F,row.names = F)

save(ladc,diff.wilcox,file = "./ladc_result/diff.Rdata")



####DEGs####
load("./ladc_result/diff.Rdata")

all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_logFC) > 1)

top500 <- all.markers%>% group_by(cluster) %>% top_n(n = 500, wt = avg_logFC)
top500_tib <- top500 
write.csv(top500_tib,file = "./ladc_result/4_top500.csv")
DEG_500 <- top500
deg_500 <- DEG_500$gene
deg_500_pt <- deg_500[1:500]
deg_500_mbt <- deg_500[501:1000]


write.table(deg_500,file = "./ladc_result/DEG_all.txt",quote = F,row.names = F)
write.table(deg_500_pt,file = "./ladc_result/DEG_pt_500.txt",quote = F,row.names = F)
write.table(deg_500_mbt,file = "./ladc_result/DEG_mbt_500.txt",quote = F,row.names = F)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_tib <- top10 
DEG <- top10
deg <- DEG$gene
deg_pt <- deg[1:10]
deg_mbt <- deg[11:20]

write.table(deg_pt,file = "DEG_pt_top10.txt",quote = F,row.names = F)
write.table(deg_mbt,file = "DEG_mbt_top10.txt",quote = F,row.names = F)

top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(ladc)) 
length(top10)
length(unique(sort(top10)))

p3_2 <- DoHeatmap(ladc, features = top10, group.by = "orig.ident")+
  labs(tag = "G")
p3_1 | p3_2

save(ladc,file="./ladc_result/seurat_pipiline.Rdata")


####molecle####
exp.matrix<-as(as.matrix(ladc@assays$RNA@data), 'sparseMatrix')
feature_ann<-data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann)<-rownames(exp.matrix)
exp_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-ladc@meta.data
rownames(sample_ann)<-colnames(exp.matrix)
exp_pd<-new("AnnotatedDataFrame", data =sample_ann)

exp.monocle<-newCellDataSet(exp.matrix,phenoData =exp_pd,featureData =exp_fd,expressionFamily=negbinomial.size())
head(pData(exp.monocle))
head(fData(exp.monocle))

exp.monocle <- estimateSizeFactors(exp.monocle)
exp.monocle <- estimateDispersions(exp.monocle)

diff_test_res<-differentialGeneTest(exp.monocle,fullModelFormulaStr = "~seurat_clusters") 
ordering_genes<-row.names (subset(diff_test_res, qval < 0.01))
exp.monocle<-setOrderingFilter(exp.monocle, ordering_genes)
plot_ordering_genes(exp.monocle)

exp.monocle<-reduceDimension(exp.monocle, max_components = 2, reduction_method = "DDRTree")
exp.monocle<-orderCells(exp.monocle)
colnames(pData(exp.monocle))

names(pData(exp.monocle))[names(pData(exp.monocle))=="seurat_clusters"]="Cluster"

plot1<-plot_cell_trajectory(exp.monocle, color_by = "orig.ident",cell_size=2)

p4_1 <- plot1 +
  labs(tag = "H")
p4_1


p <- (p1_1 | (p1_2 | p1_3) ) /
  (p2_1| p2_2|p3_2) /
  (p3_1|p4_1)
p
ggsave("./ladc_result/ladc_figure1.pdf", plot = p, width = 15, height = 18)
ggsave("./ladc_result/4_ladc_figure1.tiff", plot = p, width = 15, height = 18)

####Figure2####
#heatmap
top500 = CaseMatch(search = as.vector(DEG_500$gene), match = rownames(ladc)) 

figure2_1 <- DoHeatmap(ladc, features = top500, group.by = "orig.ident")+
  theme(axis.text.y.left = element_blank())+
  labs(tag = "A")
figure2_1

#GO
go_ladc <- read.table("./ladc_result/DEG_all.txt",sep=" ")
go_ladc <- t(go_ladc)
keytypes(org.Hs.eg.db)

go_ladc_id_trance <- bitr(go_ladc, fromType = "SYMBOL", toType = "ENSEMBL",
                          OrgDb = "org.Hs.eg.db",drop=T)
write.table(go_ladc_id_trance$ENSEMBL,"./ladc_result/go_ladc_id_trance.txt", row.names = F,col.names = F)

f <- read.table("./ladc_result/go_ladc_id_trance.txt")
f <- f[c(1)]

EG2Ensembl = toTable(org.Hs.egENSEMBL)
f = f$V1
geneLists = data.frame(ensembl_id = f)
results = merge(geneLists,EG2Ensembl,by='ensembl_id',all.x  =T)
id =na.omit(results$gene_id)

All <- enrichGO(OrgDb="org.Hs.eg.db",gene=id,ont="ALL",readable=T)

dotplot(All,showCategory=10,title="Enrichment Go")+
barplot(All, showCategory=15,title="EnrichmentGO") 
figure2_2 <- dotplot(All,split="ONTOLOGY",title ="Enrichment GO Top5",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")+
  labs(tag="B")
figure2_2


#KEGG
KEGG <- enrichKEGG(gene= id, organism  = 'hsa', pvalueCutoff = 0.05)
dotplot(KEGG,font.size=12)
barplot(KEGG,font.size=8) 
dotplot(KEGG,showCategory=5,title="Enrichment KEGG Top5")
figure2_3<-dotplot(KEGG, font.size=12, showCategory=5, title="Enrichment KEGG Top5")+
       scale_size(rang=c(5.20))+
       labs(tag="C")
figure2_3



fig2 <- figure2_1|figure2_2/figure2_3
fig2
ggsave("./ladc_result/fig_2.pdf", plot = fig2, width = 15, height = 18)
ggsave("./ladc_result/fig_2.tiff", plot = fig2, width = 15, height = 18)


####Figure4####
p4_1 <- DimPlot(ladc, reduction = "tsne",label=T,label.size = 8,pt.size = 2)  +
  theme(legend.position = "none") +
  labs(tag = "A")
p4_2 <- FeaturePlot(ladc,features = "CKAP4",reduction = "tsne")+
  labs(tag = "B")
p4_3 <- FeaturePlot(ladc,features = "SERPINA1",reduction = "tsne")+
  labs(tag = "C")
p4_4 <- FeaturePlot(ladc,features = "SDC2",reduction = "tsne")+
  labs(tag = "D")
p4_5 <- FeaturePlot(ladc,features = "GNG11",reduction = "tsne")+
  labs(tag = "E")

p4 <- CombinePlots(plots = list(p4_1,p4_2,p4_3,p4_4,p4_5))
p4
ggsave("./ladc_result/figure4_new.pdf", plot = p4_new, width = 18, height = 15,dpi=600)


