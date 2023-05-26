#Seurat官网流程#


library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggplot2)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/kabuda/Kabuda/R/SelectK/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#过滤质控
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##规范化数据
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#特征选择
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Scale
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))



#暂存一下pbmc文件
pbmc_save <- pbmc
pbmc <- pbmc_save



####clustree####
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc,algorithm = 4, resolution = c(0.1,0.5,1,2,4))
clustree(pbmc@meta.data, prefix = "RNA_snn_res.")

#把簇存起来
cluster<- cbind(pbmc@meta.data$RNA_snn_res.0.1,pbmc@meta.data$RNA_snn_res.0.5,pbmc@meta.data$RNA_snn_res.1,pbmc@meta.data$RNA_snn_res.2,pbmc@meta.data$RNA_snn_res.4,pbmc@meta.data$RNA_snn_res.8)
cluster <- data.frame(cluster)

#从Seurat中把矩阵提取出来
silhouette_data <- pbmc@reductions$pca@cell.embeddings
silhouette_data <- data.frame(silhouette_data)
silhouette_data <- cbind(silhouette_data,cluster)
save(silhouette_data,file="./silhouette_data.Rdata")
#如有必要 保存该文件
write.csv(data,file = "silhouettedata.csv",row.names = T)

pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



####读取簇的轮廓指数/以及用轮廓指数给图上色####
a<- readRDS("~/Desktop/resultspmbc/CellSilhouette_resolution_1.5_PC_10.rds")

pbmc[["Silhouette.Score"]]<-a[,3]
FeaturePlot(pbmc,features = "Silhouette.Score",reduction = "umap",pt.size=0.1,blend=FALSE,cols=c("#B4464B","#B4AF46","#4682B4")) 






