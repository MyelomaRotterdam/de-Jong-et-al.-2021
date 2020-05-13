## Analysis hematopoietic niche of CD38pos cells in 10 Multiple Myeloma samples
## Madelon de Jong, Myeloma Research Rotterdam, 2019-2020
## Part of submitted work

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

# Loading and pre-processing dataset 1
object_1 <- Read10X(data.dir = "~/PH1/filtered_feature_bc_matrix/")
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_1 <- subset(x = object_1, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "PH1"
object_1[["source"]] <- "myeloma"
object_1[["genetics"]] <- "hyperdiploid"
object_1[["kit"]] <- "kit1"

# Loading and pre-processing dataset 2. These cells were not combined with niche cells during work-up but ran as a single run.
object_2 <- Read10X(data.dir = "~/PH2/filtered_feature_bc_matrix/")
object_2 <- CreateSeuratObject(counts = object_2, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_2), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_2, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_2, slot = "counts"))
object_2[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_2 <- subset(x = object_2, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_2 <- NormalizeData(object = object_2, normalization.method = "LogNormalize", scale.factor = 1e4)
object_2 <- FindVariableFeatures(object = object_2, selection.method = "vst", nfeatures = 2000)
object_2[["state"]] <- "PH2"
object_2[["source"]] <- "myeloma"
object_2[["genetics"]] <- "hyperdiploid"
object_2[["kit"]] <- "kit1"

# Loading and pre-processing dataset 3
object_3 <- Read10X(data.dir = "~/PH3/filtered_feature_bc_matrix/")
object_3 <- CreateSeuratObject(counts = object_3, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_3), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_3, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_3, slot = "counts"))
object_3[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_3 <- subset(x = object_3, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
object_3 <- NormalizeData(object = object_3, normalization.method = "LogNormalize", scale.factor = 1e4)
object_3 <- FindVariableFeatures(object = object_3, selection.method = "vst", nfeatures = 2000)
object_3[["state"]] <- "PH3"
object_3[["source"]] <- "myeloma"
object_3[["genetics"]] <- "hyperdiploid"
object_3[["kit"]] <- "kit1"

# Loading and pre-processing dataset 4
object_4 <- Read10X(data.dir = "~/PH4/filtered_feature_bc_matrix/")
object_4 <- CreateSeuratObject(counts = object_4, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_4), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_4, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_4, slot = "counts"))
object_4[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_4 <- subset(x = object_4, subset = nFeature_RNA > 200 & nFeature_RNA <4000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_4 <- NormalizeData(object = object_4, normalization.method = "LogNormalize", scale.factor = 1e4)
object_4 <- FindVariableFeatures(object = object_4, selection.method = "vst", nfeatures = 2000)
object_4[["state"]] <- "PH4"
object_4[["source"]] <- "myeloma"
object_4[["genetics"]] <- "hyperdiploid"
object_4[["kit"]] <- "kit1"

# Loading and pre-processing dataset 5
object_5 <- Read10X(data.dir = "~/PH5/filtered_feature_bc_matrix/")
object_5 <- CreateSeuratObject(counts = object_5, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_5), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_5, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_5, slot = "counts"))
object_5[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_5, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_5 <- subset(x = object_5, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
object_5 <- NormalizeData(object = object_5, normalization.method = "LogNormalize", scale.factor = 1e4)
object_5 <- FindVariableFeatures(object = object_5, selection.method = "vst", nfeatures = 2000)
object_5[["state"]] <- "PH5"
object_5[["source"]] <- "myeloma"
object_5[["genetics"]] <- "hyperdiploid"
object_5[["kit"]] <- "kit2"

# Loading and pre-processing dataset 6
object_6 <- Read10X(data.dir = "~/PH6/filtered_feature_bc_matrix/")
object_6 <- CreateSeuratObject(counts = object_6, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_6), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_6, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_6, slot = "counts"))
object_6[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_6, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_6 <- subset(x = object_6, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
object_6 <- NormalizeData(object = object_6, normalization.method = "LogNormalize", scale.factor = 1e4)
object_6 <- FindVariableFeatures(object = object_6, selection.method = "vst", nfeatures = 2000)
object_6[["state"]] <- "PH6"
object_6[["source"]] <- "myeloma"
object_6[["genetics"]] <- "hyperdiploid"
object_6[["kit"]] <- "kit2"

# Loading and pre-processing dataset 7
object_7 <- Read10X(data.dir = "~/PT1//filtered_feature_bc_matrix/")
object_7 <- CreateSeuratObject(counts = object_7, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_7), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_7, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_7, slot = "counts"))
object_7[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_7, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_7 <- subset(x = object_7, subset = nFeature_RNA > 200 & nFeature_RNA <1000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
object_7 <- NormalizeData(object = object_7, normalization.method = "LogNormalize", scale.factor = 1e4)
object_7 <- FindVariableFeatures(object = object_7, selection.method = "vst", nfeatures = 2000)
object_7[["state"]] <- "PT1"
object_7[["source"]] <- "myeloma"
object_7[["genetics"]] <- "translocation"
object_7[["kit"]] <- "kit2"

# Loading and pre-processing dataset 8
object_8 <- Read10X(data.dir = "~/PT2/filtered_feature_bc_matrix/")
object_8 <- CreateSeuratObject(counts = object_8, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_8), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_8, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_8, slot = "counts"))
object_8[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_8, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_8 <- subset(x = object_8, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_8 <- NormalizeData(object = object_8, normalization.method = "LogNormalize", scale.factor = 1e4)
object_8 <- FindVariableFeatures(object = object_8, selection.method = "vst", nfeatures = 2000)
object_8[["state"]] <- "PT2"
object_8[["source"]] <- "myeloma"
object_8[["genetics"]] <- "translocation"
object_8[["kit"]] <- "kit2"

# Loading and pre-processing dataset 9
object_9 <- Read10X(data.dir = "~/PT3/filtered_feature_bc_matrix/")
object_9 <- CreateSeuratObject(counts = object_9, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_9), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_9, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_9, slot = "counts"))
object_9[["percent.mito"]] <- percent.mito_object
object_9 <- subset(x = object_9, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
object_9 <- NormalizeData(object = object_9, normalization.method = "LogNormalize", scale.factor = 1e4)
object_9 <- FindVariableFeatures(object = object_9, selection.method = "vst", nfeatures = 2000)
object_9[["state"]] <- "PT3"
object_9[["source"]] <- "myeloma"
object_9[["genetics"]] <- "translocation"
object_9[["kit"]] <- "kit2"

# Loading and pre-processing dataset 10
object_10 <- Read10X(data.dir = "~/PT4/filtered_feature_bc_matrix/")
object_10 <- CreateSeuratObject(counts = object_10, min.cells = 3, min.features = 200, project = "myeloma")
mito.features_object <- grep(pattern = "^MT-", x=rownames(x=object_10), value=T)
percent.mito_object <- Matrix::colSums(x = GetAssayData(object = object_10, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object_10, slot = "counts"))
object_10[["percent.mito"]] <- percent.mito_object
VlnPlot(object = object_10, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_10 <- subset(x = object_10, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <30000 & percent.mito <0.1)
object_10 <- NormalizeData(object = object_10, normalization.method = "LogNormalize", scale.factor = 1e4)
object_10 <- FindVariableFeatures(object = object_10, selection.method = "vst", nfeatures = 2000)
object_10[["state"]] <- "PT4"
object_10[["source"]] <- "myeloma"
object_10[["genetics"]] <- "del17p"
object_10[["kit"]] <- "kit2"

# Identification of integration anchors
reference.list <- c(object_1, object_2, object_3, object_4, object_5,
                    object_6, object_7, object_8, object_9, object_10)
niche.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:40)
niche.integrated <- IntegrateData(anchorset = niche.anchors, dims = 1:40)

# Post-processing merged data
DefaultAssay(object = niche.integrated) <- "integrated"
niche.integrated <- ScaleData(object = niche.integrated, verbose = F)
niche.integrated <- RunPCA(object = niche.integrated, verbose = F)
ElbowPlot(object = niche.integrated, ndims = 40)
niche.integrated <- RunTSNE(object = niche.integrated, reduction = "pca", dims = 1:30)
niche.integrated <- RunUMAP(object = niche.integrated, reduction = "pca", dims = 1:30)
niche.integrated <- FindNeighbors(object = niche.integrated, dims = 1:30)
niche.integrated <- FindClusters(object = niche.integrated, resolution = 0.3)

niche.integrated <- subset(niche.integrated, idents = c(0:2, 4:7, 9)) 
# Subset clusters that are not plasma cells based on SDC1 expression. Run clustering again
niche.integrated <- ScaleData(object = niche.integrated, verbose = F)
niche.integrated <- RunPCA(object = niche.integrated, verbose = F)
ElbowPlot(object = niche.integrated, ndims = 40)
niche.integrated <- RunTSNE(object = niche.integrated, reduction = "pca", dims = 1:30)
niche.integrated <- RunUMAP(object = niche.integrated, reduction = "pca", dims = 1:30)
niche.integrated <- FindNeighbors(object = niche.integrated, dims = 1:30)
niche.integrated <- FindClusters(object = niche.integrated, resolution = 0.3)

# Basic visualization
niche.integrated@active.ident <- factor(x=niche.integrated@active.ident, levels = c(2,4,8,10,6,0,13,5,1,3,9,14,11,7,12))

DimPlot(object = niche.integrated, reduction = "umap", pt.size = 1, label = T, 
        cols=c("#56B4E9", "#316684", "#AA811B", 
               "#2F7C6D", "#D2691E", "#E0E0DF", 
               "#DF9560", "#74BAAD", "#DAA520", 
               "#3E89A0", "#C2AC74", "#787878", 
               "#97D5CA", "#999999", "#595959"))
DimPlot(object = niche.integrated, reduction = "umap", pt.size = 1, label = T, split.by = "source", 
        cols=c("#56B4E9", "#316684", "#AA811B", 
               "#2F7C6D", "#D2691E", "#E0E0DF", 
               "#DF9560", "#74BAAD", "#DAA520", 
               "#3E89A0", "#C2AC74", "#787878", 
               "#97D5CA", "#999999", "#595959"))
DimPlot(object = niche.integrated, reduction = "umap", pt.size = 1, label = T, split.by = "state", ncol = 4, 
        cols=c("#56B4E9", "#316684", "#AA811B", 
               "#2F7C6D", "#D2691E", "#E0E0DF", 
               "#DF9560", "#74BAAD", "#DAA520", 
               "#3E89A0", "#C2AC74", "#787878", 
               "#97D5CA", "#999999", "#595959"))
FeaturePlot(object = niche.integrated, features = c("CD3D"), pt.size=1.5, reduction = "umap", 
            split.by = "source", label = F, sort.cell=T, cols=c("lightgrey", "brown"), min.cutoff = 0)
VlnPlot(object = niche.integrated, features = c("CD4"), pt.size = -1, assay="RNA",
        cols=c("#56B4E9", "#316684", "#AA811B", 
               "#2F7C6D", "#D2691E", "#E0E0DF", 
               "#DF9560", "#74BAAD", "#DAA520", 
               "#3E89A0", "#C2AC74", "#787878", 
               "#97D5CA", "#999999", "#595959")) + 
  stat_summary(fun.y = median, geom = "point", size=0.5, color="black") + NoLegend() 