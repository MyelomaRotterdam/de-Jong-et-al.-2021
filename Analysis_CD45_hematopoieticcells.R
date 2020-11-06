## Loading libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

# Loading and pre-processing dataset 1
object_1 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_1 <- subset(x = object_1, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_1 <- read.csv(file="~/Barcodes_PH1.csv", header=T, row.names = 1)
select.cells_object_1 <- as.character(select.cells_object_1$x)
select.cells_object_1_im <- object_1@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_1]
object_1 <- subset(x = object_1, cells = select.cells_object_1_im)
object_1 <- NormalizeData(object = object_1, normalization.method = "LogNormalize", scale.factor = 1e4)
object_1 <- FindVariableFeatures(object = object_1, selection.method = "vst", nfeatures = 2000)
object_1[["state"]] <- "PH1"
object_1[["source"]] <- "myeloma"
object_1[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 2
object_2 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_2 <- CreateSeuratObject(counts = object_2, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_2), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_2, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_2, slot = "counts"))
object_2[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_2 <- subset(x = object_2, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
object_2 <- NormalizeData(object = object_2, normalization.method = "LogNormalize", scale.factor = 1e4)
object_2 <- FindVariableFeatures(object = object_2, selection.method = "vst", nfeatures = 2000)
object_2[["state"]] <- "PH2"
object_2[["source"]] <- "myeloma"
object_2[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 3
object_3 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_3 <- CreateSeuratObject(counts = object_3, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_3), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_3, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_3, slot = "counts"))
object_3[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_3 <- subset(x = object_3, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_3 <- read.csv(file="~/Barcodes_PH3.csv", header=T, row.names = 1)
select.cells_object_3 <- as.character(select.cells_object_3$x)
select.cells_object_1_im <- object_3@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_3]
object_3 <- subset(x = object_3, cells = select.cells_object_1_im)
object_3 <- NormalizeData(object = object_3, normalization.method = "LogNormalize", scale.factor = 1e4)
object_3 <- FindVariableFeatures(object = object_3, selection.method = "vst", nfeatures = 2000)
object_3[["state"]] <- "PH3"
object_3[["source"]] <- "myeloma"
object_3[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 4
object_4 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_4 <- CreateSeuratObject(counts = object_4, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_4), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_4, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_4, slot = "counts"))
object_4[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_4 <- subset(x = object_4, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_4 <- read.csv(file="~/Barcodes_PH4.csv", header=T, row.names = 1)
select.cells_object_4 <- as.character(select.cells_object_4$x)
select.cells_object_1_im <- object_4@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_4]
object_4 <- subset(x = object_4, cells = select.cells_object_1_im)
object_4 <- NormalizeData(object = object_4, normalization.method = "LogNormalize", scale.factor = 1e4)
object_4 <- FindVariableFeatures(object = object_4, selection.method = "vst", nfeatures = 2000)
object_4[["state"]] <- "PH4"
object_4[["source"]] <- "myeloma"
object_4[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 5
object_5 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_5 <- CreateSeuratObject(counts = object_5, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_5), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_5, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_5, slot = "counts"))
object_5[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_5, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_5 <- subset(x = object_5, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <5000 & percent.mito <0.1)
select.cells_object_5 <- read.csv(file="~/Barcodes_PH5.csv", header=T, row.names = 1)
select.cells_object_5 <- as.character(select.cells_object_5$x)
select.cells_object_1_im <- object_5@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_5]
object_5 <- subset(x = object_5, cells = select.cells_object_1_im)
object_5 <- NormalizeData(object = object_5, normalization.method = "LogNormalize", scale.factor = 1e4)
object_5 <- FindVariableFeatures(object = object_5, selection.method = "vst", nfeatures = 2000)
object_5[["state"]] <- "PH5"
object_5[["source"]] <- "myeloma"
object_5[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 6
object_6 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_6 <- CreateSeuratObject(counts = object_6, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_6), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_6, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_6, slot = "counts"))
object_6[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_6, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_6 <- subset(x = object_6, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
select.cells_object_6 <- read.csv(file="~/Barcodes_PH6.csv", header=T, row.names = 1)
select.cells_object_6 <- as.character(select.cells_object_6$x)
select.cells_object_1_im <- object_6@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_6]
object_6 <- subset(x = object_6, cells = select.cells_object_1_im)
object_6 <- NormalizeData(object = object_6, normalization.method = "LogNormalize", scale.factor = 1e4)
object_6 <- FindVariableFeatures(object = object_6, selection.method = "vst", nfeatures = 2000)
object_6[["state"]] <- "PH6"
object_6[["source"]] <- "myeloma"
object_6[["genetics"]] <- "hyperdiploid"

# Loading and pre-processing dataset 7
object_7 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_7 <- CreateSeuratObject(counts = object_7, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_7), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_7, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_7, slot = "counts"))
object_7[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_7, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_7 <- subset(x = object_7, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_7 <- read.csv(file="~/Barcodes_PT1.csv", header=T, row.names = 1)
select.cells_object_7 <- as.character(select.cells_object_7$x)
select.cells_object_1_im <- object_7@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_7]
object_7 <- subset(x = object_7, cells = select.cells_object_1_im)
object_7 <- NormalizeData(object = object_7, normalization.method = "LogNormalize", scale.factor = 1e4)
object_7 <- FindVariableFeatures(object = object_7, selection.method = "vst", nfeatures = 2000)
object_7[["state"]] <- "PT1"
object_7[["source"]] <- "myeloma"
object_7[["genetics"]] <- "translocation"

# Loading and pre-processing dataset 8
object_8 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_8 <- CreateSeuratObject(counts = object_8, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_8), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_8, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_8, slot = "counts"))
object_8[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_8, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_8 <- subset(x = object_8, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
select.cells_object_8 <- read.csv(file="~/Barcodes_PT2.csv", header=T, row.names = 1)
select.cells_object_8 <- as.character(select.cells_object_8$x)
select.cells_object_1_im <- object_8@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_8]
object_8 <- subset(x = object_8, cells = select.cells_object_1_im)
object_8 <- NormalizeData(object = object_8, normalization.method = "LogNormalize", scale.factor = 1e4)
object_8 <- FindVariableFeatures(object = object_8, selection.method = "vst", nfeatures = 2000)
object_8[["state"]] <- "PT2"
object_8[["source"]] <- "myeloma"
object_8[["genetics"]] <- "translocation"

# Loading and pre-processing dataset 9
object_9 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_9 <- CreateSeuratObject(counts = object_9, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_9), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_9, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_9, slot = "counts"))
object_9[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_9, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_9 <- subset(x = object_9, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <11000 & percent.mito <0.1)
select.cells_object_9 <- read.csv(file="~/Barcodes_PT3.csv", header=T, row.names = 1)
select.cells_object_9 <- as.character(select.cells_object_9$x)
select.cells_object_1_im <- object_9@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_9]
object_9 <- subset(x = object_9, cells = select.cells_object_1_im)
object_9 <- NormalizeData(object = object_9, normalization.method = "LogNormalize", scale.factor = 1e4)
object_9 <- FindVariableFeatures(object = object_9, selection.method = "vst", nfeatures = 2000)
object_9[["state"]] <- "PT3"
object_9[["source"]] <- "myeloma"
object_9[["genetics"]] <- "translocation"

# Loading and pre-processing dataset 10
object_10 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_10 <- CreateSeuratObject(counts = object_10, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_10), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_10, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_10, slot = "counts"))
object_10[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_10, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_10 <- subset(x = object_10, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_10 <- read.csv(file="~/Barcodes_PT4.csv", header=T, row.names = 1)
select.cells_object_10 <- as.character(select.cells_object_10$x)
select.cells_object_1_im <- object_10@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_10]
object_10 <- subset(x = object_10, cells = select.cells_object_1_im)
object_10 <- NormalizeData(object = object_10, normalization.method = "LogNormalize", scale.factor = 1e4)
object_10 <- FindVariableFeatures(object = object_10, selection.method = "vst", nfeatures = 2000)
object_10[["state"]] <- "PT4"
object_10[["source"]] <- "myeloma"
object_10[["genetics"]] <- "del17p"

# Loading and pre-processing dataset 11
object_11 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_11 <- CreateSeuratObject(counts = object_11, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_11), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_11, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_11, slot = "counts"))
object_11[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_11, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_11 <- subset(x = object_11, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_11 <- read.csv(file="~/Barcodes_PT5.csv", header=T, row.names = 1)
select.cells_object_11 <- as.character(select.cells_object_11$x)
select.cells_object_1_im <- object_11@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_11]
object_11 <- subset(x = object_11, cells = select.cells_object_1_im)
object_11 <- NormalizeData(object = object_11, normalization.method = "LogNormalize", scale.factor = 1e4)
object_11 <- FindVariableFeatures(object = object_11, selection.method = "vst", nfeatures = 2000)
object_11[["state"]] <- "PT5"
object_11[["source"]] <- "myeloma"
object_11[["genetics"]] <- "translocation"

# Loading and pre-processing dataset 12
object_12 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_12 <- CreateSeuratObject(counts = object_12, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_12), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_12, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_12, slot = "counts"))
object_12[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_12, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_12 <- subset(x = object_12, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <11000 & percent.mito <0.1)
select.cells_object_12 <- read.csv(file="~/Barcodes_PT7.csv", header=T, row.names = 1)
select.cells_object_12 <- as.character(select.cells_object_12$x)
select.cells_object_1_im <- object_12@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_12]
object_12 <- subset(x = object_12, cells = select.cells_object_1_im)
object_12 <- NormalizeData(object = object_12, normalization.method = "LogNormalize", scale.factor = 1e4)
object_12 <- FindVariableFeatures(object = object_12, selection.method = "vst", nfeatures = 2000)
object_12[["state"]] <- "PT7"
object_12[["source"]] <- "myeloma"
object_12[["genetics"]] <- "translocation"

# Loading and pre-processing dataset 13
object_13 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_13 <- CreateSeuratObject(counts = object_13, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_13), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_13, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_13, slot = "counts"))
object_13[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_13, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_13 <- subset(x = object_13, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_13 <- read.csv(file="~/Barcodes_PT6.csv", header=T, row.names = 1)
select.cells_object_13 <- as.character(select.cells_object_13$x)
select.cells_object_1_im <- object_13@assays$RNA@counts@Dimnames[2]
select.cells_object_1_im <- unlist(select.cells_object_1_im)
select.cells_object_1_im <- select.cells_object_1_im[!select.cells_object_1_im %in% select.cells_object_13]
object_13 <- subset(x = object_13, cells = select.cells_object_1_im)
object_13 <- NormalizeData(object = object_13, normalization.method = "LogNormalize", scale.factor = 1e4)
object_13 <- FindVariableFeatures(object = object_13, selection.method = "vst", nfeatures = 2000)
object_13[["state"]] <- "PT6"
object_13[["source"]] <- "myeloma"
object_13[["genetics"]] <- "del17p"

# Loading and pre-processing dataset 14
object_14 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_14 <- CreateSeuratObject(counts = object_14, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_14), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_14, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_14, slot = "counts"))
object_14[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_14, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_14 <- subset(x = object_14, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
object_14 <- NormalizeData(object = object_14, normalization.method = "LogNormalize", scale.factor = 1e4)
object_14 <- FindVariableFeatures(object = object_14, selection.method = "vst", nfeatures = 2000)
object_14[["state"]] <- "CBM12"
object_14[["source"]] <- "control"
object_14[["genetics"]] <- "NA"

# Loading and pre-processing dataset 15
object_15 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_15 <- CreateSeuratObject(counts = object_15, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_15), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_15, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_15, slot = "counts"))
object_15[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_15, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_15 <- subset(x = object_15, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
object_15 <- NormalizeData(object = object_15, normalization.method = "LogNormalize", scale.factor = 1e4)
object_15 <- FindVariableFeatures(object = object_15, selection.method = "vst", nfeatures = 2000)
object_15[["state"]] <- "CBM13"
object_15[["source"]] <- "control"
object_15[["genetics"]] <- "NA"

# Loading and pre-processing dataset 16
object_16 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_16 <- CreateSeuratObject(counts = object_16, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_16), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_16, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_16, slot = "counts"))
object_16[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_16, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_16 <- subset(x = object_16, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <7000 & percent.mito <0.1)
object_16 <- NormalizeData(object = object_16, normalization.method = "LogNormalize", scale.factor = 1e4)
object_16 <- FindVariableFeatures(object = object_16, selection.method = "vst", nfeatures = 2000)
object_16[["state"]] <- "CBM14"
object_16[["source"]] <- "control"
object_16[["genetics"]] <- "NA"

# Loading and pre-processing dataset 17
object_17 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_17 <- CreateSeuratObject(counts = object_17, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_17), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_17, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_17, slot = "counts"))
object_17[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_17, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_17 <- subset(x = object_17, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
object_17 <- NormalizeData(object = object_17, normalization.method = "LogNormalize", scale.factor = 1e4)
object_17 <- FindVariableFeatures(object = object_17, selection.method = "vst", nfeatures = 2000)
object_17[["state"]] <- "CBM15"
object_17[["source"]] <- "control"
object_17[["genetics"]] <- "NA"

# Loading and pre-processing dataset 18
object_18 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_18 <- CreateSeuratObject(counts = object_18, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_18), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_18, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_18, slot = "counts"))
object_18[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_18, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_18 <- subset(x = object_18, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
object_18 <- NormalizeData(object = object_18, normalization.method = "LogNormalize", scale.factor = 1e4)
object_18 <- FindVariableFeatures(object = object_18, selection.method = "vst", nfeatures = 2000)
object_18[["state"]] <- "CBM16"
object_18[["source"]] <- "control"
object_18[["genetics"]] <- "NA"

# Identification of integration anchors
myeloma.list <- c(object_1, object_2, object_3, object_4, object_5, object_6, object_7, 
                  object_8, object_9, object_10, object_11, object_12, object_13)
myeloma.anchors <- FindIntegrationAnchors(object.list = myeloma.list, dims = 1:30)
myeloma.integrated <- IntegrateData(anchorset = myeloma.anchors, dims = 1:30)

# Identification of integration anchors
myeloma.list <- c(object_14, object_15, object_16, object_17, object_18)
control.anchors <- FindIntegrationAnchors(object.list = myeloma.list, dims = 1:30)
control.integrated <- IntegrateData(anchorset = control.anchors, dims = 1:30)

# Identification of integration anchors
myeloma.list <- c(myeloma.integrated, control.integrated)
total.anchors <- FindIntegrationAnchors(object.list = myeloma.list, dims = 1:30)
total.integrated <- IntegrateData(anchorset = total.anchors, dims = 1:30)

# Post-processing merged data
DefaultAssay(object = total.integrated) <- "integrated"
total.integrated <- ScaleData(object = total.integrated, verbose = F)
total.integrated <- RunPCA(object = total.integrated, verbose = F)
ElbowPlot(object = total.integrated, ndims = 40)
total.integrated <- RunUMAP(object = total.integrated, reduction = "pca", dims = 1:20)
total.integrated <- FindNeighbors(object = total.integrated, dims = 1:20)
total.integrated <- FindClusters(object = total.integrated, resolution = 0.6)

# Saving and loading RDS files
saveRDS(total.integrated, file="Analyses for myself/RDS files/CD45_immune.rds")
myeloma.integrated <- readRDS("Analyses for myself/RDS files/CD45_immune.rds")