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
select.cells_object_1 <- paste0(select.cells_object_1, "-1")
object_1 <- subset(x = object_1, cells = select.cells_object_1)
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
select.cells_object_2 <- read.csv(file="~/Barcodes_PH2.csv", header=T, row.names = 1)
select.cells_object_2 <- as.character(select.cells_object_2$x)
select.cells_object_2 <- paste0(select.cells_object_2, "-1")
object_2 <- subset(x = object_2, cells = select.cells_object_2)
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
object_3 <- subset(x = object_3, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <12000 & percent.mito <0.1)
select.cells_object_3 <- read.csv(file="~/Barcodes_PH3.csv", header=T, row.names = 1)
select.cells_object_3 <- as.character(select.cells_object_3$x)
select.cells_object_3 <- paste0(select.cells_object_3, "-1")
object_3 <- subset(x = object_3, cells = select.cells_object_3)
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
object_4 <- subset(x = object_4, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
select.cells_object_4 <- read.csv(file="~/Barcodes_PH4.csv", header=T, row.names = 1)
select.cells_object_4 <- as.character(select.cells_object_4$x)
select.cells_object_4 <- paste0(select.cells_object_4, "-1")
object_4 <- subset(x = object_4, cells = select.cells_object_4)
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
select.cells_object_5 <- paste0(select.cells_object_5, "-1")
object_5 <- subset(x = object_5, cells = select.cells_object_5)
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
object_6 <- subset(x = object_6, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <6000 & percent.mito <0.1)
select.cells_object_6 <- read.csv(file="~/Barcodes_PH6.csv", header=T, row.names = 1)
select.cells_object_6 <- as.character(select.cells_object_6$x)
select.cells_object_6 <- paste0(select.cells_object_6, "-1")
object_6 <- subset(x = object_6, cells = select.cells_object_6)
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
select.cells_object_7 <- paste0(select.cells_object_7, "-1")
object_7 <- subset(x = object_7, cells = select.cells_object_7)
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
object_8 <- subset(x = object_8, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_8 <- read.csv(file="~/Barcodes_PT2.csv", header=T, row.names = 1)
select.cells_object_8 <- as.character(select.cells_object_8$x)
select.cells_object_8 <- paste0(select.cells_object_8, "-1")
object_8 <- subset(x = object_8, cells = select.cells_object_8)
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
object_9 <- subset(x = object_9, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_9 <- read.csv(file="~/Barcodes_PT3.csv", header=T, row.names = 1)
select.cells_object_9 <- as.character(select.cells_object_9$x)
select.cells_object_9 <- paste0(select.cells_object_9, "-1")
object_9 <- subset(x = object_9, cells = select.cells_object_9)
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
object_10 <- subset(x = object_10, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
select.cells_object_10 <- read.csv(file="~/Barcodes_PT4.csv", header=T, row.names = 1)
select.cells_object_10 <- as.character(select.cells_object_10$x)
select.cells_object_10 <- paste0(select.cells_object_10, "-1")
object_10 <- subset(x = object_10, cells = select.cells_object_10)
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
object_11 <- subset(x = object_11, cells = select.cells_object_11)
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
object_12 <- subset(x = object_12, cells = select.cells_object_12)
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
object_13 <- subset(x = object_13, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_13 <- read.csv(file="~/Barcodes_PT6.csv", header=T, row.names = 1)
select.cells_object_13 <- as.character(select.cells_object_13$x)
object_13 <- subset(x = object_13, cells = select.cells_object_13)
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
object_14 <- subset(x = object_14, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <5000 & percent.mito <0.1)
select.cells_object_14 <- read.csv(file="~/Barcodes_CBM1.csv", header=T, row.names = 1)
select.cells_object_14 <- as.character(select.cells_object_14$x)
select.cells_object_14 <- paste0(select.cells_object_14, "-1")
object_14 <- subset(x = object_14, cells = select.cells_object_14)
object_14 <- NormalizeData(object = object_14, normalization.method = "LogNormalize", scale.factor = 1e4)
object_14 <- FindVariableFeatures(object = object_14, selection.method = "vst", nfeatures = 2000)
object_14[["state"]] <- "CBM1"
object_14[["source"]] <- "control"
object_14[["genetics"]] <- "NA"

# Loading and pre-processing dataset 15
object_15 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_15 <- CreateSeuratObject(counts = object_15, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_15), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_15, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_15, slot = "counts"))
object_15[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_15, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_15 <- subset(x = object_15, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <6000 & percent.mito <0.1)
select.cells_object_15 <- read.csv(file="~/Barcodes_CBM2.csv", header=T, row.names = 1)
select.cells_object_15 <- as.character(select.cells_object_15$x)
select.cells_object_15 <- paste0(select.cells_object_15, "-1")
object_15 <- subset(x = object_15, cells = select.cells_object_15)
object_15 <- NormalizeData(object = object_15, normalization.method = "LogNormalize", scale.factor = 1e4)
object_15 <- FindVariableFeatures(object = object_15, selection.method = "vst", nfeatures = 2000)
object_15[["state"]] <- "CBM2"
object_15[["source"]] <- "control"
object_15[["genetics"]] <- "NA"

# Loading and pre-processing dataset 16
object_16 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_16 <- CreateSeuratObject(counts = object_16, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_16), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_16, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_16, slot = "counts"))
object_16[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_16, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_16 <- subset(x = object_16, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <9000 & percent.mito <0.1)
select.cells_object_16 <- read.csv(file="~/Barcodes_CBM7.csv", header=T, row.names = 1)
select.cells_object_16 <- as.character(select.cells_object_16$x)
object_16 <- subset(x = object_16, cells = select.cells_object_16)
object_16 <- NormalizeData(object = object_16, normalization.method = "LogNormalize", scale.factor = 1e4)
object_16 <- FindVariableFeatures(object = object_16, selection.method = "vst", nfeatures = 2000)
object_16[["state"]] <- "CBM7"
object_16[["source"]] <- "control"
object_16[["genetics"]] <- "NA"

# Loading and pre-processing dataset 17
object_17 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_17 <- CreateSeuratObject(counts = object_17, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_17), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_17, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_17, slot = "counts"))
object_17[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_17, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_17 <- subset(x = object_17, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >200 & nCount_RNA <8000 & percent.mito <0.1)
select.cells_object_17 <- read.csv(file="~/Barcodes_CBM8.csv", header=T, row.names = 1)
select.cells_object_17 <- as.character(select.cells_object_17$x)
object_17 <- subset(x = object_17, cells = select.cells_object_17)
object_17 <- NormalizeData(object = object_17, normalization.method = "LogNormalize", scale.factor = 1e4)
object_17 <- FindVariableFeatures(object = object_17, selection.method = "vst", nfeatures = 2000)
object_17[["state"]] <- "CBM8"
object_17[["source"]] <- "control"
object_17[["genetics"]] <- "NA"

# Loading and pre-processing dataset 18
object_18 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_18 <- CreateSeuratObject(counts = object_18, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_18), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_18, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_18, slot = "counts"))
object_18[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_18, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_18 <- subset(x = object_18, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <6000 & percent.mito <0.1)
select.cells_object_18 <- read.csv(file="~/Barcodes_CBM9.csv", header=T, row.names = 1)
select.cells_object_18 <- as.character(select.cells_object_18$x)
object_18 <- subset(x = object_18, cells = select.cells_object_18)
object_18 <- NormalizeData(object = object_18, normalization.method = "LogNormalize", scale.factor = 1e4)
object_18 <- FindVariableFeatures(object = object_18, selection.method = "vst", nfeatures = 2000)
object_18[["state"]] <- "CBM9"
object_18[["source"]] <- "control"
object_18[["genetics"]] <- "NA"

# Loading and pre-processing dataset 19
object_19 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_19 <- CreateSeuratObject(counts = object_19, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_19), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_19, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_19, slot = "counts"))
object_19[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_19, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_19 <- subset(x = object_19, subset = nFeature_RNA > 200 & nFeature_RNA <3000 & nCount_RNA >200 & nCount_RNA <20000 & percent.mito <0.1)
select.cells_object_19 <- read.csv(file="~/Barcodes_CBM10.csv", header=T, row.names = 1)
select.cells_object_19 <- as.character(select.cells_object_19$x)
object_19 <- subset(x = object_19, cells = select.cells_object_19)
object_19 <- NormalizeData(object = object_19, normalization.method = "LogNormalize", scale.factor = 1e4)
object_19 <- FindVariableFeatures(object = object_19, selection.method = "vst", nfeatures = 2000)
object_19[["state"]] <- "CBM10"
object_19[["source"]] <- "control"
object_19[["genetics"]] <- "NA"

# Loading and pre-processing dataset 20
object_20 <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
object_20 <- CreateSeuratObject(counts = object_20, min.cells = 3, min.features = 200, project = "Niche")
mito.features_object2 <- grep(pattern = "^MT-", x=rownames(x=object_20), value=T)
percent.mito_object2 <- Matrix::colSums(x = GetAssayData(object = object_20, slot="counts")[mito.features_object2,]) / Matrix::colSums(x = GetAssayData(object = object_20, slot = "counts"))
object_20[["percent.mito"]] <- percent.mito_object2
VlnPlot(object = object_20, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
object_20 <- subset(x = object_20, subset = nFeature_RNA > 200 & nFeature_RNA <2000 & nCount_RNA >200 & nCount_RNA <10000 & percent.mito <0.1)
select.cells_object_20 <- read.csv(file="~/Barcodes_CBM11.csv", header=T, row.names = 1)
select.cells_object_20 <- as.character(select.cells_object_20$x)
object_20 <- subset(x = object_20, cells = select.cells_object_20)
object_20 <- NormalizeData(object = object_20, normalization.method = "LogNormalize", scale.factor = 1e4)
object_20 <- FindVariableFeatures(object = object_20, selection.method = "vst", nfeatures = 2000)
object_20[["state"]] <- "CBM11"
object_20[["source"]] <- "control"
object_20[["genetics"]] <- "NA"

# Identification of integration anchors
reference.list <- c(object_1, object_2, object_3, object_4, object_5, object_6, object_7, object_8, object_9, object_10,
                    object_11, object_12, object_13)
myeloma.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
myeloma.integrated <- IntegrateData(anchorset = myeloma.anchors, dims = 1:30)

# Identification of integration anchors
object_20 <- merge(object_20, object_19)
reference.list <- c(object_14, object_15, object_16, object_17, object_18, object_20)
control.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
control.integrated <- IntegrateData(anchorset = control.anchors, dims = 1:30)

# Identification of integration anchors
reference.list <- c(control.integrated, myeloma.integrated)
total.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
total.integrated <- IntegrateData(anchorset = total.anchors, dims = 1:30)
myeloma.integrated <- total.integrated

# Post-processing merged data
DefaultAssay(object = myeloma.integrated) <- "integrated"
myeloma.integrated <- ScaleData(object = myeloma.integrated, verbose = F)
myeloma.integrated <- RunPCA(object = myeloma.integrated, verbose = F)
ElbowPlot(object = myeloma.integrated, ndims = 40)
myeloma.integrated <- RunUMAP(object = myeloma.integrated, reduction = "pca", dims = 1:10, return.model = T)
myeloma.integrated <- FindNeighbors(object = myeloma.integrated, dims = 1:10)
myeloma.integrated <- FindClusters(object = myeloma.integrated, resolution = 0.3)
