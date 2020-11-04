## Pseudotime analysis using integrated Seurat object with Monocle 3
# Only useful for very homogeneous datasets, as it is based on the 2000 genes used for integration rather than the entire datamatrix
# Script is from user rwo012 on Github: https://github.com/satijalab/seurat/issues/1658 and user tlusardi on Github: https://github.com/cole-trapnell-lab/monocle-release/issues/388

library(Seurat)
library(monocle3)

#### Create a Monocle CDS Object
# Project PC dimensions to whole data set
seurat.object <- readRDS("Analyses for myself/RDS files/niche_merge_myway_noPTPRC.rds")
p <- DimPlot(seurat.object)
cells <- CellSelector(p)
seurat.object <- subset(seurat.object, cells = cells)
Idents(seurat.object) <- seurat.object$source
seurat.object <- subset(seurat.object, idents = "myeloma")
Idents(seurat.object) <- seurat.object$seurat_clusters

# Create an expression matrix
expression_matrix <- seurat.object@assays$RNA@counts

# Get cell metadata
cell_metadata <- seurat.object@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(seurat.object@assays$RNA), row.names = rownames(seurat.object@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

# Transfer Seurat embeddings
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
reducedDim(my.cds, type = "PCA") <- seurat.object@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- seurat.object@reductions$pca@stdev
plot_pc_variance_explained(my.cds)

# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- seurat.object@reductions$umap@cell.embeddings
plot_cells(my.cds)

# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- seurat.object@meta.data$gt_tp_cell_type_integrated_.0.9

my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-3)

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(seurat.object, reduction = "umap")
plot_cells(my.cds, color_cells_by = "partition", group_label_size = 3.5)
plot_cells(my.cds, color_cells_by = "seurat_clusters", show_trajectory_graph = F, group_label_size = 3.5)
my.cds = learn_graph(my.cds)
plot_cells(my.cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(my.cds, time_bin="6"){
  cell_ids <- which(colData(my.cds)[, "seurat_clusters"] == time_bin)
  
  closest_vertex <-
    my.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(my.cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(my.cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
my.cds <- order_cells(my.cds, root_pr_nodes=get_earliest_principal_node(my.cds))

my.cds = order_cells(my.cds, reduction_method = "UMAP")
plot_cells(my.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           trajectory_graph_color = "#000000")

plot_genes_in_pseudotime(my.cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)
