rm(list = ls())
setwd(this.path::here())

library(Biobase)
library(Seurat)
library(slingshot)

# Loading data ----

dspacker <- readRDS("./data/eset_packer.rds")
dspacker$lineage2 <- substr(dspacker$lineage, 1, 2) # make 2 char lineage

pX <- exprs(dspacker) # extract counts
table(rowSums(pX) < 10) # filter low-expressed genes (also to lighten the data)

# select annotated cells from one batch
sel_cells <- dspacker$batch == "Waterston_400_minutes" & dspacker$lineage != "unannotated"

pX <- as.matrix(pX[rowSums(pX) > 10, sel_cells])

p_fpx <- pData(dspacker)[sel_cells,]

saveRDS(pX, file = "./data/gene_count_matrix_packer.rds")
saveRDS(p_fpx, file = "./data/metadata_packer.rds")

pX <- readRDS(file = "./data/gene_count_matrix_packer.rds")
p_fpx <- readRDS(file = "./data/metadata_packer.rds")

# Preprocessing ----

seu <- CreateSeuratObject(counts = pX)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)

seu <- AddMetaData(seu, p_fpx$cell.type, 
                   col.name = "cell_type")

svg("fig/violinplot_celltypes.svg")
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "cell_type")
dev.off()

# Dimension reduction ----
seu <- RunPCA(seu)
# ElbowPlot(seu, ndims = 50)

svg("fig/dimplot_pca_celltypes.svg")
DimPlot(seu, reduction = "pca",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

seu <- RunTSNE(seu, dims = 1:50, verbose = FALSE)
svg("fig/dimplot_tsne_celltypes.svg")
DimPlot(seu, reduction = "tsne",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

seu <- RunUMAP(seu, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3, seed.use = 654)
svg("fig/dimplot_umap_celltypes.svg")
DimPlot(seu, reduction = "umap",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

# search clusters ----
seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:50)
seu <- FindClusters(seu, algorithm = 4, random.seed = 654, resolution = 1)

# Plot the clusters
svg("fig/dimplot_seurat_clusters.svg")
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)
dev.off()

# Save the objects as separate matrices for input in slingshot ---
dimred <- seu@reductions$umap@cell.embeddings
clustering <- seu$RNA_snn_res.1
counts <- as.matrix(seu@assays$RNA@counts[seu@assays$RNA@var.features, ])

saveRDS(seu, file = "./data/seu_packer.rds")
saveRDS(dimred, file = "./data/dimred_packer.rds")
saveRDS(clustering, file = "./data/clustering_packer.rds")
saveRDS(counts, file = "./data/counts_packer.rds")

seu <- readRDS(file = "./data/seu_packer.rds")
dimred <- readRDS(file = "./data/dimred_packer.rds")
clustering <- readRDS(file = "./data/clustering_packer.rds")
counts <- readRDS(file = "./data/counts_packer.rds")

# slingshot ----
set.seed(654)
sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = clustering, 
                 start.clus = 4, stretch = 0)


lineages <- getLineages(data = dimred, clusterLabels = clustering)

plot(dimred[,1:2], col = rainbow(nlevels(clustering))[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")

sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, 
                 start.clus = 4, stretch = 0)

plot(reducedDim(sds), col = rainbow(nlevels(clustering))[clustering], pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')