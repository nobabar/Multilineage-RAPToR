library(dplyr)
library(Seurat)
library(slingshot)
library(viridis)

# data ----

pX <- readRDS(file = "./data/neurons_gene_count_matrix_packer.rds")
p_fpx <- readRDS(file = "./data/neurons_metadata_packer.rds")
seu <- readRDS(file = "./data/neurons_seu_packer.rds")

# slingshot ----

set.seed(654)
sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters,
                 start.clus = 4, end.clus=c(), stretch = 0)
sds <- as.SlingshotDataSet(sds)

# create the color palette
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(seu$cell.subtype, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())


curves <- as.SlingshotDataSet(getCurves(sds, approx_points = 300, thresh = 0.01,
                                        stretch = 0.8, allow.breaks = FALSE,
                                        shrink = 0.99))

svg("fig/neurons_slingshot_lineages.svg")
plot(slingReducedDim(sds), col = cell_colors, pch = 16)
lines(curves, lwd = 3, col = "black")
dev.off()

svg("fig/neurons_slingshot_split_lineages.svg")
pt <- slingPseudotime(sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(1, 2))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(slingReducedDim(sds), col = colors, pch = 16, main = i)
  lines(curves, lwd = 2, col = 'black')
}
dev.off()

# differential expression

counts <- as.matrix(seu@assays$RNA@counts[seu@assays$RNA@var.features, ])
dim(counts)

sce <- fitGAM(counts = as.matrix(counts), sds = curves)

saveRDS(sce, file = "./data/neurons_sce.rds")
sce <- readRDS(file = "./data/neurons_sce.rds")

svg("fig/neurons_genecount.svg")
plotGeneCount(curves, counts, clusters = seu$seurat_clusters, models = sce)
dev.off()

pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

plot_differential_expression <- function(feature_id, counts, clustering) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)

svg("fig/neurons_differential_expression.svg")
plot_differential_expression(feature_id, counts, seu$seurat_clusters)
dev.off()