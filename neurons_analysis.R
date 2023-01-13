rm(list = ls())
setwd(this.path::here())

library(Biobase)
library(dplyr)
library(GEOquery, quietly = T) # to download data from GEO
library(ica)
library(limma, quietly = T) # for normalization
library(RAPToR, quietly = T) # install following instructions at www.github.com/LBMC/RAPToR
library(Seurat)
library(slingshot)
library(scales)
library(tradeSeq)
library(viridis)
library(wormRef, quietly = T) # for gene IDs, www.github.com/LBMC/wormRef

source("multige_im.R")

# Loading data ----

dspacker <- readRDS("./data/eset_packer.rds")
dspacker$lineage2 <- substr(dspacker$lineage, 1, 2) # make 2 char lineage

# select neurons lineages ASE, ASJ and AUA
csel <- grepl(pattern = "ASE|ASJ|AUA", dspacker$cell.subtype)
table(dspacker$cell.subtype[csel], dspacker$embryo.time.bin[csel])

pX <- exprs(dspacker) # extract counts
table(rowSums(pX) < 10) # filter low-expressed genes (also to lighten the data)

# select annotated cells from one batch
sel_cells <- dspacker$batch == "Waterston_400_minutes" & dspacker$lineage != "unannotated" & grepl(pattern = "ASE|ASJ|AUA", dspacker$cell.subtype)

pX <- as.matrix(pX[rowSums(pX) > 10, sel_cells])

p_fpx <- pData(dspacker)[sel_cells,]

# normalize & log
pX <- limma::normalizeBetweenArrays(pX, method = "quantile")
pX <- log1p(pX)

# compute correlation matrix
cc <- cor(pX, method="spearman")
diag(cc) <- NA # remove 1 diagonal
filt <- apply(cc, 1, quantile, probs=.99, na.rm=T)
thr <- mean(filt) - 2*sd(filt) # define threshold based on 99th percentile distribution

svg("fig/neurons_boxplot_correlation.svg")
boxplot(cc, border = 1+(filt<=thr), col = 0+7*(filt<=thr), names=NA,
        xlab = "Samples", ylab = "Spearman corr. with other samples", )
dev.off()

pX <- pX[,filt>thr]
p_fpx <- p_fpx[filt>thr,]

saveRDS(pX, file = "./data/neurons_gene_count_matrix_packer.rds")
saveRDS(p_fpx, file = "./data/neurons_metadata_packer.rds")

pX <- readRDS(file = "./data/neurons_gene_count_matrix_packer.rds")
p_fpx <- readRDS(file = "./data/neurons_metadata_packer.rds")


# Preprocessing ----

seu <- CreateSeuratObject(counts = pX)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)

seu <- AddMetaData(seu, p_fpx["cell.subtype"], 
                   col.name = "cell.subtype")

svg("fig/neurons_violinplot_celltypes.svg")
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "cell.subtype")
dev.off()

# Dimension reduction ----
seu <- RunPCA(seu)
# ElbowPlot(seu, ndims = 50)

svg("fig/neurons_dimplot_pca_celltypes.svg")
DimPlot(seu, reduction = "pca",
        group.by = "cell.subtype", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

seu <- RunTSNE(seu, dims = 1:50, verbose = FALSE)
svg("fig/neurons_dimplot_tsne_celltypes.svg")
DimPlot(seu, reduction = "tsne",
        group.by = "cell.subtype", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

seu <- RunUMAP(seu, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3, seed.use = 654)
svg("fig/neurons_dimplot_umap_celltypes.svg")
DimPlot(seu, reduction = "umap",
        group.by = "cell.subtype", pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

# search clusters ----
seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:50)
seu <- FindClusters(seu, algorithm = 1, random.seed = 654, resolution = 1.9)

# Leiden algorithm doesn't seems to work
# library(leiden)
# leiden(object = seu@graphs$RNA_snn)

# Plot the clusters
svg("fig/neurons_dimplot_seurat_clusters.svg")
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)
dev.off()

saveRDS(seu, file = "./data/neurons_seu_packer.rds")
seu <- readRDS(file = "./data/neurons_seu_packer.rds")

# slingshot ----
set.seed(654)
sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, stretch = 0)
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

# pseudotime

p_fpx$embryo.time <- p_fpx$embryo.time - min(p_fpx$embryo.time) + 1

cell_subtypes <- seu$cell.subtype

# for plotting
cols <- viridisLite::inferno(max(p_fpx$embryo.time), end=.9)[p_fpx$embryo.time]
pchs <- as.numeric(as.factor(cell_subtypes))
transp <- function(col, a=.5){
  # make color transparent 
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}

ctypes <- unique(cell_subtypes)

# define lineage cell sets
lins <- list(
  ASJ = cell_subtypes %in% ctypes[c(1, 4, 5)],
  AUA = cell_subtypes %in% ctypes[c(1, 5, 6)],
  ASE = cell_subtypes %in% ctypes[c(1, 2, 3)]
)

tpX <- scale(t(pX), scale = FALSE, center = TRUE) 
# PCA
pca_pX <- summary(prcomp(tpX, scale. = F, center = F))

ncs <- 1:8
lapply(seq_along(lins), function(li){ # for each lineage
  l <- lins[[li]]
  svg(paste0("fig/neurons_lineage_", names(lins[li]), "_pca.svg"))
  par(mfrow = c(5, max(ncs) / 2), mar=c(3,1,2,1), pty='s', bty='l')
  sapply(ncs, function(i){ # for each component
    plot(p_fpx$embryo.time, pca_pX$x[,i], 
         col = transp(cols, a=.2), pch = pchs,
         xlab = "blast stage", ylab = "PC",
         ylim=range(pca_pX$x[,i]),
         main = paste0(names(lins)[li]," lin, PC",i))
    points(p_fpx$embryo.time[l], pca_pX$x[l,i], lwd=2,
           col = cols[l], pch = pchs[l])
    
  })
  dev.off()
})

# GAM model on PCs
model_gam_pca <- function(X, p, formula, nc = ncol(X), subset=T){
  formula <- stats::as.formula(formula)
  
  tXc <- scale(t(X), center=T, scale = F)
  pca <- stats::prcomp(tXc, center=F, scale=F, rank=nc)
  pca$gcenters <- attr(tXc, "scaled:center")
  
  # subset PCs and p according to lineage
  pca$x <- pca$x[subset, ]
  p <- p[subset, ]
  
  formulas <- lapply(seq_len(nc), function(i) stats::update(formula, paste0("PC", i, " ~ .")))
  colnames(pca$x) <- paste0("PC", seq_len(nc))
  p <- cbind(p, pca$x)
  
  m <- list()
  m$model <- lapply(formulas, mgcv::gam, data = p)
  m$pca <- pca
  
  return(m)
}

# predict gam output
predict_gam_pca <- function(m, newdata, as.pc = FALSE){
  preds <- do.call(cbind, lapply(m$model, predict, newdata = newdata))
  
  if(!as.pc)
    return(apply(tcrossprod(m$pca$rotation, preds), 2, function(co) co + m$pca$gcenters))
  else return(preds)
}

ncs <- 1:8
# build models
ms <- lapply(seq_along(lins), function(li){
  m <- model_gam_pca(X = pX, p = p_fpx, # gene expr and pdata
                     nc = max(ncs), # components to use
                     formula = "X~s(embryo.time, bs='cr', k=6)", # gam formula
                     subset = lins[[li]] # lineage subset
  )
})

# build prediction data
ndat <- data.frame(embryo.time=seq(min(p_fpx$embryo.time), max(p_fpx$embryo.time), l=100))

# predict in component space
ps_c <- lapply(ms, predict_gam_pca, newdata=ndat, as.pc=T)

# build refs (=predict in gene space)
rfs <- lapply(seq_along(ms), function(i)
  list(interpGE=predict_gam_pca(ms[[i]], ndat), 
       time=ndat$embryo.time))

lapply(seq_along(lins), function(li){ # for each lineage
  l <- lins[[li]]
  svg(paste0("fig/neurons_lineage_", names(lins[li]), "_interpolation.svg"))
  par(mfrow = c(5,max(ncs) / 2), mar=c(3,1,2,1), pty='s', bty='l')
  sapply(ncs, function(i){ # for each component
    plot(p_fpx$embryo.time, pca_pX$x[,i], 
         col = transp(cols, a=.2), pch = pchs,
         xlab = "blast stage", ylab = "PC",
         ylim=range(pca_pX$x[,i]),
         main = paste0(names(lins)[li]," lin, PC",i))
    # show lineage-specific points 
    points(p_fpx$embryo.time[l], pca_pX$x[l,i], lwd=2,
           col = cols[l], pch = pchs[l])
    
    # add interpolation
    points(ndat$embryo.time, ps_c[[li]][,i], lwd=2, type = "l", col='green')
  })
  dev.off()
})

# staging

ae_s <- lapply(rfs, function(r){
  RAPToR::ae(pX, refdata = r$interpGE, ref.time_series = r$time, bootstrap.n = 30)
})

svg("fig/neurons_staging.svg")
par(mfrow = c(3,3), mar=c(3,1,2,1), pty='s', bty='l')
xl <- range(p_fpx$embryo.time)
sapply(seq_along(lins), function(i){
  l <- lins[[i]]
  sapply(seq_along(lins), function(j){
    x <- p_fpx$embryo.time[l]
    y <- ae_s[[j]]$age.estimates[l,1]
    
    plot(x, y, xlim = xl, ylim = xl,
         main = paste0(names(lins)[i], ' lin /', names(lins)[j], ' ref'), 
         col = i, xlab = "Blast. time", ylab="Estimated age")
    
    abline(a=0, b=1, lty=2) # add x=y
    # add cor.scores
    mtext(text = paste0("r = ", round(cor(x, y, method = 'pearson'), 3),
                        "\nrho =", round(cor(x, y, method = 'spearman'), 3)),
          at = 1, side = 3, line = -2, cex = .8, adj = 0)
    
  })
})
dev.off()

# Exploring lineage-inference from staging

svg("fig/neurons_lineage_inference.svg")
par(mfrow = c(2,3), mar=c(3,2,2,1), bty='l')
sapply(ctypes[c(1,5,4,6,3,2)], function(ci){
  ss <- which(cell_subtypes==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
  colnames(ccs) <- names(lins)
  boxplot(ccs, border=seq_along(lins), las=2,
          lwd=2, ylim=c(0, 1), col = 0, boxwex=.4,
          ylab = "Cor. at estimate", xlab = "reference used", main = ci)
  mtext(paste0('n=',length(ss)), at=.5, line=-1, cex=.75, side=3, adj = 0)
  
})
dev.off()

# redo with more data

ncs <- 1:40 # 90% of varexp
# build models
ms2 <- lapply(seq_along(lins), function(li){
  m <- model_gam_pca(X = pX, p = p_fpx, # gene expr and pdata
                     nc = max(ncs), # components to use
                     formula = "X~s(embryo.time, bs='cr', k=6)", # gam formula
                     subset = lins[[li]] # lineage subset
  )
})

# build refs (=predict in gene space)
rfs2 <- lapply(seq_along(ms), function(i)
  list(interpGE=predict_gam_pca(ms2[[i]], ndat), 
       time=ndat$embryo.time))

ae_s2 <- lapply(rfs2, function(r){
  RAPToR::ae(pX, refdata = r$interpGE, ref.time_series = r$time, bootstrap.n = 30)
})

svg("fig/neurons_lineage_inference_2.svg")
par(mfrow = c(2,3), mar=c(3,2,2,1), bty='l')
sapply(ctypes[c(1,5,4,6,3,2)], function(ci){
  ss <- which(cell_subtypes==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s2, function(a) a$age.estimates[ss,"cor.score"]))
  colnames(ccs) <- names(lins)
  boxplot(ccs, border=seq_along(lins), las=2,
          lwd=2, ylim=c(0, 1), col = 0, boxwex=.4,
          ylab = "Cor. at estimate", xlab = "reference used", main = ci)
  mtext(paste0('n=',length(ss)), at=.5, line=-1, cex=.75, side=3, adj = 0)
  
})
dev.off()

# function to normalize by min value in series
levelmin <- function(x){(x-min(x))}

svg("fig/neurons_lineage_inference_norm.svg")
par(mfrow = c(2,3), mar=c(3,2,2,1), bty='l', pty='s')
sapply(ctypes[c(1,5,4,6,3,2)], function(ci){
  ss <- which(cell_subtypes==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
  ccs <- t(apply(ccs,1, levelmin))
  x <- p_fpx$embryo.time[ss]
  colnames(ccs) <- names(lins)
  plot(range(p_fpx$embryo.time), range(ccs), type = 'n',
       lwd=2, ylim=c(0, .05), col = 0,
       ylab = "Correlation diff. between lineages", xlab = "Blast. time", main = ci)
  sapply(seq_along(lins), function(i){
    points(x, ccs[,i], col = transp(i, a = .9), lwd=2)
  })
  if(ci==ctypes[1])
    legend('right', bty='n', legend=names(lins), 
           lwd=3, lty=NA, pch=1, col=seq_along(lins),
           title = "Reference")
})
dev.off()

svg("fig/neurons_lineage_inference_norm_2.svg")
par(mfrow = c(2,3), mar=c(3,2,2,1), bty='l', pty='s')
sapply(ctypes[c(1,5,4,6,3,2)], function(ci){
  ss <- which(cell_subtypes==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s2, function(a) a$age.estimates[ss,"cor.score"]))
  ccs <- t(apply(ccs,1, levelmin))
  x <- p_fpx$embryo.time[ss]
  colnames(ccs) <- names(lins)
  plot(range(p_fpx$embryo.time), range(ccs), type = 'n',
       lwd=2, ylim=c(0, .11), col = 0,
       ylab = "Correlation diff. between lineages", xlab = "Blast. time", main = ci)
  sapply(seq_along(lins), function(i){
    points(x, ccs[,i], col = transp(i, a = .9), lwd=2)
  })
  if(ci==ctypes[1])
    legend('right', bty='n', legend=names(lins), 
           lwd=3, lty=NA, pch=1, col=seq_along(lins),
           title = "Reference")
})
dev.off()
