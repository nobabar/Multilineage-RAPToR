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

source("ale.R")
source("multige_im.R")
source("utils.R")

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

# build train and test
p <- 0.8
strats <- p_fpx$cell.subtype

rr <- split(1:length(strats), strats)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

pX_train <- pX[idx, ]
pX_test <- pX[-idx, ]

p_fpx_train <- p_fpx[idx, ]
p_fpx_test <- p_fpx[-idx, ]

table(p_fpx$cell.subtype) / nrow(p_fpx)
table(p_fpx_train$cell.subtype) / nrow(p_fpx_train)
table(p_fpx_test$cell.subtype) / nrow(p_fpx_test)

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
seu <- FindClusters(seu, algorithm = 1, random.seed = 654, resolution = 1.25)

# Leiden algorithm doesn't seems to work
# library(leiden)
# leiden(object = seu@graphs$RNA_snn)

# Plot the clusters
svg("fig/neurons_dimplot_seurat_clusters.svg")
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label = TRUE)
dev.off()

saveRDS(seu, file = "./data/neurons_seu_packer.rds")
seu <- readRDS(file = "./data/neurons_seu_packer.rds")

# pseudotime

p_fpx$embryo.time <- p_fpx$embryo.time - min(p_fpx$embryo.time) + 1

cell_subtypes <- p_fpx$cell.subtype

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

prepare_and_interpolate <- function(X, p, lins, dim_red="ica", nc=8){
  mt <- multige_im(X = X, p = p, lineages = lins,
                   dim_red = dim_red, nc = nc,
                   formulas = c("X~s(embryo.time, bs='cr', k=6)"))
  
  # predict new data (in comp. space)
  ndat <- data.frame(embryo.time=seq(min(p$embryo.time),
                                     max(p$embryo.time), l=100))
  nX <- predict_mgeim(mt, ndat, as.c = T)
  interpGE <- predict_mgeim(mt, ndat)
  
  return(list(mt=mt, ndat=ndat, nX=nX, interpGE=interpGE))
}

multige <- prepare_and_interpolate(pX, p_fpx, lins)

svg("fig/neurons_ica_pairs.svg",  width=8, height=6)
pairs(multige$mt$dim_red$S[,1:4], col = cols, pch=pchs, lwd=2, labels=paste0("Comp.", 1:4), oma=c(2,2,2,25))
par(xpd = TRUE)
legend("bottomright", pch = as.numeric(as.factor(unique(cell_subtypes))),
       legend = c(levels(as.factor(cell_subtypes))), title="Lineage")
legend_image <- as.raster(matrix(viridisLite::inferno(max(p_fpx$embryo.time), end=.9), ncol=1))
text(x=.85, y=.85, labels = "Time")
rasterImage(legend_image, .8, .8, .9, .5)
rect(.8, .5, .9, .8)
labs <- seq(0, max(p_fpx$embryo.time),50)
ys <- seq(.5, .78, l=length(labs))
text(x=.903, y=ys, labels = "-") 
text(x=.94, y=ys, labels = labs)
dev.off()

svg("fig/neurons_lineage.svg")
plot_lineages(multige$mt, multige$ndat, multige$nX,
              p_fpx$embryo.time, p_fpx$cell.subtype, 1:8)
dev.off()

res <- lapply(1:4, function(i){
  set.seed(i)
  p <- 0.8
  idx <- stratified_split(cell_subtypes, p)
  
  pX_train <- pX[, idx]
  pX_test <- pX[, -idx]
  
  p_fpx_train <- p_fpx[idx, ]
  p_fpx_test <- p_fpx[-idx, ]
  
  lins_train <- lapply(lins,'[', idx)
  lins_test <- lapply(lins,'[', -idx)
  
  multige <- prepare_and_interpolate(pX_train, p_fpx_train, lins_train)
  
  ae_s_train <- ale(samp = pX_train, multiref=list(interpGE=multige$interpGE, time=multige$ndat$embryo.time))
  
  ae_s_test <- ale(samp = pX_test, multiref=list(interpGE=multige$interpGE, time=multige$ndat$embryo.time))
  
  return(list(multige=multige,
              train=list(idx=idx,
                         ae_s=ae_s_train),
              test=list(idx=seq_along(cell_subtypes)[!seq_along(cell_subtypes) %in% idx],
                        ae_s=ae_s_test)))
})

saveRDS(res, file = "./data/neurons_raptor_res.rds")
res <- readRDS(file = "./data/neurons_raptor_res.rds")

svg("fig/neurons_staged_trainset.svg")
plot_staged(res[[1]]$train$ae_s$ae.lin, p_fpx[res[[1]]$train$idx, ]$embryo.time,
            lapply(lins,'[', res[[1]]$train$idx))
dev.off()

lapply(seq_along(res), function(i){
  svg(paste0("fig/neurons_staged_testset_", i,".svg"))
  plot_staged(res[[i]]$test$ae_s$ae.lin, p_fpx[res[[i]]$test$idx, ]$embryo.time,
              lapply(lins,'[', res[[i]]$test$idx), mode="diag")
  dev.off()
})

svg("fig/neurons_lineage_inference_norm_trainset.svg")
plot_lineage_inference_norm(res[[1]]$train$ae_s$ae.lin,
                            p_fpx[res[[1]]$train$idx, ]$cell.subtype,
                            ctypes[c(1,5,4,6,3,2)],
                            p_fpx[res[[1]]$train$idx, ]$embryo.time,
                            lapply(lins,'[', res[[1]]$train$idx))
dev.off()

svg("fig/neurons_lineage_inference_norm_testset.svg")
plot_lineage_inference_norm(res[[1]]$test$ae_s$ae.lin,
                            p_fpx[res[[1]]$test$idx, ]$cell.subtype,
                            ctypes[c(1,5,4,6,3,2)],
                            p_fpx[res[[1]]$test$idx, ]$embryo.time,
                            lapply(lins,'[', res[[1]]$test$idx))
dev.off()

reduce_lineage_inference(res, p_fpx, "train")

reduce_lineage_inference(res, p_fpx, "test")
