library(GEOquery, quietly = T) # to download data from GEO
library(limma, quietly = T) # for normalization

library(ica)
library(viridis)

library(RAPToR, quietly = T) # install following instructions at www.github.com/LBMC/RAPToR
library(wormRef, quietly = T) # for gene IDs, www.github.com/LBMC/wormRef

setwd(this.path::here())

# get count matrix
supfiles <- GEOquery::getGEOSuppFiles("GSE50548", fetch_files = F) # we want the full expression data (#2)
fname <- as.character(supfiles$fname[2])
furl <- as.character(supfiles$url[2])
utils::download.file(url = furl, destfile = fname)

g <- read.table(gzfile(fname), quote = '"', sep = '\t', 
                h=T, row.names = 1)

# get sample info
geo_obj <- GEOquery::getGEO("GSE50548")[[1]] 
## Found 1 file(s)
## GSE50548_series_matrix.txt.gz
p <- Biobase::pData(geo_obj)

# keep relevant fields
p <- p[, c("title", "geo_accession", 
           "time - blastomeres:ch1", "time:ch1", 
           "collection timing in reference to e lineage:ch1",
           "germ layer:ch1", "experiment set:ch1")]
colnames(p) <- c("title", "geo_accession", "blastomere_time", "time_min", 
                 "lineage_coll_time", "germ_layer", "exp_set")

# keep only single-cell data
selsamp <- !grepl("Whole embryo", p$exp_set)

p <- p[selsamp,]
g <- g[,selsamp]


file.remove(fname)
rm(supfiles, furl, fname, geo_obj, selsamp)

# format ids
g <- RAPToR::format_ids(g, wormRef::Cel_genes, from="sequence_name", to="wb_id", aggr.fun = sum)
## Kept 20168 out of 20517 - aggregated into 20168

# remove non-expr genes (at least 5 counts across all samples)
g <- g[rowSums(g)>5,]

# normalize & log
g <- limma::normalizeBetweenArrays(g, method = "quantile")
g <- log1p(g)

# compute correlation matrix
cc <- cor(g, method="spearman")
diag(cc) <- NA # remove 1 diagonal
filt <- apply(cc, 1, quantile, probs=.99, na.rm=T)
thr <- mean(filt) - 2*sd(filt) # define threshold based on 99th percentile distribution
boxplot(cc, border = 1+(filt<=thr), col = 0+7*(filt<=thr), names=NA,
        xlab = "Samples", ylab = "Spearman corr. with other samples", )

g <- g[,filt>thr]
p <- p[filt>thr,]

p$blastomere_time <- as.numeric(p$blastomere_time)

# for plotting
cols <- viridisLite::inferno(11, end=.9)[as.numeric(p$blastomere_time)]
pchs <- as.numeric(as.factor(p$germ_layer))
transp <- function(col, a=.5){
  # make color transparent 
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}

# cell types/lineages
ctypes <- unique(p$germ_layer)

# define lineage cell sets
lins <- list(
  AB=p$germ_layer %in% ctypes[2],
  P=p$germ_layer %in% ctypes[c(1, 4, 8)],
  C=p$germ_layer %in% ctypes[c(1, 4, 7)],
  E=p$germ_layer %in% ctypes[c(1, 3, 5)],
  MS=p$germ_layer %in% ctypes[c(1, 3, 6)]
)

# center on genes
tXc <- scale(t(g), scale = FALSE, center = TRUE)
# PCA
pca_x <- summary(prcomp(tXc, scale. = F, center = F))

ncs <- 1:4 # components to look at
pairs(pca_x$x[,ncs], col = cols, pch=pchs, lwd=2)

ncs <- 1:8

par(mfrow = c(5,max(ncs)), mar=c(4,4,3,1), pty='s', bty='l')
lapply(seq_along(lins), function(li){ # for each lineage
  l <- lins[[li]]
  sapply(ncs, function(i){ # for each component
    plot(p$blastomere_time, pca_x$x[,i], 
         col = transp(cols, a=.2), pch = pchs,
         xlab = "blast stage", ylab = "PC",
         ylim=range(pca_x$x[,i]),
         main = paste0(names(lins)[li]," lin, PC",i))
    points(p$blastomere_time[l], pca_x$x[l,i], lwd=2,
           col = cols[l], pch = pchs[l])
    
  })
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
  m <- model_gam_pca(X = g, p = p, # gene expr and pdata
                     nc = max(ncs), # components to use
                     formula = "X~s(blastomere_time, bs='cr', k=6)", # gam formula
                     subset = lins[[li]] # lineage subset
  )
  
})

# build prediction data
ndat <- data.frame(blastomere_time=seq(min(p$blastomere_time), max(p$blastomere_time), l=100))

# predict in component space
ps_c <- lapply(ms, predict_gam_pca, newdata=ndat, as.pc=T)

# build refs (=predict in gene space)
rfs <- lapply(seq_along(ms), function(i)
  list(interpGE=predict_gam_pca(ms[[i]], ndat), 
       time=ndat$blastomere_time))

par(mfrow = c(5,max(ncs)), mar=c(4,4,3,1), pty='s', bty='l')
lapply(seq_along(lins), function(li){ # for each lineage
  l <- lins[[li]]
  sapply(ncs, function(i){ # for each component
    plot(p$blastomere_time, pca_x$x[,i], 
         col = transp(cols, a=.2), pch = pchs,
         xlab = "blast stage", ylab = "PC",
         ylim=range(pca_x$x[,i]),
         main = paste0(names(lins)[li]," lin, PC",i))
    # show lineage-specific points 
    points(p$blastomere_time[l], pca_x$x[l,i], lwd=2,
           col = cols[l], pch = pchs[l])
    
    # add interpolation
    points(ndat$blastomere_time, ps_c[[li]][,i], lwd=2, type = "l", col='green')
  })
})

# for each reference (lineage) do age estimation of all samples
ae_s <- lapply(rfs, function(r){
  RAPToR::ae(g, refdata = r$interpGE, ref.time_series = r$time, bootstrap.n = 30)
})

par(mfrow = c(5,5), mar=c(4,4,3,1), pty='s', bty='l')

xl <- range(p$blastomere_time)
sapply(seq_along(lins), function(i){
  l <- lins[[i]]
  sapply(seq_along(lins), function(j){
    x <- p$blastomere_time[l]
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

par(mfrow = c(2,4), mar=c(4,4,3,1), bty='l')
sapply(ctypes[c(1,3,4,2,5:8)], function(ci){
  ss <- which(p$germ_layer==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
  colnames(ccs) <- names(lins)
  boxplot(ccs, border=seq_along(lins), las=2,
          lwd=2, ylim=c(.5, 1), col = 0, boxwex=.4,
          ylab = "Cor. at estimate", xlab = "reference used", main = ci)
  mtext(paste0('n=',length(ss)), at=.5, line=-1, cex=.75, side=3, adj = 0)
  
})

ncs <- 1:40 # 75% of varexp
# build models
ms2 <- lapply(seq_along(lins), function(li){
  m <- model_gam_pca(X = g, p = p, # gene expr and pdata
                     nc = max(ncs), # components to use
                     formula = "X~s(blastomere_time, bs='cr', k=6)", # gam formula
                     subset = lins[[li]] # lineage subset
  )
  
})
# build refs (=predict in gene space)
rfs2 <- lapply(seq_along(ms), function(i)
  list(interpGE=predict_gam_pca(ms2[[i]], ndat), 
       time=ndat$blastomere_time))

ae_s2 <- lapply(rfs2, function(r){
  RAPToR::ae(g, refdata = r$interpGE, ref.time_series = r$time, bootstrap.n = 30)
})

par(mfrow = c(2,4), mar=c(4,4,3,1), bty='l')
sapply(ctypes[c(1,3,4,2,5:8)], function(ci){
  ss <- which(p$germ_layer==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s2, function(a) a$age.estimates[ss,"cor.score"]))
  colnames(ccs) <- names(lins)
  boxplot(ccs, border=seq_along(lins), las=2,
          lwd=2, ylim=c(.5, 1), col = 0, boxwex=.4,
          ylab = "Cor. at estimate", xlab = "reference used", main = ci)
  mtext(paste0('n=',length(ss)), at=.5, line=-1, cex=.75, side=3, adj = 0)
  
})

par(mfrow = c(2,4), mar=c(4,4,3,1), bty='l', pty='s')

# function to normalize by min value in series
levelmin <- function(x){(x-min(x))}

yl <- c(0, .28) # plot limits

sapply(ctypes[c(1,3,4,2,5:8)], function(ci){
  ss <- which(p$germ_layer==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
  ccs <- t(apply(ccs,1, levelmin))
  x <- p$blastomere_time[ss]
  colnames(ccs) <- names(lins)
  plot(range(p$blastomere_time), range(ccs), type = 'n',
       lwd=2, ylim=yl, col = 0,
       ylab = "Correlation diff. between lineages", xlab = "Blast. time", main = ci)
  sapply(seq_along(lins), function(i){
    points(x, ccs[,i], col = transp(i, a = .9), lwd=2)
  })
  if(ci==ctypes[1])
    legend('right', bty='n', legend=names(lins), 
           lwd=3, lty=NA, pch=1, col=seq_along(lins),
           title = "Reference")
})

par(mfrow = c(2,4), mar=c(4,4,3,1), bty='l', pty='s')
sapply(ctypes[c(1,3,4,2,5:8)], function(ci){
  ss <- which(p$germ_layer==ci) # select a cell type
  
  # get cor. score at estimate for each reference
  ccs <- do.call(cbind, lapply(ae_s2, function(a) a$age.estimates[ss,"cor.score"]))
  ccs <- t(apply(ccs,1, levelmin))#t(apply(ccs,1,range01))
  x <- p$blastomere_time[ss]
  colnames(ccs) <- names(lins)
  plot(range(p$blastomere_time), range(ccs), type = 'n',
       lwd=2, ylim=yl, col = 0,
       ylab = "Correlation diff. between lineages", xlab = "Blast. time", main = ci)
  sapply(seq_along(lins), function(i){
    points(x, ccs[,i], col = transp(i, a = .9), lwd=2)
  })
  
  if(ci==ctypes[1])
    legend('right', bty='n', legend=names(lins), 
           lwd=3, lty=NA, pch=1, col=seq_along(lins),
           title = "Reference")
  
})
