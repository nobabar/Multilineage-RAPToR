library(GEOquery, quietly = T) # to download data from GEO
library(limma, quietly = T) # for normalization

library(ica)
library(viridis)

library(RAPToR, quietly = T) # install following instructions at www.github.com/LBMC/RAPToR
library(wormRef, quietly = T) # for gene IDs, www.github.com/LBMC/wormRef

library(slingshot)

setwd(this.path::here())
rm(list = ls())

library(VisCello.celegans)

load("D:/Software/R/R-4.1.3/library/VisCello.celegans/app/data/g_all.rda")
