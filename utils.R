stratified_split <- function(strats, prop){
  rr <- split(1:length(strats), strats)
  idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))
  return(idx)
}

# function to normalize by min value in series
levelmin <- function(x) (x-min(x))

transp <- function(col, a=.5){
  # make color transparent 
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}

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

reduce_cor <- function(res, lins, m, set, method="pearson"){
  cor_df <- lapply(seq_along(res), function(i){
    cor_mat <- sapply(seq_along(lins), function(j){
      l <- lins[[j]][res[[i]][[set]]$idx]
      sapply(seq_along(lins), function(k){
        r <- round(cor(m[res[[i]][[set]]$idx, ]$embryo.time[l],
                       res[[i]][[set]]$ae_s[[k]]$age.estimates[l,1],
                       method = method), 3)
      })
    })
    as.data.frame(cor_mat)
  })
  
  cor_df <- Reduce(`+`, cor_df)/length(cor_df)
  colnames(cor_df) <- names(lins)
  rownames(cor_df) <- names(lins)
  return(cor_df)
}

reduce_summary <- function(res, lins, ctypes, set){
  summary_df <- lapply(ctypes, function(ci){
    summary_df <- lapply(seq_along(res), function(i){
      ss <- which(cell_subtypes[res[[i]][[set]]$idx]==ci)
      
      # get cor. score at estimate for each reference
      ccs <- do.call(cbind, lapply(res[[i]][[set]]$ae_s, function(a) a$age.estimates[ss,"cor.score"]))
      colnames(ccs) <- names(lins)
      as.data.frame(apply(ccs, 2, summary))
    })
    Reduce(`+`, summary_df)/length(summary_df)
  })
  names(summary_df) <- ctypes
  return(summary_df)
}

reduce_norm_summary <- function(res, lins, ctypes, set){
  summary_df <- lapply(ctypes, function(ci){
    summary_df <- lapply(seq_along(res), function(i){
      ss <- which(cell_subtypes[res[[i]][[set]]$idx]==ci)
      
      # get cor. score at estimate for each reference
      ccs <- do.call(cbind, lapply(res[[i]][[set]]$ae_s, function(a) a$age.estimates[ss,"cor.score"]))
      ccs <- t(apply(ccs, 1, levelmin))
      colnames(ccs) <- names(lins)
      as.data.frame(apply(ccs, 2, summary))
    })
    Reduce(`+`, summary_df)/length(summary_df)
  })
  names(summary_df) <- ctypes
  return(summary_df)
}

plot_lineages <- function(mt, ndat, nX, time, ctypes, ncs){
  cols <- viridisLite::inferno(max(time), end=.9)[time]
  pchs <- as.numeric(as.factor(ctypes))
  
  par(mfrow = c(3,8), mar=c(3,1,2,1), bty='l', pty='s')
  sapply(seq_along(mt$lineages), function(li){
    l <- mt$lineages[[li]]
    sapply(ncs, function(ci){
      plot(time[l], mt$dim_red$S[l,ci], xlab = "time", ylab="IC",
           col = transp(cols, a=.2), pch = pchs,
           main = paste0(names(mt$lineages)[li], "lin, IC", mt$dim_red$ncs[ci]),
           ylim = range(mt$dim_red$S[,ci]), xlim = range(time))
      points(time[-l], mt$dim_red$S[-l,ci], cex = .1, col = cols[l], pch = pchs[l], lw=2)
      points(ndat$embryo.time, nX[[li]][,ci], type='l', col = 'red')
    })
  })
}

plot_staged <- function(ae_s, time, lins){
  par(mfrow = c(length(lins), length(lins)), mar=c(3,1,2,1), pty='s', bty='l')
  xl <- range(time)
  sapply(seq_along(lins), function(i){
    l <- lins[[i]]
    sapply(seq_along(lins), function(j){
      x <- time[l]
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
}

plot_lineages_inference <- function(ae_s, cell_subtypes, ctypes, lins){
  par(mfrow = c(2,length(ctypes) / 2), mar=c(3,3,2,1), bty='l')
  sapply(ctypes, function(ci){
    ss <- which(cell_subtypes==ci) # select a cell type
    
    # get cor. score at estimate for each reference
    ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
    colnames(ccs) <- names(lins)
    boxplot(ccs, border=seq_along(lins), las=2,
            lwd=2, ylim=c(0, 1), col = 0, boxwex=.4,
            ylab = "Cor. at estimate", xlab = "reference used", main = ci)
    mtext(paste0('n=',length(ss)), at=.5, line=-1, cex=.75, side=3, adj = 0)
  })
}

plot_lineage_inference_norm <- function(ae_s, cell_subtypes, ctypes, time, lins){
  par(mfrow = c(2,length(ctypes) / 2), mar=c(3,2,2,1), bty='l', pty='s')
  sapply(ctypes, function(ci){
    ss <- which(cell_subtypes==ci) # select a cell type
    
    # get cor. score at estimate for each reference
    ccs <- do.call(cbind, lapply(ae_s, function(a) a$age.estimates[ss,"cor.score"]))
    ccs <- t(apply(ccs,1, levelmin))
    x <- time[ss]
    colnames(ccs) <- names(lins)
    plot(range(time), range(ccs), type = 'n',
         lwd=2, ylim=range(ccs), col = 0,
         ylab = "Correlation diff. between lineages", xlab = "Blast. time", main = ci)
    sapply(seq_along(lins), function(i){
      points(x, ccs[,i], col = transp(i, a = .9), lwd=2)
    })
    if(ci==ctypes[1])
      legend('right', bty='n', legend=names(lins), 
             lwd=3, lty=NA, pch=1, col=seq_along(lins),
             title = "Reference")
  })
}
