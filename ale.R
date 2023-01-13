ale <- function(samp, multiref, cor.method=c("spearman", "pearson"), nb.cores=2){
  # samp: a gene expression matrix (genes by samples).
  # multiref: a list with
  #  - $interpGE, a named list of gene expression matrices (one per lineage), as returned by predict_mgeim()
  #  - $time, a vector of time values (length of ncol(interpGE[[1]]))
  
  lins <- names(multiref$interpGE)
  gl <- rownames(multiref$interpGE[[1]])
  cor.method <- match.arg(cor.method)
  
  if (!identical(gl, rownames(samp))) {
    overlap <- RAPToR::format_to_ref(samp, multiref$interpGE[[1]], verbose = T)
    samp <- overlap$samp
    multiref$interpGE <- lapply(multiref$interpGE, function(ge) ge[overlap$inter.genes, , drop=F])
    rm(overlap)
    gc(verbose = F)
  }
  samp.seq <- 1:ncol(samp)
  
  cl <- parallel::makeCluster(nb.cores, type = ifelse(.Platform$OS.type == 
                                                        "windows", "PSOCK", "FORK"))

  mres <- parallel::parLapply(cl = cl, X = multiref$interpGE, function(refdata){
    cors <- RAPToR::cor.gene_expr(samp, refdata, cor.method = cor.method)
    cmax.i <- apply(cors, 2, which.max)
    cmax <- cors[cbind(cmax.i,samp.seq)]
    
    age.est <- multiref$time[cmax.i]

    res <- data.frame(age.est=age.est, cor.score=cmax, 
                      row.names = colnames(samp))
    return(list(age.estimates=res, cors=cors))
  })
  
  
  # get age & cor. score per lineage
  cor_lin <- do.call(cbind, lapply(mres, function(res) res$age.estimates$cor.score))
  age_lin <- do.call(cbind, lapply(mres, function(res) res$age.estimates$age.est))
  
  l.rank <- t(apply(cor_lin, 1, order, decreasing=T))
  lmax <- cor_lin[cbind(samp.seq, l.rank[,1])]
  lmin <- cor_lin[cbind(samp.seq, l.rank[,ncol(l.rank)])]
  cf <- -1 + (lmax / lmin)
  

  ale_res <- data.frame(age.est=age_lin[cbind(samp.seq, l.rank[,1])], 
                        lin.est=lins[l.rank[,1]],
                        lin.conf=cf,
                        cor.score=lmax,
                        row.names = colnames(samp))
  
  parallel::stopCluster(cl)
  gc(verbose=F)
  
  return(list(ale=ale_res, ae.lin=mres))
}


# ale_t <- ale(dstinto$X, list(interpGE=interpGE, time=ndat$blastomere_time), 
#              nb.cores = 4)

