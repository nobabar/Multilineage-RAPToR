rm(list = ls())
setwd(this.path::here())

multige_im <- function(X, p, 
                       lineages, formulas,
                       dim_red=c('pca','ica'), nc=min(c(10, ncol(X))), ncs=NULL){
  # X: Expression matrix (genes by samples), expecting log1p(tpm)
  # p: Pheno data, samples as rows. Must at the minimum have the variables used in formulas.
  # lineages: named list of vectors of indices or booleans defining each lineage
  # formulas: the formulas to use for each model. can be given as 
  #     - a single formula to use the same for all lineages/components
  #     - a vector of size length(lineages) to use lineage-specific formulas
  #     - a vector of size length(ncs) to use component-specific formulas
  #     - a matrix of dim [length(lineages), length(ncs)] to use lineage- and component-specific formulas
  #
  # dim_red: whether to use pca/ica
  # nc: the number of components to extract
  # ncs (optional): which components to keeo. defaults to 1:nc.
  
  nli <- length(lineages)
  if(typeof(lineages[[1]])==typeof(T)){lineages <- lapply(lineages, which)}
  if(is.null(names(lineages))){names(lineages) <- paste0('lin',1L:nli)}
  dim_red <- match.arg(dim_red)
  if(is.null(ncs)){ncs <- 1L:nc}
  lnc <- length(ncs)
  
  if(is.null(dim(formulas)) | any(dim(formulas)!=c(nli, lnc))){
    lf <- length(formulas)
    if(lf==nli){
      message("Using lineage-specific formulas.")
      formulas <- matrix(formulas, nrow = nli, ncol=lnc, byrow = F)
    } else if(lf==lnc){
      message("Using component-specific formulas.")
      formulas <- matrix(formulas, nrow = nli, ncol=lnc, byrow = T)
    } else {
      if (lf != 1)
        warning("formulas length doesn't match lineages or components.\n  --> Using the first formula for all models.")
      formulas <- matrix(formulas[1], nrow = nli, ncol=lnc)
    }
  }
  # ncX <- ncol(X)
  tXc <- scale(t(X), center=T, scale=F)
  
  if('pca' == dim_red){
    pca <- stats::prcomp(tXc, center=F, scale=F)
    dr <- list(
      G = pca$rotation[, ncs],
      S = pca$x[, ncs]
    )
  } else if('ica' == dim_red){
    ica_ <- ica::icafast(t(tXc), center=T, nc=nc)
    dr <- list(
      G = ica_$S[, ncs],
      S = ica_$M[, ncs]
    )
  }
  dr$ncs <- ncs
  dr$gcenters <- attr(tXc, "scaled:center")
  dr$backtransf <- function(G,S, gcent) apply(tcrossprod(G,S), 2, function(co) co+gcent)
  dr$dim_red <- dim_red
  
  models <- lapply(seq_len(nli), function(li){
    Sli <- dr$S[lineages[[li]], , drop=F]
    colnames(Sli) <- paste0(".__C", ncs)
    pli <- cbind(p[lineages[[li]], , drop=F], Sli)
    
    fs <- lapply(seq_along(ncs), 
                 function(i) 
                   stats::update(stats::as.formula(formulas[li, i]), 
                                 paste0(".__C", ncs[i], " ~ ."))
                 )
    return(lapply(fs, mgcv::gam, data = pli))
  })
  names(models) <- names(lineages)
  
  m <- list(model = models,
            lineages = lineages,
            dim_red = dr)

}

predict_mgeim <- function(m, newdata=NULL, as.c=F){
  preds <- lapply(seq_along(m$lineages), 
                  function(li)
                    do.call(cbind, lapply(m$model[[li]], predict, newdata = newdata))
                  )
  names(preds) <- names(m$lineages)
  
  if(!as.c)
    return(lapply(preds, function(predli) m$dim_red$backtransf(m$dim_red$G, predli, m$dim_red$gcenters)))
  else return(preds)
}


# Test/ example
# Uses the g & p objects as defined in the multi_traj_tests.Rmd Report
# 
# mt <- multige_im(X = g, p = p, lineages = lins,
#                  dim_red = 'ica', nc = 8,
#                  formulas = c("X~s(blastomere_time, bs='cr', k=6)",
#                               "X~s(blastomere_time, bs='cr', k=6)",
#                               "X~s(blastomere_time, bs='cs', k=6)",
#                               "X~s(blastomere_time, bs='cr', k=5)",
#                               "X~s(blastomere_time, bs='cr', k=5)",
#                               "X~s(blastomere_time, bs='cr', k=5)",
#                               "X~s(blastomere_time, bs='cr', k=4)",
#                               "X~s(blastomere_time, bs='cr', k=4)"))
# 
# 
# 
# # predict new data (in comp. space)
# newdat <- data.frame(blastomere_time=seq(1,11, l=100))
# nX <- predict_mgeim(mt, newdat, as.c = T)
# 
# # plot results
# par(mfrow = c(5,8), bty='l', pty='s')
# sapply(seq_along(mt$lineages), function(li){
#   sl <- mt$lineages[[li]]
#   sapply(1:8, function(ci){
#     plot(p$blastomere_time[sl], mt$dim_red$S[sl,ci], xlab = "time", ylab="IC",
#          main = paste0(names(mt$lineages)[li], "lin, IC", mt$dim_red$ncs[ci]),
#          ylim = range(mt$dim_red$S[,ci]), xlim = range(p$blastomere_time))
#     points(p$blastomere_time[-sl], mt$dim_red$S[-sl,ci], cex = .1, col = 'grey')
#     points(newdat$blastomere_time, nX[[li]][,ci], type  ='l', col = 'red')
#   })
# })

## predict new data in gene expr. space:
#interpGE <- predict_mgeim(mt, newdat)


# then apply standard RAPToR ae using the different refs:
# ae_multi <- lapply(interpGE, function(rdat){
#   RAPToR::ae(samp = X, refdata = rdat, ref.time_series = ndat$blastomere_time,
#      bootstrap.n = 1) 
# })