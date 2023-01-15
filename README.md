# Adding support to multiple branching trajectories on RAPToR

This repository contains a demonstration of the RAPToR tool with multiple branching trajectories. The algorithm is described in [this paper]{https://www.nature.com/articles/s41592-022-01540-0}. The code is based on the [original RAPToR implementation]{https://github.com/LBMC/RAPToR}.

# Implementation

The main function is `multige_im ` in the `multige_im.R` file.

Parameters:
X: Expression matrix (genes by samples), expecting log1p(tpm)
p: Pheno data, samples as rows. Must at the minimum have the variables used in formulas.
lineages: named list of vectors of indices or booleans defining each lineage
formulas: the formulas to use for each model. can be given as 
    - a single formula to use the same for all lineages/components
    - a vector of size length(lineages) to use lineage-specific formulas
    - a vector of size length(ncs) to use component-specific formulas
    - a matrix of dim [length(lineages), length(ncs)] to use lineage- and component-specific formulas

Returns:
dim_red: whether to use pca/ica
nc: the number of components to extract
ncs (optional): which components to keeo. defaults to 1:nc.

We also developed a function to infer lineages on top of the usual RAPToR output. This function is the `ale` function in the `ale.R` file. It is based on the [original `ae` function]{https://github.com/LBMC/RAPToR/blob/master/R/ae.R}

Parameters:
samp: a gene expression matrix (genes by samples).
multiref: a list with
 - $interpGE, a named list of gene expression matrices (one per lineage), as returned by predict_mgeim()
 - $time, a vector of time values (length of ncol(interpGE[[1]]))

In addition to these two you will find several utility functions in the `utils.R` file.

# Data

For this analysis we worked with a single-cell C.elegans dataset [(Packer et al. 2019)]{https://www.science.org/doi/10.1126/science.aax1971}. The data is accessible on [GEO: GSE126954]{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126954}. The data is also accessible via the VisCello R package (Zhu, Murray, Tan, and Kim 2019) and the expression matrix is directly downloadable from their GitHub repository [here]{https://media.githubusercontent.com/media/qinzhu/VisCello.celegans/master/inst/app/data/eset.rds}.

# Usage

The whole analysis can be found in the `neurons_analysis.R` file as well as the `neurons_analysis.Rms` notebook, compiled with knitr into the `neurons_analysis.html` file.

Several figures were generated in the `fig` folder.