#sample R script to do module-trait correlation analysis
#by Minghui Wang
#

#Simulate datasets
#simulate gene expression data
dat=matrix(rnorm(1000),ncol=50,dimnames=list(paste0('S',1:20),paste0('G',1:50))) # dat is a matrix of gene expression data (say a module) with genes in the columns and samples in the rows
#simulate trait
trait=rnorm(nrow(dat))

#PCA
library('pcaMethods')
tdat <- prep(dat, scale='uv', center=TRUE)
pcs <- pca(tdat, method='ppca', nPcs=3)@scores #compute the top 3 PCs
row.names(pcs) <- row.names(tdat)
# check PC1 for sign to be consistent with at least half of the genes
if(sum(cor(pcs[,1], dat, use="pairwise.complete.obs") > 0) < ncol(dat)/2) pcs[,1] =  0 - pcs[,1]

#Spearman correlation between the top PC and a trait
fit=cor.test(pcs[,1],trait,method='spearman',exact=FALSE)
cat('correlation =',fit$estimate,'\n')
cat('p.value =',fit$p.value,'\n')
