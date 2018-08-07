library(heatmaply)
library(dpf)
load("extras/mazurka1results.Rdata")
source("extras/djm_optimization.R")
Rmat = apply(pvec_ml,1,ContoR)
Rmat = t(Rmat)
Rmat = Rmat[,-5]
Rmatsc = scale(Rmat)
heatmaply(as.matrix(dist(Rmatsc)),k_row=8,k_col=8,symm=TRUE,
          labCol=rep(NA,nrow(Rmat)), file='extras/dirty_heatmap.pdf', width = 1280, height = 800)
