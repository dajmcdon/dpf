library(heatmaply)
library(dpf)
load("extras/mazurka1results.Rdata")
Rmat = apply(mat,1,ContoR)
Rmat = t(Rmat)
Rmat = Rmat[,-5]
Rmatsc = scale(Rmat)
heatmaply(as.matrix(dist(Rmatsc)),k_row=8,k_col=8,symm=TRUE,
          labCol=rep(NA,nrow(Rmat)), file='extras/dirty_heatmap.pdf')
