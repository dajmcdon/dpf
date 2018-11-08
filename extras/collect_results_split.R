library(batchtools)
loadRegistry('mazurka_split', writeable=TRUE)

getStatus()

outL = reduceResultsList()
fun <- function(listy){
  fun2 <- function(x) x[which.min(x$value),]
  s = sapply(listy, fun2)
  s = matrix(s,nrow=1)
  s
}

pvec_ml = t(sapply(outL,fun))
rownames(pvec_ml) = names(tempos[-c(1:3)])
pvec_ml = as.data.frame(pvec_ml)
nn = names(outL[[1]][[1]])
nn = c(nn,paste0(nn,'_2'))
names(pvec_ml) = nn

save(pvec_ml, file='mazurka_code/mazurka_split_results.Rdata')