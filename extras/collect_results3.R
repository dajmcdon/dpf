library(batchtools)
loadRegistry('mazurka1', writeable=TRUE)

getStatus()

outL = reduceResultsList()
fun = function(x) x[which.min(x$value),]
pvec_ml = t(sapply(outL,fun))
rownames(pvec_ml) = names(tempos[-c(1:3)])

save(pvec_ml, file='mazurka_code/mazurka3results.Rdata')