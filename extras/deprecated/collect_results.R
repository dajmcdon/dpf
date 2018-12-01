library(batchtools)
loadRegistry('mazurka1', writeable=TRUE)

getStatus()

outL = reduceResultsList()
pvec_ml = as.data.frame(do.call(rbind,outL))
rownames(pvec_ml) =  names(tempos[-c(1:3)])
colnames(pvec_ml) = c("sig2eps", "mu1", "mu2", "mu3", 
                 "sig1", "sig2", "sig4", "sig3", "p1", 
                 "p2", "p3", "p4")

save(pvec_ml, file='mazurka_code/mazurka1results.Rdata')