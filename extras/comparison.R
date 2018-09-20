load("manuscript/mazurka2results.Rdata")
source('extras/djm_optimization.R')
library(dpf)
fliere = ContoR(unlist(pvec_ml[substr(rownames(pvec_ml), 1, 6) == "Fliere",]))[-c(2,5)]
tomsic = ContoR(unlist(pvec_ml[substr(rownames(pvec_ml), 1, 6) == "Tomsic",]))[-c(2,5)]
data("tempos")
lt = diff(c(tempos$note_onset,61))
yf = matrix(tempos$Fliere_1977,1)
yt = matrix(tempos$Tomsic_1995, 1)
out = array(dim = c(2,2))#each row is different parameters, each column is different performance
out[1,1] = toOptimize(fliere, yf, lt, 200)
out[2,1] = toOptimize(tomsic, yf, lt, 200)
out[1,2] = toOptimize(fliere, yt, lt, 200)
out[2,2] = toOptimize(tomsic, yt, lt, 200)
