library(dpf)
library(batchtools)
source("my_model-rev1.R")
temps = dpf::tempos[,-c(1:3)]

load("mazurkaResults-update.Rdata")
lt = diff(c(tempos$note_onset,61))

performers = gsub("_"," ",rownames(pvec_ml))
lab_lookup = c("sigma[epsilon]^2", "mu[tempo]",
               "mu[acc]", "mu[stress]", "sigma[tempo]^2",
               "p[1*','*1]", "p[1*','*2]", "p[3*','*1]","p[1*','*3]",
               "p[2*','*1]","p[3*','*2]","p[2*','*2]")


getci <- function(i, pvec_ml, lt) {
  temps = dpf::tempos[,-c(1:3)]
  lab_lookup = c("sigma[epsilon]^2", "mu[tempo]",
                 "mu[acc]", "mu[stress]", "sigma[tempo]^2",
                 "p[1*','*1]", "p[1*','*2]", "p[3*','*1]","p[1*','*3]",
                 "p[2*','*1]","p[3*','*2]","p[2*','*2]")
  yt = matrix(temps[,i],nrow=1)
  theta = unlist(pvec_ml[i,1:12])
  hess = optimHess(theta, toOptimize, lt=lt, Npart=100, samp_mean=mean(yt),
                   priorfun=logprior, yt=yt)
  
  FLAG = tryCatch(solve(hess), error = function(e) "bad")
  if(is.character(FLAG)){ 
    B = 1/diag(hess)
  } else{ 
    B = diag(FLAG)
    FLAG = "good"
  }
  
  confidence_intervals = tibble::tibble(
    performer = gsub("_"," ",names(temps)[i]),
    ests = theta, 
    se2 = sqrt(B)*2, 
    lb = c(0,0,-Inf,-Inf,rep(0,8)),
    ub = c(Inf,Inf,0,0,Inf,rep(1,7)),
    lci = pmax(ests-se2, lb),
    uci = pmin(ests+se2, ub),
    param = lab_lookup,
    flag = FLAG)
  return(list(confidence_intervals=confidence_intervals, hess=hess))
}

makeRegistry("mazurka-cis", packages=c('dpf'), 
             source = "my_model-rev1.R",
             seed = 20181109)

ids = batchMap(getci, i=1:length(performers),  more.args = list(pvec_ml=pvec_ml, lt=lt))
ids$chunk = chunk(ids$job.id, 8)