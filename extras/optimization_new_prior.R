logprior <- function(pvec, samp_mean){
  theta = pvec
  p1s = c(theta[8:9], 1-sum(theta[8:9]))
  sig2eps = dnorm(theta[1], 400.0, sd=100.0, log = TRUE)
  #mu1 = dnorm(theta[2], samp_mean, sd=20.0, log = TRUE)
  mu2 = dnorm(theta[3], -10.0, sd=10.0, log = TRUE)
  mu3 = dnorm(theta[4], -30.0, sd=10.0, log = TRUE)
  sig2acc = dnorm(theta[5], 400, sd=100, log=TRUE)
  sig2dec = dnorm(theta[6], 400, sd=100, log=TRUE)
  sig2stress = dnorm(theta[7], 900, sd=100, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(20,5,5))
  p22 = dbeta(theta[10], 9, 6, log=TRUE)
  p31 = dbeta(theta[11], 10, 10, log=TRUE)
  lp = sum(sig2eps, #mu1, 
           mu2, mu3, #sig2obs, 
           sig2acc, sig2dec, sig2stress, p1, p22, p31)
  lp
}

ddirichlet <- function(theta, alpha) sum(alpha*log(theta)) # on the log scale, no constant

init <- function(samp_mean, noise = 0){
  means = c(400, samp_mean, -10, -30, .0001, 400, 400, 900, .85, .1, 9/15, 1/2)
  m = ContoR(means)
  ini = rnorm(length(m), m, sd = noise)
  ini = ini[-c(2,5)]
  ini
}

toab <- function(x, a, b) x*(b-a) + a # maps [0,1] to [a,b]
toabInv <- function(y, a, b) (y - a)/(b - a)
logistic <- function(x) 1/(1+exp(-x)) # maps R to [0,1]
invlogistic <- function(x) log(x/(1-x))
nonneg <- function(x) log(x)
invnonneg <- function(x) exp(x)
RtoCon <- function(p){
  tp = logistic(p[9:12])
  tp[2] = toab(tp[2], 0, 1-tp[1])
  return(c(invnonneg(p[1]), p[2], -invnonneg(p[3:4]), 0, invnonneg(p[6:8]), tp))
}
ContoR <- function(p){
  #deal with the simplex
  p[10] = toabInv(p[10], 0, 1-p[9])
  return(c(nonneg(p[1]), p[2], nonneg(-p[3:4]), 0, nonneg(p[6:8]), 
           invlogistic(p[9:12])))
}

logStatesGivenParams <- function(states,transProbs){
  ind = cbind(states[1:(length(states)-1)], states[2:length(states)]) + 1
  return(sum(log(transProbs[ind])))
}

toOptimize <- function(theta, yt, lt, Npart, badvals=Inf){
  samp_mean = mean(yt)
  theta = c(theta[1], samp_mean, theta[2:3], 0, theta[4:10])
  pvec = RtoCon(theta)
  pmats = yupengMats(lt, pvec[1], pvec[2:4], pvec[5:8], pvec[9:12],
                     initialMean = c(132,0), # 132 is marked tempo, 0 is unused
                     initialVariance = c(400,10)) # sd of 20, 10 is unused
  beam = beamSearch(pmats$a0, pmats$P0, c(1,0,0,0,0,0,0,0), 
                    pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                    pmats$Rt, pmats$Qt, pmats$GGt, yt, pmats$transMat, Npart)
  if(beam$LastStep < length(lt)){
    cat('beam$LastStep < length(lt)\n')
    return(badvals)
  }
  if(all(is.na(beam$weights))){
    cat('all weights are NA\n')
    return(badvals)
  }
  states = beam$paths[which.max(beam$weights),]
  negllike = getloglike(pmats, states, yt)# -log(P(y|params, states))
  sgp = -1 * logStatesGivenParams(states, pmats$transMat)
  logp = -1 * logprior(pvec, samp_mean)
  obj = negllike + logp + sgp
  obj
}

## add transition from 2->1??

## Code to run on cluster



# Cluster funs ------------------------------------------------------------

optimizer <- function(perf, lt, Npart=200, ntries = 5, spread_init=3,
                      badvals=1e8){
  yt = matrix(perf, nrow=1)
  samp_mean = mean(yt)
  randos = NULL
  if(ntries > 1) randos = t(replicate(ntries-1, init(samp_mean, spread_init)))
  init_vals = rbind(init(samp_mean,0), randos)
  out1 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                   badvals=badvals,
                   method='Nelder-Mead',
                   control=list(trace=0, maxit=5000,badval=badvals))
  out2 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                    badvals=badvals,
                    method='SANN',
                    control=list(trace=0, maxit=5000,badval=badvals))
  out = rbind.data.frame(out1, out2)
  if(all(is.na(out$value))) return(rep(NA, 12))
  theta = unlist(out[which.min(out$value), 1:10])
  theta = c(theta[1], samp_mean, theta[2:3], 0, theta[4:10])
  pvec = RtoCon(theta)
  pvec
}
