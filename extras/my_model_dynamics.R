# Source file for prior and likelihood functions. Not needed separately.

logprior <- function(theta,yt){
  p2s = c(theta[9],1-sum(theta[c(9,10,11)]),theta[c(10,11)])
  p2 = ddirichlet(p2s, alpha=c(5,85,5,5),log=TRUE)
  p4 = dbeta(theta[12],2,8, log=TRUE)
  muc = dgamma(theta[1], shape=100, scale=.1, log = TRUE)
  sig2eps = dgamma(theta[2], shape=10, scale=.5, log = TRUE)
  mu0 = dnorm(theta[3], mean(yt), sd=10, log = TRUE)
  mu1 = dnorm(theta[4], 0, sd=.25, log = TRUE)
  mu2 = dnorm(theta[5], 0, sd=.25, log = TRUE)
  sig0 = dgamma(theta[6], shape=10, scale=1, log=TRUE)
  sig1 = dgamma(theta[7], shape=3, scale=1, log=TRUE)
  sig2 = dgamma(theta[8], shape=3, scale=1, log=TRUE)
  mue = dgamma(-theta[13], shape=100, scale=.1, log = TRUE)
  lp = sum(muc, sig2eps, mu0, mu1, mu2, sig0, sig1, sig2, 
           p2, p4, mue)
  lp
}

prior_means <- function(yt){
  c(10, 5,
    mean(yt), 0, 0,
    10, 3, 3,
    0.05,0.05,.05,0.2,
    -10)
}

rprior <- function(n,yt){
  muc = rgamma(n, shape=100, scale=.1)
  sig2eps = rgamma(n, shape=10, scale=.5)
  mu0 = rnorm(n, mean=mean(yt), sd = 10)
  mu1 = rnorm(n, mean=0, sd = .25)
  mu2 = rnorm(n, mean=0, sd = .25)
  sig0 = rgamma(n, shape=10, scale=1)
  sig1 = rgamma(n, shape=3, scale=1)
  sig2 = rgamma(n, shape=3, scale=1)
  p2 = rdirichlet(n, alpha=c(5,85,5,5))
  p21 = p2[1]
  p23 = p2[3]
  p24 = p2[4]
  p41 = rbeta(n,2,8)
  mue = -rgamma(n, shape=100, scale=.1)
  cbind(muc, sig2eps, mu0, mu1, mu2, sig0, sig1, sig2, 
        p21, p23, p24, p41, mue)
}



logStatesGivenParams <- function(states,transProbs){
  ind = cbind(states[1:(length(states)-1)], states[2:length(states)]) + 1
  return(sum(log(transProbs[ind])))
}


toOptimize <- function(theta, yt, lt, Npart, badvals=Inf){
  if(any(theta[c(2,6:8)]<0)){
    cat('do not allow negative variance')
    return(badvals)
  }
  pmats = musicModeldynamics(lt, theta[1], theta[2], theta[3:5], theta[6:8], theta[9:12], theta[13],
                             c(18,0,0), 
                             c(10,2,2))
  beam = beamSearch(pmats$a0, pmats$P0, c(1,0,0,0), 
                    pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                    pmats$HHt, pmats$GGt, yt, pmats$transMat, Npart, samplemethod = 1)
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
  logp = -1 * logprior(theta,yt)
  obj = negllike + logp + sgp
  obj
}


# Cluster funs ------------------------------------------------------------

optimizer <- function(perf, lt, Npart=500, ntries = 10, badvals=1e8){
  yt = matrix(perf, nrow=1)
  randos = NULL
  if(ntries > 1) randos = rprior(ntries-1,yt)
  init_vals = rbind(prior_means(yt), randos)
  out1 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                    badvals=badvals,
                    method='Nelder-Mead',
                    control=list(trace=0, maxit=5000, badval=badvals))
  out2 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                    badvals=badvals,
                    method='SANN',
                    control=list(trace=0, maxit=5000,badval=badvals))
  out = rbind.data.frame(out1, out2)
}
