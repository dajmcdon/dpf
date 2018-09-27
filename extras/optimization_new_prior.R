logprior <- function(theta, samp_mean=132){
  p1s = c(theta[8:9], 1-sum(theta[8:9]))
  sig2eps = dgamma(theta[1], shape=40, scale=10, log = TRUE)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 30, scale=1/2, log = TRUE)
  sig2tempo = dgamma(theta[5], shape=40, scale=10, log=TRUE)
  sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(34,2,4))
  p22 = dbeta(theta[10], 10, 5, log=TRUE)
  p31 = dbeta(theta[11], 5, 10, log=TRUE)
  lp = sum(sig2eps, mu1, 
           mu2, mu3, #sig2obs, 
           sig2tempo, sig2acc, sig2stress, p1, p22, p31)
  lp
}

prior_means <- function(samp_mean=132){
  c(400, samp_mean, -10, -15, 400, 1, 1, .85, 1/20, 10/15, 5/15)
}

rprior <- function(n, samp_mean=132){
  sig2eps = rgamma(n, shape=40, scale=10)
  mu1 = rgamma(n, samp_mean^2/100, scale=100/samp_mean)
  mu2 = -1*rgamma(n, 15, scale=2/3)
  mu3 = -1*rgamma(n, 30, scale=1/2)
  sig2tempo = rgamma(n, shape=40, scale=10)
  sig2acc = rgamma(n, shape=1, scale=1)
  sig2stress = rgamma(n, shape=1, scale=1)
  p1 = rdirichlet(n, alpha=c(34,2,4))[,-3,drop=FALSE]
  colnames(p1) = c('p11', 'p12')
  p22 = rbeta(n, 10, 5)
  p31 = rbeta(n, 5, 10)
  cbind(sig2eps, mu1, mu2, mu3, sig2tempo, sig2acc, sig2stress, p1, p22, p31)
}

# ddirichlet <- function(theta, alpha) sum(alpha*log(theta)) # on the log scale, no constant

init <- function(samp_mean=132, noise = 0){
  if(noise > 0){
    x = rprior(1, samp_mean)
  } else {
    x = prior_means(samp_mean)
  }
  x
}


logStatesGivenParams <- function(states,transProbs){
  ind = cbind(states[1:(length(states)-1)], states[2:length(states)]) + 1
  return(sum(log(transProbs[ind])))
}

toOptimize <- function(theta, yt, lt, Npart, samp_mean = 132, badvals=Inf){
  pmats = yupengMats(lt, theta[1], theta[2:4], theta[5:7], theta[8:11],
                     initialMean = c(132,0), # 132 is marked tempo, 0 is unused
                     initialVariance = c(400,10)) # sd of 20, 10 is unused
  beam = beamSearch(pmats$a0, pmats$P0, c(1,0,0,0,0,0,0,0), 
                    pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                    pmats$HHt, pmats$GGt, yt, pmats$transMat, Npart)
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
  logp = -1 * logprior(theta, samp_mean)
  obj = negllike + logp + sgp
  obj
}

## add transition from 2->1??

## Code to run on cluster



# Cluster funs ------------------------------------------------------------

optimizer <- function(perf, lt, Npart=200, ntries = 5, samp_mean=132, badvals=1e8){
  yt = matrix(perf, nrow=1)
  if(is.null(samp_mean)) samp_mean = mean(yt)
  randos = NULL
  if(ntries > 1) randos = rprior(ntries-1, samp_mean)
  init_vals = rbind(prior_means(samp_mean), randos)
  out1 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                   badvals=badvals,samp_mean=samp_mean,
                   method='Nelder-Mead',
                   control=list(trace=0, maxit=5000, badval=badvals))
  out2 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                    badvals=badvals,samp_mean=samp_mean,
                    method='SANN',
                    control=list(trace=0, maxit=5000,badval=badvals))
  out = rbind.data.frame(out1, out2)
  out
}
