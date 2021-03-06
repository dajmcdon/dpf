# Source file for prior and likelihood functions. Not needed separately.

logprior <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9)]))
  p2s = c(theta[10], 1-theta[10])
  p3s = c(theta[11], 1-theta[11])
  sig2eps = dgamma(theta[1], shape=40, scale=10, log = TRUE)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=2, log = TRUE)
  sig2tempo = dgamma(theta[5], shape=40, scale=10, log=TRUE)
  #sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  #sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,10))
  p22 = ddirichlet(p2s, alpha=c(10,5))
  p31 = ddirichlet(p3s, alpha=c(5,10))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, #sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}

prior_means <- function(samp_mean=132){
  c(400, samp_mean, -10, -25, 400, #1, 1, 
    .85, 1/20, 10/15, 5/15)
}

rprior <- function(n, samp_mean=132){
  sig2eps = rgamma(n, shape=40, scale=10)
  mu1 = rgamma(n, samp_mean^2/100, scale=100/samp_mean)
  mu2 = -1*rgamma(n, 15, scale=2/3)
  mu3 = -1*rgamma(n, 20, scale=2)
  sig2tempo = rgamma(n, shape=40, scale=10)
  sig2acc = rgamma(n, shape=1, scale=1)
  sig2stress = rgamma(n, shape=1, scale=1)
  p1 = rdirichlet(n, alpha=c(85,5,10))
  p11 = p1[,1]
  p12 = p1[,2]
  p2 = rdirichlet(n, alpha=c(10,5))
  p22 = p2[,1]
  p3 = rdirichlet(n, alpha=c(5,10))
  p31 = p3[,1]
  cbind(sig2eps, mu1, mu2, mu3, sig2tempo, #sig2acc, sig2stress, 
        p11, p12, p22, p31)
}


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

toOptimize <- function(theta, yt, lt, Npart, samp_mean = 132, badvals=Inf, 
                       priorfun = logprior){
  theta = c(theta[1:5],1,1,theta[6:9])
  pmats = yupengMats(lt, theta[1], theta[2:4], theta[5:7], theta[8:11],
                     initialMean = c(samp_mean,0), # 132 is marked tempo, 0 is unused
                     initialVariance = c(400,10)) # sd of 20, 10 is unused
  beam = with(pmats, 
              beamSearch(a0, P0, c(1,0,0,0,0,0,0,0), 
                         dt, ct, Tt, Zt,
                         HHt, GGt, yt, transMat, Npart))
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
  logp = -1 * priorfun(theta, samp_mean)
  obj = negllike + logp + sgp
  obj
}


# Cluster funs ------------------------------------------------------------

optimizer <- function(perf, lt, Npart=200, ntries = 5, samp_mean=132, badvals=1e8, 
                      priorfun=logprior){
  yt = matrix(perf, nrow=1)
  if(is.null(samp_mean)) samp_mean = mean(yt)
  randos = NULL
  if(ntries > 1) randos = rprior(ntries-1, samp_mean)
  init_vals = rbind(prior_means(samp_mean), randos)
  out1 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart, 
                   badvals=badvals,samp_mean=samp_mean,
                   priorfun=priorfun,
                   method='Nelder-Mead',
                   control=list(trace=0, maxit=5000, badval=badvals))
  # out2 = multistart(init_vals, toOptimize, yt=yt, lt=lt, Npart=Npart,
  #                   badvals=badvals,samp_mean=samp_mean,
  #                   priorfun=priorfun,
  #                   method='SANN',
  #                   control=list(trace=0, maxit=5000,badval=badvals))
  # out = rbind.data.frame(out1, out2)
  out1
}
