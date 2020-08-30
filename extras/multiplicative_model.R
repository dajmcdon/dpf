# Source file for prior and likelihood functions. Not needed separately.

logprior <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = dgamma(theta[1], shape=40, scale=10, log = TRUE)
  # mu1 on the log scale is, evaluate at exp(theta) and multiply by Jacobian
  mu1 = dgamma(exp(theta[2]), samp_mean^2/100, 
               scale=100/samp_mean, log = TRUE) + theta[2]
  mu2 = dbeta(-theta[3], 15, 185, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=2, log = TRUE)
  sig2tempo = dgamma(theta[5], 3.6, rate=30, log=TRUE)
  sig2acc = 0 # dgamma(theta[6], shape=1, scale=1, log=TRUE)
  sig2stress = 0 # dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,8,2))
  p22 = ddirichlet(p2s, alpha=c(10,1,4))
  p31 = ddirichlet(p3s, alpha=c(5,7,3))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}

prior_means <- function(samp_mean=132){
  c(400, digamma(samp_mean^2/100) + log(100/samp_mean), 
    -.075, -40, .12, #1, 1, 
    .85, 1/20, 10/15, 5/15, 1/50, 4/15, 3/15)
}

rprior <- function(n, samp_mean=132){
  sig2eps = rgamma(n, shape=40, scale=10)
  mu1 = log(rgamma(n, samp_mean^2/100, scale=100/samp_mean)) # changed
  mu2 = -1*rbeta(n, 15, 185)
  mu3 = -1*rgamma(n, 20, scale=2)
  sig2tempo = rgamma(n, shape=3.6, rate=30)
  sig2acc = rgamma(n, shape=1, scale=1)
  sig2stress = rgamma(n, shape=1, scale=1)
  p1 = rdirichlet(n, alpha=c(85,5,8,2))
  p13 = p1[,4]
  p11 = p1[,1]
  p12 = p1[,2]
  p2 = rdirichlet(n, alpha=c(10,1,4))
  p22 = p2[,1]
  p21 = p2[,3]
  p3 = rdirichlet(n, alpha=c(5,7,3))
  p31 = p3[,1]
  p32 = p3[,3]
  cbind(sig2eps, mu1, mu2, mu3, sig2tempo, #sig2acc, sig2stress, 
        p11, p12, p22, p31, p13, p21, p32)
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

toOptimize <- function(theta, yt, lt, Npart, samp_mean = 132, badvals=Inf){
  theta = c(theta[1:5], .05^2, 1, theta[6:12])
  trans_par = theta[8:14]
  transProbs = dpf:::ecreateTransMat(trans_par)
  sig2eps = theta[1]
  mus = theta[2:4]
  sig2eta = theta[5:7]
  a0 = c(log(samp_mean), 0)
  P0 = c(.12, 0)
  logp = -1 * logprior(theta, samp_mean)
  if(logp > 1e8) return(badvals)
  if(mus[1] > 7) return(badvals)
  if(mus[1] < 0) return(badvals)
  if(is.na(logp)) return(badvals)
  # print(theta)
  beam = dpf:::ebeamSearch(
    lt, c(1,rep(0,10)), sig2eps, mus, sig2eta, a0, P0, yt, transProbs, Npart)
  if(beam$LastStep < length(lt)){
    cat('beam$LastStep < length(lt)\n')
    return(badvals)
  }
  if(all(is.na(beam$weights))){
    cat('all weights are NA\n')
    return(badvals)
  }
  states = beam$paths[which.max(beam$weights),]
  negllike = dpf:::egetloglike(lt, states, yt, sig2eps, mus, sig2eta, a0, P0)
  sgp = -1 * logStatesGivenParams(states, transProbs)
  obj = negllike + logp + sgp
  obj
}


# Cluster funs ------------------------------------------------------------

optimizer <- function(perf, lt, Npart=200, ntries = 5, samp_mean=NULL, badvals=1e8){
  yt = as.vector(perf)
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
