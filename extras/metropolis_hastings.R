library(dpf)
logprior <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = dgamma(theta[1], shape=40, scale=10, log = TRUE)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=2, log = TRUE)
  sig2tempo = dgamma(theta[5], shape=40, scale=10, log=TRUE)
  sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,8,2))
  p22 = ddirichlet(p2s, alpha=c(10,1,4))
  p31 = ddirichlet(p3s, alpha=c(5,7,3))
  lp = sum(sig2eps, mu1, 
           mu2, mu3, #sig2obs, 
           sig2tempo, sig2acc, sig2stress, p1, p22, p31)
  lp
}

prior_means <- function(samp_mean=132){
  c(400, samp_mean, -10, -40, 400, 1, 1, .85, 1/20, 10/15, 5/15, 1/50, 4/15, 3/15)
}

rprior <- function(n, samp_mean=132){
  sig2eps = rgamma(n, shape=40, scale=10)
  mu1 = rgamma(n, samp_mean^2/100, scale=100/samp_mean)
  mu2 = -1*rgamma(n, 15, scale=2/3)
  mu3 = -1*rgamma(n, 20, scale=2)
  sig2tempo = rgamma(n, shape=40, scale=10)
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
  cbind(sig2eps, mu1, mu2, mu3, sig2tempo, sig2acc, sig2stress, 
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

mlogit <- function(p, tol = 1e-16){
    # Transforms probabilities on the simplex to real values
    # checks and clean up
    stopifnot(all(p <= 1), all(p >= 0))
    # add a column if we're using only K-1 classes
    sumto1 = sum(p)
    allp = TRUE
    if(1-sumto1 > tol){
        allp = FALSE # assumes you put in fewer probs than classes
        p <- c(p, 1-sumto1)
    }
    K = length(p)
    xy = log(p[-K] / p[K])
    if(allp) xy = c(xy, 0)
    xy
}

invmlogit <- function(xy, tol=1e-16){
    # transforms real values to probabilities on the simplex
    K = length(xy)
    allp = TRUE
    # check if all
    if(abs(1-xy[K]) > tol){
        allp = FALSE
        K = K + 1
    }
    xy = xy[1:(K-1)]
    exy = exp(xy)
    denom = 1+sum(exy)
    p = exy / denom
    if(allp) p <- c(p, 1 - sum(p))
    p
}


RtoCon <- function(p){
    p1s = invmlogit(p[c(8,9,12)])
    p2s = invmlogit(p[c(10,13)])
    p3s = invmlogit(p[c(11, 14)])
    return(c(exp(p[1]), p[2:4], exp(p[5:7]), p1s[1], p1s[2], p2s[1], p3s[1], p1s[3], p2s[2], p3s[2]))
}
ContoR <- function(p){
    p1s = mlogit(p[c(8,9,12)])
    p2s = mlogit(p[c(10, 13)])
    p3s = mlogit(p[c(11, 14)])
    return(c(log(p[1]), p[2:4], log(p[5:7]), p1s[1], p1s[2], p2s[1], p3s[1], p1s[3], p2s[2], p3s[2]))
}

Det <- function(P){
    a = P[1,1]
    b = P[1,2]
    c = P[1,3]
    d = P[2,1]
    e = P[2,2]
    f = P[2,3]
    g = P[3,1]
    h = P[3,2]
    k = P[3,3]
    return(a*(e*k - f*h) - b*(d*k - f*g) + c*(d*h - e*g))
}

invmlogit_prime <- function(p){
    S = 1 + sum(exp(p))
    mat = (-exp(p) %o% exp(p) + diag(exp(p)*S))/S^2
    return(abs(Det(mat)))
}

Jacob <- function(p){
    j1 = -log(p[1])
    # no transformation for 2:4
    j57 = -log(p[5:7])
    jp1s = log(invmlogit_prime(p[c(8,9,12)]))
    jp2s = log(invmlogit_prime(p[c(10,13)]))
    jp3s = log(invmlogit_prime(p[c(11, 14)]))
    det = j1+sum(j57)+sum(jp1s)+sum(jp2s)+sum(jp3s)
    return(log(abs(det)))
}

propSample <- function(old, propsd){
    o = ContoR(old)
    out = rnorm(length(old), o, propsd)
    return(RtoCon(out))
}

posteriorSample <- function(y,lt,propsd, initialMean = c(mean(y),0), # 132 is marked tempo, 0 is unused
                            initialVariance = c(400,10), burn = 10000,
                            n=50000, keep_every = 5, Npart = 100){
    if(is.vector(y)) dim(y) = c(1, length(y))
    w0 = rep(0, 8)
    w0[1] = 1
    samp_mean = mean(y)
    initParams = prior_means(mean(y))
    stopifnot((p <- length(initParams)) == 14, (Yn <- ncol(y))==length(lt))
    Params = array(dim = c(n+burn,p))
    States = array(dim = c(n+burn,Yn))
    sig2eps = initParams[1]
    mus = initParams[2:4]
    sig2etas = initParams[5:7]
    transprobs = initParams[8:14]
    Params[1,] = initParams
    pmats = musicModel(lt, sig2eps, mus, sig2etas, transprobs, initialMean, initialVariance)
    beam = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                      pmats$HHt, pmats$GGt, y, pmats$transMat, Npart)
    States[1,] = beam$paths[sample(1:nrow(beam$paths), size = 1, prob = beam$weights),]
    LogLik = -1*getloglike(pmats, States[1,], y)#log(P(y|params, states))
    sgp = logStatesGivenParams(States[1,],pmats$transMat)
    pri = logprior(Params[1,], mean(y))
    oldLogProb = LogLik + sgp + pri + Jacob(Params[1,])
    accept = 0
    for(i in 2:(n+burn)){
        print(i)
        proposal = propSample(Params[i-1,], propsd)
        sig2eps = proposal[1]
        mus = proposal[2:4]
        sig2etas = proposal[5:7]
        transprobs = proposal[8:14]
        pmats = musicModel(lt, sig2eps, mus, sig2etas, transprobs, initialMean, initialVariance)
        ## P(params|y,states) = P(y|params,states)P(states|params)p(params)  this is what I want. Will log for computation
        LogLik = -1*getloglike(pmats, States[i-1,], y)#log(P(y|params, states))
        sgp = logStatesGivenParams(States[i-1,], pmats$transMat)
        pri = logprior(proposal, mean(y))
        p = min(c(1, exp(LogLik + sgp + pri + Jacob(proposal) - oldLogProb)))
        if(is.nan(p)){
            p = 0
        }
        if(rbinom(1, 1, p)){
            Params[i,] = proposal
            accept = accept + 1
            oldLogProb = LogLik + sgp + pri + Jacob(proposal)
        }
        else{
            Params[i,] = Params[i-1,]
            sig2eps = Params[i,1]
            mus = Params[i,2:4]
            sig2etas = Params[i,5:7]
            transprobs = Params[i,8:14]
            pmats = musicModel(lt, sig2eps, mus, sig2etas, transprobs, initialMean, initialVariance)
            #P(params|y,states) = P(y|params,states)P(states|params)p(params)  this is what I want. Will log for computation
            LogLik = -1*getloglike(pmats, States[i-1,], y)#log(P(y|params, states))
            sgp = logStatesGivenParams(States[i-1,], pmats$transMat)
            pri = logprior(Params[i,], mean(y))
            oldLogProb = LogLik + sgp + pri + Jacob(Params[i,])
        }
        beam = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                          pmats$HHt, pmats$GGt, y, pmats$transMat, Npart)
        States[i,] = beam$paths[sample(1:nrow(beam$paths), size = 1, prob = beam$weights),]
    }
    return(list(params = Params[seq(1+burn, n+burn, keep_every),], states = States[seq(1+burn, n+burn, keep_every),], acceptance = accept/n))
}

optimal_sample <- function(y, lt, name, targetb = 0.28, targets = 0.22, initial_step = 0.5, ntries = 10){
    propstep = initial_step
    print(paste('trying ', propstep, sep = ''))
    out = posteriorSample(y, lt, propstep, n=1000)
    oldRatio = out$acceptance
    newRatio = oldRatio
    iter = 0
    while(iter < ntries){
        if(newRatio < targetb & newRatio > targets) break
        oldRatio = newRatio
        oldstep = propstep
        if(oldRatio > targetb){
            propstep = 1.1*propstep
        }
        else{
            propstep = 0.9*propstep
        }
        print(paste('trying ', propstep, sep = ''))
        out = posteriorSample(y,lt,propstep, n=1000)
        newRatio = out$acceptance
        print(newRatio)
        iter = iter + 1
    }
    print("OUT OF LOOP")
    print(paste("acceptance ratio:",out$acceptance,"ntries:",iter,"/",ntries))
    bestStep = propstep
    out = posteriorSample(y,lt,bestStep)
    print(paste("acceptance ratio: ",out$acceptance,sep = ''))
    #save(bestStep, out, file = 'music_samp.Rdata')
    #write.table(out$params, file = paste(name, '_params.csv'), row.names = FALSE)
    #write.csv(out$states, file = paste(name, '_states.csv'), row.names = FALSE)
    #write.table(c(bestStep, out$acceptance), file = paste(name, '_results.csv'), row.names = FALSE)
    return(list(params = out$params,
                states = out$states,
                sampling = c(bestStep, out$acceptance)))
    # pushover("Done!",
    #          user = "u1i3udnicjcosuaxx6615zhpivu73j", 
    #          app = "arzw55n8kkf18voqamjyer7vrk1dwq")
}

toOptimize <- function(theta, yt, lt, Npart, samp_mean = 132, badvals=Inf){
  pmats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14],
                     initialMean = c(samp_mean,0), # 132 is marked tempo, 0 is unused
                     initialVariance = c(400,10)) # sd of 20, 10 is unused
  beam = beamSearch(pmats$a0, pmats$P0, c(1,0,0,0,0,0,0,0,0,0), 
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
