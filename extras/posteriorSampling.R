library(dpf)
data(tempos)
#library(pushoverr)
#a = load("~/Documents/dpf/extras/music.Rdata")
lt = diff(c(tempos$note_onset,61))
y = matrix(tempos$Rubinstein_1961,1)
w0 = rep(0, 8)
w0[1] = 1
toab <- function(x, a, b) x*(b-a) + a # maps [0,1] to [a,b]
toabInv <- function(y, a, b) (y - a)/(b - a)
logistic <- function(x) 1/(1+exp(-x)) # maps R to [0,1]
invlogistic <- function(x) log(x/(1-x))
RtoCon <- function(p){
  tp = logistic(p[9:12])
  tp[2] = toab(tp[2], 0, 1-tp[1])
  return(c(exp(p[1]), p[2:4], exp(p[5:8]), tp))
}
ContoR <- function(p){
  tp2 = toabInv(p[10], 0, 1 - p[9])
  return(c(log(p[1]), p[2:4], log(p[5:8]), invlogistic(p[9]), invlogistic(tp2), invlogistic(p[11:12])))
}

#prior distributions
sig2epsMean = log(25)
sig2epsSd = log(25)/2 #(log(25) - log(1))/2 
mu1Mean = 116
mu1Sd = 20
mu2Mean = -40
mu2Sd = 20
mu3Mean = -40
mu3Sd = 20
sig1Mean = log(.01^2) #log(1)
sig1Sd = .1 #log(4)/2
sig2Mean = log(400)
sig2Sd = (log(400) - log(225))/2
sig4Mean = log(20)
sig4Sd = (log(400) - log(225))/2
sig3Mean = log(25)
sig3Sd = log(25)/2
p1Mean = invlogistic(0.85)
p1Sd = (invlogistic(0.95) - invlogistic(0.85))/2
p2Mean = invlogistic(0.1)
p2Sd = (invlogistic(0.2) - invlogistic(0.1))/2
p3Mean = invlogistic(0.8)
p3Sd = (invlogistic(0.2) - invlogistic(0.1))/2
p4Mean = invlogistic(0.4)
p4Sd = (invlogistic(0.2) - invlogistic(0.1))/2
means = c(sig2epsMean, mu1Mean, mu2Mean, mu3Mean, sig1Mean, sig2Mean, sig4Mean, sig3Mean, p1Mean, p2Mean, p3Mean, p4Mean)
sds = c(sig2epsSd, mu1Sd, mu2Sd, mu3Sd, sig1Sd, sig2Sd, sig4Sd, sig3Sd, p1Sd, p2Sd, p3Sd, p4Sd)

propSample <- function(old, propsd){
  o = ContoR(old)
  out = rnorm(length(old), o, propsd)
  out[5] = o[5] # don't sample this one
  return(RtoCon(out))
}

logprior <- function(params, means, sds){
  p = ContoR(params)
  out = dnorm(p, means, sds, log=TRUE)
  out = out[-5] # not sampling this one
  return(sum(out))
}

logStatesGivenParams <- function(states,transProbs){
  ind = cbind(states[1:(length(states)-1)], states[2:length(states)]) + 1
  return(sum(log(transProbs[ind])))
}

posteriorSample <- function(y,lt,propsd, priormeans, priorsds, initParams = RtoCon(priormeans),
                            n=10000, Npart = 100){
  stopifnot((p <- length(initParams)) == 12, length(priormeans) == p, p ==length(priorsds),
            (Yn <- ncol(y))==length(lt)
            )
  Params = array(dim = c(n,p))
  States = array(dim = c(n,Yn))
  sig2eps = initParams[1]
  mus = initParams[2:4]
  sig2etas = initParams[5:8]
  transprobs = initParams[9:12]
  Params[1,] = initParams
  pmats = yupengMats(lt, sig2eps, mus, sig2etas, transprobs)
  beam = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                    pmats$Rt, pmats$Qt, pmats$GGt, y, pmats$transMat, Npart)
  States[1,] = beam$paths[sample(1:nrow(beam$paths), size = 1, prob = beam$weights),]
  LogLik = -1*getloglike(pmats, States[1,], y)#log(P(y|params, states))
  sgp = logStatesGivenParams(States[1,],pmats$transMat)
  pri = logprior(Params[1,], priormeans, priorsds)
  oldLogProb = LogLik + sgp + pri
  accept = 0
  for(i in 2:n){
    #print(i)
    proposal = propSample(Params[i-1,], propsd)
    sig2eps = proposal[1]
    mus = proposal[2:4]
    sig2etas = proposal[5:8]
    transprobs = proposal[9:12]
    pmats = yupengMats(lt, sig2eps, mus, sig2etas, transprobs)
    ## P(params|y,states) = P(y|params,states)P(states|params)p(params)  this is what I want. Will log for computation
    LogLik = -1*getloglike(pmats, States[i-1,], y)#log(P(y|params, states))
    sgp = logStatesGivenParams(States[i-1,], pmats$transMat)
    pri = logprior(proposal, priormeans, priorsds)
    p = min(c(1, exp(LogLik + sgp + pri - oldLogProb)))
    if(rbinom(1, 1, p)){
      Params[i,] = proposal
      accept = accept + 1
      oldLogProb = LogLik + sgp + pri
    }
    else{
      Params[i,] = Params[i-1,]
      sig2eps = Params[i,1]
      mus = Params[i,2:4]
      sig2etas = Params[i,5:8]
      transprobs = Params[i,9:12]
      pmats = yupengMats(lt, sig2eps, mus, sig2etas, transprobs)
      #P(params|y,states) = P(y|params,states)P(states|params)p(params)  this is what I want. Will log for computation
      LogLik = -1*getloglike(pmats, States[i-1,], y)#log(P(y|params, states))
      sgp = logStatesGivenParams(States[i-1,], pmats$transMat)
      pri = logprior(proposal, priormeans, priorsds)
      oldLogProb = LogLik + sgp + pri
    }
    beam = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                      pmats$Rt, pmats$Qt, pmats$GGt, y, pmats$transMat, Npart)
    States[i,] = beam$paths[sample(1:nrow(beam$paths), size = 1, prob = beam$weights),]
  }
  return(list(params = Params, states = States, acceptance = accept/n))
}

targetb = 0.28
targets = 0.22
propstep = 0.5
print(paste('trying ', propstep, sep = ''))
out = posteriorSample(y, lt, propstep, means, sds, n=1000)
oldRatio = out$acceptance
newRatio = oldRatio
ntries = 10
iter = 0

tryCatch({
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
    out = posteriorSample(y,lt,propstep, means, sds,n=1000)
    newRatio = out$acceptance
    print(newRatio)
    iter = iter + 1
  }
  print("OUT OF LOOP")
  print(paste("acceptance ratio:",out$acceptance,"ntries:",iter,"/",ntries))
  bestStep = propstep
  out = posteriorSample(y,lt,bestStep,means,sds,n=25000)
  print(paste("acceptance ratio: ",out$acceptance,sep = ''))
  save(bestStep, out, file = 'music_samp.Rdata')
  #write.table(out$params, file = 'params.csv', row.names = FALSE)
  #write.csv(out$states, file = 'states.csv', row.names = FALSE)
  #write.table(c(bestStep, out$acceptance), file = 'results.csv', row.names = FALSE)
  # pushover("Done!",
  #          user = "u1i3udnicjcosuaxx6615zhpivu73j", 
  #          app = "arzw55n8kkf18voqamjyer7vrk1dwq")
}, 
error = function(e){
  print(e)
  # pushover(message = "ERROR", 
  #         user = "u1i3udnicjcosuaxx6615zhpivu73j", app = "arzw55n8kkf18voqamjyer7vrk1dwq")
}
)