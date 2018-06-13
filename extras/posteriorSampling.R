library(dpf)
library(pushoverr)
#a = load("~/Documents/dpf/extras/music.Rdata")
lt = diff(c(note.onset,61))
y = matrix(tempo.mat[,2],1)
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
sig1Mean = 0 #log(1)
sig1Sd = log(4)/2
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
p3Mean = invlogistic(0.2)
p3Sd = (invlogistic(0.2) - invlogistic(0.1))/2
p4Mean = invlogistic(0.2)
p4Sd = (invlogistic(0.2) - invlogistic(0.1))/2
means = c(sig2epsMean, mu1Mean, mu2Mean, mu3Mean, sig1Mean, sig2Mean, sig4Mean, sig3Mean, p1Mean, p2Mean, p3Mean, p4Mean)

propSample <- function(old, propsd){
    o = ContoR(old)
    out = array(dim = length(old))
    out[1] = rnorm(1, mean = o[1], sd = propsd)
    out[2] = rnorm(1, mean = o[2], sd = propsd)
    out[3] = rnorm(1, mean = o[3], sd = propsd)
    out[4] = rnorm(1, mean = o[4], sd = propsd)
    out[5] = rnorm(1, mean = o[5], sd = propsd)
    out[6] = rnorm(1, mean = o[6], sd = propsd)
    out[7] = rnorm(1, mean = o[7], sd = propsd)
    out[8] = rnorm(1, mean = o[8], sd = propsd)
    out[9] = rnorm(1, mean = o[9], sd = propsd)
    out[10] = rnorm(1, mean = o[10], sd = propsd)
    out[11] = rnorm(1, mean = o[11], sd = propsd)
    out[12] = rnorm(1, mean = o[12], sd = propsd)
    return(RtoCon(out))
}

logprior <- function(params){
    p = ContoR(params)
    out = log(dnorm(p[1], mean = sig2epsMean, sd = sig2epsSd))+
        log(dnorm(p[2], mean = mu1Mean, sd = mu1Sd))+
        log(dnorm(p[3], mean = mu2Mean, sd = mu2Sd))+
        log(dnorm(p[4], mean = mu3Mean, sd = mu3Sd))+
        log(dnorm(p[5], mean = sig1Mean, sd = sig1Sd))+
        log(dnorm(p[6], mean = sig2Mean, sd = sig2Sd))+
        log(dnorm(p[7], mean = sig4Mean, sd = sig4Sd))+
        log(dnorm(p[8], mean = sig3Mean, sd = sig3Sd))+
        log(dnorm(p[9], mean = p1Mean, sd = p1Sd))+
        log(dnorm(p[10], mean = p2Mean, sd = p2Sd))+
        log(dnorm(p[11], mean = p3Mean, sd = p3Sd))+
        log(dnorm(p[12], mean = p4Mean, sd = p4Sd))
    #print(which(is.infinite(c(log(dnorm(p[1], mean = sig2epsMean, sd = sig2epsSd)),
    #                        log(dnorm(p[2], mean = mu1Mean, sd = mu1Sd)),
    #                        log(dnorm(p[3], mean = mu2Mean, sd = mu2Sd)),
    #                        log(dnorm(p[4], mean = mu3Mean, sd = mu3Sd)),
    #                        log(dnorm(p[5], mean = sig1Mean, sd = sig1Sd)),
    #                        log(dnorm(p[6], mean = sig2Mean, sd = sig2Sd)),
    #                        log(dnorm(p[7], mean = sig4Mean, sd = sig4Sd)),
    #                        log(dnorm(p[8], mean = sig3Mean, sd = sig3Sd)),
    #                        log(dnorm(p[9], mean = p1Mean, sd = p1Sd)),
    #                        log(dnorm(p[10], mean = p2Mean, sd = p2Sd)),
    #                        log(dnorm(p[11], mean = p3Mean, sd = p3Sd)),
    #                        log(dnorm(p[12], mean = p4Mean, sd = p4Sd))))))
    return(out)
}

logStatesGivenParams <- function(states,transProbs){
    ind = cbind(states[1:(length(states)-1)], states[2:length(states)]) + 1
    return(sum(log(transProbs[ind])))
}

posteriorSample <- function(y,lt,propsd,initParams = RtoCon(means),#initParams = c(400,mean(y[1:100]),-40,-40,100,400,100,900,.8,.1,.8,.25),
                            n=10000, Npart = 100){
    p = length(initParams)
    Yn = ncol(y)
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
    pri = logprior(Params[1,])
    oldLogProb = LogLik + sgp + pri
    accept = 0
    for(i in 2:n){
        #print(i)
        proposal = propSample(Params[i-1,],propsd)
        sig2eps = proposal[1]
        mus = proposal[2:4]
        sig2etas = proposal[5:8]
        transprobs = proposal[9:12]
        pmats = yupengMats(lt, sig2eps, mus, sig2etas, transprobs)
        #P(params|y,states) = P(y|params,states)P(states|params)p(params)  this is what I want. Will log for computation
        LogLik = -1*getloglike(pmats, States[i-1,], y)#log(P(y|params, states))
        sgp = logStatesGivenParams(States[i-1,],pmats$transMat)
        pri = logprior(proposal)
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
            sgp = logStatesGivenParams(States[i-1,],pmats$transMat)
            pri = logprior(proposal)
            oldLogProb = LogLik + sgp + pri
        }
        beam = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct, pmats$Tt, pmats$Zt,
                          pmats$Rt, pmats$Qt, pmats$GGt, y, pmats$transMat, Npart)
        States[i,] = beam$paths[sample(1:nrow(beam$paths), size = 1, prob = beam$weights),]
    }
    return(list(params = Params, states = States, acceptance = accept/n))
}

target = 0.25
propstep = 0.5
print(paste('trying ', propstep, sep = ''))
out = posteriorSample(y,lt,propstep,n=100)
oldRatio = out$acceptance
newRatio = oldRatio

tryCatch({
    while(sign(target - newRatio) == sign(target - oldRatio)){
        oldRatio = newRatio
        oldstep = propstep
        if(oldRatio > target){
            propstep = propstep + 0.1*propstep
        }
        else{
            propstep = propstep - 0.1*propstep
        }
        print(paste('trying ', propstep, sep = ''))
        out = posteriorSample(y,lt,propstep,n=100)
        newRatio = out$acceptance
        print(newRatio)
    }
    print("OUT OF LOOP")
    bestStep = ifelse(abs(target - newRatio) < abs(target - oldRatio), propstep, oldstep)
    out = posteriorSample(y,lt,bestStep)
    print(paste("acceptance ratio: ",out$acceptance,sep = ''))
    write.table(out$params, file = 'params.csv', row.names = FALSE)
    write.csv(out$states, file = 'states.csv', row.names = FALSE)
    write.table(c(bestStep, out$acceptance), file = 'results.csv', row.names = FALSE)
    pushover("Done!",
             user = "u1i3udnicjcosuaxx6615zhpivu73j", app = "arzw55n8kkf18voqamjyer7vrk1dwq")
}, error = function(e){
    print(e)
    pushover(message = "ERROR", 
             user = "u1i3udnicjcosuaxx6615zhpivu73j", app = "arzw55n8kkf18voqamjyer7vrk1dwq")
})