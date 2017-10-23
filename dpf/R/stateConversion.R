convert4to8 <- function(path){
    cpps = 0:7
    tm1 = c(1,1,1,2,2,3,3,4)
    t1 = c(1,2,4,2,3,1,3,1)
    n = length(path)
    npath = c(1,path[-n])
    path4 = integer(n)
    for(i in 1:n){
        path4[i] = cpps[tm1==npath[i] & t1==path[i]]
    }
    return(path4)
}

convert8to4 <- function(path){
    t1 = c(1,2,4,2,3,1,3,1)
    path8 = t1[path+1]
    return(path8)
}

beamSearchWrap <- function(pmats, w0, y, Npart){
    if(!is.matrix(y)) dim(y) = c(1,length(y))
    test = beamSearch(pmats$a0, pmats$P0, w0, pmats$dt, pmats$ct,
                      pmats$Tt, pmats$Zt,
                      pmats$Rt, pmats$Qt, pmats$GGt, y, pmats$transMat, Npart)
    return(test)
}

kf11 <- function(pmats, y, iter, s){
    if(!is.matrix(y)) dim(y) = c(1,length(y))
    if(dim(pmats$a0)[2] !=1){
        a0 = pmats$a0[,s,drop=FALSE]
        P0 = pmats$P0[,s,drop=FALSE]
    }else{
        a0 = pmats$a0
        P0 = pmats$P0
    }
    dt = matrix(pmats$dt[,s,iter],ncol=1)
    ct = matrix(pmats$ct[,s,],ncol=1)
    Tt = matrix(pmats$Tt[,s,iter],ncol=1)
    Zt = matrix(pmats$Zt[,s,],ncol=1)
    Rt = matrix(pmats$Rt[,s,iter],ncol=1)
    Qt = matrix(pmats$Qt[,s,],ncol=1)
    HHt = HHcreate(Rt, Qt, sqrt(nrow(Rt)), sqrt(nrow(Qt)))
    GGt = matrix(pmats$GGt[,s,],ncol=1)
    out = kf1step(a0, P0, dt, ct, Tt, Zt, HHt, GGt, matrix(y[,iter],ncol=1))
    return(out)
}

kf1wrap <- function(ssMod, y, a0=matrix(0,nrow=nrow(ssMod$dt)), 
                    P0=diag(nrow(ssMod$dt))){
    n = ncol(y)
    m = nrow(P0)
    P0 = matrix(P0,ncol=1)
    at = matrix(0, m, n)
    Pt = matrix(0, m^2, n)
    lik = double(n)
    print(a0)
    print(P0)
    Tt = matrix(ssMod$Tt,ncol=1)
    Zt = matrix(ssMod$Zt,ncol=1)
    HHt = matrix(ssMod$HHt,ncol=1)
    GGt = matrix(ssMod$GGt,ncol=1)
    for(i in 1:n){
        step = kf1step(a0, P0, ssMod$dt, ssMod$ct, Tt, Zt, HHt, GGt, y[,i,drop=FALSE])
        at[,i] = step$a1
        Pt[,i] = step$P1
        lik[i] = step$lik
        a0 = step$a1
        P0 = step$P1
    }
    return(list(at=at,Pt=Pt,lik=lik))
}

logtransitions <- function(path, transmat){
    n = length(path)
    tm1 = c(1,path[-n])
    return(log(transmat[cbind(tm1,path)]))
}
