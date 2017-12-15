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

kf11 <- function(pmats, y, iter, s, a0=NULL, P0=NULL){
    if(!is.matrix(y)) dim(y) = c(1,length(y))
    if(is.null(a0) & is.null(P0)){
        if(dim(pmats$a0)[2] !=1){
            a0 = pmats$a0[,s,drop=FALSE]
            P0 = pmats$P0[,s,drop=FALSE]
        }else{
            a0 = pmats$a0
            P0 = pmats$P0
        }
    }
    dt = matrix(pmats$dt[,s,iter],ncol=1)
    ct = matrix(pmats$ct[,s,],ncol=1)
    Tt = matrix(pmats$Tt[,s,iter],ncol=1)
    Zt = matrix(pmats$Zt[,s,],ncol=1)
    Rt = matrix(pmats$Rt[,s,iter],ncol=1)
    Qt = matrix(pmats$Qt[,s,],ncol=1)
    HHt = HHcreate(Rt, Qt, sqrt(nrow(Rt)), sqrt(nrow(Qt)))
    GGt = matrix(pmats$GGt[,s,],ncol=1)
    out = kf1stepR(a0, P0, dt, ct, Tt, Zt, HHt, GGt, matrix(y[,iter],ncol=1))
    return(out)
}

predictiveDensities <- function(timepoint, pmats, path, y){
    if(!is.matrix(y)) dim(y) = c(1,length(y))
    n = length(y)
    d = 1 # nrow(y)
    m = nrow(pmats$a0)
    K = ncol(pmats$a0)
    if(timepoint>n+1) stop("Can't examine predictive densities after last obs.")
    pathout = pathStuff(pmats, path[1:(timepoint-1)], 
                        y[,1:(timepoint-1),drop=FALSE])
    cs = path[timepoint-1]
    HHt = HHcreate(pmats$Rt[,,timepoint-1], pmats$Qt[,,,drop=TRUE], m, 
                   sqrt(nrow(pmats$Qt)))
    dim(pmats$ct) = c(d,K)
    dim(pmats$Zt) = c(m,K)
    dim(pmats$GGt) = c(d,K)
    stepahead = dpf(cs, matrix(1, 1,1), K+1, pmats$transMat, 
                    a0 = pathout$at[,timepoint-1,drop=FALSE],
                    P0 = pathout$Pt[,,timepoint-1],
                    dt = pmats$dt[,,timepoint-1],
                    ct = pmats$ct, Tt=pmats$Tt[,,timepoint-1],
                    Zt = pmats$Zt, GGt=pmats$GGt, HHt=HHt, y[,timepoint-1])
    pos = stepahead$newstates+1
    ss = c(1,2,4,2,3,1,3,1)
    preds = double(length(pos)) # depends on d, assumed 1
    predsd = double(length(pos))
    for(i in 1:length(pos)){
        cti = pmats$ct[,pos[i]]
        Zti = pmats$Zt[,pos[i]]
        GGti = pmats$GGt[,pos[i]]
        preds[i] = cti + Zti %*% stepahead$a1[,i]
        predsd[i] = sqrt(t(Zti) %*% matrix(stepahead$P1[,i], m, m) %*% Zti + GGti)
    }
    ## Prepare for plotting
    xlims = c(0,max(y,preds+3*predsd))
    ylims = c(0, max(dnorm(0,0,predsd)*stepahead$newW))
    par(mar=c(4,3,1,1))
    plot(NA,NA, xlim = xlims, ylim = ylims,
         main = paste0('timepoint = ',timepoint),ylab='',xlab='')
    for(i in 1:length(pos)){
        curve(dnorm(x,preds[i],predsd[i])*stepahead$newW[i], 
              from=min(xlims), to=max(xlims), add = TRUE, col=ss[pos[i]])
    }
    abline(v=y[,timepoint],col=ss[path[timepoint]+1])
}
    
        
        

kf1wrap <- function(ssMod, y, a0=matrix(0,nrow=nrow(ssMod$dt)), 
                    P0=diag(nrow(ssMod$dt))){
    n = ncol(y)
    m = nrow(P0)
    P0 = matrix(P0,ncol=1)
    at = matrix(0, m, n+1)
    Pt = matrix(0, m^2, n+1)
    at[,1] = a0
    Pt[,1] = P0
    lik = double(n)
    # print(a0)
    # print(P0)
    Tt = matrix(ssMod$Tt,ncol=1)
    Zt = matrix(ssMod$Zt,ncol=1)
    HHt = matrix(ssMod$HHt,ncol=1)
    GGt = matrix(ssMod$GGt,ncol=1)
    for(i in 1:n){
        step = kf1step(a0, P0, ssMod$dt, ssMod$ct, Tt, Zt, HHt, GGt, y[,i,drop=FALSE])
        at[,i+1] = step$a1
        Pt[,i+1] = step$P1
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

# kf1stepR(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt) {
#     d = nrow(a0)
#     m = nrow(yt)
#     dim(P0) = c(d,d)
#     dim(Tt) = c(d,d)
#     dim(Zt) = c(m,d)
#     dim(HHt) = c(d,d)
#     dim(GGt) = c(m,m)
#     
#     ## KF
#     pred = ct + Zt %*% a0
#     vt = yt - pred
#     Ft = GGt + Zt %*% P0 %*% t(Zt.t)
#     Ftinv = solve(Ft)
#     Kt = P0 %*% t(Zt.t) %*% Ftinv
#     a1 += Kt * vt;
#     arma::mat P1 = P0 - Kt * Zt * P0;
#     a1 = dt + Tt * a1;
#     P1 = HHt + Tt * P1 * Tt.t();
#     
#     // Calculate likelihood
#     double ftdet = arma::det(Ftinv);
#     ftdet = 1 / std::abs(ftdet);
#     double mahalanobis = arma::as_scalar(vt.t() * Ftinv * vt);
#     mahalanobis += yt.size()* log(2*M_PI) + log(ftdet);
#     double lik = exp(-1.0*mahalanobis/2);
#     P1.reshape(d*d, 1);
#     return List::create(Named("a1") = a1,
#                         Named("P1") = P1,
#                         Named("lik") = lik,
#                         Named("pred") = pred);
# }