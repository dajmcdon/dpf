KLdirichlet <- function(a, b, symmetric = TRUE){
  a0 = sum(a)
  b0 = sum(b)
  f = sum( (a-b)*(digamma(a)-digamma(a0)) )
  if(symmetric){
    g = sum( (b-a)*(digamma(b)-digamma(b0)) ) 
    f = (f + g)/2
  }else{
    back = sum(lgamma(b)-lgamma(a)) + lgamma(a0)-lgamma(b0)
    f = f + back
  }
  f
}

KLbeta <- function(a, b, symmetric=TRUE){
  KLdirichlet(a,b,symmetric)
}

KLgaussian <- function(m1, m2, s1, s2, symmetric=TRUE){
  f = log(s2) - log(s1)
  g = (s1 + (m1-m2)^2)/(2*s2)
  sc = -1/2
  if(symmetric) g = (g + (s2 + (m2-m1)^2)/(2*s1))/2 + sc
  else g = f + g + sc
  g
}

KLmultinom <- function(a, b, symmetric=TRUE){
  one = sum(a*log(a)) - sum(a*log(b))
  if(symmetric){
    two = sum(b*log(b)) - sum(b*log(a))
    one = (one + two)/2
  }
  one
}


#' Compute the squared Hellinger distance between two gaussian distributions
#'
#' @param m1 the mean of the first distribution
#' @param m2 the mean of the second distribution
#' @param s1 the standard deviation of the first distribution
#' @param s2 the standard deviation of the second distribution
#' @param subtract1 a boolean indicating whether the squared distance should be computed instead of one minus the squared distance
#'
#' @return the squared Hellinger distance
#' @export
H2gauss <- function(m1, m2, s1, s2, subtract1=FALSE){
  sig1 = sqrt(s1)
  sig2 = sqrt(s2)
  const = sqrt(2*sig1*sig2/(s1+s2))
  ee = (m1-m2)^2/(s1+s2)
  ret = const * exp(-ee/4)
  if(subtract1) ret = 1 - ret
  ret
}


#' Compute the squared Hellinger distance between two multinomial distributions
#'
#' @param p the probability vector for the first distribution
#' @param q the probability vector for the second distribution
#' @param subtract1 a boolean indicating whether the squared distance should be computed instead of one minus the squared distance
#'
#' @return the squared Hellinger distance
#' @export
H2multinom <- function(p, q, subtract1=FALSE){
  ret = sum(sqrt(p*q))
  if(subtract1) ret = 1 - ret
  ret
}

full_multinom <- function(p) c(p, 1-sum(p))


#' Compute the squared Hellinger distance between two models specified by passing the parameters to musicModel
#'
#' @param p1 parameters for the first model
#' @param p2 parameters for the second model
#'
#' @return the squared Hellinger distance
#' @export
H2total <- function(p1, p2){
  p1 = as.vector(unlist(p1))
  p2 = as.vector(unlist(p2))
  observed = H2gauss(0,0,p1[1],p2[1])#sig2eps
  constant = H2gauss(p1[2],p2[2],p1[5],p2[5])#mu1 and sig2tempo
  accel = H2gauss(p1[3],p2[3],p1[6],p2[6])#mu2 and sig2acc
  stress = H2gauss(p1[4],p2[4],p1[7],p2[7])#mu3 and sig2stress
  state1 = H2multinom(full_multinom(p1[c(8,9,12)]),#p's out of state 1
                      full_multinom(p2[c(8,9,12)]))
  state2 = H2multinom(full_multinom(p1[c(10,13)]),#p's out of state 2
                      full_multinom(p2[c(10,13)]))
  state3 = H2multinom(full_multinom(p1[11]),#p's out of state 3
                      full_multinom(p2[11]))
  return(1 - prod(
    c(observed,constant,accel,stress,state1,state2,state3)
    ))
    
}


#' create a distance matrix
#'
#' @param x data matrix
#' @param dist_fun distance function
#' @param returnmat boolean indicating whether to ensure that the return type is a matrix
#'
#' @return a distance matrix
#' @export
Dist <- function(x, dist_fun=H2total, returnmat = TRUE){
  d = usedist::dist_make(x, dist_fun)
  if(returnmat) d = as.matrix(d)
  d
}

full_samp_kde <- function(x, y, indices){
  outL = list()
  outL$fx = ks::kde(x[,indices], binned=FALSE, eval.points = y[,indices])
  outL$fy = ks::kde(y[,indices], binned=FALSE, eval.points = x[,indices])
  outL
}

loo_kde <- function(x, fx, indices){
  d = length(indices)
  n = nrow(x)
  fminus = double(n)
  if(d==1){
    h = fx$h
    for(i in 1:n){
      fminus[i] = ks::kde(x[-i,indices],
                          eval.points = x[i,indices], h=h)$estimate
    }
  } else {
    H = fx$H
    for(i in n){
      fminus[i] = ks::kde(x[-i,indices],
                          eval.points = x[i,indices], H=H)$estimate
    } 
  }
  fminus
}

stat_kern <- function(x, y, indices, lb = 1e-12){
  fxfy = full_samp_kde(x, y, indices)
  looX = loo_kde(x, fxfy$fx, indices)
  looY = loo_kde(y, fxfy$fy, indices)
  
  kern = list()
  kern$x = sqrt( pmax(fxfy$fx$estimate, lb) / pmax(looX, lb) )
  kern$y = sqrt( pmax(fxfy$fy$estimate, lb) / pmax(looY, lb) )

  kern
}

perc_overlap <- function(a, b){
  if(a[1] > b[2] || b[1] > a[2]) return(0)
  da = a[2] - a[1]
  db = b[2] - b[1]
  os = max(a[1], b[1])
  ol = min(a[2], b[2])
  (ol-os) / max(da,db)
}

kill_check <- function(x, y, ind, kill = 1e-6){
  # determine if support overlap is too small
  xmat = x[,ind,drop=FALSE]
  ymat = y[,ind,drop=FALSE]
  xr = apply(xmat, 2, range)
  yr = apply(ymat, 2, range)
  if(any(xr[2,]-xr[1,] < kill , yr[2,]-yr[1,] < kill)) return(TRUE)
  xr = split(xr, rep(ind, each=2))
  yr = split(yr, rep(ind, each=2))
  if(any(mapply(FUN=perc_overlap, xr, yr) < kill)) return(TRUE)
  FALSE
}

clamp <- function(x, lb=0, ub=1) min(ub, max(lb, x))
  

#' Estimate the Hellinger divergence between two distributions
#'
#' @param x a sample from the first distribution (n x d)
#' @param y a sample from the second distribution (n x d)
#' @param index.list list of columns to be used (default 1), see details
#' @param lb lower-bound for the density truncation
#' @param all return the kde for x and y rather than the divergence estimate
#' 
#' @details 
#' 
#' This estimates the Hellinger divergence between samples from two 
#' distributions. The methodology implements the estimator in 
#' "Influence Functions for Machine Learning: Nonparametric Estimators for 
#' Entropies, Divergences and Mutual Informations" by Kandasamy, Krishnamurthy,
#' Poczos, Wasserman, and Robins. See their Matlab functions. Here, we use only
#' the Hellinger divergence.
#' 
#' Each sample (x and y) is assumed to be a matrix with n samples in d columns.
#' The estimator can only use up to 6 dimensions jointly. The list breaks up the
#' d dimensions into smaller collections. Generically,
#' 
#' \deqn{H = 2 - \int \sqrt(p_x p_y).}
#' 
#' If \eqn{p_x = p_{x1} p_{x2}} and \eqn{p_y = p_{y1}p_{y2}}, then 
#' \deqn{H = 2 - \int \sqrt(p_{x1} p_{y1}) \int \sqrt(p_{x1} p_{y1}).} So, for example, to 
#' treat the first column as independent from columns 2 and 3, use
#' \code{index.list = list(1,2:3)}.
#' 
#'
#' @return the Hellinger divergence
#' 
#' @examples
#' x = matrix(rnorm(3*1000)) %*% matrix(c(1,0,0,0,2,.5,0,.5,4),ncol=3)
#' y = matrix(rnorm(3*1000)) %*% matrix(c(2,0,0,0,1,.5,0,.5,4),ncol=3)
#' hellinger_kde(x, y, list(1, 2:3))
#' @export
hellinger_kde <- function(x, y, index.list = 1, lb = 1e-6, kill = 1e-6){
  
  if(!is.list(index.list)) index.list = list(index.list)
  folds = length(index.list)
  
  out = 2*folds
  for(f in 1:folds){
    if(kill_check(x,y,index.list[[f]])){
      kern = list(x=0, y=0)  
    } else {
      kern = stat_kern(x, y, index.list[[f]])
    }
    out = out - clamp(mean(kern$x)) - clamp(mean(kern$y))
  }
  out
}
