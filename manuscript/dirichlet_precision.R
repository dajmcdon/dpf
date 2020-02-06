# Source functions for clustering dirichlet parameters

dirichlet_precision <- function(alpha) {
  a0 = sum(alpha)
  A = - tcrossprod(alpha)
  D = diag(alpha)*a0
  kern = D + A
  mult = a0^2*(a0+1)
  MASS::ginv(kern) * mult
}

dist_calc <- function(x, y, prec = diag(ncol(pvec))){
  t((x-y)) %*% prec %*% (x-y)
}

full_multinom <- function(p) c(p, 1-sum(p))

Dmats <- function(pvec, precisions = list(
  1, 1, 1, 1, 1, 
  #1/4000, 1/100, 1/(15*(2/3)^2), 1/(20*4), 1/4000,
  dirichlet_precision(c(85,5,2,8)),
  dirichlet_precision(c(4,10,1)),
  dirichlet_precision(c(5,3,7))),
  params = list(1,2,3,4,5,c(6,7,10),c(8,11),c(9,12)),
  dirs = c(F,F,F,F,F,T,T,T)
  ){
  stopifnot((L <-length(precisions))==length(params), L==length(dirs))
  Ds = list()
  for(l in 1:L){
    if(!dirs[l]){
      M = pvec[, params[[l]], drop=FALSE]
    } else {
      M = t(apply(pvec[,params[[l]]], 1, full_multinom))
    }
    Ds[[l]] = usedist::dist_make(
      M,
      distance_fcn = function(x,y) dist_calc(x,y, prec = precisions[[l]])
      )
  }
  Ds = lapply(Ds, as.matrix)
  Ds = lapply(Ds, function(x) x/max(x))
  Ds
}