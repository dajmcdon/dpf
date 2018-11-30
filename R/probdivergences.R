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


#' Title
#'
#' @param m1 
#' @param m2 
#' @param s1 
#' @param s2 
#' @param subtract1 
#'
#' @return
#' @export
#'
#' @examples
H2gauss <- function(m1, m2, s1, s2, subtract1=FALSE){
  sig1 = sqrt(s1)
  sig2 = sqrt(s2)
  const = sqrt(2*sig1*sig2/(s1+s2))
  ee = (m1-m2)^2/(s1+s2)
  ret = const * exp(-ee/4)
  if(subtract1) ret = 1 - ret
  ret
}


#' Title
#'
#' @param p 
#' @param q 
#' @param subtract1 
#'
#' @return
#' @export
#'
#' @examples
H2multinom <- function(p, q, subtract1=FALSE){
  ret = sum(sqrt(p*q))
  if(subtract1) ret = 1 - ret
  ret
}

full_multinom <- function(p) c(p, 1-sum(p))


#' Title
#'
#' @param p1 
#' @param p2 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param x 
#' @param dist_fun 
#' @param returnmat 
#'
#' @return
#' @export
#'
#' @examples
Dist <- function(x, dist_fun=H2total, returnmat = TRUE){
  d = usedist::dist_make(x, dist_fun)
  if(returnmat) d = as.matrix(d)
  d
}

