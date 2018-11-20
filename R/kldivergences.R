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
