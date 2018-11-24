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


full_multinom <- function(p){
    return(c(p, 1-sum(p)))
}

KLtotal <- function(p1, p2, symmetric=TRUE){
    p1 = as.vector(unlist(p1))
    p2 = as.vector(unlist(p2))
    observed = KLgaussian(0,0,p1[1],p2[1],symmetric)#sig2eps
    constant = KLgaussian(p1[2],p2[2],p1[5],p2[5],symmetric)#mu1 and sig2tempo
    accel = KLgaussian(p1[3],p2[3],p1[6],p2[6],symmetric)#mu2 and sig2acc
    stress = KLgaussian(p1[4],p2[4],p1[7],p2[7],symmetric)#mu3 and sig2stress
    state1 = KLmultinom(full_multinom(p1[c(8,9,12)]),#p's out of state 1
                        full_multinom(p2[c(8,9,12)]),
                        symmetric)
    state2 = KLmultinom(full_multinom(p1[c(10,13)]),#p's out of state 2
                        full_multinom(p2[c(10,13)]),
                        symmetric)
    state3 = KLmultinom(full_multinom(p1[11]),#p's out of state 3
                        full_multinom(p2[11]))
    return(observed+constant+accel+stress+state1+state2+state3)
}

Dist <- function(pvec_ml){
    n = nrow(pvec_ml)
    out = array(dim = c(n,n))
    for(i in 1:n){
        for(j in 1:i){
            out[i,j] = out[j,i] = KLtotal(pvec_ml[i,],pvec_ml[j,])
            if(out[i,j] > 1000){
                out[i,j] = out[j,i] = 1000
            }
        }
    }
    out = dist(out)
    return(out)
}
