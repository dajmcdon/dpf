

# more likely to use states 2,3,4?
lessvar <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = dgamma(theta[1], shape=20, scale=10, log = TRUE)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=1, log = TRUE)
  sig2tempo = dgamma(theta[5], shape=40, scale=10, log=TRUE)
  #sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  #sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,8,2))
  p22 = ddirichlet(p2s, alpha=c(10,1,4))
  p31 = ddirichlet(p3s, alpha=c(5,7,3))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, #sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}



## alternative priors
# more state 2

# invgaussian_pars <- function(m, v) c(shape=m^2/v +2, rate=m*(m^2/v+1))
library(invgamma)
#curve(dgamma(x, 40, scale=10,log=TRUE), 0, 1000)
#curve(dinvgamma(x, 42, rate=16400, log=TRUE), 0, 1000, col=2, add=TRUE)

logprior_igs <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = dinvgamma(theta[1], 42, scale=16400, log = TRUE)
  mu1 = dinvgamma(theta[2], samp_mean^2/100+2, rate=samp_mean*(samp_mean^2/100+1), log = TRUE)
  mu2 = dinvgamma(-theta[3], 17, 160, log = TRUE)
  mu3 = dinvgamma(-theta[4], 22, 840, log = TRUE)
  sig2tempo = dinvgamma(theta[5], 42, 16400, log=TRUE)
  #sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  #sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,8,2))
  p22 = ddirichlet(p2s, alpha=c(10,1,4))
  p31 = ddirichlet(p3s, alpha=c(5,7,3))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, #sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}


## uniform
logprior_uvar <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = ifelse(theta[1] > 0, 0, -Inf)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=2, log = TRUE)
  sig2tempo = ifelse(theta[5] > 0, 0, -Inf)
  #sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  #sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(85,5,8,2))
  p22 = ddirichlet(p2s, alpha=c(10,1,4))
  p31 = ddirichlet(p3s, alpha=c(5,7,3))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, #sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}


## Uniform transitions
logprior_uprobs <- function(theta, samp_mean=132){
  p1s = c(theta[c(8,9)], 1-sum(theta[c(8,9,12)]), theta[12])
  p2s = c(theta[10], 1-sum(theta[c(10,13)]), theta[13])
  p3s = c(theta[11], 1-sum(theta[c(11,14)]), theta[14])
  sig2eps = dgamma(theta[1], shape=40, scale=10, log = TRUE)
  mu1 = dgamma(theta[2], samp_mean^2/100, scale=100/samp_mean, log = TRUE)
  mu2 = dgamma(-theta[3], 15, scale=2/3, log = TRUE)
  mu3 = dgamma(-theta[4], 20, scale=2, log = TRUE)
  sig2tempo = dgamma(theta[5], shape=40, scale=10, log=TRUE)
  #sig2acc = dgamma(theta[6], shape=1, scale=1, log=TRUE)
  #sig2stress = dgamma(theta[7], shape=1, scale=1, log=TRUE)
  p1 = ddirichlet(p1s, alpha=c(1,1,1,1))
  p22 = ddirichlet(p2s, alpha=c(1,1,1))
  p31 = ddirichlet(p3s, alpha=c(1,1,1))
  lp = sum(sig2eps, mu1, mu2, mu3, 
           sig2tempo, #sig2acc, sig2stress, 
           p1, p22, p31)
  lp
}
