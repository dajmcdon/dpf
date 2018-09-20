

#' Dirichlet random variables
#'
#' @param x vector of probabilities supported on the K-1 simplex. ddirichlet is vectorized over x.
#' @param n number of observations
#' @param K number of components
#' @param alpha vector of parameters. Must be positive. These functions are not vectorized over alpha.
#'
#' @return \code{ddirichlet} gives the density on the log scale, unnormalized by default. \code{rdirichlet} generates random variates.
#' 
#' @export
#'
#' @examples
#' 
#' randos = rdirichlet(10, 3, c(7,2,1))
#' density = ddirichlet(c(.5,.5,.5), c(7,2,1))
rdirichlet <- function(n, K=length(alpha), alpha = rep(1,K)){
  stopifnot(n > 0, length(alpha)==K, all(alpha>0))
  mat = matrix(rgamma(n*K, alpha), nrow = n, byrow = TRUE)
  ys = rowSums(mat)
  mat = mat/ys
  mat
}

#' @rdname rdirichlet
#' @export
ddirichlet <- function(x, alpha, log=TRUE, normalize=FALSE){
  p = length(x)
  stopifnot(length(alpha)==p, all(alpha>0))
  if(any(x>1) || any(x<0) || !isTRUE(all.equal(1,sum(x)))){
    return(ifelse(log,-Inf,0))
  }
  kern = sum(log(x) * alpha) # on the log scale, no constant
  if(normalize){
    const = -sum(lgamma(alpha)) + lgamma(sum(alpha))
    kern = kern + const
  }
  if(!log) kern = exp(kern)
  kern
}