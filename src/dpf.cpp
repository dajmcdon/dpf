#include <RcppArmadillo.h>
#include "samplers.h"
#include <set>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


struct KFOUT{
  arma::mat att; // a_t|t
  arma::mat at;  // a_t+1 (but this is done at the beginning of the KF
                 // instead of the end)
  arma::mat Ptt; // same deal
  arma::mat Pt;
  double lik;
  arma::mat pred; // f(at) to estimate yt
};

    

KFOUT kf1step(arma::mat a0, arma::mat P0, arma::mat dt,
             arma::mat ct, arma::mat Tt,
             arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
  // Reshape matrices
  int m = a0.size();//number of elements.
  int d = yt.size();
  P0.reshape(m,m);
  Tt.reshape(m,m);
  Zt.reshape(d,m);
  HHt.reshape(m,m); // State var (see fkf vs Durbin)
  GGt.reshape(d,d); // Obs var
  arma::mat Ptemp(m,m);
  arma::mat atemp(m,1);

  // KF
  // Update a_t-1|t-1 to a_t using CURRENT state
  atemp = dt + Tt * a0;
  Ptemp = HHt + Tt * P0 * Tt.t();
  
  // Make predictions using a_t
  arma::mat pred = ct + Zt * atemp;
  arma::mat vt = yt - pred;
  arma::mat Ft = GGt + Zt * Ptemp * Zt.t();
  arma::mat Ftinv = arma::inv(Ft);
  arma::mat Kt = Ptemp * Zt.t() * Ftinv;
  
  // Incorporate current information into a_t|t
  arma::mat a1 = atemp + Kt * vt;
  arma::mat P1 = Ptemp - Ptemp * Zt.t() * Kt.t();

  // Calculate likelihood
  double ftdet = arma::det(Ftinv);
  ftdet = 1 / std::abs(ftdet);
  double mahalanobis = arma::as_scalar(vt.t() * Ftinv * vt);
  mahalanobis += yt.size()* log(2*M_PI) + log(ftdet);
  double lik = exp(-1.0*mahalanobis/2);
  P1.reshape(m*m, 1);
  KFOUT output = {a1, atemp, P1, Ptemp, lik, pred};
  return output;
}



  

KFOUT ks1step(arma::mat r1, arma::mat N1, 
              arma::mat pred, arma::mat Pt, 
              arma::mat ct, arma::mat Tt,
              arma::mat Zt, arma::mat GGt, arma::mat yt) {
    
    // Reshape matrices
    int m = r1.size();//number of elements.
    int d = yt.size();
    Pt.reshape(m,m);
    Tt.reshape(m,m);
    Zt.reshape(d,m);
    GGt.reshape(d,d); // Obs var
    
    // KS
    arma::mat vt = yt - pred;
    arma::mat Ft = GGt + Zt * Pt * Zt.t();
    arma::mat Ftinv = arma::inv(Ft);
    arma::mat Lt = Tt - Tt * Pt * Zt.t() * Ftinv * Zt;
    arma::mat r0 = Zt.t() * Ftinv * vt + Lt.t() * r1;
    arma::mat N0 = Zt.t() * Ftinv * Zt + Lt.t() * N1 * Lt;
    
    KFOUT output = {r0, r0, N0, N0, 1, r0}; // Note: we don't need anything but r0 and N0
    return output;
}


List dpf(arma::uvec currentStates, arma::colvec w, int N,
         arma::mat transProbs,
         arma::mat a0, arma::mat P0,
         arma::mat dt, arma::mat ct, arma::mat Tt, arma::mat Zt, //Matrices stored in columns of these variables. Re turned into matrices in other functions.
         arma::mat HHt, arma::mat GGt, arma::vec yt){
  int npart = currentStates.size();
  int m = a0.n_rows;                            //dimension of each particle (the continuous part)
  int K = dt.n_cols;                            //number of possible discrete states
  int mm = m*m;
  arma::cube a1(m, K, npart);
  arma::cube P1(mm, K, npart);                   
  arma::mat lik(npart, K);
  arma::mat a11(m, K);
  arma::mat P11(mm, K);
  for(int part=0; part < npart; part++){
    for(int k=0;  k < K; k++){   //both loops together have effect of for each particle
      KFOUT kfout = kf1step(a0.col(part), P0.col(part), dt.col(k), ct.col(k),
                                       Tt.col(k), Zt.col(k), HHt.col(k), GGt.col(k),
                                       yt);
      arma::mat atmp = kfout.att;
      arma::mat Ptmp = kfout.Ptt;
      a11.col(k) = atmp;
      P11.col(k) = Ptmp;
      lik(part,k) = kfout.lik;
    }
    a1.slice(part) = a11;
    P1.slice(part) = P11;
  }
  arma::mat testLik = lik;
  lik %= transProbs.rows(currentStates);
  lik.each_col() %= w;
  w = arma::vectorise(lik);
  w = arma::normalise(w,1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec positive = arma::find(w);
  int nkeep = positive.size();
  arma::colvec newW = w(positive);
  arma::mat aout(m, nkeep);
  arma::mat Pout(mm, nkeep);
  arma::uvec newstates(nkeep);
  arma::uvec oldstates(nkeep);
  int i=0;
  for(int k=0; k<K; k++){
    for(int part=0; part < npart; part++){
      if(w(part + k*npart)>0){
        newstates(i) = k;
        oldstates(i) = part;
        arma::mat a1tmp = a1.subcube(0,k,part,m-1,k,part);
        arma::mat P1tmp = P1.subcube(0,k,part,mm-1,k,part);
        aout.col(i) = a1tmp;
        Pout.col(i) = P1tmp;
        i++;
      }
    }
  }
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("oldstates") = oldstates,
                      Named("newstates") = newstates,
                      Named("newW") = newW,
                      Named("testLik") = testLik);
}


//' Fast likelihood evaluation given parameters and discrete states
//' 
//' This function quickly computes the (negative) log likelihood for a possibly
//' time varying state-space model. The goal is to enable ML (or Bayesian) estimation
//' of parameters without extra computational overhead
//' 
//' @param pmats a list of parameter matrices for the kalman filter. 
//' This can either be the output of \code{\link{musicModel}}, or a list with 
//' the same names as the output of \code{\link{musicModel}}. This list must contain
//' the following elements: \code{a0, P0, dt, ct, Tt, Zt, HHt, GGt}. See the details below.
//' @param path vector giving the path for hidden discrete states
//' @param y observations, each time point in a column
//' 
//' @return the negative log-likelihood
//'
//' @details 
//' \describe{
//' \item{a0}{a pxd matrix of the initial means of the hidden state. The j'th column corresponds to the initial mean when starting in the j'th discrete state.}
//' \item{P0}{a (p^2)xd matrix of the initial covariances of the hidden state. The j'th column corresponds to the initial covariances stored columnwise when starting in the j'th discrete state.}
//' \item{dt}{a pxdxn cube of state intercepts. The j'th column of the i'th slice corresponds to the intercept specified by the j'th discrete state at time i.}
//' \item{ct}{a kxdx1 cube of observation intercepts. The j'th column corresponds to the intercept specified by the j'th discrete state.}
//' \item{Tt}{a (p^2)xdxn cube of state slopes. The j'th column of the i'th slice corresponds to the slope matrix stored columnwise of the j'th discrete state at time i.}
//' \item{Zt}{a pkxdx1 cube of obvervation slopes. The j'th column corresponds to the slope matrix stored columnwise of the j'th discrete state.}
//' \item{HHt}{a (p^2)xdxn cube of state covariances. The j'th column of the i'th slice corresponds to the covariance matrix stored columnwise of the j'th discrete state at time i.}
//' \item{GGt}{a (k^2)xdx1 cube of observation covariances. The j'th column corresponds to the covariance matrix stored columnwise of the j'th discrete state.}
//' }
//' 
//' @examples
//' data(tempos)
//' theta = c(426.69980736, 136.33213703, -11.84256691, -34.82234559, 
//'           439.37886221, 1, 1, 0.84916635, 0.04611644, 0.74119571, 
//'           0.43966082, 0.02116317, 0.24513563, 0.17253254)
//' y = matrix(tempos[,'Richter_1976'], 1)
//' lt = diff(c(tempos$note_onset, 61))
//' pmats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14], 
//'                   c(132,0), c(400,10)) # prior means and variances on X_1
//' beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0,0,0,0,0,0,0), dt, ct, Tt, Zt,
//'             HHt, GGt, y, transMat, 200))
//' bestpath = with(beam, paths[which.max(weights),])
//' getloglike(pmats, bestpath, y)
//' 
//' @export 
// [[Rcpp::export]]
double getloglike(List pmats, arma::uvec path, arma::mat y){
  arma::mat a0 = pmats["a0"];
  arma::mat P0 = pmats["P0"];
  arma::cube dt = pmats["dt"];
  arma::cube ct = pmats["ct"];
  arma::cube Tt = pmats["Tt"];
  arma::cube Zt = pmats["Zt"];
  arma::cube HHt = pmats["HHt"];
  arma::cube GGt = pmats["GGt"];
  
  arma::uword n = y.n_cols;
  arma::uword m = a0.n_rows;
  arma::uword mm = m*m;
  arma::uword d = y.n_rows;
  arma::uword dm = d*m;
  arma::uword dd = d*d;
  
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword HHtvar = HHt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  
  arma::mat aa0 = a0.col(path(0));
  arma::mat PP0 = reshape(P0.col(path(0)), m, m);
  
  double llik = 0;
  double liktmp;
  
  for(arma::uword iter=0; iter<n; iter++){
    arma::uword s = path(iter);
    KFOUT step = kf1step(aa0, PP0, dt.subcube(0,s,iter*dtvar,arma::size(m,1,1)),
                                    ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)), 
                                    Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                                    Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)), 
                                    HHt.subcube(0,s,iter*HHtvar,arma::size(mm,1,1)), 
                                    GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                                    y.col(iter));
    
    aa0 = step.att;
    PP0 = step.Ptt;
    liktmp = step.lik;
    llik -= log(liktmp);
  }
  return llik;
}
    


List initializeParticles(arma::vec w0, int N, arma::mat a0, arma::mat P0,
                         arma::mat dt, arma::mat ct, arma::mat Tt, arma::mat Zt,
                         arma::mat HHt, arma::mat GGt, arma::vec yt){
  arma::uvec particles = arma::find(w0);
  int npart = particles.n_elem;
  if(npart == 0){return List::create(Named("BadPars") = 1);}
  arma::vec weights = w0(particles);
  arma::mat a11(a0.n_rows, npart);
  arma::mat P11(P0.n_rows, npart);
  arma::vec lik(npart);
  for(int i = 0; i < npart; i++){
    int p = particles(i);
    KFOUT kfout = kf1step(a0.col(p), P0.col(p), dt.col(p), ct.col(p),
                          Tt.col(p), Zt.col(p), HHt.col(p), GGt.col(p), yt);
    arma::mat atmp = kfout.att;
    arma::mat Ptmp = kfout.Ptt;
    a11.col(i) = atmp;
    P11.col(i) = Ptmp;
    lik(i) = kfout.lik;
  }
  lik %= weights;
  arma::vec w = arma::normalise(lik, 1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec nz = arma::find(w);
  arma::vec newW = w(nz);
  arma::uvec newstates = particles(nz);
  arma::mat aout = a11.cols(nz);
  arma::mat Pout = P11.cols(nz);
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("newstates") = newstates,
                      Named("newW") = newW);
}


//' Greedy HMM estimation given continuous hidden states
//' 
//' This function maintains a "beam" of size \code{N} representing the highest likelihood
//' states up to time \code{t}. For a typical switching state space model, ML evaluation
//' requires computing each of the \code{d^N} possible paths through the discrete
//' state space. Beam Search greedily approximates these computations be propagating
//' the most likely states at each time point.
//' 
//' @param a0 a pxd matrix of state prior means, where the i'th column is the mean for the i'th state
//' @param P0 a (p^2)xd state prior covariance matrix, where the i'th column is the covariance matrix stored columnwise for the i'th state
//' @param w0 a vector specifying the prior probability of starting in each of the d discrete states
//' @param dt a pxdxn (or pxdx1 if all n slices are the same) cube of state intercepts. The j'th column of the i'th slice corresponds to the intercept specified by the j'th discrete state at time i.
//' @param ct a kxdxn (or kxdx1 if all n slices are the same) cube of observation intercepts. The j'th column of the i'th slice corresponds to the intercept specified by the j'th discrete state at time i.
//' @param Tt a (p^2)xdxn (or (p^2)xdx1 if all n slices are the same) cube of state slopes. The j'th column of the i'th slice corresponds to the slope matrix stored columnwise of the j'th discrete state at time i.
//' @param Zt a pkxdxn (or pkxdx1 if all n slices are the same) cube of obvervation slopes. The j'th column of the i'th slice corresponds to the slope matrix stored columnwise of the j'th discrete state at time i.
//' @param HHt a (p^2)xdxn (or (p^2)xdx1 if all n slices are the same) cube of state covariances. The j'th column of the i'th slice corresponds to the covariance matrix stored columnwise of the j'th discrete state at time i.
//' @param GGt a (k^2)xdxn (or kxdx1 if all n slices are the same) cube of observation covariances. The j'th column of the i'th slice corresponds to the covariance matrix stored columnwise of the j'th discrete state at time i.
//' @param yt a kxn matrix of obervations
//' @param transProbs a dxd matrix of transition probabilities for the discrete states
//' @param N the maximum particle number
//' 
//' @return List with components "paths", "weights", and "LastStep". 
//' \describe{
//' \item{"paths"}{is an Nxn matrix specifying the paths of each particle,}
//' \item{"weights"}{is a vector giving the final sampling weight of each path,}
//' \item{"LastStep"}{is the timestep computed. "LastStep" will alway equal n unless all of the sampling weights vanish.}
//' }
//' 
//' @examples
//' data(tempos)
//' theta = c(426.69980736, 136.33213703, -11.84256691, -34.82234559, 
//'           439.37886221, 1, 1, 0.84916635, 0.04611644, 0.74119571, 
//'           0.43966082, 0.02116317, 0.24513563, 0.17253254)
//' y = matrix(tempos[,'Richter_1976'], 1)
//' lt = diff(c(tempos$note_onset, 61))
//' pmats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14], 
//'                   c(132,0), c(400,10)) # prior means and variances on X_1
//' beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0,0,0,0,0,0,0), dt, ct, Tt, Zt,
//'             HHt, GGt, y, transMat, 200))
//' 
//' @export 
// [[Rcpp::export]]
List beamSearch(arma::mat a0, arma::mat P0, arma::vec w0,
                arma::cube dt, arma::cube ct, arma::cube Tt, arma::cube Zt,
                arma::cube HHt, arma::cube GGt, arma::mat yt,
                arma::mat transProbs, int N){
  arma::uword n = yt.n_cols;
  arma::uword K = transProbs.n_cols;
  double maxpart = pow(K,n);
  if(maxpart > N) maxpart = N;
  int m = a0.n_rows;
  a0.resize(m, maxpart);
  P0.resize(m*m, maxpart);
  // Determine if any parameter matrices are time varying. Assumes that each has either:
  // (1) 1 column per state or (2) n columns. Does not perform any checks.
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword HHtvar = HHt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  arma::uvec particles(maxpart, arma::fill::zeros);
  arma::uvec filler = arma::regspace<arma::uvec>(0,K-1);
  arma::umat paths(maxpart, n, arma::fill::zeros);
  arma::colvec weights(maxpart,arma::fill::zeros);
  List step = initializeParticles(w0, N, a0, P0, dt.slice(0), ct.slice(0),
                                  Tt.slice(0), Zt.slice(0), 
                                  HHt.slice(0), 
                                  GGt.slice(0), yt.col(0));
  int BP = step["BadPars"];
  if(BP) return List::create(Named("LastStep") = 1);
  arma::vec newW = step["newW"];
  arma::uword CurrentPartNum = newW.n_elem;
  weights.head(CurrentPartNum) = newW;
  arma::mat a1tmp = step["a1"];
  a0.head_cols(CurrentPartNum) = a1tmp;
  arma::mat P1tmp = step["P1"];
  P0.head_cols(CurrentPartNum) = P1tmp;
  arma::uvec newS = step["newstates"];
  particles.head(CurrentPartNum) = newS;
  paths(arma::span(0, CurrentPartNum - 1), 0) = newS;
  arma::uword iter = 1;
  while(iter < n){
    step = dpf(particles.head(CurrentPartNum), weights.head(CurrentPartNum), 
                              maxpart, transProbs, 
                              a0.head_cols(CurrentPartNum), P0.head_cols(CurrentPartNum),
                              dt.slice(iter * dtvar), ct.slice(iter * ctvar), 
                              Tt.slice(iter * Ttvar), Zt.slice(iter * Ztvar), 
                              HHt.slice(iter * HHtvar), 
                              GGt.slice(iter * GGtvar), yt.col(iter));
    int BP = step["BadPars"];
    if(BP) break;
    arma::vec newW = step["newW"];
    CurrentPartNum = newW.n_elem;
    weights.head(CurrentPartNum) = newW;
    if(CurrentPartNum < maxpart) weights.tail(maxpart-CurrentPartNum) *= 0; // fix decreasing CurrentPartNum
    arma::mat a1tmp = step["a1"];
    a0.head_cols(CurrentPartNum) = a1tmp;
    arma::mat P1tmp = step["P1"];
    P0.head_cols(CurrentPartNum) = P1tmp;
    arma::uvec newS = step["newstates"];
    particles.head(CurrentPartNum) = newS;
    arma::uvec old = step["oldstates"];
    arma::umat tmp = paths.rows(old);
    paths.head_rows(CurrentPartNum) = tmp;   //could create issues if CurrentPartNum decreases over time?
    paths(arma::span(0,CurrentPartNum-1), iter) = newS;
    iter++;
  }
  return List::create(Named("paths") = paths,
                      Named("weights") = weights,
                      Named("LastStep") = iter++);
}




//' Estimate continuous states given parameters and discrete hidden states
//' 
//' Performs Kalman filtering and smoothing for the switching state space model.
//' Also works on any (potentially time varying) state space model.
//' 
//' @param pmats a list with componente \code{a0, P0, dt, ct, Tt, Zt, HHt,} and \code{GGt} 
//' e.g., as output from musicModel. See the details below.
//' @param path path of discrete hidden states
//' @param y observations
//' 
//' @return List with components from Kalman filter and Smoother. See details.
//' 
//' @details
//' 
//' \strong{State space form:}
//' 
//' The following notation is closest to the one of Koopman et al. The state
//' space model is represented by the transition equation and the measurement
//' equation. Let m be the dimension of the state variable, d be the dimension of
//' the observations, and n the number of observations. The transition equation
//' and the measurement equation are given by \deqn{a(t + 1) = d(t) + T(t) a(t) +
//' H(t) \eta(t)} \deqn{y(t) = c(t) + Z(t) a(t) + G(t) \epsilon(t),} where
//' \eqn{\eta(t)} and \eqn{\epsilon(t)} are iid \eqn{N(0, I_m)} and iid \eqn{N(0,
//' I_d)}, respectively, and \eqn{\alpha(t)} denotes the state variable. The
//' parameters admit the following dimensions:
//' 
//' \tabular{lll}{ \eqn{a_t \in R^m}{a[t] in R^m} \tab \eqn{d_t \in R^m}{d[t] in
//' R^m} \tab \eqn{\eta_t \in R^m}{eta[t] in R^m} \cr \eqn{T_t \in R^{m \times
//' m}}{d[t] in R^{m \times m}} \tab \eqn{H_t \in R^{m \times m}}{d[t] in R^{m \times m}}
//' \tab \cr \eqn{y_t \in R^d}{y[t] in R^d} \tab \eqn{c_t \in R^d}{c[t] in R^d}
//' \tab \eqn{\epsilon_t \in R^d}{epsilon[t] in R^d}. \cr \eqn{Z_t \in R^{d
//' \times m}}{Z[t] in R^{d \times m}} \tab \eqn{G_t \in R^{d \times d}}{G[t] in R^{d
//' \times d}} \tab}
//' 
//' Note that \code{kalman} takes as input \code{HHt} and \code{GGt} which
//' corresponds to \eqn{H_t H_t'}{H[t] \%*\% t(H[t])} and \eqn{G_t G_t'}{G[t]
//' \%*\% t(G[t])}.
//' 
//' \strong{Iteration:}
//' 
//' Let \code{i} be the loop variable. The filter iterations are implemented the
//' following way (in case of no NA's): Initialization:\cr \code{if(i == 1)\{}\cr
//' \code{  at[, i] = a0}\cr \code{  Pt[,, i] = P0}\cr \code{\}}
//' 
//' Updating equations:\cr \code{vt[, i] = yt[, i] - ct[, i] - Zt[,,i] \%*\% at[,
//' i]}\cr \code{Ft[,, i] = Zt[,, i] \%*\% Pt[,, i] \%*\% t(Zt[,, i]) + GGt[,,
//' i]}\cr \code{Kt[,, i] = Pt[,, i] \%*\% t(Zt[,, i]) \%*\% solve(Ft[,, i])}\cr
//' \code{att[, i] = at[, i] + Kt[,, i] \%*\% vt[, i]}\cr \code{Ptt[, i] = Pt[,,
//' i] - Pt[,, i] \%*\% t(Zt[,, i]) \%*\% t(Kt[,, i])}
//' 
//' Prediction equations:\cr \code{at[, i + 1] = dt[, i] + Tt[,, i] \%*\% att[,
//' i]}\cr \code{Pt[,, i + 1] = Tt[,, i] \%*\% Ptt[,, i] \%*\% t(Tt[,, i]) +
//' HHt[,, i]}
//' 
//' Next iteration:\cr \code{i <- i + 1}\cr goto \dQuote{Updating equations}.
//' 
//' 
//' \strong{Parameters:}
//' 
//' The parameters can either be constant or deterministic time-varying. Assume
//' the number of observations is \eqn{n} (i.e. \eqn{y = (y_t)_{t = 1, \ldots,
//' n}, y_t = (y_{t1}, \ldots, y_{td})}{y = y[,1:n]}). Then, the parameters admit
//' the following classes and dimensions: \tabular{ll}{ \code{dt} \tab either a
//' \eqn{m \times n}{m * n} (time-varying) or a \eqn{m \times 1}{m * 1}
//' (constant) matrix. \cr \code{Tt} \tab either a \eqn{m \times m \times n}{m *
//' m * n} or a \eqn{m \times m \times 1}{m * m * 1} array. \cr \code{HHt} \tab
//' either a \eqn{m \times m \times n}{m * m * n} or a \eqn{m \times m \times
//' 1}{m * m * 1} array. \cr \code{ct} \tab either a \eqn{d \times n}{d * n} or a
//' \eqn{d \times 1}{d * 1} matrix. \cr \code{Zt} \tab either a \eqn{d \times m
//' \times n}{d * m * n} or a \eqn{d \times m \times 1}{d * m * 1} array. \cr
//' \code{GGt} \tab either a \eqn{d \times d \times n}{d * d * n} or a \eqn{d
//' \times d \times 1}{d * d * 1} array. \cr \code{yt} \tab a \eqn{d \times n}{d
//' * n} matrix. }
//' 
//' \strong{Output:}
//' 
//' The output values are the filter estimates, smoother estimates, 
//' associated forecasts, and the negative log likelihood.
//' 
//' \describe{
//' \item {\code{at}}{One step ahead expectation of the hidden state}
//' \item {\code{Pt}{One step ahead variance of the hidden state}
//' \item {\code{preds}}{One step ahead forecast of the observation}
//' \item {\code{ahat}}{Smoothed state mean (conditional on all data)}
//' \item {\code{att}}{State expectation given information up to time t}
//' \item {\code{Ptt}}{State variance given information up to time t}
//' \item {\code{ests}}{Smoothed forecast (conditionall on all data)}
//' \item {\code{llik}}{log likelihood (not negative)}
//' }
//' 
//' @examples
//' data(tempos)
//' theta = c(426.69980736, 136.33213703, -11.84256691, -34.82234559, 
//'           439.37886221, 1, 1, 0.84916635, 0.04611644, 0.74119571, 
//'           0.43966082, 0.02116317, 0.24513563, 0.17253254)
//' y = matrix(tempos[,'Richter_1976'], 1)
//' lt = diff(c(tempos$note_onset, 61))
//' pmats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14], 
//'                   c(132,0), c(400,10)) # prior means and variances on X_1
//' beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0,0,0,0,0,0,0), dt, ct, Tt, Zt,
//'             HHt, GGt, y, transMat, 200))
//' bestpath = with(beam, paths[which.max(weights),])
//' kal = kalman(pmats, bestpath, y)
//' 
//' @export 
// [[Rcpp::export]]
List kalman(List pmats, arma::uvec path, arma::mat y){
  // What if I want different initial state (instead of 0)?
  
  arma::mat a0 = pmats["a0"];
  arma::mat P0 = pmats["P0"];
  arma::cube dt = pmats["dt"];
  arma::cube ct = pmats["ct"];
  arma::cube Tt = pmats["Tt"];
  arma::cube Zt = pmats["Zt"];
  arma::cube HHt = pmats["HHt"];
  arma::cube GGt = pmats["GGt"];
  
  arma::uword n = y.n_cols;
  arma::uword m = a0.n_rows;
  arma::uword mm = m*m;
  arma::uword d = y.n_rows;
  arma::uword dm = d*m;
  arma::uword dd = d*d;
  
  arma::uword dtvar = dt.n_slices > 1;
  arma::uword ctvar = ct.n_slices > 1;
  arma::uword Ttvar = Tt.n_slices > 1;
  arma::uword Ztvar = Zt.n_slices > 1;
  arma::uword HHtvar = HHt.n_slices > 1;
  arma::uword GGtvar = GGt.n_slices > 1;
  
  // output storage
  arma::colvec llik = arma::zeros(n);
  arma::mat at(m,n,arma::fill::zeros); // predicted mean E[a_t | y_1:t-1]
  arma::mat att(m,n,arma::fill::zeros); // updated mean E[a_t | y_1:t]
  arma::cube Pt(m,m,n,arma::fill::zeros); // predicted variance
  arma::cube Ptt(m,m,n,arma::fill::zeros); // predicted variance
  arma::mat preds(d,n,arma::fill::zeros); // predictions  
  arma::mat ahat(m,n,arma::fill::zeros); // smoothed mean E[a_t | y_1:n]
  arma::mat ests(d,n,arma::fill::zeros); //smoothed estimates
  double liktmp = 0.0;
  
  // initialization
  arma::uword s = path(0); // See comment above, replace 0 with appropriate initial s.
  arma::mat a00 =  a0.col(s);
  arma::mat P00 = reshape(P0.col(s), m, m);
  
  
  // Kalman filtering
  arma::uword iter=0;
  while(iter < n){
    s = path(iter);
    KFOUT step = kf1step(a00, P00, 
                         dt.subcube(0,s,iter*dtvar,arma::size(m,1,1)),
                         ct.subcube(0,s,iter*ctvar,arma::size(d,1,1)), 
                         Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1)),
                         Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1)), 
                         HHt.subcube(0,s,iter*HHtvar,arma::size(mm,1,1)), 
                         GGt.subcube(0,s,iter*GGtvar,arma::size(dd,1,1)), 
                         y.col(iter));
    at.col(iter) = step.at;
    att.col(iter) = step.att;
    Pt.slice(iter) = arma::reshape(step.Pt, m, m);
    Ptt.slice(iter) = arma::reshape(step.Ptt, m, m);
    preds.col(iter) = step.pred;
    a00 = step.att;
    P00 = arma::reshape(step.Ptt, m, m);
    liktmp = step.lik;
    llik(iter) += log(liktmp);
    iter++;
  }
  
  // Kalman smoothing. Based on Durbin and Koopman (4.70)
  iter--;
  ahat.col(iter) = att.col(iter);
  while(iter > 0){
    arma::mat cc = ct.subcube(0,s,iter*ctvar,arma::size(d,1,1));
    arma::mat zz = Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1));
    zz.reshape(d,m);
    ests.col(iter) = cc + zz * ahat.col(iter);
    P00 = arma::pinv(Pt.slice(iter));
    a00 = ahat.col(iter) - at.col(iter);
    
    iter--; // This is important, need T after the increment rather than before.
    arma::mat Tmat = Tt.subcube(0,s,iter*Ttvar,arma::size(mm,1,1));
    s = path(iter);
    Tmat.reshape(m,m);
    ahat.col(iter) = att.col(iter) + Ptt.slice(iter) * Tmat.t() * P00 * a00;
  }
  arma::mat cc = ct.subcube(0,s,iter*ctvar,arma::size(d,1,1));
  arma::mat zz = Zt.subcube(0,s,iter*Ztvar,arma::size(dm,1,1));
  zz.reshape(d,m);
  ests.col(iter) = cc + zz * ahat.col(iter);
    
  return List::create(Named("at") = at, Named("Pt") = Pt, Named("preds") = preds,
                            Named("ahat") = ahat,
                            Named("att") = att, Named("Ptt") = Ptt,
                            Named("ests") = ests, Named("llik") = llik);
}


/* ------------------------------------------------- */
/*                                                   */
/* Extended KF code begins below                     */
/*                                                   */
/* ------------------------------------------------- */

struct eKFOUT{
  double att; // a_t|t
  double at;  // a_t+1 (but this is done at the beginning of the KF
  // instead of the end)
  double Ptt; // same deal
  double Pt;
  double lik;
  double pred; // f(at) to estimate yt
};


eKFOUT ekfmm1step(
    arma::vec llt, int state, int t, arma::vec yyt,
    double sig2eps, arma::vec mus, arma::vec sig2eta,
    double a0, double P0){
  
  // Creates parameter matrices and does 1 step ekf
  // Specific to the simplified multiplicative model
  // 
  // States are: (1,1) (1,2) (1,4) (2,2) (2,3) (3,1) (3,3) (4,1), 
  //            (1,3), (2,1), (3,2) (same as with usual model)
  
  double ct;
  double Zt;
  double HHt = 0.0;
  double dt = 0.0;
  double Tt = 1.0;
  double GGt = sig2eps;
  
  double lt = llt(t);
  double yt = yyt(t);
  
  std::set<int> nonlin = {1, 3, 4, 6, 8, 10};
  std::set<int> acc = {4, 6, 8};
  
  if(state == 5 || state == 9){
    dt = mus(0);
    Tt = 0.0;
    HHt = sig2eta(0);
  }
  if(nonlin.find(state) != nonlin.end()){
    dt = lt;
    dt *= acc.find(state) == acc.end() ? log(1+mus(1)) : log(1-mus(1));
    HHt = sig2eta(1);
  }
  
  // KF
  // Update a_t-1|t-1 to a_t using CURRENT state
  double atemp = dt + Tt * a0;
  double Ptemp = HHt + Tt * P0 * Tt;
  // Evaluate Taylor approximation of exp(x)
  Zt = exp(atemp);
  ct = Zt * (1.0 - atemp);
  if(state == 2){
    ct += mus(2);
    GGt += sig2eta(2);
  }
  // Make predictions using a_t
  double pred = ct + Zt * atemp;
  double vt = yt - pred;
  double Ft = GGt + Zt * Ptemp * Zt;
  double Ftinv = 1/Ft;
  double Kt = Ptemp * Zt * Ftinv;
  // Incorporate current information into a_t|t
  double a1 = atemp + Kt * vt;
  double P1 = Ptemp - Ptemp * Zt * Kt;
  
  // Calculate likelihood
  double ftdet = std::abs(Ftinv);
  ftdet = 1 / std::abs(ftdet);
  double mahalanobis = vt * Ftinv * vt;
  mahalanobis += log(2*M_PI) + log(ftdet);
  double lik = exp(-1.0*mahalanobis/2);
  eKFOUT output = {a1, atemp, P1, Ptemp, lik, pred};
  return output;
}







// [[Rcpp::export]]
List edpf(arma::uvec currentStates, arma::colvec w, int N,
          arma::mat transProbs, arma::vec lt, int step,
          arma::vec a0, arma::vec P0,
          double sig2eps, arma::vec mus, arma::vec sig2eta, arma::vec yt){
  int npart = currentStates.size();
  int K = transProbs.n_cols;  //number of possible discrete states
  
  arma::mat a1(K, npart);
  arma::mat P1(K, npart);                   
  arma::mat lik(npart, K);
  arma::vec a11(K);
  arma::vec P11(K);
  for(int part=0; part < npart; part++){
    for(int k=0;  k < K; k++){   //both loops together have effect of for each particle
      eKFOUT ekfout = ekfmm1step(lt, k, step, yt, sig2eps, mus, sig2eta,
                                 a0(part), P0(part));
      double atmp = ekfout.att;
      double Ptmp = ekfout.Ptt;
      a11(k) = atmp;
      P11(k) = Ptmp;
      lik(part,k) = ekfout.lik;
    }
    a1.col(part) = a11;
    P1.col(part) = P11;
  }
  arma::mat testLik = lik;
  lik %= transProbs.rows(currentStates);
  lik.each_col() %= w;
  w = arma::vectorise(lik);
  w = arma::normalise(w,1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec positive = arma::find(w);
  int nkeep = positive.size();
  arma::colvec newW = w(positive);
  arma::vec aout(nkeep);
  arma::vec Pout(nkeep);
  arma::uvec newstates(nkeep);
  arma::uvec oldstates(nkeep);
  int i=0;
  for(int k=0; k<K; k++){
    for(int part=0; part < npart; part++){
      if(w(part + k*npart)>0){
        newstates(i) = k;
        oldstates(i) = part;
        aout(i) = a1(k,part);
        Pout(i) = P1(k,part);
        i++;
      }
    }
  }
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("oldstates") = oldstates,
                      Named("newstates") = newstates,
                      Named("newW") = newW,
                      Named("testLik") = testLik);
}


// [[Rcpp::export]]
double egetloglike(arma::vec lt, arma::uvec path, arma::vec y,
                   double sig2eps, arma::vec mus, arma::vec sig2eta,
                   arma::vec a0, arma::vec P0){
  
  arma::uword n = y.n_elem;
  
  double a0_ = a0(path(0));
  double P0_ = P0(path(0));
  
  double llik = 0;
  double liktmp;
  int s;
  
  for(int iter=0; iter<n; iter++){
    s = path(iter);
    eKFOUT step = ekfmm1step(lt, s, iter, y, sig2eps, mus, sig2eta,
                             a0_, P0_);
    a0 = step.att;
    P0 = step.Ptt;
    liktmp = step.lik;
    llik -= log(liktmp);
  }
  return llik;
}



// [[Rcpp::export]]
List einitializeParticles(arma::vec w0, int N,  arma::vec lt, arma::vec yt,
                          double sig2eps, arma::vec mus, arma::vec sig2eta,
                          arma::vec a0, arma::vec P0){
  arma::uvec particles = arma::find(w0);
  int npart = particles.n_elem;
  if(npart == 0){return List::create(Named("BadPars") = 1);}
  arma::vec weights = w0(particles);
  arma::vec a11(npart);
  arma::vec P11(npart);
  arma::vec lik(npart);
  for(int i = 0; i < npart; i++){
    int p = particles(i);
    eKFOUT kfout = ekfmm1step(lt, p, 0, yt, sig2eps, mus, sig2eta,
                              a0(p), P0(p));
    a11(i) = kfout.att;
    P11(i) = kfout.Ptt;
    lik(i) = kfout.lik;
  }
  lik %= weights;
  arma::vec w = arma::normalise(lik, 1);
  w = resampleSubOptimal(w, N);
  if( !arma::any(w)) return List::create(Named("BadPars") = 1);
  arma::uvec nz = arma::find(w);
  arma::vec newW = w(nz);
  arma::uvec newstates = particles(nz);
  arma::mat aout = a11(nz);
  arma::mat Pout = P11(nz);
  return List::create(Named("BadPars") = 0,
                      Named("a1") = aout,
                      Named("P1") = Pout,
                      Named("newstates") = newstates,
                      Named("newW") = newW);
}


// [[Rcpp::export]]
List ebeamSearch(arma::vec lt, arma::vec w0,
                 double sig2eps, arma::vec mus, arma::vec sig2eta,
                 arma::vec a0, arma::vec P0, arma::vec yt,
                 arma::mat transProbs, int N){
  
  arma::uword n = yt.n_elem;
  
  arma::uword K = transProbs.n_cols;
  double maxpart = pow(K,n);
  if(maxpart > N) maxpart = N;
  a0.resize(maxpart);
  P0.resize(maxpart);
  
  arma::uvec particles(maxpart, arma::fill::zeros);
  arma::uvec filler = arma::regspace<arma::uvec>(0,K-1);
  arma::umat paths(maxpart, n, arma::fill::zeros);
  arma::colvec weights(maxpart,arma::fill::zeros);
  List step = einitializeParticles(w0, N, lt, yt, sig2eps, mus, sig2eta, a0, P0);
  
  
  int BP = step["BadPars"];
  if(BP) return List::create(Named("LastStep") = 1);
  arma::vec newW = step["newW"];
  arma::uword CurrentPartNum = newW.n_elem;
  weights.head(CurrentPartNum) = newW;
  arma::vec a1tmp = step["a1"];
  a0.head(CurrentPartNum) = a1tmp;
  arma::vec P1tmp = step["P1"];
  P0.head(CurrentPartNum) = P1tmp;
  arma::uvec newS = step["newstates"];
  particles.head(CurrentPartNum) = newS;
  paths(arma::span(0, CurrentPartNum - 1), 0) = newS;
  arma::uword iter = 1;
  
  while(iter < n){
    step = edpf(particles.head(CurrentPartNum), weights.head(CurrentPartNum), 
                maxpart, transProbs, lt, iter,
                a0.head(CurrentPartNum), P0.head(CurrentPartNum),
                sig2eps, mus, sig2eta, yt);
    int BP = step["BadPars"];
    if(BP) break;
    arma::vec newW = step["newW"];
    CurrentPartNum = newW.n_elem;
    weights.head(CurrentPartNum) = newW;
    if(CurrentPartNum < maxpart) weights.tail(maxpart-CurrentPartNum) *= 0; // fix decreasing CurrentPartNum
    arma::vec a1tmp = step["a1"];
    a0.head(CurrentPartNum) = a1tmp;
    arma::vec P1tmp = step["P1"];
    P0.head(CurrentPartNum) = P1tmp;
    arma::uvec newS = step["newstates"];
    particles.head(CurrentPartNum) = newS;
    arma::uvec old = step["oldstates"];
    arma::umat tmp = paths.rows(old);
    paths.head_rows(CurrentPartNum) = tmp;   //could create issues if CurrentPartNum decreases over time?
    paths(arma::span(0,CurrentPartNum-1), iter) = newS;
    iter++;
  }
  return List::create(Named("paths") = paths,
                      Named("weights") = weights,
                      Named("LastStep") = iter++);
}

// [[Rcpp::export]]
arma::mat ecreateTransMat(arma::vec transprobs){
  
  int nstates = 11;
  arma::mat transMat(nstates, nstates, arma::fill::zeros);
  if(transprobs.n_elem<6) transprobs.resize(6);
  transMat(0,0) = transprobs(0); // (1,1) -> (1,1)
  transMat(0,1) = transprobs(1); // (1,1) -> (1,2)
  transMat(0,2) = 1-transprobs(0)-transprobs(1)-transprobs(4); // (1,1) -> (1,4)
  transMat(0,8) = transprobs(4); // (1,1) -> (1,3)
  transMat(1,3) = 1; // (1,2) -> (2,2)
  transMat(2,7) = 1; // (1,4) -> (4,1)
  transMat(3,3) = transprobs(2); // (2,2) -> (2,2)
  transMat(3,4) = 1 - transprobs(2) - transprobs(5); // (2,2) -> (2,3)
  transMat(3,9) = transprobs(5); // (2,2) -> (2,1)
  transMat(4,6) = 1; // (2,3) -> (3,3)
  transMat(5,0) = 1; // (3,1) -> (1,1)
  transMat(6,5) = transprobs(3); // (3,3) -> (3,1)
  transMat(6,10) = transprobs(6); // (3,3) -> (3,2)
  transMat(6,6) = 1 - transprobs(3) - transprobs(6); // (3,3) -> (3,3)
  transMat(7,0) = 1; // (4,1) -> (1,1)
  transMat(8,6) = 1; // (1,3) -> (3,3)
  transMat(9,0) = 1; // (2,1) -> (1,1)
  transMat(10,3) = 1; // (3,2) -> (2,2)
  
  return transMat;
}





// [[Rcpp::export]]
List ekalman(arma::vec lt, arma::uvec path, arma::vec y,
             double sig2eps, arma::vec mus, arma::vec sig2eta,
             double a0, double P0){
  // What if I want different initial state (instead of 0)?
  
  
  arma::uword n = y.n_elem;
  
  // output storage
  arma::vec llik = arma::zeros(n);
  arma::vec at(n, arma::fill::zeros); // predicted mean E[a_t | y_1:t-1]
  arma::vec att(n, arma::fill::zeros); // updated mean E[a_t | y_1:t]
  arma::vec Pt(n,arma::fill::zeros); // predicted variance
  arma::vec Ptt(n,arma::fill::zeros); // predicted variance
  arma::vec preds(n,arma::fill::zeros); // predictions  
  arma::vec ahat(n,arma::fill::zeros); // smoothed mean E[a_t | y_1:n]
  arma::vec ests(n,arma::fill::zeros); //smoothed estimates
  double liktmp = 0.0;
  
  
  // Kalman filtering
  arma::uword iter=0;
  double a0_ = a0;
  double P0_ = P0;
  int s;
  while(iter < n){
    s = path(iter);
    eKFOUT step = ekfmm1step(lt, s, iter, y, sig2eps, mus, sig2eta, a0_, P0_);
    
    at(iter) = step.at;
    att(iter) = step.att;
    Pt(iter) = step.Pt;
    Ptt(iter) = step.Ptt;
    preds(iter) = step.pred;
    a0_ = step.att;
    P0_ = step.Ptt;
    liktmp = step.lik;
    llik(iter) += log(liktmp);
    iter++;
  }
  
  // Kalman smoothing. Based on Durbin and Koopman (4.70)
  iter--;
  ahat(iter) = att(iter);
  double Tt;
  double P00;
  double a00;
  
  while(iter > 0){
    s = path(iter);
    
    ests(iter) = exp(ahat(iter));
    if(s == 2) ests(iter) += mus(2);
    P00 = 1/Pt(iter);
    a00 = ahat(iter) - at(iter);
    
    iter--; // This is important, need T after the increment rather than before.
    s = path(iter);
    Tt = (s==5 || s==9) ? 0.0 : 1.0;
    ahat(iter) = att(iter) + Ptt(iter) * Tt * P00 * a00;
  }
  
  ests(0) = exp(ahat(0));
  if(s==2) ests(0) += mus(2);
  
  return List::create(Named("at") = at, Named("Pt") = Pt, Named("preds") = preds,
                            Named("ahat") = ahat,
                            Named("att") = att, Named("Ptt") = Ptt,
                            Named("ests") = ests, Named("llik") = llik);
}



