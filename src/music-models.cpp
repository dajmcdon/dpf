#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

//' Parameter matrices for our music model
//' 
//' This function accepts a number of parameters and creates a list of matrices
//' for Kalman filter evaluation. See the paper for the particular form of the model.
//' 
//' @param lt vector of durations between successive notes in the score
//' @param sig2eps variance of the observation noise
//' @param mus vector of 3 mean parameters (\eqn{\mu, \tau, and \varphi})
//' @param sig2eta vector of 3 state variance parameters (\eqn{\sigma_3^2, \sigma_2^2,and \sigma_4^2})
//' @param transprobs vector of 7 transition probabilities
//' @param initialMean a vector of length 2 giving the prior tempo and the prior acceleration for when state 1 or 3 is the starting state
//' @param initialVariance a vector of length 2 giving the prior variance for the tempo and the prior variance for the acceleration for when state 1 or 3 is the starting state 
//' 
//' @return List with components as appropriate for Kalman filtering or Beam Search. These include: \describe{
//' \item{a0}{a pxd matrix of the initial means of the hidden state. The j'th column corresponds to the initial mean when starting in the j'th discrete state.}
//' \item{P0}{a (p^2)xd matrix of the initial covariances of the hidden state. The j'th column corresponds to the initial covariances stored columnwise when starting in the j'th discrete state.}
//' \item{dt}{a pxdxn cube of state intercepts. The j'th column of the i'th slice corresponds to the intercept specified by the j'th discrete state at time i.}
//' \item{ct}{a kxdx1 cube of observation intercepts. The j'th column corresponds to the intercept specified by the j'th discrete state.}
//' \item{Tt}{a (p^2)xdxn cube of state slopes. The j'th column of the i'th slice corresponds to the slope matrix stored columnwise of the j'th discrete state at time i.}
//' \item{Zt}{a pkxdx1 cube of obvervation slopes. The j'th column corresponds to the slope matrix stored columnwise of the j'th discrete state.}
//' \item{HHt}{a (p^2)xdxn cube of state covariances. The j'th column of the i'th slice corresponds to the covariance matrix stored columnwise of the j'th discrete state at time i.}
//' \item{GGt}{a (k^2)xdx1 cube of observation covariances. The j'th column corresponds to the covariance matrix stored columnwise of the j'th discrete state.}
//' \item{transMat}{a dxd matrix of transition probabilities for the discrete states}
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
//'                   c(132,0), c(400,10))
//'                   
//' @export    
// [[Rcpp::export]]
List musicModel(arma::vec lt, double sig2eps, arma::vec mus,
                arma::vec sig2eta, arma::vec transprobs,
                arma::vec initialMean, arma::vec initialVariance){ 
  //confirm that t's stay in same order for each matrix
  // in each section, we have:
  //   3 state means (tempo, accel, stress), 3 state variances (same), 1 obs variance
  //   7 transition matrix parameters
  // up to 22 parameters per performance
  // for now, everything constant, except mu_tempo
  // States are: (1,1) (1,2) (1,4) (2,2) (2,3) (3,1) (3,3) (4,1) (1,3) (2,1), (3,2)
  int nstates = 11;
  int d = 1;
  int m = 2;
  int mm = m*m;
  int n = lt.n_elem;
  
  arma::mat a0(m, nstates, arma::fill::zeros);
  a0.row(0) += initialMean(0);
  a0(1,4) += initialMean(1);
  a0(1,6) += initialMean(1);
  
  arma::mat P0(mm, nstates, arma::fill::zeros);
  P0.row(0) += initialVariance(0);
  P0(3,4) += initialVariance(1);
  P0(3,6) += initialVariance(1);
  
  arma::cube dt(m, nstates, n, arma::fill::zeros);
  dt.tube(0,1) = lt*mus(1);
  dt.tube(0,4) = -lt*mus(1);
  dt.tube(0,8) = -lt*mus(1);
  dt.tube(0,9) = lt*mus(1);
  dt.tube(1,1) += mus(1);
  dt.tube(1,4) -= mus(1);
  dt.tube(1,8) -= mus(1);
  dt.tube(1,2) += mus(2);
  dt.tube(0,5) += mus(0);
  dt.tube(0,9) += mus(0);
  dt.tube(0,10) = lt*mus(1);
  dt.tube(1,10) += mus(1);
  
  arma::cube ct(d, nstates, 1, arma::fill::zeros);
  
  arma::cube Tt(mm, nstates, n, arma::fill::zeros);
  Tt.tube(0,0,0,4) += 1;
  Tt.tube(0,6,0,8) += 1;
  Tt.tube(3,3) += 1;
  Tt.tube(3,6) += 1;
  Tt.tube(2,6) = lt;
  Tt.tube(2,3) = lt;
  Tt.tube(0,10) += 1;
  
  arma::cube Zt(m, nstates, 1, arma::fill::zeros);
  Zt.tube(0,0,0,9) += 1;
  Zt(1,2,0) += 1;
  
  arma::cube HHt(mm, nstates, n, arma::fill::zeros);
  //(1,2)
  HHt.tube(0,1) = sig2eta(1)*lt%lt;
  for(int iter=1; iter<3; iter++) HHt.tube(iter,1) = sig2eta(1)*lt;
  HHt.tube(3,1) += sig2eta(1);
  //(1,4)
  HHt.tube(3,2) += sig2eta(2);
  //(2,3)
  HHt.tube(0,4) = sig2eta(1)*lt%lt;
  for(int iter=1; iter<3; iter++) HHt.tube(iter,4) = sig2eta(1)*lt;
  HHt.tube(3,4) += sig2eta(1);
  //(3,1)
  HHt.tube(0,5) += sig2eta(0);
  //(1,3)
  HHt.tube(0,8) = sig2eta(1)*lt%lt;
  for(int iter=1; iter<3; iter++) HHt.tube(iter,8) = sig2eta(1)*lt;
  HHt.tube(3,8) += sig2eta(1);
  //(2,1)
  HHt.tube(0,9) += sig2eta(0);
  //(3,2)
  HHt.tube(0,10) = sig2eta(1)*lt%lt;
  for(int iter=1; iter<3; iter++) HHt.tube(iter,10) = sig2eta(1)*lt;
  HHt.tube(3,10) += sig2eta(1);
  // Extra noise
  if(sig2eta.n_elem>3) HHt.tube(0,0,0,9) += sig2eta(3); // set to zero, generally
  
  arma::cube GGt(d, nstates, 1, arma::fill::ones);
  GGt *= sig2eps;
  
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
  
  return List::create(Named("a0") = a0, Named("P0") = P0,
                      Named("dt") = dt, Named("ct") = ct,
                      Named("Tt") = Tt, Named("Zt") = Zt,
                      Named("HHt") = HHt,
                      // Named("Rt") = Rt, Named("Qt") = Qt,
                      Named("GGt") = GGt,
                      Named("transMat") = transMat);
}


// [[Rcpp::export]]
List yupengMats(arma::vec lt, double sig2eps, arma::vec mus,
                arma::vec sig2eta, arma::vec transprobs,
                arma::vec initialMean, arma::vec initialVariance){ 
  //   3 state means (tempo, accel, stress), 3 state variances (same), 1 obs variance
  //   4 transition matrix parameters
  int nstates = 8;
  int d = 1;
  int m = 2;
  int mm = m*m;
  int n = lt.n_elem;
  arma::mat a0(m, nstates, arma::fill::zeros);
  a0.row(0) += initialMean(0);
  a0(1,4) += initialMean(1);
  a0(1,6) += initialMean(1);
  
  arma::mat P0(mm, nstates, arma::fill::zeros);
  P0.row(0) += initialVariance(0);
  P0(3,4) += initialVariance(1);
  P0(3,6) += initialVariance(1);
  
  arma::cube dt(m, nstates, n, arma::fill::zeros);
  dt.tube(0,1) = lt*mus(1);
  dt.tube(0,4) = lt*mus(1);
  dt.tube(1,1) += mus(1);
  dt.tube(1,4) += mus(1);
  dt.tube(1,2) += mus(2);
  dt.tube(0,5) += mus(0);
  // Rcout << "dt done" << std::endl;
  arma::cube ct(d, nstates, 1, arma::fill::zeros);
  // Rcout << "ct done" << std::endl;
  arma::cube Tt(mm, nstates, n, arma::fill::zeros);
  Tt.tube(0,0,0,4) += 1;
  Tt.tube(0,6,0,7) += 1;
  Tt.tube(3,3) += 1;
  Tt.tube(3,6) += 1;
  Tt.tube(2,4) = lt;
  Tt.tube(2,6) = lt;
  Tt.tube(2,3) = lt;
  // Rcout << "Tt done" << std::endl;
  arma::cube Zt(m, nstates, 1, arma::fill::zeros);
  Zt.tube(0,0,0,7) += 1;
  Zt(1,2,0) += 1;
  // Rcout << "Zt done" << std::endl;
  arma::cube Rt(mm, nstates, n, arma::fill::zeros);
  Rt.tube(0,5) += 1;
  Rt.tube(3,1,3,2) += 1;
  Rt.tube(3,4) += 1;
  Rt.tube(2,1) = lt;
  Rt.tube(2,4) = lt;
  // Rcout << "Rt done" << std::endl;
  arma::mat Qt(mm, nstates, arma::fill::zeros);
  Qt.submat(0,0,0,7) += sig2eta(0);
  Qt.submat(3,0,3,7) += sig2eta(1);
  Qt(3,2) -= sig2eta(1) - sig2eta(2);
  // Rcout << "Qt done" << std::endl;
  arma::cube HHt(mm, nstates, n, arma::fill::zeros);
  arma::mat R(mm,1);
  arma::mat Q(mm,1);
  arma::mat HH(m,m);
  for (int i=0; i < n; i++){
    for (int j=0; j<nstates; j++){
      R = Rt.subcube(0,j,i, arma::size(mm,1,1));
      Q = Qt.col(j);
      R.reshape(m,m);
      Q.reshape(m,m);
      HH = R * Q * R.t();
      HH.reshape(mm,1);
      HHt.subcube(0, j, i, arma::size(mm,1,1)) += HH;
    }
  }
  arma::cube GGt(d, nstates, 1, arma::fill::ones);
  GGt *= sig2eps;
  // Rcout << "GGt done" << std::endl;
  arma::mat transMat(nstates, nstates, arma::fill::zeros);
  transMat(0,0) = transprobs(0); // (1,1) -> (1,1)
  transMat(0,1) = transprobs(1); // (1,1) -> (1,2)
  transMat(0,2) = 1-transprobs(0)-transprobs(1); // (1,1) -> (1,4)
  transMat(1,3) = transprobs(2); // (1,2) -> (2,2)
  transMat(1,4) = 1-transprobs(2); // (1,2) -> (2,3)
  transMat(2,7) = 1; // (1,4) -> (4,1)
  transMat(3,3) = transprobs(2); // (2,2) -> (2,2)
  transMat(3,4) = 1-transprobs(2); // (2,2) -> (2,3)
  transMat(4,5) = transprobs(3); // (2,3) -> (3,1)
  transMat(4,6) = 1-transprobs(3); // (2,3) -> (3,3)
  transMat(5,0) = transprobs(0); // (3,1) -> (1,1)
  transMat(5,1) = transprobs(1); // (3,1) -> (1,2)
  transMat(5,2) = 1-transprobs(0)-transprobs(1); // (3,1) -> (1,4)
  transMat(6,5) = transprobs(3); // (3,3) -> (3,1)
  transMat(6,6) = 1-transprobs(3); // (3,3) -> (3,3)
  transMat(7,0) = transprobs(0); // (4,1) -> (1,1)
  transMat(7,1) = transprobs(1); // (4,1) -> (1,2)
  transMat(7,2) = 1-transprobs(0)-transprobs(1); // (4,1) -> (1,4)
  // Rcout << "transMat done" << std::endl;
  return List::create(Named("a0") = a0, Named("P0") = P0,
                      Named("dt") = dt, Named("ct") = ct,
                      Named("Tt") = Tt, Named("Zt") = Zt,
                      Named("HHt") = HHt,
                      Named("GGt") = GGt,
                      Named("transMat") = transMat);
}
