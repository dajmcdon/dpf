#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


arma::uvec SampleNoReplace(arma::uvec x, int size) {
  int nOrig = x.size();
  arma::uvec index(size);
  arma::uvec sub(nOrig);
  for (int ii = 0; ii < nOrig; ii++) sub(ii) = ii;
  RNGScope scope4;
  for (int ii = 0; ii < size; ii++) {
    int jj = floor(nOrig * runif(1)[0]);
    index(ii) = sub(jj);
    // replace sampled element with last, decrement
    sub(jj) = sub(--nOrig);
  }
  arma::uvec ret(size);
  ret = x(index);
  return(ret);
}

//' @export 
// [[Rcpp::export]]
arma::vec resampleSubOptimal(arma::vec w, int N){
  //Rcout << "The value of w is:" << std::endl << w << std::endl;
  int M = w.size();
  double tol = 1e-10;
  arma::vec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  arma::vec nzz = arma::nonzeros(ws);
  arma::uword nz = nzz.n_elem;
  if(M <= N || nz <= N){
    return ws;
  }
  //Rcout << "The value of N is:" << std::endl << N << std::endl;
  //Rcout << "The value of M is:" << std::endl << M << std::endl;
  //Rcout << "The value of tol is:" << std::endl << tol << std::endl;
  //Rcout << "The value of ws is:" << std::endl << ws << std::endl;
  //Rcout << "The value of nzz is:" << std::endl << nzz << std::endl;
  //Rcout << "The value of nz is:" << std::endl << nz << std::endl;
  
  arma::vec w1 = ws;
  typedef std::vector<double> stdvec;
  stdvec z = arma::conv_to<stdvec>::from(w1);
  std::nth_element(z.begin(), z.end()-N, z.end());
  arma::vec z1 = arma::conv_to<arma::vec>::from(z);
  double minprob = z1(M-N);                              //This is the Nth largest element of z1.
  
  //Rcout << "The value of minprob is:" << std::endl << minprob << std::endl;
  
  ws.elem(arma::find(ws < minprob)).zeros();
  arma::uvec ties = arma::find(ws==minprob);
  arma::uvec keep = arma::find(ws > minprob);
  
  //Rcout << "The value of ties is:" << std::endl << ties << std::endl;
  //Rcout << "The value of keep is:" << std::endl << keep << std::endl;
  
  int nkeep = keep.size();
  int tt = ties.size();
  
  //Rcout << "The value of nkeep is:" << std::endl << nkeep << std::endl;
  //Rcout << "The value of tt is:" << std::endl << tt << std::endl;
  
  if(tt > 1){
    arma::vec probs(tt); //This simply declares a vector probs of size tt
    double p = 1.0 / tt;
    probs.fill(p);  //This fills the vector probs with p  //MICHAEL: probs is never used?
    int dontneed = tt - (N - nkeep);
    
    //Rcout << "The value of probs is:" << std::endl << probs << std::endl;
    //Rcout << "The value of dontneed is:" << std::endl << dontneed << std::endl;
    
    arma::uvec idx = SampleNoReplace(ties, dontneed);
    //Rcout << "The value of idx is:" << std::endl << idx << std::endl;
    ws.elem(idx).zeros();
  }
  ws = arma::normalise(ws,1);
  //Rcout << "The value of ws is:" << std::endl << ws << std::endl;
  
  
  return ws;
}

//' @export 
// [[Rcpp::export]]
arma::colvec resampleOptimal(arma::colvec w, int N){

  //Rcout << "The value of w is:" << std::endl << w << std::endl;
  
  int M = w.size();
  double tol = 1e-10;
  arma::colvec ws = w;
  ws.elem( arma::find( ws < tol) ).zeros();
  arma::vec nzz = arma::nonzeros(ws);
  
  //Rcout << "The value of ws is:" << std::endl << ws << std::endl;
  
  arma::uword nz = nzz.n_elem;
  if(M <= N || nz <= N){ 
    return ws;
  }
  
  // Hard case.
  arma::colvec wssortnorm = ws;
  wssortnorm = arma::sort(wssortnorm);
  wssortnorm = arma::normalise(wssortnorm,1);
  //Rcout << "The value of wssortnorm is:" << std::endl << wssortnorm << std::endl;
  int Ak = 0;
  double Bk = 1.0;
  int i = M;
  
  while(i > M - N){
    i--;
    if(wssortnorm(i) < tol || Bk <= tol){
      ws.elem( arma::find( ws < tol) ) .zeros();
      return ws;
    }
    if(Bk/wssortnorm(i) + Ak >= N) break;
    Ak++;
    Bk -= wssortnorm(i);
  }
  double cinv = Bk / (N-Ak);
  
  //Rcout << "The value of cinv is:" << std::endl << cinv << std::endl;
  
  // Set 1
  arma::vec NewW;
  NewW.zeros(M);
  arma::colvec wnormalized = ws;
  wnormalized = arma::normalise(ws,1);
  
  //Rcout << "The value of wnormalized is:" << std::endl << wnormalized << std::endl;
  
  int L=0;
  double K=0;
  for(i=0; i<M; i++){
    if(cinv < wnormalized(i)){
      //Rcout << "The value of i is:" << std::endl << i << std::endl;
      NewW(i) += ws(i);
      L++;
      //Rcout << "The value of L is:" << std::endl << L << std::endl;
    }else{
      K += ws(i);  
    }
  }
  K /= (N-L);
  //Rcout << "The value of N is:" << std::endl << N << std::endl;
  //Rcout << "The value of L is:" << std::endl << L << std::endl;
  //Rcout << "The value of K is:" << std::endl << K << std::endl;
  
  // Set 2
  RNGScope scope;
  double U1 = runif(1)[0] * K;
  //Rcout << "The value of U1 is:" << std::endl << U1 << std::endl;
  for(int i=0; i<M; i++){
    if(NewW(i)==0){
      U1 -= ws(i);
      if(U1 < 0){
        NewW(i) += cinv;
        U1 += K;
      }
    }
  }
  NewW = arma::normalise(NewW,1);
  //Rcout << "The value of NewW is:" << std::endl << NewW << std::endl;
  return NewW;
}

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


List kf1stepR(arma::mat a0, arma::mat P0, arma::mat dt,
              arma::mat ct, arma::mat Tt,
              arma::mat Zt, arma::mat HHt, arma::mat GGt, arma::mat yt) {
    KFOUT output = kf1step(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt);
    
    return List::create(Named("att") = output.att,
                        Named("at") = output.at,
                        Named("Ptt") = output.Ptt,
                        Named("Pt") = output.Pt,
                        Named("lik") = output.lik,
                        Named("pred") = output.pred);
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
  
  //Rcout << "The value of lik is:" << std::endl << lik << std::endl;
  
  w = arma::vectorise(lik);
  w = arma::normalise(w,1);
  
  //Rcout << "The value of w is:" << std::endl << w << std::endl;
  
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
  //   4 transition matrix parameters
  // up to 22 parameters per performance
  // for now, everything constant, except mu_tempo
  // States (in order) are: (1,1) (1,2) (1,4) (2,2) (2,3) (3,1) (3,3) (4,1) (1,3) (2,1)  (3,2)  #Rob
  int nstates = 11;
  int d = 1;
  int m = 2;
  int mm = m*m;
  int n = lt.n_elem; //n is length of onset vector (231)
  
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
  
  arma::cube ct(d, nstates, 1, arma::fill::zeros);
  
  arma::cube Tt(mm, nstates, n, arma::fill::zeros);
  Tt.tube(0,0,0,4) += 1; // .tube(first row, first col, law row, last col)
  Tt.tube(0,6,0,8) += 1;
  Tt.tube(3,3) += 1;
  Tt.tube(3,6) += 1;
  Tt.tube(2,6) = lt;
  Tt.tube(2,3) = lt;
  
  //Rcout << "The value of Tt is:" << std::endl << Tt << std::endl;
  
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


//' Parameter matrices for our music model on dynamics
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
List musicModeldynamics(arma::vec lt, double mueps, double sig2eps, arma::vec mus,
                arma::vec sig2eta, arma::vec transprobs,
                arma::vec initialMean, arma::vec initialVariance){ 
  // States (si,si-1) are: (1,1) (1,2) (2,1) (2,2)
  int nstates = 4; 
  int d = 1; //length of observed state
  int m = 3; //length of hidden state
  int mm = m*m;
  int n = lt.n_elem;
  
  arma::mat a0(m, nstates, arma::fill::zeros);
  a0.row(0) += initialMean(0);
  a0.row(1) += initialMean(1);
  a0.row(2) += initialMean(2);
  
  arma::mat P0(mm, nstates, arma::fill::zeros);
  P0(0,0) += initialVariance(0);
  P0(4,0) += initialVariance(1);
  P0(8,0) += initialVariance(2);
  P0(0,1) += initialVariance(0);
  P0(4,1) += initialVariance(1);
  P0(8,1) += initialVariance(2);
  P0(0,2) += initialVariance(0);
  P0(4,2) += initialVariance(1);
  P0(8,2) += initialVariance(2);
  P0(0,3) += initialVariance(0);
  P0(4,3) += initialVariance(1);
  P0(8,3) += initialVariance(2);
  
  arma::cube dt(m, nstates, n, arma::fill::zeros);
  dt.tube(0,0) += mus(0);
  dt.tube(0,2) += mus(0);
  dt.tube(1,2) += mus(1);
  dt.tube(2,2) += mus(2);
  
  arma::cube ct(d, nstates, 1, arma::fill::zeros);
  ct.tube(0,0) += mueps;
  ct.tube(0,1) += mueps;
  ct.tube(0,2) += mueps;
  ct.tube(0,3) += mueps;
  
  arma::cube Tt(mm, nstates, n, arma::fill::zeros);
  Tt.tube(0,1) += 1;
  Tt.tube(3,1) += 1;
  Tt.tube(4,1) += 1;
  Tt.tube(7,1) += 1;
  Tt.tube(8,1) += 1;
  Tt.tube(0,3) += 1;
  Tt.tube(3,3) += 1;
  Tt.tube(4,3) += 1;
  Tt.tube(7,3) += 1;
  Tt.tube(8,3) += 1;
  
  arma::cube Zt(m, nstates, 1, arma::fill::ones);
  
  arma::cube HHt(mm, nstates, n, arma::fill::zeros);
  HHt.tube(0,0) += sig2eta(0);
  HHt.tube(0,2) += sig2eta(0);
  HHt.tube(4,2) += sig2eta(1);
  HHt.tube(8,2) += sig2eta(2);
  
  arma::cube GGt(d, nstates, 1, arma::fill::ones);
  GGt *= sig2eps;
  
  arma::mat transMat(nstates, nstates, arma::fill::zeros);
  transMat(0,0) = transprobs(0); // (1,1) -> (1,1)
  transMat(0,1) = 1 - transprobs(0); // (1,1) -> (1,2)
  transMat(1,2) = transprobs(1); // (1,2) -> (2,1)
  transMat(1,3) = 1 - transprobs(1); // (1,2) -> (2,2)
  transMat(2,0) = transprobs(2); // (2,1) -> (1,1)
  transMat(2,1) = 1 - transprobs(2); // (2,1) -> (1,2)
  transMat(3,2) = transprobs(3); // (2,2) -> (1,1)
  transMat(3,3) = 1 - transprobs(3); // (2,2) -> (2,2)
  
  return List::create(Named("a0") = a0, Named("P0") = P0,
                      Named("dt") = dt, Named("ct") = ct,
                      Named("Tt") = Tt, Named("Zt") = Zt,
                      Named("HHt") = HHt,
                      // Named("Rt") = Rt, Named("Qt") = Qt,
                      Named("GGt") = GGt,
                      Named("transMat") = transMat);
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
    ahat.col(iter) = att.col(iter) + Ptt.slice(iter) * Tmat * P00 * a00;
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
