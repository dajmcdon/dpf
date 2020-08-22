// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// resampleSubOptimal
arma::vec resampleSubOptimal(arma::vec w, int N);
RcppExport SEXP _dpf_resampleSubOptimal(SEXP wSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(resampleSubOptimal(w, N));
    return rcpp_result_gen;
END_RCPP
}
// resampleOptimal
List resampleOptimal(arma::colvec w, int N);
RcppExport SEXP _dpf_resampleOptimal(SEXP wSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(resampleOptimal(w, N));
    return rcpp_result_gen;
END_RCPP
}
// getloglike
double getloglike(List pmats, arma::uvec path, arma::mat y);
RcppExport SEXP _dpf_getloglike(SEXP pmatsSEXP, SEXP pathSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pmats(pmatsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type path(pathSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(getloglike(pmats, path, y));
    return rcpp_result_gen;
END_RCPP
}
// musicModel
List musicModel(arma::vec lt, double sig2eps, arma::vec mus, arma::vec sig2eta, arma::vec transprobs, arma::vec initialMean, arma::vec initialVariance);
RcppExport SEXP _dpf_musicModel(SEXP ltSEXP, SEXP sig2epsSEXP, SEXP musSEXP, SEXP sig2etaSEXP, SEXP transprobsSEXP, SEXP initialMeanSEXP, SEXP initialVarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lt(ltSEXP);
    Rcpp::traits::input_parameter< double >::type sig2eps(sig2epsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig2eta(sig2etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type transprobs(transprobsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialMean(initialMeanSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialVariance(initialVarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(musicModel(lt, sig2eps, mus, sig2eta, transprobs, initialMean, initialVariance));
    return rcpp_result_gen;
END_RCPP
}
// musicModeldynamics
List musicModeldynamics(arma::vec lt, double mueps, double sig2eps, arma::vec mus, arma::vec sig2eta, arma::vec transprobs, double muerror, arma::vec initialMean, arma::vec initialVariance);
RcppExport SEXP _dpf_musicModeldynamics(SEXP ltSEXP, SEXP muepsSEXP, SEXP sig2epsSEXP, SEXP musSEXP, SEXP sig2etaSEXP, SEXP transprobsSEXP, SEXP muerrorSEXP, SEXP initialMeanSEXP, SEXP initialVarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lt(ltSEXP);
    Rcpp::traits::input_parameter< double >::type mueps(muepsSEXP);
    Rcpp::traits::input_parameter< double >::type sig2eps(sig2epsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sig2eta(sig2etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type transprobs(transprobsSEXP);
    Rcpp::traits::input_parameter< double >::type muerror(muerrorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialMean(initialMeanSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialVariance(initialVarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(musicModeldynamics(lt, mueps, sig2eps, mus, sig2eta, transprobs, muerror, initialMean, initialVariance));
    return rcpp_result_gen;
END_RCPP
}
// beamSearch
List beamSearch(arma::mat a0, arma::mat P0, arma::vec w0, arma::cube dt, arma::cube ct, arma::cube Tt, arma::cube Zt, arma::cube HHt, arma::cube GGt, arma::mat yt, arma::mat transProbs, int N, int samplemethod);
RcppExport SEXP _dpf_beamSearch(SEXP a0SEXP, SEXP P0SEXP, SEXP w0SEXP, SEXP dtSEXP, SEXP ctSEXP, SEXP TtSEXP, SEXP ZtSEXP, SEXP HHtSEXP, SEXP GGtSEXP, SEXP ytSEXP, SEXP transProbsSEXP, SEXP NSEXP, SEXP samplemethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< arma::cube >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ct(ctSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Tt(TtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type HHt(HHtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type GGt(GGtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type yt(ytSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type transProbs(transProbsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type samplemethod(samplemethodSEXP);
    rcpp_result_gen = Rcpp::wrap(beamSearch(a0, P0, w0, dt, ct, Tt, Zt, HHt, GGt, yt, transProbs, N, samplemethod));
    return rcpp_result_gen;
END_RCPP
}
// kalman
List kalman(List pmats, arma::uvec path, arma::mat y);
RcppExport SEXP _dpf_kalman(SEXP pmatsSEXP, SEXP pathSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pmats(pmatsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type path(pathSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(kalman(pmats, path, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dpf_resampleSubOptimal", (DL_FUNC) &_dpf_resampleSubOptimal, 2},
    {"_dpf_resampleOptimal", (DL_FUNC) &_dpf_resampleOptimal, 2},
    {"_dpf_getloglike", (DL_FUNC) &_dpf_getloglike, 3},
    {"_dpf_musicModel", (DL_FUNC) &_dpf_musicModel, 7},
    {"_dpf_musicModeldynamics", (DL_FUNC) &_dpf_musicModeldynamics, 9},
    {"_dpf_beamSearch", (DL_FUNC) &_dpf_beamSearch, 13},
    {"_dpf_kalman", (DL_FUNC) &_dpf_kalman, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dpf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
