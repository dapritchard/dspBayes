#ifndef DSP_BAYES_SRC_U_PROD_TAU_H
#define DSP_BAYES_SRC_U_PROD_TAU_H

#include "Rcpp.h"


class UProdTau {

 public:

    Rcpp::NumericVector& m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    const Rcpp::NumericVector& m_coefs;

    UProdTau(Rcpp::NumericVector& utau, Rcpp::List& tau_coefs);

    double* vals() { return m_vals; }
    const double* vals() const { return m_vals; }
    int n_days() const { return m_vals_rcpp.size(); }
};


#endif
