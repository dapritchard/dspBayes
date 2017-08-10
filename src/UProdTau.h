#ifndef DSP_BAYES_SRC_U_PROD_TAU_H
#define DSP_BAYES_SRC_U_PROD_TAU_H

#include "Rcpp.h"


class UProdTau {

 public:

    Rcpp::NumericVector& m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    const Rcpp::NumericVector& m_coefs;

    const double m_sex_coef;

    UProdTau(Rcpp::NumericVector& utau_rcpp, Rcpp::NumericVector& tau_coefs, double tau_prev);
    ~UProdTau();

    double* vals() { return m_vals; }
    const double* vals() const { return m_vals; }
    double sex_coef() const { return m_sex_coef; }
    int n_days() const { return m_vals_rcpp.size(); }
};


#endif
