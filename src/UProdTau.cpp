#include "Rcpp.h"
#include "UProdTau.h"


UProdTau::UProdTau(Rcpp::NumericVector& utau, Rcpp::NumericVector& tau_coefs, double tau_prev) :
    // initialization list
    m_vals_rcpp(utau),
    m_vals(m_vals_rcpp.begin()),
    m_coefs(tau_coefs),
    m_sex_coef(tau_prev) {
}
