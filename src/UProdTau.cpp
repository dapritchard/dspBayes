#include "Rcpp.h"
#include "UProdTau.h"

using Rcpp::as;


UProdTau::UProdTau(Rcpp::NumericVector& utau, Rcpp::List& tau_coefs) :
    // initialization list
    m_vals_rcpp(utau),
    m_vals(m_vals_rcpp.begin()),
    m_coefs(as<Rcpp::NumericVector>(tau_coefs["u_coefs"])) {
}
