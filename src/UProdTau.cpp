#include "Rcpp.h"
#include "UProdTau.h"

using Rcpp::as;


UProdTau::UProdTau(Rcpp::NumericVector& utau, Rcpp::List& tau_coefs) :
    m_vals_rcpp(utau),
    m_vals(m_vals_rcpp.begin()),
    m_coefs_rcpp(as<Rcpp::NumericVector>(tau_coefs["u_coefs"])),
    m_coefs(m_coefs_rcpp.begin()) {
}
