#include "Rcpp.h"
#include "UProdTau.h"


UProdTau::UProdTau(Rcpp::NumericVector& utau_rcpp) :
    // initialization list
    m_vals_rcpp(utau_rcpp),
    m_vals(m_vals_rcpp.begin()) {
}
