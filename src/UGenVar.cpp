#include "Rcpp.h"
#include "UGenVar.h"


UGenVar::UGenVar() :
    m_u_var_col(0),
    m_n_days(0),
    m_w_idx(0),
    m_x_idx(0) {
}


UGenVar::UGenVar(Rcpp::NumericMatrix& u_rcpp,
		 Rcpp::IntegerVector& preg_map,
		 Rcpp::IntegerVector& sex_map,
		 int u_col) :
    m_u_var_col(u_rcpp.begin() + ((int) u_rcpp.nrow()) * u_col),
    m_n_days(u_rcpp.nrow()),
    m_w_idx(preg_map.begin()),
    m_x_idx(sex_map.begin()) {
}
