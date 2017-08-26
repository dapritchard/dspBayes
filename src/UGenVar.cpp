#include "Rcpp.h"
#include "UGenVar.h"


UGenVar::UGenVar(Rcpp::IntegerVector& preg_map,
		 Rcpp::IntegerVector& sex_map) :
    m_w_idx(preg_map.begin()),
    m_x_idx(sex_map.begin()) {
}
