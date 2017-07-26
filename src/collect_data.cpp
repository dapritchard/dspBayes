#include "Rcpp.h"
#include "CoefGen.h"
#include "PhiGen.h"
#include "WGen.h"
#include "XiGen.h"


Rcpp::List collect_output(const CoefGen& regr_coefs,
			  const XiGen& xi,
			  const PhiGen& phi) {

    Rcpp::NumericVector coefs_data(regr_coefs.output_start(), regr_coefs.output_end());

    return Rcpp::List::create(Rcpp::Named("coefs") = coefs_data,
			      Rcpp::Named("xi")    = xi.m_vals_rcpp,
			      Rcpp::Named("phi")   = phi.m_vals_rcpp);
}
