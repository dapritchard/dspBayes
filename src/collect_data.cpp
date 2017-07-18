#include "Rcpp.h"
#include "CoefGen.h"
#include "PhiGen.h"
#include "WGen.h"
#include "XiGen.h"


Rcpp::List collect_output(const CoefGen& regr_coefs,
			  const XiGen& xi,
			  const PhiGen& phi) {

    Rcpp::NumericVector coefs_data(regr_coefs.output_start(), regr_coefs.output_end());
    Rcpp::NumericVector xi_data(xi.output_start(), xi.output_end());
    Rcpp::NumericVector phi_data(phi.output_start(), phi.output_end());

    return Rcpp::List::create(Rcpp::Named("coefs") = coefs_data,
			      Rcpp::Named("xi")    = xi_data,
			      Rcpp::Named("phi")   = phi_data);
}
