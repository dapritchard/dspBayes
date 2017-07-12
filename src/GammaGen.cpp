#include "GammaGen.h"

GammaGen::GammaGen(Rcpp::NumericMatrix& U, int h) :
    Uh_begin(U.begin() + U.nrow() * h),
    Uh_end(Uh_begin + U.nrow()),
    val(h),
    hyp_p(0.5),
    hyp_a(1),
    hyp_b(1),
    bnd_l(0),
    bnd_u(9999) {

    for (int k = 0; k < 5; k++) {
	Rcpp::Rcout << *(Uh_begin + k) << "  ";
    }
    Rcpp::Rcout << "    " << val << std::endl;
}
