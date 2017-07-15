#include "Rcpp.h"
#include "GammaGen.h"

// [[Rcpp::export]]
void dsp_sampler(Rcpp::NumericMatrix U) {

    WGen W();
    XiGen xi();
    CoefGen coefs();
    PhiGen phi();

    for (int s = 0; s < nSamp; s++) {

	// update the latent day-specific pregnancy variables W
	W.sample();

	// update the woman-specific fecundability multipliers xi
	xi.sample();

	// update the regression coefficients gamma and psi
	coefs.sample();

	// update exp_u_prod_beta

	// update phi, the variance parameter for xi
	phi.sample();
    }



}
