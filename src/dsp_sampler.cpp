#include "Rcpp.h"
#include "GammaGen.h"

// [[Rcpp::export]]
void dsp_sampler(Rcpp::NumericMatrix U,
		 // W
		 Rcpp::List preg_cyc,
		 Rcpp::IntegerVector w_days_idx,
		 Rcpp::IntegerVector w_cyc_idx,
		 int fw_len,
		 // xi
		 Rcpp::NumericVector xi_initial,
		 Rcpp::List subj_days,
		 // phi
		 Rcpp::NumericVector phi_hyper) {

    WGen W(preg_cyc, w_days_idx, w_cyc_idx, fw_len);
    XiGen xi(xi_initial, subj_days);
    CoefGen coefs();
    PhiGen phi();
    UProdBeta u_prod_beta(n_days);

    for (int s = 0; s < nSamp; s++) {

	// update the latent day-specific pregnancy variables W
	W.sample(xi, u_prod_beta);

	// update the woman-specific fecundability multipliers xi
	xi.sample(W, phi, u_prod_beta);

	// update the regression coefficients gamma and psi
	coefs.sample();
	u_prod_beta.update_exp();

	// update phi, the variance parameter for xi
	phi.sample(xi);
    }



}
