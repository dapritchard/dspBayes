#include "Rcpp.h"
#include "CoefGen.h"
#include "PhiGen.h"
#include "WGen.h"
#include "XiGen.h"

int* d2s;

#define DSP_BAYES_N_INTER_CHECK 1000

Rcpp::List collect_output(const CoefGen& regr_coefs,
			  const XiGen& xi,
			  const PhiGen& phi);

// preg_cycle     used when sampling W
// w_days_idx     categorical gamma: a_tilde
// w_cyc_idx      used when sampling xi (first term)
// fw_len         how much memory to set aside when sampling W in a cycle
// subj_days      used when sampling xi (second term)
// gamma_specs    gamma hyperparameters
// phi_hyper      phi hyperparameters




// [[Rcpp::export]]
Rcpp::List dsp_sampler(Rcpp::NumericMatrix U,
		       Rcpp::IntegerVector X_rcpp,
		       Rcpp::List w_day_blocks,
		       Rcpp::IntegerVector w_to_days_idx,
		       Rcpp::IntegerVector w_cyc_to_cyc_idx,
		       Rcpp::List subj_day_blocks,
		       Rcpp::IntegerVector day_to_subj_idx,
		       Rcpp::List gamma_specs,
		       Rcpp::NumericVector phi_hyper,
		       int fw_len,
		       int n_burn,
		       int n_samp) {

    // create data objects
    WGen W(preg_cyc, w_days_idx, w_cyc_idx, fw_len);
    XiGen xi(subj_days, n_samp);
    CoefGen regr_coefs(U, gamma_specs, n_samp);
    PhiGen phi(phi_hyper, n_samp);
    UProdBeta u_prod_beta(U.size());
    int* X = X_rcpp.begin();
    d2s = subj_idx.begin();

    // begin sampler loop
    for (int s = 0; s < n_samp; s++) {

    	// update the latent day-specific pregnancy variables W
    	W.sample(xi, u_prod_beta);

    	// // update the woman-specific fecundability multipliers xi
    	xi.sample(W, phi, u_prod_beta);

    	// update the regression coefficients gamma and psi
    	regr_coefs.sample(W, xi, u_prod_beta, X);
    	u_prod_beta.update_exp(X);

    	// update phi, the variance parameter for xi
    	phi.sample(xi);

	// call the `record` family of functions to inform the various functions
	// to begin recording their MCMC samples
	if (s == 0) {
	    // xi.record();
	    // regr_coefs.record();
	    phi.record();
	}

	// check for user interrupt every `DSP_BAYES_N_INTER_CHECK` iterations
	if ((s % DSP_BAYES_N_INTER_CHECK) == 0) Rcpp::checkUserInterrupt();
    }

    return collect_output(regr_coefs, xi, phi);
    // Rcpp::List ret;
    // return ret;
}
