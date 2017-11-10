#include "Rcpp.h"
#include "CoefGen.h"
#include "GammaGen.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"

extern bool g_record_status;




CoefGen::CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp) :
    m_gamma(GammaGen::create_arr(U, gamma_specs)),
    m_vals_rcpp(Rcpp::NumericVector(Rcpp::no_init(gamma_specs.size() * n_samp))),
    m_vals(m_vals_rcpp.begin()),
    m_n_psi(0),
    m_n_gamma(gamma_specs.size()),
    m_dsp_ar_vals(),
    m_n_dsp_ar_vals(),
    m_mu_0(0.1),
    m_nu_0(1.0),
    m_mu_accept_ctr(),
    m_nu_accept_ctr(),
    m_mu_hyp_a(),
    m_mu_hyp_b(),
    m_nu_hyp_a(),
    m_nu_hyp_b(),
    m_delta(),
    m_proposal_fcn() {
}




CoefGen::~CoefGen() {
    for (int h = 0; h < m_n_gamma; ++h) {
    	delete m_gamma[h];
    }
    delete[] m_gamma;
}




void CoefGen::sample(const WGen& W, const XiGen& xi, UProdBeta& ubeta, const int* X) {

    bool is_first_ar = true;

    // if we're past the burn-in phase then update `m_vals` so that we don't
    // overwrite the previous samples in the current scan
    if (g_record_status) {
	m_vals += m_n_gamma;
    }

    // each iteration updates one gamma_h term and correspondingly udjusts
    // the value of `ubeta`.
    for (int j = 0; j < m_n_gamma; ++j) {

	// case: coefficient is a DSP autregressive gamma
	if (m_gamma[j]->is_dsp_ar()) {

	    // set rate parameter of current gamma variable
	    m_gamma[j]->set_hyp_b(m_nu_0);

	    // note that the following if / then statements assume that the
	    // coefficients corresponding to the DSPs are in relative order

	    // case: first day of the fertile window
	    if (is_first_ar) {
		m_gamma[j]->set_hyp_a(m_mu_0 * m_nu_0);
		is_first_ar = false;
	    }
	    // case: not the first day of the fertile window
	    else {
		m_gamma[j]->set_hyp_a(m_vals[j - 1] * m_nu_0);
	    }
	}

	// sample new gamma coefficient
	m_vals[j] = m_gamma[j]->sample(W, xi, ubeta, X);
    }

    // case: the coefficients corresponding to the DSP have AR priors
    if (m_n_dsp_ar_vals != 0) {
	sample_mu();
	sample_nu();
    }
}




void CoefGen::sample_mu() {

    // assume that the coefficients corresponding to the DSPs are in relative
    // order
    double gamma_1_val = *m_dsp_ar_vals;

    // sample the proposal value for Metropolis step
    double proposal_mu = m_proposal_fcn(m_mu, m_delta);
    double mu_diff = proposal_val - m_mu_0;

    // calculate the log proposal ratio r
    double log_r = (m_nu_0 * mu_diff * log(m_nu_0 * gamma_1_val)
		    - R::lgammafn(proposal_mu * m_nu_0)
		    + log(m_mu_0 * m_nu_0)
		    - (m_mu_hyp_a - 1) * log(proposal_mu / m_mu_0)
		    - m_mu_hyp_b * mu_diff);

    // sample new ratio
    if ((log_r < 0) && (log(R::unif_rand()) < log_r)) {
	*m_mu_0 = proposal_mu;
	++m_mu_accept_ctr;
    }
}




void CoefGen::sample_nu() {

    // sample the proposal value for Metropolis step
    double proposal_nu = m_proposal_fcn(m_nu, m_delta);
    double nu_diff = proposal_val - m_nu_0;

    // calculate calculate each of the terms of the log proposal ratio r
    double log_nu_term_gamma_1      = calc_log_nu_term_gamma_1();
    double log_nu_term_gamma_remain = calc_log_nu_term_gamma_remain();
    double log_nu_term_nu_prior     = calc_log_nu_term_gamma_prior();

    // calculate the log proposal ratio r
    double log_r = log_nu_term_gamma_1 + log_nu_term_gamma_remain + log_nu_term_nu_prior;

    // sample new ratio
    if ((log_r < 0) && (log(R::unif_rand()) < log_r)) {
	*m_nu_0 = proposal_nu;
	++m_mu_accept_ctr;
    }
}




static double CoefGen::calc_nu_term_gamma_1(double proposal_nu, double nu_diff) {

    // assume that the coefficients corresponding to the DSPs are in relative
    // order
    double gamma_1_val = *m_dsp_ar_vals;

    return (m_mu_0 * proposal_nu * log(proposal_nu)
	    - R::lgammafn(m_mu_0 * proposal_nu)
	    - m_mu_0 * m_nu_0 * log(m_nu_0)
	    + R::lgamma(m_mu_0 * m_nu_0)
	    + m_mu_0 * nu_diff * log(gamma_1_val)
	    - gamma_1_val * nu_diff);
}




static double CoefGen::calc_log_nu_term_gamma_remain(double proposal_nu, double nu_diff) {

    // start gamma_k at the second entry.  Assume that the coefficients
    // corresponding to the DSPs are in relative order.
    double* gamma_k   = m_dsp_ar_vals + 1;
    double* gamma_end = m_dsp_ar_vals + m_n_dsp_ar_vals;

    // each iteration calculates the log ratio term contributed by gamma_k and
    // adds it to `sum_log_terms`
    double gamma_k_minus_1_val = *(gamma_k - 1);
    double sum_log_terms = 0.0;
    for ( ; gamma_k < gamma_end; ++gamma_k) {

	// update gamma_k value
	double gamma_k_val = *gamma_k;

	// calculate k-th term and add to running total
	sum_log_terms += (gamma_k_minus_1_val * proposal_nu * log(proposal_nu)
			  - R::lgammafn(gamma_k_minus_1_val * proposal_nu)
			  - gamma_k_minus_1_val * m_nu_0 * log(m_nu_0)
			  + R::lgammafn(gamma_k_minus_1_val * m_nu_0)
			  + gamma_k_minus_1_val * nu_diff * log(gamma_k_val)
			  - gamma_k_val * nu_diff);

	// update gamma_{k-1} value
	gamma_k_minus_1_val = gamma_k_val;
    }

    return sum_log_terms;
}




static double CoefGen::calc_log_nu_term_gamma_prior(double proposal_nu, double nu_diff) {

    return ((m_nu_hyp_a - 1) * log(proposal_nu / m_nu_0)
	    - m_nu_hyp_b * nu_diff);
}
