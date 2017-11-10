#include "Rcpp.h"
#include "CoefGen.h"
#include "GammaGen.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"

extern bool g_record_status;




CoefGen::CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp) :
    // initialization list
    m_gamma(GammaGen::create_arr(U, gamma_specs)),
    m_vals_rcpp(Rcpp::NumericVector(Rcpp::no_init(gamma_specs.size() * n_samp))),
    m_vals(m_vals_rcpp.begin()),
    m_n_psi(0),
    m_n_gamma(gamma_specs.size()),
    m_mu_0(0.1),
    m_nu_0(1.0) {
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

	if (m_gamma[j]->is_dsp_ar()) {

	    m_gamma[j]->set_hyp_b(m_nu_0);

	    // note that this assumes that the coefficients corresponding to the
	    // DSPs are in order relative to each other
	    if (is_first_ar) {
		m_gamma[j]->set_hyp_a(m_mu_0 * m_nu_0);
		is_first_ar = false;
	    }
	    else {
		m_gamma[j]->set_hyp_a(m_vals[j - 1] * m_nu_0);
	    }
	}

	m_vals[j] = m_gamma[j]->sample(W, xi, ubeta, X);
    }


}
