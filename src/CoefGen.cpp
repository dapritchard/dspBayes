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
    m_n_gamma(gamma_specs.size()) {
}




CoefGen::~CoefGen() {
    for (int h = 0; h < m_n_gamma; ++h) {
    	delete m_gamma[h];
    }
    delete[] m_gamma;
}




void CoefGen::sample(const WGen& W, const XiGen& xi, UProdBeta& ubeta, const int* X) {

    // each iteration updates one gamma_h term and correspondingly udjusts
    // the value of `ubeta`.
    GammaGen** end = m_gamma + m_n_gamma;
    for (GammaGen** curr = m_gamma; curr != end; ++curr) {

	// update the regression coefficient gamma_h.  If we are past the
	// burn-in phase then store the samples in `m_vals`.
	if (g_record_status) {
	    *m_vals++ = (*curr)->sample(W, xi, ubeta, X);
	} else {
	    // same function call as the first case, but now we don't keep the
	    // result
	    (*curr)->sample(W, xi, ubeta, X);
	}
    }
}
