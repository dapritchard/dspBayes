#include "Rcpp.h"
#include "CoefGen.h"
#include "GammaGen.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"




CoefGen::CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp) :
    // initialization list
    m_gamma(GammaGen::create_arr(U, gamma_specs)),
    m_vals(new double[gamma_specs.size() * n_samp]),
    m_output_start(m_vals),
    m_output_end(m_output_start + (gamma_specs.size() * n_samp)),
    m_n_psi(0),
    m_n_gamma(gamma_specs.size()) {
}




CoefGen::~CoefGen() {
    for (int h = 0; h < m_n_gamma; ++h) {
    	delete m_gamma[h];
    }
    delete[] m_gamma;
    delete[] m_output_start;
}




void CoefGen::sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta, const int* X) {

    GammaGen** end = m_gamma + m_n_gamma;

    // each iteration updates one gamma_h term and correspondingly udjusts
    // the value of `u_prod_beta`.
    for (GammaGen** curr = m_gamma; curr != end; ++curr) {

	// update the regression coefficient gamma_h
	*m_vals++ = (*curr)->sample(W, xi, u_prod_beta, X);

	// if (++curr != end) {
	//     curr->set_ar_param(gam_val);
	// }
    }
}
