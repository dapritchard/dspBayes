#include "Rcpp.h"
#include "CoefGen.h"
#include "GammaGen.h"


CoefGen::CoefGen(Rcpp::List& gamma_specs) :
    // initialization list
    m_gamma(GammaGen::list_to_arr(gamma_specs)),
    m_n_psi(0),
    m_n_gamma(gamma_specs.size()) {
}




CoefGen::~CoefGen() {
    for (int h = 0; h < m_n_gamma; ++h) {
    	delete m_gamma[h];
    }
    delete[] m_gamma;
}




void CoefGen::sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta) {

    GammaGen** end = m_gamma + m_n_gamma;

    // each iteration updates one gamma_h term and correspondingly udjusts
    // the value of `u_prod_beta`.
    for (GammaGen** curr_gamma = m_gamma; curr_gamma != end; ++curr_gamma) {

	// update the regression coefficient gamma_h
	(*curr_gamma)->sample(W, xi, u_prod_beta);

	++curr_gamma;

	// if (++curr != end) {
	//     curr->set_ar_param(gam_val);
	// }
    }
}
