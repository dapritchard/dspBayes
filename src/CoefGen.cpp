#include "Rcpp.h"


CoefGen::CoefGen(Rcpp::List& gamma_list) :
    // initialization list
    gamma(init_gamma(gamma_list)) {

    m_coef = new GammaGen*[p];

    for (int h = 0; h < p; ++h) {
    	*(m_coef + h) = new GammaGen(U, h);
    }
}




CoefGen::~CoefGen() {


    // for (int h = 0; h < p; ++h) {
    // 	delete *(m_coef + h);
    // }
    // delete m_coef;
}




CoefGen::sample() {


    // each iteration updates one gamma_h term and correspondingly udjusts
    // the value of `u_prod_beta`.
    for (GammaGen** curr = gamma; curr != end; ++curr) {

	// update the regression coefficient gamma_h
	gam_val = curr->sample();
	// if (++curr != end) {
	//     curr->set_ar_param(gam_val);
	// }
    }
}
