#include "Rcpp.h"
#include "GammaGen.h"

// [[Rcpp::export]]
void dsp_sampler(Rcpp::NumericMatrix U) {

    // int n = U.nrow;
    // int p = U.ncol;
    int p = 5;

    GammaGen** gamma = new GammaGen*[p];
    GammaGen** gamma_end = gamma + p;

    for (int h = 0; h < p; ++h) {
    	*(gamma + h) = new GammaGen(U, h);
    }




    // int s;

    // bool is_burn;

    // vector<GammaGen>::const_iterator var_end = var_gen.end();


    // for (int s = 0; s < nSamp; s++) {

	// samp_W();

	// // sample regression coefficients gamma
	// for (vector<VarGen>::const_iterator curr = var_gen.begin(); curr != var_end; ) {
	//     gam_val = curr->sample_gam();
	//     if (++curr != end) {
	// 	curr->set_ar_param(gam_val);
	//     }
	// }

	// // sample woman-specific fecundability multiplier xi
	// xi = samp_xi();

	// // sample phi, the variance parameter for xi
	// phi = samp_phi();


	// samp_X(X);

	// for (int i = 0; i < var.size(); i++) {
	//     var.samp_miss(U);
	// }
    // }


    for (int h = 0; h < p; ++h) {
    	delete *(gamma + h);
    }
    delete gamma;

}
