#include "Rcpp.h"

// [[Rcpp::export]]
void dsp_sampler() {

    int s;

    bool is_burn;

    vector<VarGen>::const_iterator var_end = var_gen.end();


    for (int s = 0; s < nSamp; s++) {

	samp_W();

	for (vector<VarGen>::const_iterator curr = var_gen.begin(); curr != var_end; curr++) {
	    curr->sample_gam();
	}

	xi = samp_xi();

	phi = samp_phi();

	samp_X(X);

	for (int i = 0; i < var.size(); i++) {
	    var.samp_miss(U);
	}
    }




}
