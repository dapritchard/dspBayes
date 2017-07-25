#include "Rcpp.h"
#include "XiGen.h"




XiGen::XiGen(Rcpp::List subj_day_blocks, int n_samp) :
    // initialization list
    m_xi_vals(new double[subj_day_blocks.size() * n_samp]),
    m_output_start(m_xi_vals),
    m_output_end(m_output_start + (subj_day_blocks.size() * n_samp)),
    m_subj(DayBlock::list_to_arr(subj_day_blocks)),
    m_n_subj(subj_day_blocks.size()),
    m_record_status(false) {

    // initialize values for all subjects to 1 (i.e. no fecundability effect)
    for (int i = 0; i < m_n_subj; ++i) {
    	m_xi_vals[i] = 1;
    }
}




XiGen::~XiGen() {
    delete[] m_output_start;
    delete[] m_subj;
}




void XiGen::sample(const WGen& W, const PhiGen& phi, const UProdBeta& u_prod_beta) {

    int curr_idx, curr_end;
    double curr_w_sum, curr_sum_exp_ubeta;

    const int* w_cyc_idx = W.cyc_idx();
    const int* w_sum_vals = W.sum_vals();
    const double* exp_uprod_beta = u_prod_beta.exp_vals();
    double phi_val = phi.val();

    // each iteration samples the i-th value of `xi_i` and stores it `m_xi_vals`
    for (int i = 0; i < m_n_subj; ++i) {

    	// obtain `sum_jk W_ijk`
    	if (i == *w_cyc_idx) {
    	    curr_w_sum = *w_sum_vals++;
    	    ++w_cyc_idx;
    	}
    	else {
    	    curr_w_sum = 0;
    	}

    	// index in the day-specific data of the first day and one past the last
    	// day for the current subject
    	curr_idx = m_subj[i].beg_idx;
    	curr_end = curr_idx + m_subj[i].n_days;

    	// obtain `sum_jk { X_ijk * exp( u_{ijk}^T beta ) }`
    	curr_sum_exp_ubeta = 0;
    	for ( ; curr_idx < curr_end; ++curr_idx) {
    	    // exp_uprod_beta is already 0 if corresponding `X_ijk` is 0
    	    curr_sum_exp_ubeta += exp_uprod_beta[curr_idx];
    	}

    	// sample new value of `xi_i`
    	m_xi_vals[i] = R::rgamma(phi_val + curr_w_sum, 1 / (phi_val + curr_sum_exp_ubeta));
    }

    // if we are past the burn phase then move the pointer past the samples so
    // that we don't overwrite them
    if (m_record_status) {
	m_xi_vals += m_n_subj;
    }
}
