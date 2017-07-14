#include "XiGen.h"
#include "global_vars.h"


const double* XiGen::sample(WGen W, PhiGen phi, double* exp_uprod_beta) {

    int curr_idx, curr_end;
    int* preg_cyc_idx, w_sum;
    double curr_w_sum, curr_exp_ubeta_sum;
    DayBlock curr_subj;

    preg_cyc_idx = W.preg_cyc_idx();
    w_sum = W.w_sum();
    phi_val = phi.val();

    // each iteration samples the i-th value of `xi_i` and stores it `m_xi`
    for (int i = 0; i < m_n_subj; ++i) {

	// obtain `sum_jk W_ijk`
	if (i == *preg_cyc_idx) {
	    curr_w_sum = *w_sum++;
	    ++preg_cyc_idx;
	}
	else {
	    curr_w_sum = 0;
	}

	// index in the day-specific data of the first day and one past the last
	// day for the current subject
	curr_subj = m_subj[i];
	curr_idx = curr_subj.day_beg_idx;
	curr_end = curr_idx + curr_subj.n_days;

	// obtain `sum_jk { exp( u_{ijk}^T beta ) }`
	curr_exp_ubeta_sum = 0;
	for ( ; curr_idx < curr_end; ++curr_idx) {
	    if (X[curr_idx]) {
		curr_exp_ubeta += exp_uprod_beta[curr_idx]);
	}

	// sample new value of xi_i
	m_xi[i] = rgamma(phi_val + curr_w_sum, 1 / (phi_val + curr_exp_ubeta));
    }

    return m_xi;
}
