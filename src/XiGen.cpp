#include "Rcpp.h"
#include "XiGen.h"
#include "DayBlock.h"


XiGen::XiGen(Rcpp::NumericVector xi_initial, Rcpp::List subj_days) :
    // initialization list
    m_xi_vals(xi_initial.begin()),
    m_subj(new DayBlock[subj_days.size()]),
    m_n_subj(subj_days.size()) {

    // create `DayBlock` structs that provide the indices in the day-specific
    // data that correspond to the t-th individual
    for (int t = 0; t < m_n_subj; ++t) {
	m_subj[t] = DayBock(subj_days[t]["beg_idx"], preg_cyc[t]["n_days"]);
    }
}


XiGen::~XiGen() {
    delete m_subj;
}


const double* XiGen::sample(WGen W, PhiGen phi, UProdBeta& u_prod_beta) {

    int curr_idx, curr_end;
    int* preg_cyc_idx, w_sum;
    double curr_w_sum, curr_exp_ubeta_sum;
    double* exp_uprod_beta;
    DayBlock curr_subj;

    preg_cyc_idx = W.preg_cyc_idx();
    w_sum = W.sum_vals();
    phi_val = phi.val();
    exp_uprod_beta = u_prod_beta.exp_vals();

    // each iteration samples the i-th value of `xi_i` and stores it `m_xi_vals`
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
	curr_idx = curr_subj.beg_idx;
	curr_end = curr_idx + curr_subj.n_days;

	// obtain `sum_jk { exp( u_{ijk}^T beta ) }`
	curr_exp_ubeta_sum = 0;
	for ( ; curr_idx < curr_end; ++curr_idx) {
	    if (X[curr_idx]) {
		curr_exp_ubeta += exp_uprod_beta[curr_idx]);
	}

	// sample new value of `xi_i`
	m_xi_vals[i] = rgamma(phi_val + curr_w_sum, 1 / (phi_val + curr_exp_ubeta));
    }

    return m_xi_vals;
}
