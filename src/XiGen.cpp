class XiGen {

    double* m_xi;

    const DayBlock* m_subj;

    const int m_n_days;
    const int m_n_subj;
};


void XiGen::sample(WGen W, PhiGen phi) {

    int curr_idx, curr_end;
    int* preg_cyc_idx, w_sum;
    double curr_w_sum, curr_exp_ubeta;
    DayBlock curr_subj;

    preg_cyc_idx = W.preg_cyc_idx();
    w_sum = W.w_sum();

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

	// ************  redo using established values for exp_ubeta *************
	// obtain `sum_jk { exp( u_{ijk}^T beta ) }`
	curr_exp_ubeta = 0;
	for ( ; curr_idx < curr_end; ++curr_idx) {

	    if (X[curr_idx]) {
		curr_w_sum += W[curr_idx];
		curr_exp_ubeta += exp(U_prod_beta[curr_idx]);
	    }
	}

	// sample new value of xi_i
	double phi_val = phi.val();
	m_xi[i] = rgamma(phi_val + curr_w_sum, 1 / (phi_val + curr_exp_ubeta));
    }
}
