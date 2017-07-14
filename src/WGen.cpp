#include "WGen.h"

WGen::WGen(Rcpp::List& cycs_w_preg) {}




void WGen::samp() {

    // point to the beginning of the arrays storing the `W_ijk` and `sum_k
    // W_ijk`
    int* curr_W = m_W;
    int* curr_w_sum = m_w_sum;

    // each iteration samples new values for the `W_ijk` that were both in
    // cycles that resulted in a pregnancy and were also also days in which
    // intercourse occurred (or at least there was a missing value for
    // intercourse)
    for (int q = 0; q < m_n_preg_cyc; ++q) {`

	// the day-specific index and number of days in the current cycle
	PregCyc curr_cyc = m_preg_cyc[q];
	int curr_beg_idx = curr_cyc.day_beg_idx;
	int curr_n_days = curr_cyc.n_days;

	// variable to store the value of sum_k W_ijk for the current cycle
	double curr_w_sum = 0;

	// each iteration calculates `X_ijk * exp( u_{ijk}^T beta )` for the
	// v-th day in the current cycle with a random term in the cycle,
	// (i.e. a day with intercourse or at least a missing value for
	// intercourse), and adds it to `curr_w_sum`
	for (int v = 0; v < curr_n_days; ++v) {

	    // day-specific index of the v-th day in the current cycle
	    r = curr_beg_idx + v;

	    // calculate `X_ijk * exp( u_{ijk}^T beta )`
	    m_mult_probs[v] = X[r] ? exp(U_prod_beta[r]) : 0;

	    // add in the `X_ijk * exp( u_{ijk}^T beta )` term to the running
	    // total for `sum_k W_ijk`
	    curr_w_sum += m_mult_probs[v];
	}

	// normalize the multinomial probabilities
	for (int v = 0; v < n_days; ++v) {
	    m_mult_probs[v] /= curr_w_sum;
	}

	// calculate `xi_i * sum_k { X_ijk * exp( u_{ijk}^T beta ) }`
	pois_mean = xi[ curr_cyc.subj_idx ] * curr_w_sum;

	// sample new `sum_k W_ijk`
	*curr_w_sum = sample_w_sum();

	// sample new `W_ij | { sum_k W_ijk }`
	rmultinom(*curr_w_sum, m_mult_probs, n_days, curr_W);

	// update pointers to storage for `W_ijk` and `sum_k W_ijk`
	curr_W += n_days;
	++curr_w_sum;
    }

}
