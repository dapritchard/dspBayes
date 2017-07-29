#include "Rcpp.h"
#include "WGen.h"
#include "DayBlock.h"




WGen::WGen(Rcpp::List& preg_cyc,
	   Rcpp::IntegerVector& w_to_days_idx,
	   Rcpp::IntegerVector& w_cyc_to_subj_idx,
	   int fw_len) :
    // initialization list
    m_vals(new int[w_to_days_idx.size()]),
    m_sums(new int[preg_cyc.size()]),
    m_days_idx(w_to_days_idx),
    m_subj_idx(w_cyc_to_subj_idx),
    m_preg_cyc(PregCyc::list_to_arr(preg_cyc)),
    m_n_preg_days(w_to_days_idx.size()),
    m_n_preg_cyc(preg_cyc.size()),
    m_fw_len(fw_len) {
}




WGen::~WGen() {
    delete[] m_vals;
    delete[] m_sums;
    delete[] m_preg_cyc;
}




void WGen::sample(XiGen& xi, UProdBeta& ubeta) {

    // point to the beginning of the arrays storing the `W_ijk` and `sum_k
    // W_ijk`
    int* curr_w = m_vals;
    int* curr_w_sum = m_sums;
    // point to the beginning of the array storing the current values of `xi`
    const double* xi_vals = xi.vals();
    // point to the beginning of the array storing the current values of `X_ijk
    // * exp( u_{ijk}^T beta )`
    const double* ubeta_exp_vals = ubeta.exp_vals();

    // scratch storage for multinomial probabilities
    double mult_probs[m_fw_len];

    // each iteration samples new values for the `W_ijk` that were both (i) in
    // cycles that resulted in a pregnancy and were also (ii) days in which
    // intercourse occurred (or at least there was a missing value for
    // intercourse).  Values of `sum_k W_ijk` are also stored for these cycles.
    for (int q = 0; q < m_n_preg_cyc; ++q) {

	// the day-specific index and number of days in the current cycle
	PregCyc curr_cyc = m_preg_cyc[q];
	int curr_beg_idx = curr_cyc.beg_idx;
	int curr_n_days = curr_cyc.n_days;
	int curr_subj_idx = curr_cyc.subj_idx;

	// variable to store the value of `sum_k W_ijk` for the current cycle
	double curr_sum_val = 0;

	// each iteration calculates `X_ijk * exp( u_{ijk}^T beta )` for the
	// v-th day in the current cycle with a random `W_ijk` in the cycle,
	// (i.e. a day with intercourse or at least a missing value for
	// intercourse), and adds it to `sum_val`
	for (int v = 0; v < curr_n_days; ++v) {

	    // day-specific index of the v-th day in the current cycle
	    int r = curr_beg_idx + v;

	    // copy and add in the `X_ijk * exp( u_{ijk}^T beta )` term to the
	    // running total for `sum_k W_ijk`
	    curr_sum_val += mult_probs[v] = ubeta_exp_vals[r];
	}

	// normalize the multinomial probabilities
	for (int v = 0; v < curr_n_days; ++v) {
	    mult_probs[v] /= curr_sum_val;
	}

	// calculate `xi_i * sum_k { X_ijk * exp( u_{ijk}^T beta ) }`
	double pois_mean = xi_vals[ curr_subj_idx ] * curr_sum_val;

	// sample new `sum_k W_ijk`   **********  TODO: must be truncated  *********************
	*curr_w_sum = rpois_zero_tr(pois_mean);

	// sample new `W_ij | { sum_k W_ijk }`
	rmultinom(*curr_w_sum, mult_probs, curr_n_days, curr_w);

	// update pointers to point to the next elements of `W_ijk` and `sum_k
	// W_ijk`
	curr_w += curr_n_days;
	++curr_w_sum;
    }
}




double WGen::rpois_zero_tr(double lambda) {
    double u = R::runif(exp(-lambda), 1);
    return R::qpois(u, lambda, 1, 0);
}
