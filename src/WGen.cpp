#include "Rcpp.h"
#include "WGen.h"


WGen::WGen(Rcpp::List& preg_cyc,
	   Rcpp::IntegerVector& w_days_idx,
	   Rcpp::IntegerVector& w_cyc_idx,
	   int fw_len) :
    // initialization list
    m_w_vals(new int[w_len]),
    m_w_sums(new int[preg_cyc.size()]),
    m_w_days_idx(w_days_idx.begin()),
    m_w_cyc_idx(w_cyc_idx.begin()),
    m_preg_cyc(new DayBlock[preg_cyc.size()]),
    m_n_preg_cyc(n_preg_cyc),
    m_mult_probs(new double[fw_len]) {

    // create `DayBlock` structs that provide the indices in the day-specific
    // data that correspond to the t-th pregnancy cycle
    for (int t = 0; t < m_n_preg_cyc; ++t) {
	m_preg_cyc[t] = DayBock(preg_cyc[t]["beg_idx"], preg_cyc[t]["n_days"]);
    }
}




WGen::~WGen() {
    delete m_w_vals;
    delete m_w_sums;
    delete m_mult_probs;
    delete m_preg_cyc;
}




const int* WGen::sample(XiGen& xi, UProdBeta& u_prod_beta) {

    // point to the beginning of the arrays storing the `W_ijk` and `sum_k
    // W_ijk`
    int* curr_w = m_w_vals;
    int* curr_w_sum = m_w_sums;
    // point to the beginning of the array storing the current values of `xi`
    const double* xi_vals = xi.vals();
    // point to the beginning of the array storing the current values of `X_ijk
    // * exp( u_{ijk}^T beta )`
    const double* exp_uprod_beta = u_prod_beta.exp_vals();

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
	double sum_val = 0;

	// each iteration calculates `X_ijk * exp( u_{ijk}^T beta )` for the
	// v-th day in the current cycle with a random `W_ijk` in the cycle,
	// (i.e. a day with intercourse or at least a missing value for
	// intercourse), and adds it to `sum_val`
	for (int v = 0; v < curr_n_days; ++v) {

	    // day-specific index of the v-th day in the current cycle
	    r = curr_beg_idx + v;

	    // copy and add in the `X_ijk * exp( u_{ijk}^T beta )` term to the running
	    // total for `sum_k W_ijk`
	    m_mult_probs[v] = exp_uprod_beta[r];
	    curr_w_sum += m_mult_probs[v];
	}

	// normalize the multinomial probabilities
	for (int v = 0; v < curr_n_days; ++v) {
	    m_mult_probs[v] /= curr_w_sum;
	}

	// calculate `xi_i * sum_k { X_ijk * exp( u_{ijk}^T beta ) }`
	double pois_mean = xi_vals[ curr_subj_idx ] * exp(curr_w_sum);

	// sample new `sum_k W_ijk`
	*curr_w_sum = R::rpois(pois_mean);

	// sample new `W_ij | { sum_k W_ijk }`
	rmultinom(*curr_w_sum, m_mult_probs, curr_n_days, curr_w);

	// update pointers to point to the next elements of `W_ijk` and `sum_k
	// W_ijk`
	curr_w += n_days;
	++curr_w_sum;
    }

    return m_w_vals;
}
