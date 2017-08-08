#include <cmath>
#include "Rcpp.h"
#include "XGen.h"

#define SEX_NO    0
#define SEX_YES   1
#define SEX_MISS  2

#using std::exp;




void XGen::sample() {

    double prior_probs[m_n_max_miss][2];
    double posterior_probs[m_n_max_miss_pow2];

    XMiss* curr_miss_day = m_miss_day;

    for (int i = 0; i < m_n_xmiss; ++i) {

	day_before_fw_sex = calc_prior_probs(m_xmiss[i], prior_probs);

	n_perms = calc_model_probs();

	// sample the t-th permutation
	t = sample_x_perm(probs, n_perms);

	// update X with the values corresponding to the k-th permutation
	update_x(curr_miss_day, t);

	// update `curr_miss_day` to point to the next missing day
	curr_miss_day += n_miss;
    }

}




int XGen::calc_prior_probs(double* prior_probs, XMissCyc* curr_miss_cyc, XMissDay* curr_miss_day) {

    int curr_beg_idx = curr_miss_cyc->beg_idx;

    // calculate the prior probability of intercourse for the first missing element.  The reason
    // that we handle this seperately is because how we deal with a missing
    // intercourse value from the day before the fertile window differs than
    // for the remaining cases.  If this intercourse value is indeed missing
    // then we sample it before sampling the rest of the data.  Thus
    // whatever the cycle day that the first missing observation may be, the
    // day before it always have a nonmissing intercourse value.


    // Perform the calculations for the first missing element.  The reason
    // that we handle this seperately is because how we deal with a missing
    // intercourse value from the day before the fertile window differs than
    // for the remaining cases.  If this intercourse value is indeed missing
    // then we sample it before sampling the rest of the data.
    //
    // Note that the first missing day in the cycle need not be the first
    // day in the cycle.  But if that's the case then it won't have a
    // missing previous day of intercourse, so we don't need to consider a
    // special case for circumstance. Thus whatever the cycle day that the
    // first missing observation may be, after this procedure, that day
    // before it always have a nonmissing previous intercourse value.

    if ((curr_miss_day->prev == SEX_NO) ||
	((curr_miss_day->prev == SEX_MISS) && (R::unif_rand() > m_cohort_sex_prob))) {
	day_before_fw_sex = 0;
    }
    else {
	day_before_fw_sex = 1;
    }

    switch (day_before_fw_sex) {
    case (0):
	prior_probs[0][0] = 1 / (1 + exp(utau[curr_beg_idx]));
	day_before_fw_sex = 0;
	break;
    case (1):
	prior_probs[0][1] = 1 / (1 + exp(utau[curr_beg_idx] + tau.sex_coef));
	day_before_fw_sex = 1;
	break;
    }

    // calculate the remaining cases
    for (int k = 0; k < xmiss_cyc->n_miss; ++k) {

	switch (prev_day_sex[k]) {
	case SEX_NO:
	case SEX_MISS_X0_NO:
	    prior_probs[k][0] = 1 / (1 + exp(utau[xmiss_cyc->beg_idx + k]));
	    break;
	case SEX_YES:
	case SEX_MISS_X0_YES:
	    prior_probs[k][1] = 1 / (1 + exp(utau[xmiss_cyc->beg_idx + k] + tau.sex_prev_coef));
	    break;
	case SEX_MISS_DURING:
	    prior_probs[k][0] = 1 / (1 + exp(utau[xmiss_cyc->beg_idx + k]));
	    prior_probs[k][1] = 1 / (1 + exp(utau[xmiss_cyc->beg_idx + k] + tau.sex_prev_coef));
	    break;
	}
    }

    return day_before_fw_sex;
}




void XGen::calc_posterior_probs(double* posterior_probs,
				const XMissCyc* curr_miss_cyc,
				const XMissDay* curr_miss_day,
				const double* prior_probs,
				const UProdBeta& ubeta) {

    // the bit shift operator calculates `2^{ curr_miss_cyc->n_miss_days }`.
    // Thus we set aside enough memory to store the unnormalized posterior
    // probability of every possible realization of the missing intercourse
    // terms in the current cycle.
    int n_perms = 1 << curr_miss_cyc->n_miss_days;
    double probs[n_perms];

    // calculate the `sum_{k: nonrand} x_ijk * exp(u_{ijk}^T beta)`, i.e. the
    // nonrandom part of the sum that stays constant over the various
    // permutations of the random parts of `X`.
    double* ubeta_exp_vals = ubeta.exp_vals();
    double sum_nonrand_exp_ubeta; = calc_sum_nonrand_exp_ubeta(curr_miss_cyc,
							       curr_miss_day,
							       ubeta);

    // Each iteration `t` indexes a permutation of a vector of 0's and 1's from
    // among the possible realizations of missing itercourse data for the
    // `ij`-th cycle.  Then for each of the `x_{ijk_r}` elements in the vector
    // of missing intercourse data `(x_{ijk_1}, ..., x_{ijk_d})`, the
    // calculations
    //
    //     x_{ijk_r} * exp(u_{ijk_r}^T beta)
    //
    // and
    //
    //     P(X_{ijk_r} = x_{ijk_r} | X_{ij{k_r-1}})
    //
    // are performed, from which we can calculate the unnormalized posterior
    // probability of the `t`-th permutation.
    //
    // The permutation is determined by filling in the 0-th missing day with the
    // value of the 0-th bit of the binary representation of `t`, the 1-th
    // missing day with the value of the of the 1-th bit of the binary
    // representation of `t`, and so on.

    for (int t = 0; t < n_perms; ++t) {

	double sum_exp_ubeta = nonrandom_sum_exp_ubeta;
	double prod_priors;
	int curr_x;

	// perform the calculations for the first missing element.  The reason
	// that we handle this seperately is because for the first element we
	// are guaranteed to have a nonmissing value of intercourse for the
	// previous day, and its status is stored in the variable
	// `day_before_fw_sex`.

	// case: x_{ijk_1} = 1 for the `t`-th permutation
	if (t & 1) {
	    sum_exp_ubeta += ubeta_exp_vals[curr_miss_day[0]];
	    prod_priors = priors_probs[0][day_before_fw_sex];
	    curr_x = 1;
	}
	// case: x_{ijk_1} = 0 for the `t`-th permutation.  Note that the
	// probability is `1 - P(X_{ijk_1} = 1 | X_{ij{k_1-1}})`.
	else {
	    prod_priors = 1 - prior_probs[0][day_before_fw_sex];
	    curr_x = 0;
	}

	// perform the calculations for the remaining missing elements.  Now it
	// is possible that the previous day of itercourse was missing.  If it
	// was not missing, then we find its value using `curr_miss_day[r].prev`
	// and look up the corresponding prior probability in the `prior_probs`
	// table.  If it was missing, then it must have been the immediately
	// preceeding value of intercourse under the `t`-th permutation, which
	// is stored as `curr_x`.

	for (int r = 1; r < curr_miss_cyc->n_miss; ++r) {

	    // case: x_{ijk_r} = 1 for the current permutation
	    if (t & (1<<r)) {
		sum_exp_ubeta += ubeta_exp_vals[miss_day[r].idx];
		prod_priors *= (curr_miss_day[r].prev != SEX_MISS) ?
		    prior_probs[r][curr_miss_day[r].prev] :
		    prior_probs[r][curr_x];
		curr_x = 1;
	    }
	    // case: x_{ijk_r} = 0 for the current permutation.  Note that the
	    // probability is `1 - P(X_{ijk_r} = 1 | X_{ij{k_r-1}})`.
	    else {
		prod_priors *= (curr_miss_day[r].prev != SEX_MISS) ?
		    1 - prior_probs[r][curr_miss_day[r].prev] :
		    1 - prior_probs[r][curr_x];
		curr_x = 0;
	    }
	}

	// Now calculate the unnormalized posterior probability for the `t`-th
	// permutation.  The first step in the random variable case calculates
	//
	//     a := xi_i * sum_r x_{ijk_r} * exp(u_{ijk_r}^T beta)
	//
	// The second step calculates `P(Y_ij = y_ij | data)`, where
	//
	//     P(Y_ij = 1 | data) = 1 - exp(-a)
	//
	// Then in the third step, the unnormalized posterior probability is
	// given by
	//
	//     P(Y_ij = y_ij | data) P(X_ij | data)
	//
	//         = P(Y_ij = y_ij | data) prod_k { P(X_ijk | X_{ij{k-1}}) }

	// case: at least one day with intercourse that occurred during the
	// cycle, and hence `Y_ij` is a random variable.
	if (sum_exp_ubeta > 0) {
	    sum_exp_ubeta *= *(xi.vals() + curr_miss_cyc.subj_idx);
	    double prob_y_ij = (curr_miss_cyc->preg) ?
		1 - exp(-sum_exp_ubeta) :
		exp(-sum_exp_ubeta);
	    probs[t] = prob_y_ij * prod_priors;
	}
	// case: no days of intercourse occurred during the cycle under the
	// `t`-th permutation of missing intercourse values, so no pregnancy can
	// possibly occur.
	else {
	    probs[t] = curr_miss_cyc->preg ? 0 : 1;
	}
    }

    return n_perms;
}




// Sampling from a multinomial distribution with n bins and 1 trial.  The idea
// is to sample a uniform random variable and find which bin it falls into,
// where the bins are as shown below.
//
//    probs[0]    probs[1]   probs[2]                  probs[n-1]
//  /         \ /         \ /        \                 /        \
// |-----------|-----------|----------|----  ...  ----|----------|

inline int XGen::sample_x_perm(double* probs, int n_perms) {

    double sum_probs, u;
    int k;

    // sum of the unnormalized probabilities
    sum_probs = std::accumulate(probs, probs + n_miss, 0.0);

    // sample from a `unif(0, sum_probs)` distribution
    u = R::unif_rand() * sum_probs;

    // find the bin that `u` falls into
    sum_probs = *probs++;
    k = 0;
    while (u > sum_probs) {
	sum_probs += *probs;
	++probs;
	++k;
    }

    return k;
}




// fill in the updated values of X for the current cycle.  Fills in the 0-th
// missing day with the value of the 0-th bit of the binary representation of
// `t`, the 1-th missing day with the value of the of the 1-th bit of the binary
// representation of `t`, and so on.

inline void XGen::update_cyc_x(XMissDay* curr_miss_day, int n_miss, int t) {
    for (int r = 0; r < curr_miss_day->n_miss; ++r) {
	m_vals[curr_miss_day[r].idx] = t & (1<<r) ? 1 : 0;
    }
}




inline double XGen::calc_nonrandom_sum_exp_ubeta(XMissCyc* curr_miss_cyc,
						 XMissDay* curr_miss_day,
						 UProdBeta& ubeta) {

    double* ubeta_exp_vals = ubeta.exp_vals();
    double sum_val = 0;
    int curr_idx = curr_miss_cyc->beg_idx;
    int end_idx = curr_idx + curr_miss_cyc->n_days;

    for ( ; curr_idx < end_idx; ++curr_idx) {

	// case: the current index isn't a missing value for `X_ijk`.  Note that
	// guarantees that it has a value of 1 or else it wouldn't be in the
	// data.
	if (curr_idx != curr_miss_day.idx) {
	    sum_val += ubeta_exp_vals[curr_idx];
	}
	// case: the current index is a missing value for `X_ijk`.  Don't do
	// anything with `sum_val`, and increment `curr_miss_day` so as to look
	// for the next day with missing intercourse data.
	else {
	    ++curr_miss_day;
	}
    }

    return sum_val;
}
