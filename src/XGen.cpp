#include <cmath>
#include "Rcpp.h"
#include "DayBlock.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"

// #define SEX_YES            1
// #define SEX_IMPUTE_SHIFT   2
// #define NON_PREG_CYC      -1

#define IS_FW_ZERO_MISS -99
#define SEX_SHIFT         2



// constructor
XGen::XGen(Rcpp::IntegerVector& X_rcpp,
	   Rcpp::IntegerVector& x_miss,
	   Rcpp::IntegerVector& sex_miss_to_w,
	   int fw_len,
	   double cohort_sex_prob,
	   double sex_coef) :
    m_x_rcpp(X_rcpp),
    m_vals(X_rcpp.begin()),
    m_x_miss_rcpp(x_miss),
    m_x_miss(x_miss.begin()),
    m_x_miss_end(x_miss.end()),
    m_sex_miss_to_w_rcpp(sex_miss_to_w_rcpp),
    m_sex_miss_to_w(sex_miss_to_w.begin()),
    m_ext_fw_len(fw_len + 1),
    m_cohort_sex_prob(cohort_sex_prob),
    m_sex_coef(sex_coef) {
}




// // destructor
// XGen::~XGen() {
// }




// update the missing values of X
void XGen::sample(const WGen& W,
		  const XiGen& xi,
		  const UProdBeta& ubeta,
		  const UProdTau& utau) {

    const int* w_vals      = W.vals();
    const int* utau_vals   = utau.vals();
    const int* curr_x_to_w = m_sex_miss_to_w;

    // each iteration samples the missing intercourse values for the current
    // cycle
    for (curr_x_miss = m_x_miss;
	 curr_x_miss != m_x_miss_end;
	 curr_x_miss += m_ext_fw_len) {

	sample_cycle(curr_x_miss, &curr_x_to_w, w_vals, xi, ubeta, utau);
    }
}




void XGen::sample_cycle(const int* x_cyc,
			const int** miss_to_w_map,
			const int* w_vals,
			const XiGen& xi,
			const UProdBeta& ubeta,
			const UProdTau& utau) {

    int sex_prev_day, sex_curr_day, sex_next_day, x_curr_idx, x_next_idx;
    boolean is_sex_miss, need_to_samp_bool;

    // value of xi for the subject that `miss_cycl` corresponds to
    const double xi_i = xi.vals()[miss_cyc->subj_idx];

    // obtain intercourse for the 0-th fertile window day
    if (x_cyc[0] == IS_FW_ZERO_MISS) {
	sex_prev_day = sample_day_before_fw_sex();
    }
    else {
	sex_prev_day = x_cyc[0] + SEX_SHIFT;
    }

    // obtain intercourse status for the 1-th fertile window day, and track
    // whether or not it is missing
    x_curr_idx = x_cyc[1];
    if (check_if_sex_miss(x_curr_idx)) {
	sex_curr_day = m_vals[x_curr_idx];
	is_sex_miss = true;
    }
    else {
	sex_curr_day = x_curr_idx + SEX_SHIFT;
	is_sex_miss = false;
    }

    for (int r = 1; r < m_ext_fw_len; ++r) {

	// set initial state of `need_to_samp_bool` to the value of
	// `is_sex_miss`, which was obtained when looking at the (r + 1)-th day
	// in the previous iteration
	need_to_samp_bool = is_sex_miss;

	// get value of W_ijk and update `miss_to_w_map`
	curr_w = (*miss_to_w_map == NON_PREG_CYC) ? 0 : w_vals[*miss_to_w_map];
	miss_to_w_map++;
	// case: W_ijk > 0, so X_ijk must have a value of 1
	if (curr_w > 0) {
	    m_vals[x_curr_idx] = sex_prev_day = 1;
	    need_to_samp_bool = false;
	}

	// obtain intercourse for the (r + 1)-th fertile window day
	if (r + 1 < m_ext_fw_len) {

	    // note that `is_sex_miss` informs the next iteration whether we
	    // need to sample for intercourse
	    x_next_idx = x_cyc[r + 1];
	    is_sex_miss = check_if_sex_miss(x_next_idx);

	    sex_next_day = is_sex_miss ?
		m_vals[x_next_idx] :
		x_next_idx + SEX_SHIFT;
	}

	// case: need to sample because (i) intercourse status was missing today
	// and (ii) `W_ijk = 0`
	if (need_to_samp_bool) {

	    // obtain `P(W = 0 | X_ijk = 1)`.  Note that we don't need to
	    // calculate this when `X_ijk = 0`, since then the expression is
	    // simply 1.
	    p_w_xIsOne = calc_p_w_zero(ubeta, xi_i, curr_day_idx);

	    // obtain `P(X_ijk = 1 | X_{ij,k-1})`
	    p_xTodayIsOne = calc_p_xTodayIsOne(utau, r, sex_prev_day);

	    // obtain `p(X_{ij,k+1} | X_ijk = 0)` and `p(X_{ij,k+1} | X_ijk = 1)`
	    //
	    // case: this is the last day of the fertile window, so we can
	    // effectively drop the `p(X_{ij,k+1} | X_ijk = t)` terms by making
	    // them both equal to 1
	    if (r + 1 == m_ext_fw_len) {
		p_xTom_xTodayIsZero = p_xTom_xTodayIsOne = 1;
	    }
	    // case: not the last day of the fertile window, so have to
	    // calculate both terms
	    else {
		p_xTom_xTodayIsZero = calc_p_xTom(sex_next_day, 0, utau_idx);
		p_xTom_xTodayIsOne  = calc_p_xTom(sex_next_day, 1, utau_idx);
	    }

	    // sample missing intercourse for today
	    m_vals[x_idx] = sex_curr_day = sample_x_ijk(p_w_xIsOne,
							p_xTodayIsOne,
							p_xTom_xTodayIsZero,
							p_xTom_xTodayIsOne);
	}

	// slide the 3-day sequence over 1 in preparation for the next iteration
	sex_prev_day = sex_curr_day;
	sex_curr_day = sex_next_day;
	x_curr_idx   = x_next_idx;
    }
}




double XGen::calc_p_xTom(int x_tomorrow,
			 int x_today,
			 int utau_idx) {

    const double* utau_vals = utau.vals();

    // calculate for `X_ij{k+1} = 1` and then adjust at the end for the 0 case
    // if necessary
    double p_xTomIsOne = x_today ?
	1 / (1 + exp(-utau_vals[utau_idx] - m_sex_coef)) :
	1 / (1 + exp(-utau_vals[utau_idx]));

    return x_tommorow ?
	p_xTomIsOne :
	1 - p_xTomIsOne;
}




inline double XGen::calc_p_xTodayIsOne(const UProdTau& utau,
				       const int miss_day_idx,
				       const int sex_prev_day) const {

    const double* utau_vals = utau.vals();
    return sex_prev_day ?
	1 / (1 + exp(-utau_vals[miss_day_idx] - m_sex_coef)) :
	1 / (1 + exp(-utau_vals[miss_day_idx]));
}




inline double XGen::calc_p_w_zero(const double* exp_ubeta_vals,
				  const double xi_i,
				  const int day_idx) {

    return exp(-xi_i * ubeta_exp_vals[day_idx]);
}




inline int XGen::sample_x_ijk(const double p_w_xIsOne,
			      const double p_xTodayIsOne,
			      const double p_xTom_xTodayIsZero,
			      const double p_xTom_xTodayIsOne) {

    double unnorm_prob_xIsZero, unnorm_prob_xIsOne, u;

    // when `X_ijk = 0`.  Note that `P(W_ijk = 0 | X_ijk = 0) = 1` so we are
    // ignoring that term in this expression.
    unnorm_prob_xIsZero = p_xTom_xTodayIsZero * (1 - p_xTodayIsOne);

    // when `X_ijk = 1`
    unnorm_prob_xIsOne = p_w_xIsOne * p_xTom_xTodayIsOne * p_xTodayIsOne;

    // sample a uniform r.v. and scale by sum of the unnormalized probabilities
    u = R::unif_rand() * (unnorm_prob_xIsZero + unnorm_prob_xIsOne);

    // in conjunction with the last statement, this has the effect of sampling
    // from `p(U * (p_0 + p_1) < p_0)`
    return (u < unnorm_prob_xIsZero) ? 0 : 1;
}




// intercourse is coded as -2 and -1 for days in which intercourse was observed.
// so that we can distinguish them from array indices which are nonnegative.
// Thus if the value is nonnegative then it was missing.

inline bool XGen::check_if_sex_miss(int miss_day) const {
    return (miss_day >= 0);
}
