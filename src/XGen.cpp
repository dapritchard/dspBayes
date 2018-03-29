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

#define IS_FW_ZERO_MISS  -99
#define SEX_SHIFT          2
#define NON_PREG_CYC      -1



// constructor
XGen::XGen(Rcpp::IntegerVector& X_rcpp,
	   Rcpp::IntegerVector& x_miss,
	   Rcpp::IntegerVector& utau_miss,
	   Rcpp::IntegerVector& sex_miss_to_w,
	   Rcpp::IntegerVector& sex_miss_to_xi,
	   int fw_len,
	   double cohort_sex_prob,
	   double sex_coef) :
    m_x_rcpp(X_rcpp),
    m_vals(X_rcpp.begin()),
    m_x_miss_rcpp(x_miss),
    m_x_miss(x_miss.begin()),
    m_x_miss_end(x_miss.end()),
    m_utau_miss_rcpp(utau_miss),
    m_utau_miss(utau_miss.begin()),
    m_sex_miss_to_w_rcpp(sex_miss_to_w),
    m_sex_miss_to_w(sex_miss_to_w.begin()),
    m_cyc_to_subj(sex_miss_to_xi.begin()),
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

    const int* x_cyc         = m_x_miss;
    const int* utau_cyc      = m_utau_miss;
    const int* miss_to_w_map = m_sex_miss_to_w;
    const int* subj_idx      = m_cyc_to_subj;

    // each iteration samples the missing intercourse values for the current
    // cycle
    while (x_cyc < m_x_miss_end) {

	sample_cycle(x_cyc, utau_cyc, &miss_to_w_map, *subj_idx, W, xi, ubeta, utau);

	x_cyc += m_ext_fw_len;
	utau_cyc += m_ext_fw_len;
	++subj_idx;
    }
}




void XGen::sample_cycle(const int* x_cyc,
			const int* utau_cyc,
			const int** miss_to_w_map,
			const int subj_idx,
			const WGen& W,
			const XiGen& xi,
			const UProdBeta& ubeta,
			const UProdTau& utau) {

    int sex_prev_day, sex_curr_day, sex_next_day;
    int x_curr_idx, x_next_idx;
    int utau_curr_idx, utau_next_idx;
    int curr_w;
    bool is_sex_miss, need_to_samp_bool;
    double p_w_xIsOne, p_xTodayIsOne, p_xTom_xTodayIsZero, p_xTom_xTodayIsOne;

    // value of xi for the subject that `miss_cycl` corresponds to
    const double xi_i = xi.vals()[subj_idx];
    // const double* ubeta_exp_vals = ubeta.exp_vals();
    // const double* utau_vals = utau.vals();
    const int* w_vals = W.vals();

    // obtain intercourse for the 0-th fertile window day
    if (x_cyc[0] == IS_FW_ZERO_MISS) {
	sex_prev_day = sample_day_before_fw_sex();
    }
    else {
	sex_prev_day = x_cyc[0] + SEX_SHIFT;
    }

    // for `x_curr_idx`, obtain index for missing or yes/no sex for nonmissing
    // intercourse day.  Similarly for `utau_curr_idx`, obtain index for needed
    // days and a dummy flag value for unneeded days.
    x_curr_idx = x_cyc[1];
    utau_curr_idx = utau_cyc[1];

    // obtain intercourse status for the 1-th fertile window day, and track
    // whether or not it is missing
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
	curr_w = (**miss_to_w_map == NON_PREG_CYC) ? 0 : w_vals[**miss_to_w_map];
	++*miss_to_w_map;
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
	    utau_next_idx = utau_cyc[r + 1];
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
	    p_w_xIsOne = calc_p_w_zero(ubeta, x_curr_idx, xi_i);

	    // obtain `P(X_ijk = 1 | X_{ij,k-1})`
	    p_xTodayIsOne = calc_p_xTodayIsOne(utau, utau_curr_idx, sex_prev_day);

	    // obtain `p(X_{ij,k+1} | X_ijk = 0)` and `p(X_{ij,k+1} | X_ijk = 1)`
	    //
	    // case: this is the last day of the fertile window, so we can have
	    // the effect of dropping the `p(X_{ij,k+1} | X_ijk = t)` terms by
	    // making them both equal to 1
	    if (r + 1 == m_ext_fw_len) {
		p_xTom_xTodayIsZero = p_xTom_xTodayIsOne = 1;
	    }
	    // case: not the last day of the fertile window, so have to
	    // calculate both terms
	    else {
		p_xTom_xTodayIsZero = calc_p_xTom(utau, utau_next_idx, sex_next_day, 0);
		p_xTom_xTodayIsOne  = calc_p_xTom(utau, utau_next_idx, sex_next_day, 1);
	    }

	    // sample missing intercourse for today
	    sex_curr_day = sample_x_ijk(p_w_xIsOne,
					p_xTodayIsOne,
					p_xTom_xTodayIsZero,
					p_xTom_xTodayIsOne);
	    m_vals[x_curr_idx] = sex_curr_day;
	}

	// slide the 3-day sequence over 1 in preparation for the next iteration
	sex_prev_day  = sex_curr_day;
	sex_curr_day  = sex_next_day;
	x_curr_idx    = x_next_idx;
	utau_curr_idx = utau_next_idx;
    }
}




inline double XGen::calc_p_w_zero(const UProdBeta& ubeta,
				  const int day_idx,
				  const double xi_i) {

    const double* ubeta_exp_vals = ubeta.exp_vals();
    return exp(-xi_i * ubeta_exp_vals[day_idx]);
}




double XGen::calc_p_xTom(const UProdTau& utau,
			 int utau_tomorrow_idx,
			 int x_tomorrow,
			 int x_today) const {

    const double* utau_vals = utau.vals();

    // calculate for `X_ij{k+1} = 1` and then adjust at the end for the 0 case
    // if necessary
    double p_xTomIsOne = x_today ?
	1 / (1 + exp(-utau_vals[utau_tomorrow_idx] - m_sex_coef)) :
	1 / (1 + exp(-utau_vals[utau_tomorrow_idx]));

    return x_tomorrow ?
	p_xTomIsOne :
	1 - p_xTomIsOne;
}




inline double XGen::calc_p_xTodayIsOne(const UProdTau& utau,
				       const int utau_today_idx,
				       const int sex_prev_day) const {

    const double* utau_vals = utau.vals();
    return sex_prev_day ?
	1 / (1 + exp(-utau_vals[utau_today_idx] - m_sex_coef)) :
	1 / (1 + exp(-utau_vals[utau_today_idx]));
}




int XGen::sample_x_ijk(double p_w_xIsOne,
		       double p_xTodayIsOne,
		       double p_xTom_xTodayIsZero,
		       double p_xTom_xTodayIsOne) {

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

inline bool XGen::check_if_sex_miss(int miss_day) {
    return (miss_day >= 0);
}
