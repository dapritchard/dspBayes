#include <cmath>
#include "Rcpp.h"
#include "DayBlock.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"


#define IS_FW_ZERO_MISS  -99
#define NON_PREG_CYC      -1



// constructor
XGen::XGen(Rcpp::IntegerVector& X_rcpp,
	   Rcpp::IntegerVector& x_miss,
	   Rcpp::List& sex_miss_info,
	   // Rcpp::IntegerVector& utau_miss,
	   // Rcpp::IntegerVector& sex_miss_to_w,
	   // Rcpp::IntegerVector& sex_miss_to_xi,
	   int fw_len,
	   double cohort_sex_prob,
	   double sex_coef) :
    m_x_rcpp(X_rcpp),
    m_vals(X_rcpp.begin()),
    m_x_miss_rcpp(x_miss),
    m_x_miss(x_miss.begin()),
    m_x_miss_end(x_miss.end()),
    m_miss_info(list_to_arr(sex_miss_info)),
    m_miss_info_end(m_miss_info + ((int) sex_miss_info["n_cyc"])),
    // m_utau_miss_rcpp(utau_miss),
    // m_utau_miss(utau_miss.begin()),
    // m_sex_miss_to_w_rcpp(sex_miss_to_w),
    // m_sex_miss_to_w(sex_miss_to_w.begin()),
    // m_cyc_to_subj(sex_miss_to_xi.begin()),
    m_fw_len(fw_len),
    m_cohort_sex_prob(cohort_sex_prob),
    m_sex_coef(sex_coef) {
}




XGen::~XGen() {
    delete[] m_miss_info;
}




// construct an array of `DayBlock`s based upon an `Rcpp::List` that provides
// the specifications for each block.  In more detail, an array of `DayBlock`s
// is created with length given by the length of `block_list`.  Furthermore, it
// is assumed that each element in `block_list` stores integer values that can
// be accessed using names `beg_idx` and `n_days` and which are used to
// initialize the struct's member values of the same name.  The return value is
// a pointer to the beginning of the array.

XGen::SexMissCycInfo* XGen::list_to_arr(Rcpp::List& sex_info_list) {

    // save one of the Rcpp vectors so that we can ask for its length later
    Rcpp::IntegerVector x_cyc_rcpp = Rcpp::as<Rcpp::IntegerVector>(sex_info_list["x_cyc"]);

    // get a pointer to the data for each of the vectors
    int* x_cyc    = x_cyc_rcpp.begin();
    int* w_cyc    = (Rcpp::as<Rcpp::IntegerVector>(sex_info_list["w_cyc"])).begin();
    int* xi_idx   = (Rcpp::as<Rcpp::IntegerVector>(sex_info_list["xi_idx"])).begin();
    int* sex_prev = (Rcpp::as<Rcpp::IntegerVector>(sex_info_list["sex_prev"])).begin();

    // allocate memory for our structs
    SexMissCycInfo* block_arr = new SexMissCycInfo[sex_info_list.size()];

    // each iteration constructs a new struct based upon the t-th elements of
    // the various vectors
    for (int t = 0; t < x_cyc_rcpp.size(); ++t) {

	block_arr[t] = SexMissCycInfo(x_cyc[t], w_cyc[t], xi_idx[t], sex_prev[t]);
    }

    return block_arr;
}









// update the missing values of X
void XGen::sample(const WGen& W,
		  const XiGen& xi,
		  const UProdBeta& ubeta,
		  const UProdTau& utau) {


    // each iteration samples the missing intercourse values for the current
    // cycle
    for (const SexMissCycInfo* curr_miss_info = m_miss_info;
	 curr_miss_info < m_miss_info_end;
	 ++curr_miss_info) {

	sample_cycle(curr_miss_info, W, xi, ubeta, utau);
    }
}




void XGen::sample_cycle(const SexMissCycInfo* curr_miss_info,
			const WGen& W,
			const XiGen& xi,
			const UProdBeta& ubeta,
			const UProdTau& utau) {

    int sex_prev_day, sex_curr_day, sex_next_day, utau_curr_val, utau_next_val;
    double p_w_xIsOne, p_xTodayIsOne, p_xTom_xTodayIsZero, p_xTom_xTodayIsOne;

    // point to the first day of the cycle in the (i) X data, (ii) the missing
    // intercourse status data, and (iii) the `U * tau` data
    int* x_cyc                  = m_vals + curr_miss_info->x_cyc;
    const int* x_miss_cyc       = m_x_miss + curr_miss_info->x_cyc;
    const double* ubeta_exp_cyc = ubeta.exp_vals() + curr_miss_info->x_cyc;
    const double* utau_cyc      = utau.vals() + curr_miss_info->x_cyc;
    const int *w_cyc;

    // record whether pregnancy occurred this cycle, and if it did then set
    // `w_cyc` to point to the first day of the cycle in the W data
    bool is_preg_cyc = (curr_miss_info->w_cyc != NON_PREG_CYC);
    if (is_preg_cyc) {
	w_cyc = W.vals() + curr_miss_info->w_cyc;
    }

    // the value in xi for the subject that the cycle corresponds to
    const double xi_i = xi.vals()[curr_miss_info->xi_idx];

    // obtain intercourse for the 0-th and the 1-th fertile window day
    sex_prev_day = curr_miss_info->sex_prev;
    if (sex_prev_day == IS_FW_ZERO_MISS) {
    	sex_prev_day = sample_day_before_fw_sex();
    }
    sex_curr_day = x_cyc[0];
    // obtain `U_ijk^T tau` for the 1-th fertile window day
    utau_curr_val = utau_cyc[0];


    for (int r = 0; r < m_fw_len; ++r) {

	// obtain intercourse and `U_ijk^T tau` for the (r + 1)-th fertile
	// window day
	if (r + 1 < m_fw_len) {
	    sex_next_day  = x_cyc[r + 1];
	    utau_next_val = utau_cyc[r + 1];
	}

	// the next if / else if / else block is in charge of conditionally
	// updating X_ijk, if needed

	// case: sex not missing today, so do nothing
	if (! x_miss_cyc[r]) {
	    // noop
	}

	// case: sex was missing and `W_ijk > 0`, so consequently X_ijk is
	// required to have a value of 1
	else if (is_preg_cyc && (w_cyc[r] > 0)) {
	    x_cyc[r] = 1;
	}

	// case: intercourse status was missing today and `W_ijk = 0`, so need
	// to sample X_ijk
	else {

	    // obtain `P(W = 0 | X_ijk = 1)`.  Note that we don't need to
	    // calculate this when `X_ijk = 0`, since then the expression is
	    // simply 1.
	    p_w_xIsOne = calc_p_w_zero(ubeta_exp_cyc[r], xi_i);

	    // obtain `P(X_ijk = 1 | X_{ij,k-1})`
	    p_xTodayIsOne = calc_p_xTodayIsOne(utau_curr_val, sex_prev_day);

	    // obtain `p(X_{ij,k+1} | X_ijk = 0)` and `p(X_{ij,k+1} | X_ijk = 1)`
	    //
	    // case: this is the last day of the fertile window, so we can have
	    // the effect of dropping the `p(X_{ij,k+1} | X_ijk = t)` terms by
	    // making them both equal to 1
	    if (r + 1 == m_fw_len) {
		p_xTom_xTodayIsZero = p_xTom_xTodayIsOne = 1;
	    }
	    // case: not the last day of the fertile window, so have to
	    // calculate both terms
	    else {
		p_xTom_xTodayIsZero = calc_p_xTom(utau_next_val, sex_next_day, 0);
		p_xTom_xTodayIsOne  = calc_p_xTom(utau_next_val, sex_next_day, 1);
	    }

	    // sample missing intercourse for today
	    sex_curr_day = sample_x_ijk(p_w_xIsOne,
					p_xTodayIsOne,
					p_xTom_xTodayIsZero,
					p_xTom_xTodayIsOne);
	    x_cyc[r] = sex_curr_day;
	}

	// slide the 3-day sequence over 1 in preparation for the next iteration
	sex_prev_day  = sex_curr_day;
	sex_curr_day  = sex_next_day;
	utau_curr_val = utau_next_val;
    }
}




inline double XGen::calc_p_w_zero(double ubeta_exp_val,
				  double xi_i) {
    return exp(-xi_i * ubeta_exp_val);
}




inline double XGen::calc_p_xTom(double utau_val_tom,
				int x_tomorrow,
				int x_today) const {

    // calculate for `X_ij{k+1} = 1` and then adjust at the end for the 0 case
    // if necessary
    double p_xTomIsOne = x_today ?
	1 / (1 + exp(-utau_val_tom - m_sex_coef)) :
	1 / (1 + exp(-utau_val_tom));

    return x_tomorrow ?
	p_xTomIsOne :
	1 - p_xTomIsOne;
}




inline double XGen::calc_p_xTodayIsOne(double utau_val_today,
				       int sex_prev_day) const {

    return sex_prev_day ?
	1 / (1 + exp(-utau_val_today - m_sex_coef)) :
	1 / (1 + exp(-utau_val_today));
}




int XGen::sample_x_ijk(double p_w_xIsOne,
		       double p_xTodayIsOne,
		       double p_xTom_xTodayIsZero,
		       double p_xTom_xTodayIsOne) {

    double unnorm_prob_xIsZero, unnorm_prob_xIsOne, u;

    // the term for when `X_ijk = 0`.  Note that `P(W_ijk = 0 | X_ijk = 0) = 1`
    // so we are ignoring that term in this expression.
    unnorm_prob_xIsZero = p_xTom_xTodayIsZero * (1 - p_xTodayIsOne);

    // the term for when `X_ijk = 1`
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




inline int XGen::sample_day_before_fw_sex() const {
    return (R::unif_rand() < m_cohort_sex_prob) ? 1 : 0;
}
