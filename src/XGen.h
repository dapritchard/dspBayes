#ifndef DSP_BAYES_SRC_X_GEN_H
#define DSP_BAYES_SRC_X_GEN_H

#include "Rcpp.h"
#include "DayBlock.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XiGen.h"




class XGen {

public:

    // storage for the X values
    Rcpp::IntegerVector& m_x_rcpp;
    Rcpp::IntegerVector::iterator m_vals;

    // tracks whether intercourse occured on a given day if known, or the index
    // of the corresponding day in X if missing
    Rcpp::IntegerVector& m_x_miss_rcpp;
    Rcpp::IntegerVector::iterator m_x_miss;
    Rcpp::IntegerVector::iterator m_x_miss_end;

    // tracks the index in the `U * tau` data` for the a given day if it is in
    // the data, or a gives a signal flag otherwise
    Rcpp::IntegerVector& m_utau_miss_rcpp;
    Rcpp::IntegerVector::iterator m_utau_miss;

    // one entry for each missing value in intercourse, providing either the
    // index in W that the day corresponds to, or -1 if the day was part of a
    // non-pregnancy cycle
    Rcpp::IntegerVector& m_sex_miss_to_w_rcpp;
    Rcpp::IntegerVector::iterator  m_sex_miss_to_w;

    // one entry per cycle that has missing sex, with each element providing the
    // index in Xi that the cycle corresponds to.  Thus it has the same length
    // as `m_x_miss / m_ext_fw_len`.
    const int* m_cyc_to_subj;

    // number of days in the fertile window and extended fertile window
    // (i.e. includes the day before)
    const int m_ext_fw_len;

    // global probability used to sample missing values of intercourse for the
    // day before the fertile window
    const double m_cohort_sex_prob;

    // the regression coefficient in the missing intercourse prior probabilities
    // model coresponding to the previous day of intercourse
    const double m_sex_coef;


    XGen(Rcpp::IntegerVector& X_rcpp,
	 Rcpp::IntegerVector& x_miss,
	 Rcpp::IntegerVector& utau_miss,
	 Rcpp::IntegerVector& sex_miss_to_w,
	 Rcpp::IntegerVector& sex_miss_to_xi,
	 int fw_len,
	 double cohort_sex_prob,
	 double sex_coef);

    ~XGen() {}

    void sample(const WGen& W,
		const XiGen& xi,
		const UProdBeta& ubeta,
		const UProdTau& utau);

    void sample_cycle(const int* x_cyc,
		      const int* utau_cyc,
		      const int** miss_to_w_map,
		      const int subj_idx,
		      const WGen& W,
		      const XiGen& xi,
		      const UProdBeta& ubeta,
		      const UProdTau& utau);

    static double calc_p_w_zero(const UProdBeta& ubeta,
				const int day_idx,
				const double xi_i);

    double calc_p_xTom(const UProdTau& utau,
		       int utau_tomorrow_idx,
		       int x_tomorrow,
		       int x_today) const;

    double calc_p_xTodayIsOne(const UProdTau& utau,
			      const int utau_today_idx,
			      const int sex_prev_day) const;

    static int sample_x_ijk(double p_w_xIsOne,
			    double p_xTodayIsOne,
			    double p_xTom_xTodayIsZero,
			    double p_xTom_xTodayIsOne);

    int sample_day_before_fw_sex() const;
    static bool check_if_sex_miss(int miss_day);

    int* vals() { return m_vals; }
    const int* vals() const { return m_vals; }
    // const XMissDay* miss_day() const { return m_miss_day; }
    int n_days() const { return m_x_rcpp.size(); }
    // double sex_coef() const { return m_sex_coef; }

};


#endif
