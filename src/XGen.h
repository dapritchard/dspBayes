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

    class SexMissCycInfo;

    // storage for the X values
    Rcpp::IntegerVector& m_x_rcpp;
    Rcpp::IntegerVector::iterator m_vals;

    // tracks whether intercourse occured on a given day if known, or the index
    // of the corresponding day in X if missing
    Rcpp::IntegerVector& m_x_miss_rcpp;
    Rcpp::IntegerVector::iterator m_x_miss;
    Rcpp::IntegerVector::iterator m_x_miss_end;

    SexMissCycInfo* m_miss_info;
    SexMissCycInfo* m_miss_info_end;

    // // tracks the index in the `U * tau` data` for the a given day if it is in
    // // the data, or a gives a signal flag otherwise
    // Rcpp::IntegerVector& m_utau_miss_rcpp;
    // Rcpp::IntegerVector::iterator m_utau_miss;

    // // one entry for each missing value in intercourse, providing either the
    // // index in W that the day corresponds to, or -1 if the day was part of a
    // // non-pregnancy cycle
    // Rcpp::IntegerVector& m_sex_miss_to_w_rcpp;
    // Rcpp::IntegerVector::iterator  m_sex_miss_to_w;

    // one entry per cycle that has missing sex, with each element providing the
    // index in Xi that the cycle corresponds to.  Thus it has the same length
    // as `m_x_miss / m_ext_fw_len`.
    const int* m_cyc_to_subj;

    // number of days in the fertile window
    const int m_fw_len;

    // // number of cycles with missing intercourse
    // const int m_n_miss_cyc;

    // global probability used to sample missing values of intercourse for the
    // day before the fertile window
    const double m_cohort_sex_prob;

    // the regression coefficient in the missing intercourse prior probabilities
    // model coresponding to the previous day of intercourse
    const double m_sex_coef;


    XGen(Rcpp::IntegerVector& X_rcpp,
	 Rcpp::IntegerVector& x_miss,
	 Rcpp::List& sex_miss_info,
	 int fw_len,
	 double cohort_sex_prob,
	 double sex_coef);

    ~XGen();

    SexMissCycInfo* list_to_arr(Rcpp::List& sex_info_list);

    void sample(const WGen& W,
		const XiGen& xi,
		const UProdBeta& ubeta,
		const UProdTau& utau);

    void sample_cycle(const SexMissCycInfo* curr_miss_info,
		      const WGen& W,
		      const XiGen& xi,
		      const UProdBeta& ubeta,
		      const UProdTau& utau);

    static double calc_p_w_zero(double ubeta_exp_val,
				double xi_i);

    double calc_p_xTom(double utau_val_tom,
		       int x_tomorrow,
		       int x_today) const;

    double calc_p_xTodayIsOne(double utau_val_today,
			      int sex_prev_day) const;

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




class XGen::SexMissCycInfo {

public:

    int x_cyc;
    int w_cyc;
    int xi_idx;
    int sex_prev;

    SexMissCycInfo() : x_cyc(0), w_cyc(0), xi_idx(0), sex_prev(0) {}

    SexMissCycInfo(int x_cyc, int w_cyc, int xi_idx, int sex_prev) :
	x_cyc(x_cyc),
	w_cyc(w_cyc),
	xi_idx(xi_idx),
	sex_prev(sex_prev) {
    }

    static SexMissCycInfo* list_to_arr(Rcpp::List& block_list);
};


#endif
