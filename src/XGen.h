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

    // one entry for each missing value in intercourse, providing either the
    // index in W that the day corresponds to, or -1 if the day was part of a
    // non-pregnancy cycle
    Rcpp::IntegerVector& m_sex_miss_to_w_rcpp;
    Rcpp::IntegerVector::iterator  m_sex_miss_to_w;

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
	 Rcpp::List& miss_cyc,
	 Rcpp::List& miss_day,
	 double cohort_sex_prob,
	 double sex_coef);
    ~XGen() {}

    void sample(const WGen& W,
		const XiGen& xi,
		const UProdBeta& ubeta,
		const UProdTau& utau);

    void sample_cycle(const XMissCyc* miss_cyc,
		      const int* W,
		      const XiGen& xi,
		      const UProdBeta& ubeta,
		      const UProdTau& utau);

    double calc_prior_prob(const UProdTau& utau,
			   const int miss_day_idx,
			   const int prev_day_sex) const;

    static double calc_posterior_prob(const UProdBeta& ubeta,
				      const double xi_i,
				      const int day_idx);

    static int sample_x_ijk(const double prior_prob_yes,
			    const double posterior_prob_yes);

    int sample_day_before_fw_sex() const;
    bool check_if_prev_sex_miss(int miss_day_idx) const;

    int* vals() { return m_vals; }
    const int* vals() const { return m_vals; }
    const XMissDay* miss_day() const { return m_miss_day; }
    int n_days() const { return m_x_rcpp.size(); }
    double sex_coef() const { return m_sex_coef; }

};


#endif
