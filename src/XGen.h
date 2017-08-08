#ifndef DSP_BAYES_SRC_X_GEN_H
#define DSP_BAYES_SRC_X_GEN_H

#include "Rcpp.h"
#include "DayBlock.h"
#include "UProdBeta.h"
#include "UProdTau.h"


class XGen {

public:

    class XMissDay;

    // storage for the X values
    Rcpp::IntegerVector& m_x_rcpp;
    Rcpp::IntegerVector::iterator m_vals;

    // information about the number of X missing for a given cycle
    const MissCyc* m_miss_cyc;
    const int m_n_miss_cyc;

    // tracks the indices of the missing values in X as well as whether
    // intercourse occured on the previous day
    const XMissDay* m_miss_day;

    // the number of missing in the cycle with the largest number of missing,
    // and 2 to the power of this number
    const int m_n_max_miss;
    const int m_n_max_miss_pow2;

    XGen(Rcpp::IntegerVector& X_rcpp, Rcpp::List& miss_cyc, Rcpp::List& miss_day, int n_max_miss);

    void sample(UProdBeta ubeta, UProdTau utau);

    int calc_prior_probs(double* prior_probs,
			 MissCyc* curr_miss_cyc,
			 XMissDay* curr_miss_day) const;

    static int calc_posterior_probs(double* posterior_probs,
				    const MissCyc* curr_miss_cyc,
				    const XMissDay* curr_miss_day,
				    const double* prior_probs,
				    const UProdBeta& ubeta);

    static int sample_x_perm(double* probs, int n_perms);

    static void XGen::update_cyc_x(MissCyc* curr_miss_cyc, XMissDay* curr_miss_day, int t);

    static double XGen::calc_nonrandom_sum_exp_ubeta(MissCyc* curr_miss_cyc,
						     XMissDay* curr_miss_day,
						     UProdBeta& ubeta);

    const int* vals() const { m_vals; }
};




class XGen::XMissDay {

    int idx;
    int prev;

    XMissDay() : idx(0), prev(0) {}

    XMissDay(int idx, int prev) :
	idx(idx),
	prev(prev) {
    }

    static XMissDay* list_to_arr(Rcpp::List& block_list);
};


#endif
