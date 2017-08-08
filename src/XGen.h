#ifndef DSP_BAYES_SRC_X_GEN_H
#define DSP_BAYES_SRC_X_GEN_H

#include "Rcpp.h"
#include "DayBlock.h"


class XGen {

public:

    const Rcpp::IntegerVector& m_x_rcpp;
    Rcpp::IntegerVector::iterator m_vals;

    const XMiss* m_miss_cyc;
    const int m_n_miss_cyc;

    const int* m_miss_day;

    const int m_n_max_miss;
    const int m_n_max_miss_pow2;


    class XMissDay {

	int idx;
	int prev;
    };
};


#endif
