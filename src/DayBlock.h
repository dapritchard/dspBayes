#ifndef DSP_BAYES_SRC_DAY_BLOCK_H
#define DSP_BAYES_SRC_DAY_BLOCK_H

#include "Rcpp.h"


struct DayBlock {

    int beg_idx;
    int n_days;

    DayBlock() : beg_idx(0), n_days(0) {}

    DayBlock(int beg_idx, int n_days) :
	beg_idx(beg_idx),
	n_days(n_days) {
    }

    static DayBlock* list_to_arr(Rcpp::List& block_list);
};




struct PregCyc : public DayBlock {

    int subj_idx;

    PregCyc() : DayBlock(), subj_idx(0) {}

    PregCyc(int beg_idx, int n_days, int subj_idx) :
	DayBlock(beg_idx, n_days),
	subj_idx(subj_idx) {
    }

    static PregCyc* list_to_arr(Rcpp::List& block_list);
};




struct MissCyc : public PregCyc {

    int n_miss;
    int preg;

    MissCyc() : PregCyc(), n_miss(0) {}

    MissCyc(int beg_idx, int n_days, int subj_idx, int n_miss, int preg) :
	PregCyc(beg_idx, n_days, subj_idx),
	n_miss(n_miss),
	preg(preg) {
    }

    static MissCyc* list_to_arr(Rcpp::List& block_list);
};


#endif
