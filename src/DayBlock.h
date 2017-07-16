#ifndef DSP_BAYES_SRC_DAY_BLOCK_H
#define DSP_BAYES_SRC_DAY_BLOCK_H


struct DayBlock {

    const int beg_idx;
    const int n_days;

    DayBlock() : beg_idx(0), n_days(0) {}

    DayBlock(int beg_idx, int n_days) :
	beg_idx(beg_idx),
	n_days(n_days) {
    }

    static DayBlock* list_to_arr(Rcpp::List block_list);
};


struct PregCyc : public DayBlock {

    const int subj_idx;

    PregCyc(int beg_idx, int n_days, int subj_idx) :
	DayBlock(beg_idx, n_days),
	subj_idx(subj_idx) {
    }

    static PregCyc* list_to_arr(Rcpp::List block_list);
};


#endif
