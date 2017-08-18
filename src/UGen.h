#ifndef DSP_BAYES_SRC_U_GEN_H
#define DSP_BAYES_SRC_U_GEN_H

#include "Rcpp.h"


class UGen {

    class UMissBlock;
    class UMissDay;

    const int m_col_start;
    const int col_end;
    /* const int m_n_cols; */
    const int m_n_categs;

    const bool is_ref_cell_coding;

    UMissBlock* m_miss_block;
    UMissBlock* m_curr_block;
    const UMissBlock* const m_end_block;
    /* const int m_n_miss_block; */
    /* int m_block_idx; */

    const UMissDay* miss_day;
    const int m_n_miss_day;
    int m_day_idx;


};




class UMissBlock {

    int beg_idx;
    int n_days;
    int u_col;
    int subj_idx;
};




class UMissDay {

    int preg_day_idx;
    int day_idx;
};


#endif
