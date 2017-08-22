#ifndef DSP_BAYES_SRC_U_GEN_H
#define DSP_BAYES_SRC_U_GEN_H

#include "Rcpp.h"


class UGen {

public:

    struct UMissBlock;
    struct UMissDay;

    const int m_col_start;
    const int m_col_end;
    const int m_ref_col;
    /* const int m_n_cols; */
    // const int m_n_categs;

    // const bool m_is_ref_cell_coding;

    UMissBlock* m_miss_block;
    /* UMissBlock* m_curr_block; */
    /* const UMissBlock* const m_end_block; */
    /* const int m_n_miss_block; */
    /* int m_block_idx; */

    const UMissDay* m_miss_day;
    /* const int m_n_miss_day; */
    /* int m_day_idx; */

    int* m_w_idx;
    int* m_x_idx;

    UGen();

    // void calc_posterior_w_probs(double* posterior_w_probs,
    // 				const WGen& W,
    // 				const XiGen& xi,
    // 				const CoefGen& coefs,
    // 				const UProdBeta& ubeta,
    // 				const UMissBlock* const miss_block) const;

    // void calc_ubeta_possibs(double* alt_exp_ubeta_vals,
    // 			    const double** categ_exp_ubeta_arrs,
    // 			    const CoefGen& coefs,
    // 			    const UProdBeta& ubeta,
    // 			    const UMissBlock* const miss_block) const;

    void calc_posterior_w(double* posterior_w_probs,
			  const WGen& W,
			  const XiGen& xi,
			  const CoefGen& coefs,
			  const UProdBeta& ubeta,
			  const UMissBlock* const miss_block) const;

    void calc_posterior_x(double* posterior_x_probs,
			  const XGen& X,
			  const UProdTau& utau,
			  const UMissBlock* const miss_block) const;
};




struct UGen::UMissBlock {

    int beg_day_idx;
    int n_days;

    int beg_w_idx;

    int beg_sex_idx;
    int n_sex_days;

    int u_col;
    int subj_idx;
};




// struct UGen::UMissDay {

//     int day_idx;
//     int w_day_idx;
//     // int x_day_idx;
// };


#endif
