#ifndef DSP_BAYES_SRC_U_GEN_H
#define DSP_BAYES_SRC_U_GEN_H

#include "Rcpp.h"

#include "CoefGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"


class UGen {

public:

    struct UMissBlock;

    const int m_col_start;
    const int m_col_end;
    const int m_ref_col;
    const int m_n_categs;

    const int m_max_alt_exp_ubeta_size;
    const int m_max_alt_utau_size;

    const double* m_u_prior_probs;

    UMissBlock* m_miss_block;
    const UMissBlock* const m_end_block;

    const int* m_w_idx;
    const int* m_x_idx;

    UGen(Rcpp::IntegerVector& var_info,
	 Rcpp::NumericVector& u_prior_probs,
	 Rcpp::List& var_block_list,
	 Rcpp::IntegerVector& w_idx,
	 Rcpp::IntegerVector& u_idx);
    ~UGen();

    void sample(const WGen& W,
		const XiGen& xi,
		const CoefGen& coefs,
		const XGen& X,
		UProdBeta& ubeta,
		UProdTau& utau);

    void calc_posterior_w(double* posterior_w_probs,
			  double* alt_exp_ubeta_vals,
			  const WGen& W,
			  const XiGen& xi,
			  const CoefGen& coefs,
			  const UProdBeta& ubeta,
			  const UMissBlock* const miss_block) const;

    void calc_posterior_x(double* posterior_x_probs,
			  double* alt_utau_vals,
			  const XGen& X,
			  const UProdTau& utau,
			  const UMissBlock* const miss_block) const;

    int sample_covariate(const double* posterior_w_probs,
			 const double* posterior_x_probs) const;

    static void update_ubeta(UProdBeta& ubeta,
			     const int u_categ,
			     const double* alt_exp_ubeta_vals,
			     const UMissBlock* const miss_block);

    static void update_utau(UProdTau& utau,
			    const int u_categ,
			    const double* alt_utau_vals,
			    const UMissBlock* const miss_block);

};




struct UGen::UMissBlock {

    int beg_day_idx;
    int n_days;

    int beg_w_idx;

    int beg_sex_idx;
    int n_sex_days;

    int u_col;
    int subj_idx;

    UMissBlock() :
	beg_day_idx(0),
	n_days(0),
	beg_w_idx(0),
	beg_sex_idx(0),
	n_sex_days(0),
	u_col(0),
	subj_idx(0) {
    }

    UMissBlock(int beg_day_idx,
	       int n_days,
	       int beg_w_idx,
	       int beg_sex_idx,
	       int n_sex_days,
	       int u_col,
	       int subj_idx) :
	beg_day_idx(beg_day_idx),
	n_days(n_days),
	beg_w_idx(beg_w_idx),
	beg_sex_idx(beg_sex_idx),
	n_sex_days(n_sex_days),
	u_col(u_col),
	subj_idx(subj_idx) {
    }

    static UMissBlock* list_to_arr(Rcpp::List& block_list);
};


#endif
