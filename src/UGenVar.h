#ifndef DSP_BAYES_SRC_U_GEN_VAR_H
#define DSP_BAYES_SRC_U_GEN_VAR_H

#include "Rcpp.h"

#include "CoefGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"




class UGenVar {

public:

    struct UMissBlock;

    double* const m_u_var_col;
    const int m_n_days;

    const int* const m_w_idx;
    const int* const m_x_idx;

    bool m_record_status;
    Rcpp::NumericVector m_vals_rcpp;

    UGenVar();
    virtual ~UGenVar() {}

    UGenVar(Rcpp::NumericMatrix& u_rcpp,
            Rcpp::IntegerVector& preg_map,
            Rcpp::IntegerVector& sex_map,
            int u_col,
            bool record_status);

    virtual void sample(const WGen& W,
                        const XiGen& xi,
                        const CoefGen& coefs,
                        const XGen& X,
                        UProdBeta& ubeta,
                        UProdTau& utau) = 0;
};




class UGenVarCateg : public UGenVar {

public:

    struct UMissBlockCateg;

    const int m_col_start;
    const int m_col_end;
    const int m_ref_col;
    const int m_n_categs;

    const int m_max_n_days_miss;
    const int m_max_n_sex_days_miss;

    const double* m_log_u_prior_probs;

    UMissBlockCateg* m_miss_block;
    const UMissBlockCateg* const m_end_block;

    // UGenVarCateg(Rcpp::List& var_info_list);

    UGenVarCateg(Rcpp::NumericMatrix& u_rcpp,
                 Rcpp::IntegerVector& var_info,
                 Rcpp::NumericVector& log_u_prior_probs,
                 Rcpp::List& var_block_list,
                 Rcpp::IntegerVector& preg_map,
                 Rcpp::IntegerVector& sex_map,
                 bool record_status);
    ~UGenVarCateg();

    void sample(const WGen& W,
                const XiGen& xi,
                const CoefGen& coefs,
                const XGen& X,
                UProdBeta& ubeta,
                UProdTau& utau);

    void calc_log_condit_w(double* posterior_w_probs,
                           double* alt_exp_ubeta_vals,
                           const WGen& W,
                           const XiGen& xi,
                           const CoefGen& coefs,
                           const XGen& X,
                           const UProdBeta& ubeta,
                           const UMissBlockCateg* const miss_block) const;

    void calc_log_condit_x(double* posterior_x_probs,
                           double* alt_utau_vals,
                           const XGen& X,
                           const UProdTau& utau,
                           const UMissBlockCateg* const miss_block) const;

    int sample_covariate(const double* posterior_w_probs,
                         const double* posterior_x_probs) const;

    void update_u(const UMissBlockCateg* const miss_block);

    static void update_ubeta(UProdBeta& ubeta,
                             const int u_categ,
                             const double* alt_exp_ubeta_vals,
                             const UMissBlockCateg* const miss_block);

    void update_utau(UProdTau& utau,
                     const int u_categ,
                     const double* alt_utau_vals,
                     const UMissBlockCateg* const miss_block) const;
};




struct UGenVar::UMissBlock {

    int beg_day_idx;
    int n_days;

    int beg_w_idx;

    int beg_sex_idx;
    int n_sex_days;

    int subj_idx;

    UMissBlock() :
        beg_day_idx(0),
        n_days(0),
        beg_w_idx(0),
        beg_sex_idx(0),
        n_sex_days(0),
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
        subj_idx(subj_idx) {
    }
};




struct UGenVarCateg::UMissBlockCateg : public UGenVar::UMissBlock {

    int u_col;

    UMissBlockCateg() :
        UMissBlock(),
        u_col(0) {
    }

    UMissBlockCateg(int beg_day_idx,
                    int n_days,
                    int beg_w_idx,
                    int beg_sex_idx,
                    int n_sex_days,
                    int u_col,
                    int subj_idx) :
        UMissBlock(beg_day_idx, n_days, beg_w_idx, beg_sex_idx, n_sex_days, u_col, subj_idx),
        u_col(u_col) {
    }

    static UMissBlockCateg* list_to_arr(Rcpp::List& block_list);
};


#endif
