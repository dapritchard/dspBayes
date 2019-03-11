#ifndef DSP_BAYES_SRC_W_GEN_H
#define DSP_BAYES_SRC_W_GEN_H

class XiGen;
#include "DayBlock.h"
#include "UProdBeta.h"
#include "XGen.h"


class WGen {

public:

    // storage for the `W_ijk` and `sum_k W_ijk`.  Note that we only provide
    // storage for cycles that a pregnancy occurs, because these are the only
    // cycles for which W is random.
    //
    // `m_vals` has storage for each day in the day-specific data that occurs
    // during a cycle in which a pregnancy occurs, and `m_sums` has storage for
    // each cycle in which a pregnancy occurs.
    int* m_vals;
    int* m_sums;

    // maps the r-th element of `m_vals` to the t-th index in the day-specific
    // data.  In other words, if `m_days_idx[r]` has a value of `t`, then
    // `m_vals[r]` is the value of `m_vals` for the `t`-th day.
    const Rcpp::IntegerVector& m_days_idx;

    // maps the r-th element of `m_sums` to the t-th index in the
    // subject-specific data.  In other words, if `m_subj_idx[r]` has a value of
    // `t`, then `m_sums[r]` is the value of `m_sums` for the `t`-th subject.
    const Rcpp::IntegerVector& m_subj_idx;

    // the elements of `m_preg_cyc` each map a pregnancy cycle to a block of
    // days in the day-specific data
    const PregCyc* m_preg_cyc;

    // the number of days for which intercourse occured during a cycle that
    // resulted in a pregnancy.  This value provides the amount of storage that
    // is associated with `m_vals` and `m_days_idx`.
    const int m_n_preg_days;

    // the number of cycles in the data in which a pregnancy occurred.  This
    // value provides the amount of storage that is associated with `m_sums`,
    // `m_subj_idx`, and `m_preg_cyc`.
    const int m_n_preg_cyc;

    // the number of days in the fertile window under the model
    int m_fw_len;


    WGen(Rcpp::List& preg_cyc,
         Rcpp::IntegerVector& w_to_days_idx,
         Rcpp::IntegerVector& w_cyc_to_subj_idx,
         int fw_len);
    ~WGen();

    void sample(XiGen& xi, UProdBeta& ubeta, XGen& X);
    const int* vals() const { return m_vals; }
    const int* sum_vals() const { return m_sums; }
    const int* days_idx() const { return m_days_idx.begin(); }
    const int* subj_idx() const { return m_subj_idx.begin(); }

    int n_preg_days() const { return m_n_preg_days; }
    int n_preg_cyc() const { return m_n_preg_cyc; }

    double rpois_zero_tr(double lambda);
};


#include "XiGen.h"

#endif
