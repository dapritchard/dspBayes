#ifndef DSP_BAYES_SRC_W_GEN_H
#define DSP_BAYES_SRC_W_GEN_H

class XiGen;
#include "DayBlock.h"
#include "UProdBeta.h"


class WGen {

public:

    // storage for the `W_ijk` and `sum_k W_ijk`.  Note that we only provide
    // storage for cycles that a pregnancy occurs, because these are the only
    // cycles for which W is random.
    //
    // `m_w_vals` has storage for each day in the day-specific data that occurs
    // during a cycle in which a pregnancy occurs, and `m_w_sums` has storage
    // for each cycle in which a pregnancy occurs.
    int* m_w_vals;
    int* m_w_sums;

    // maps the r-th element of `m_w_vals` to the t-th index in the day-specific
    // data.  In other words, if `m_w_days_idx[r]` has a value of `t`, then
    // `m_w_vals[r]` is the value of `m_w_vals` for the `t`-th day.
    const int* m_w_days_idx;

    // maps the r-th element of `m_w_sums` to the t-th index in the
    // cycle-specific data.  In other words, if `m_w_cyc_idx[r]` has a value of
    // `t`, then `m_w_sums[r]` is the value of `m_w_sums` for the `t`-th cycle.
    const int* m_w_cyc_idx;

    // the elements of `m_preg_cyc` each map a cycle to a block of days in the
    // day-specific data
    const PregCyc* m_preg_cyc;

    // the number of days for which intercourse occured during a cycle that
    // resulted in a pregnancy.  This value provides the amount of storage that
    // is associated with `m_w_vals`.
    const int m_n_days;

    // the number of cycles in the data in which a pregnancy occurred.  This
    // value provides the amount of storage that is associated with `m_w_sums`,
    // `m_w_cyc_idx`, and `m_preg_cyc`.
    const int m_n_preg_cyc;

    // some scratch storage that we use to place multinomial probabilities into
    double* m_mult_probs;


    WGen(Rcpp::List& preg_cyc,
	 Rcpp::IntegerVector& w_days_idx,
	 Rcpp::IntegerVector& w_cyc_idx,
	 int fw_len);
    ~WGen();

    void sample(XiGen& xi, UProdBeta& u_prod_beta);
    const int* vals() const { return m_w_vals; }
    const int* sum_vals() const { return m_w_sums; }
    const int* days_idx() const { return m_w_days_idx; }
    const int* cyc_idx() const { return m_w_cyc_idx; }
    int n_days() const { return m_n_days; }
    int n_preg_cyc() const { return m_n_preg_cyc; }

};


#include "XiGen.h"

#endif
