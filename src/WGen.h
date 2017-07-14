#ifndef DSP_BAYES_W_GEN_H_
#define DSP_BAYES_W_GEN_H_

class WGen {

public:

    // storage for the W_ijk and sum_k W_ijk.  Note that we only provide storage
    // for cycles that a pregnancy occurs, because these are the only cycles for
    // which W is random.
    int* m_W;
    int* m_w_sum;

    const int* m_preg_cyc_idx;

    const DayBlock* m_preg_cyc;
    const DayBlock* m_preg_end;

    m_n_preg_cyc;

    double* m_pois_mean;

    WGen(Rcpp::List& preg_cyc_list);

    int* preg_cyc_idx() { return m_preg_cyc_idx; }
    int* w_sum() { return m_w_sum; }
};

#endef
