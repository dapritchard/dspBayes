#ifndef DSP_BAYES_W_GEN_H_
#define DSP_BAYES_W_GEN_H_

class WGen {

public:

    // storage for the W_ijk and sum_k W_ijk.  Note that we only provide storage
    // for cycles that a pregnancy occurs, because these are the only cycles for
    // which W is random.
    int* m_W;
    int* m_w_sum;

    const DayBlock* m_preg_cyc;
    const DayBlock* m_preg_end;

    m_n_preg_cyc;

    double* m_pois_mean;

    WGen(Rcpp::List& preg_cyc_list);


};

#endef
