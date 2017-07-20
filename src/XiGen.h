#ifndef DSP_BAYES_SRC_XI_GEN_H
#define DSP_BAYES_SRC_XI_GEN_H

#include "Rcpp.h"
class WGen;
class PhiGen;
#include "DayBlock.h"
#include "UProdBeta.h"


class XiGen {

public:

    // storage for the `xi_i` values
    double* m_xi_vals;
    double* const m_output_start;
    double* const m_output_end;

    // the elements of `m_subj` each map an individual to a block of days from
    // the day-specific data
    const DayBlock* m_subj;

    // the number of subjects in the data.  This is the amount of data that is
    // sampled during one MCMC scan.
    const int m_n_subj;

    // tracks whether we are past the burn-in phase
    bool m_record_status;

    XiGen(Rcpp::List subj_days, int n_samp);
    ~XiGen();

    void sample(const WGen& W, const PhiGen& phi, const UProdBeta& u_prod_beta);

    const double* vals() const { return m_xi_vals; }
    const int n_subj() const { return m_n_subj; }

    double* output_start() const { return m_output_start; }
    double* output_end() const { return m_output_end; }
};


#include "WGen.h"
#include "PhiGen.h"

#endif
