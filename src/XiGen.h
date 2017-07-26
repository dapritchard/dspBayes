#ifndef DSP_BAYES_SRC_XI_GEN_H
#define DSP_BAYES_SRC_XI_GEN_H

#include "Rcpp.h"
class WGen;
class PhiGen;
#include "DayBlock.h"
#include "UProdBeta.h"

extern bool g_record_status;


class XiGen {

public:

    // storage for the `xi_i` values
    Rcpp::NumericMatrix m_vals_rcpp;
    Rcpp::NumericMatrix::iterator m_vals;

    // the elements of `m_subj` each map an individual to a block of days from
    // the day-specific data
    const DayBlock* m_subj;

    // the number of subjects in the data.  This is the amount of data that is
    // sampled during one MCMC scan, and is also the amount of storage
    // associated with `m_subj`.
    const int m_n_subj;

    // tracks whether we wish to save the samples of xi to return to the user
    bool m_record_status;

    XiGen(Rcpp::List subj_day_blocks, int n_samp, bool record_status);
    ~XiGen();

    void sample(const WGen& W, const PhiGen& phi, const UProdBeta& u_prod_beta);

    const double* vals() const { return m_vals; }
    const int n_subj() const { return m_n_subj; }
};


#include "WGen.h"
#include "PhiGen.h"

#endif
