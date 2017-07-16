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

    // the elements of `m_subj` each map an individual to a block of days from
    // the day-specific data
    const DayBlock* m_subj;

    // the number of subjects in the data.  This value provides the amount of
    // storage that is associated with `m_n_subj`.
    const int m_n_subj;

    XiGen(Rcpp::NumericVector xi_initial, Rcpp::List subj_days);
    ~XiGen();

    void sample(WGen W, PhiGen phi, UProdBeta& u_prod_beta);
    const double* vals() const { return m_xi_vals; }
    const int n_subj() const { return m_n_subj; }

};


#include "WGen.h"
#include "PhiGen.h"

#endif
