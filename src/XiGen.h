#ifndef DSP_BAYES_SRC_XI_GEN_H
#define DSP_BAYES_SRC_XI_GEN_H

#include "global_vars.h"


class XiGen {

    // storage for the `xi_i` values
    double* m_xi;

    // the elements of `m_subj` each map an individual to a block of days from
    // the day-specific data
    const DayBlock* m_subj;

    // the number of subjects in the data.  This value provides the amount of
    // storage that is associated with `m_n_subj`.
    const int m_n_subj;

    XiGen();
    ~XiGen();

    const double* XiGen::sample(WGen W, PhiGen phi, double* exp_uprod_beta);
};


#endif
