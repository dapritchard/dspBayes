#ifndef DSP_BAYES_SRC_PHI_GEN_H
#define DSP_BAYES_SRC_PHI_GEN_H

#include "global_vars.h"


class PhiGen {

    // gamma distribution hyperparameters c1 and c2 and tuning parameter for the
    // proposal distribution delta
    const double m_hyp_c1;
    const double m_hyp_c2;
    const double m_delta;

    // current value of phi
    double m_val;

    // number of subjects in the data
    const int m_n_subj;

    // tracks the number of times that the proposal distribution was accepted
    int m_accept_ctr;

    // whether the proposal distribution was not accepted so that the value of
    // phi is unchanged from the last scan.  When this is the case some
    // calculations can be reused.
    bool m_is_same_as_prev;
    //
    double m_prev_log_norm_const;


    double sample();
    double calc_log_r(const XiGen& xi, double proposal_val);
    double update_phi(double log_r, double proposal_val);
    double calc_log_proportion_dgamma_xi(const XiGen& xi, double proposal_val);
    double calc_log_proportion_dgamma_phi(const XiGen& xi, double proposal_val);
    double log_dgamma_norm_const(double a);
};


#endif
