#ifndef DSP_BAYES_SRC_PHI_GEN_H
#define DSP_BAYES_SRC_PHI_GEN_H

#include "Rcpp.h"
class XiGen;


class PhiGen {

public:

    // gamma distribution hyperparameters c1 and c2 and tuning parameter for the
    // proposal distribution delta
    const double m_hyp_c1;
    const double m_hyp_c2;
    const double m_delta;

    // current value of phi
    double m_phi_val;

    // tracks the number of times that the proposal distribution was accepted
    int m_accept_ctr;

    // whether the proposal distribution was not accepted so that the value of
    // phi is unchanged from the last scan.  When this is the case then the
    // calculation for `m_log_norm_const` can be reused
    bool m_is_same_as_prev;
    // the log of the normalizing constant of a gamma distribution where both
    // parameters are given by the current value of phi and which is given by
    // `log(phi^phi / gamma(phi))`
    double m_log_norm_const;


    PhiGen(Rcpp::NumericVector phi_hyper);
    double val() const { return m_phi_val; }
    int n_accept() const { return m_accept_ctr; }
    void sample(const XiGen& xi);

    double calc_log_r(const XiGen& xi, double proposal_val);
    double update_phi(double log_r, double proposal_val);
    double calc_log_proportion_dgamma_xi(const XiGen& xi, double proposal_val);
    double calc_log_proportion_dgamma_phi(const XiGen& xi, double proposal_val);
    double log_dgamma_norm_const(double a);

};


#include "XiGen.h"

#endif
