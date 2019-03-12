#ifndef DSP_BAYES_SRC_PHI_GEN_H
#define DSP_BAYES_SRC_PHI_GEN_H

#include "Rcpp.h"
class XiGen;

extern bool g_record_status;


class PhiGen {

public:

    // gamma distribution hyperparameters c1 and c2 and
    const double m_hyp_c1;
    const double m_hyp_c2;

    //tuning parameter for the proposal distribution
    const double m_delta;

    // current value of phi and storage for previous values
    Rcpp::NumericVector m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    // tracks the number of times that the proposal distribution was accepted
    int m_accept_ctr;

    // tracks whether we wish to save the samples of phi to return to the user
    const bool m_record_status;

    // whether the proposal distribution was not accepted so that the value of
    // phi is unchanged from the last scan.  When this is the case then the
    // calculation for `m_log_norm_const` can be reused
    bool m_is_same_as_prev;
    // the log of the normalizing constant of a gamma distribution where both
    // parameters are given by the current value of phi and which is given by
    // `log(phi^phi / gamma(phi))`
    double m_log_norm_const;


    PhiGen(Rcpp::NumericVector phi_hyper, int n_samp, bool record_status);

    void sample(const XiGen& xi);
    double val() const { return *m_vals; }
    int n_accept() const { return m_accept_ctr; }

    double calc_log_r(const XiGen& xi, double proposal_val);
    double update_phi(double log_r, double proposal_val);
    double calc_log_proportion_dgamma_xi(const XiGen& xi, double proposal_val);
    double calc_log_proportion_dgamma_phi(double proposal_val) const;
    double log_dgamma_norm_const(double a) const;
};


#include "XiGen.h"

#endif
