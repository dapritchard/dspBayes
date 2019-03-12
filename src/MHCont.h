#ifndef DSP_BAYES_MH_CONT_H_
#define DSP_BAYES_MH_CONT_H_

#include "Rcpp.h"


class MHCont {
public:

    // current value of phi and storage for previous values
    Rcpp::NumericVector m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    // proposal distribution dispersion parameter
    const double m_prp_disp;

    // tracks the number of times that the proposal distribution was accepted
    int m_accept_ctr;

    // tracks whether we wish to save the samples of phi to return to the user
    const bool m_record_status;

    // constructor and base functions
    MHCont(double proposal_dispersion, int n_samp, bool record_status);
    double update(double log_r, double proposal_val);
    double val() const { return *m_vals; }
};


#endif
