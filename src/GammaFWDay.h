#ifndef DSP_BAYES_GAMMA_FWDAY_H_
#define DSP_BAYES_GAMMA_FWDAY_H_

#include "Rcpp.h"
#include "FWPriors.h"


class GammaFWDay : public GammaContMH {
public:

    // // the current value of gamma_h
    // double m_beta_val;
    // double m_gam_val;

    // the FW day index
    int m_day_idx;

    // the decay values and the index of the midpoint of the array
    std::vector<double> m_decay_vals;
    int m_midpoint_idx;

    // the proposal distribution
    double (*m_proposal_fcn)(double cond, double delta);

    GammaFWDay(const Rcpp::NumericMatrix& U,
               const Rcpp::NumericVector& gamma_specs);
    double sample(const WGen& W,
                  const XiGen& xi,
                  UProdBeta& ubeta,
                  const int* X,
                  const FWPriors& fw_priors);
    void inject_decay_vals(const Rcpp::NumericVector& decay_vals);
    double get_log_r(const WGen& W,
                     const XiGen& xi,
                     const UProdBeta& ubeta,
                     const int* X,
                     double proposal_beta,
                     double proposal_gam,
                     const FWPriors& fw_priors);
    // double get_gam_log_lik(double proposal_beta,
    //                        double proposal_gam,
    //                        FWPriors fw_priors);
    double get_gam_log_lik(double proposal_beta,
                           double proposal_gam,
                           const FWPriors& fw_priors) const;
    double decay(int peak_intensity_idx) const { return m_decay_vals[m_midpoint_idx + m_day_idx - peak_intensity_idx]; }
    bool is_fw_day() const { return true; }
};


#endif
