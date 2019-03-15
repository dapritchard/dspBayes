#ifndef DSP_BAYES_FW_PRIORS_H_
#define DSP_BAYES_FW_PRIORS_H_

#include "MHCont.h"

// note that CoefGen.h included at the bottom of the file
class CoefGen;
class MDay;
class Mu;
class Nu;
class Delta;




class MDay {
public:

    double m_peak_idx;

    MDay() : m_peak_idx {3} {}
    double val() const { return m_peak_idx; }

    // void sample();
};




class Mu : public MHCont {
public:

    // hyperparameters for mu
    const double m_alpha_0_minus_1;
    const double m_beta_0;

    // current value and log-value
    double m_mu_val;
    double m_log_mu_val;

    Mu(double proposal_dispersion, int n_samp, bool record_status);
    void sample(const CoefGen& coefs, const MDay& mday, const Nu& nu, const Delta& delta);

    double calc_log_r(const CoefGen& coefs,
                      const MDay& mday,
                      const Nu& nu,
                      const Delta& delta,
                      double proposal_val,
                      double log_proposal_val) const;

    double calc_log_lik_gamma_term(const CoefGen& coefs,
                                   const MDay& mday,
                                   const Nu& nu,
                                   const Delta& delta,
                                   double proposal_val,
                                   double log_proposal_val) const;

    double calc_log_lik_mu_term(double proposal_val, double log_proposal_val) const;
};




class Nu {
public:

    double m_nu;

    Nu() : m_nu {100} {}
    double val() const { return m_nu; }
};




class Delta {
public:

    double m_delta;

    Delta() : m_delta {0.5} {}
    double val() const { return m_delta; }
};




// FIXME: need to call destructors for these classes?
class FWPriors {
public:

    MDay m_mday;
    Mu m_mu;
    Nu m_nu;
    Delta m_delta;

    FWPriors(const Rcpp::List& fw_prior_specs) :
        m_mday  {MDay()},
        m_mu    {build_mu(fw_prior_specs)},
        m_nu    {Nu()},
        m_delta {Delta()}
    {}

    void sample(const CoefGen& coefs) { // FIXME
        // m_mday.sample();
        m_mu.sample(coefs, m_mday, m_nu, m_delta);
        // m_nu.sample();
        // m_delta.sample();
    }

    Mu build_mu(const Rcpp::List& fw_prior_specs) { // FIXME
        // Rcpp::List mu_specs {fw_prior_specs["mu_specs"]};
        Mu out {Mu(1.0, 1.0, true)};
        return out;
    }
};


#include "CoefGen.h"

#endif
