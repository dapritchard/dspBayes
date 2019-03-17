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

    Mu(int n_samp, bool record_status, double proposal_dispersion);
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




class Nu : public MHCont {
public:

    // hyperparameters for nu
    const double m_alpha_0_minus_1;
    const double m_beta_0;

    // current value and log-value
    double m_nu_val;
    double m_log_nu_val;

    Nu(int n_samp, bool record_status, double proposal_dispersion);
    void sample(const CoefGen& coefs, const MDay& mday, const Mu& mu, const Delta& delta);

    double calc_log_r(const CoefGen& coefs,
                      const MDay& mday,
                      const Mu& mu,
                      const Delta& delta,
                      double proposal_val,
                      double log_proposal_val) const;

    double calc_log_lik_gamma_term(const CoefGen& coefs,
                                   const MDay& mday,
                                   const Mu& mu,
                                   const Delta& delta,
                                   double proposal_val,
                                   double log_proposal_val) const;

    double calc_log_lik_nu_term(double proposal_val, double log_proposal_val) const;
};




class Delta : public MHCont {
public:

    // hyperparameters for nu
    const double m_alpha_0_minus_1;
    const double m_beta_0_minus_1;

    // current value and log-value
    double m_delta_val;
    double m_log_delta_val;

    Delta(int n_samp, bool record_status, double proposal_dispersion);
    void sample(const CoefGen& coefs, const MDay& mday, const Mu& mu, const Nu& nu);

    double calc_log_r(const CoefGen& coefs,
                      const MDay& mday,
                      const Mu& mu,
                      const Nu& nu,
                      double proposal_val,
                      double log_proposal_val) const;

    double calc_log_lik_gamma_term(const CoefGen& coefs,
                                   const MDay& mday,
                                   const Mu& mu,
                                   const Nu& nu,
                                   double proposal_val,
                                   double log_proposal_val) const;

    double calc_log_lik_nu_term(double proposal_val, double log_proposal_val) const;
};





class FWPriors {
public:

    MDay m_mday;
    Mu m_mu;
    Nu m_nu;
    Delta m_delta;

    FWPriors(const Rcpp::List& fw_prior_specs) :
        m_mday  {MDay()},
        m_mu    {build_mu(fw_prior_specs)},
        m_nu    {build_nu(fw_prior_specs)},
        m_delta {build_delta(fw_prior_specs)}
    {}

    void sample(const CoefGen& coefs) { // FIXME
        // m_mday.sample();
        m_mu.sample(coefs, m_mday, m_nu, m_delta);
        m_nu.sample(coefs, m_mday, m_mu, m_delta);
        m_delta.sample(coefs, m_mday, m_mu, m_nu);
    }

    // FIXME
    Delta build_delta(const Rcpp::List& fw_prior_specs) { // FIXME
        // Rcpp::List mu_specs {fw_prior_specs["mu_specs"]};
        Delta out {Delta(50000, true, 0.1)};
        return out;
    }

    // FIXME
    Mu build_mu(const Rcpp::List& fw_prior_specs) { // FIXME
        // Rcpp::List mu_specs {fw_prior_specs["mu_specs"]};
        Mu out {Mu(50000, true, 0.1)};
        return out;
    }

    // FIXME
    Nu build_nu(const Rcpp::List& fw_prior_specs) { // FIXME
        // Rcpp::List mu_specs {fw_prior_specs["mu_specs"]};
        Nu out {Nu(50000, true, 0.1)};
        return out;
    }
};


#include "CoefGen.h"

#endif
