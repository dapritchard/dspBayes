#include <cmath>
#include "Rcpp.h"
#include "ProposalFcns.h"


Mu::Mu(double proposal_dispersion, int n_samp, bool record_status) :
    MHCont(proposal_dispersion, n_samp, record_status),
    m_alpha_0    {1.0},
    m_beta_0     {1.0},
    m_mu_val     {0.44},
    m_log_mu_val {std::log(0.44)}
{}




void Mu::sample(const CoefGen& coefs,
                const MDay& mday,
                const Nu& nu,
                const Delta& delta) {

    // sample the proposal value for Metropolis step.  TODO: generalize
    const double proposal_val = ProposalFcns::abs_unif(*m_vals, m_prp_disp);

    // calculate `log(r)` where `r` is the acceptance ratio for the Metropolis
    // step
    double log_r = calc_log_r(coefs, proposal_val);

    // sample the updated value by either accepting the proposal value or by
    // keeping the current value
    double new_val = update(log_r, proposal_val);

    // save the value of the new sample.  If we are recording samples then
    // increment `m_vals` so that we don't overwrite the previous sample.
    if (m_record_status && g_record_status) {
        ++m_vals;
    }
    *m_vals = new_val;
}




double Mu::calc_log_r(const CoefGen& coefs,
                      double proposal_val,
                      double log_proposal_val) {

    return get_w_log_lik(W, xi, ubeta, X, proposal_beta)
        + get_mu_prior_log_lik(proposal_val, log_proposal_val);
}



double Mu::calc_mu_log_lik_gamma_term(const CoefGen& coefs,
                                      const MDay& mday,
                                      const Nu& nu,
                                      const Delta& delta) {

    // extract priors for convenience
    double mday_val  = fw_priors.mday.val();
    double nu_val    = fw_priors.nu.val();
    double delta_val = fw_priors.delta.val();

    //
    double sum_log_lik = 0;

    for (int j = 0; j < coefs.m_n_gamma; ++j) {
        if (m_gamma[j]->is_fw_day()) {

            int curr_day_idx    = m_gamma[j]->m_day_idx;
            double curr_gam_val = m_gamma[j]->m_gam_val;

            double term1 = nu_val * (log_proposal_val - m_log_mu_val);

            double day_dist         = abs(curr_day_idx - mday_val);
            double ar_delta_pow     = pow(delta_val, day_dist);
            double term2_multiplier = -nu_val * curr_gam_val / ar_delta_pow;
            double term2_diff       = (1 / proposal_val) - (1 / m_mu_val);

            sum_log_lik += -term1 - term2;
        }
    }
}




double Mu::calc_mu_log_lik_prior_term(double proposal_val, double log_proposal_val) {

    double term1 = m_alpha_0_minus_1 * (log_proposal_val - m_log_mu_val);
    double term2 = m_beta_0 * (proposal_val - m_mu_val);

    return term1 - term2;
}
