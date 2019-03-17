#include <cmath>
#include "Rcpp.h"

#include "CoefGen.h"
#include "GammaFWDay.h"
#include "FWPriors.h"
#include "ProposalFcns.h"


Nu::Nu(int n_samp, bool record_status, double proposal_dispersion) :
    MHCont(n_samp, record_status, proposal_dispersion),
    m_alpha_0_minus_1 {1.0},
    m_beta_0          {1.0},
    m_nu_val          {2.0},
    m_log_nu_val      {std::log(2.0)}
{
    // initialize with an arbitrary value
    *m_vals = m_nu_val;
}




// update the value of `m_nu_val` using a Metropolis step
void Nu::sample(const CoefGen& coefs,
                const MDay& mday,
                const Mu& mu,
                const Delta& delta) {

    // sample the proposal value for Metropolis step.  TODO: generalize proposal function
    const double proposal_val     = ProposalFcns::abs_unif(*m_vals, m_prp_disp);
    const double log_proposal_val = std::log(proposal_val);

    // calculate `log(r)` where `r` is the acceptance ratio for the Metropolis
    // step
    double log_r = calc_log_r(coefs, mday, mu, delta, proposal_val, log_proposal_val);

    // sample the updated value by either accepting the proposal value or by
    // keeping the current value
    double new_val = update(log_r, proposal_val);
    if (new_val != m_nu_val) {
        m_nu_val = new_val;
        m_log_nu_val = std::log(m_nu_val);
    }

    // save the value of the new sample.  If we are recording samples then
    // increment `m_vals` so that we don't overwrite the previous sample.
    if (m_record_status && g_record_status) {
        ++m_vals;
    }
    *m_vals = new_val;
}




// calculate the value of `log(r)`, where `r` is the value such that the
// proposal value is accepted with probability `min(1,r)`
double Nu::calc_log_r(const CoefGen& coefs,
                      const MDay& mday,
                      const Mu& mu,
                      const Delta& delta,
                      double proposal_val,
                      double log_proposal_val) const {

    return calc_log_lik_gamma_term(coefs, mday, mu, delta, proposal_val, log_proposal_val)
        + calc_log_lik_nu_term(proposal_val, log_proposal_val);
}




// calculate
//
//         p(gamma_1, ..., gamma_K | m, mu, proposal_val, delta)
//     log ----------------------------------------------------
//         p(gamma_1, ..., gamma_K | m, mu, m_nu_val, delta)
//
//                           p(gamma_k | m, mu, proposal_val, delta)
//         = sum_{k=1}^K log --------------------------------------.
//                           p(gamma_k | m, mu, m_nu_val, delta)
//
// Furthermore, the k-th term of the sum simplifies to
//
//     proposal_val * log(proposal_val)
//
//         - m_nu_val * log(m_nu_val)
//
//         - log(GammaFn(proposal_val))
//
//         + log(GammaFn(m_nu_val))
//
//         + (proposal_val - m_nu_val) * (log(gamma_k) - log(delta^{|k-m|}) - (gamma_k / (delta^{|k-m|}))).
//
// Additionally, we note that terms not indexed by `k` can be factored out of
// each term in the sum for efficiency.

double Nu::calc_log_lik_gamma_term(const CoefGen& coefs,
                                   const MDay& mday,
                                   const Mu& mu,
                                   const Delta& delta,
                                   double proposal_val,
                                   double log_proposal_val) const {

    // extract priors for convenience
    const double mday_val  = mday.val();
    const double mu_val    = mu.val();
    const double delta_val = delta.val();
    const int K            = coefs.m_n_gamma;

    // used to collect the sum log-likelihood ratio of the individual gamma
    // terms (not including the terms that aren't indexed by the sum)
    double sum_log_lik = 0;

    // each iteration conditionally calculates the log-likelihood ratio of the
    // j-th gamma term and adds it to `sum_log-lik`
    for (int k = 0; k < coefs.m_n_gamma; ++k) {

        // only consider the gamma terms that are part of the fertile window
        // TODO: can we make an interface for gamma params?
        if (coefs.m_gamma[k]->is_fw_day()) {

            // current FW day index and coefficient value
            int curr_day_idx     = dynamic_cast<GammaFWDay*>(coefs.m_gamma[k])->m_day_idx;
            double curr_beta_val = coefs.m_gamma[k]->m_beta_val;
            double curr_gam_val  = coefs.m_gamma[k]->m_gam_val;

            // calculate `delta^{|k-m|} * mu`
            double day_dist          = abs(curr_day_idx - mday_val);
            double delta_pow         = pow(delta_val, day_dist);
            double delta_pow_prod_mu = delta_pow * mu_val;

            // calculate terms in the sum
            double term_a       = curr_beta_val;
            double term_b       = std::log(delta_pow_prod_mu);
            double term_c       = curr_gam_val / delta_pow_prod_mu;

            sum_log_lik += term_a - term_b - term_c;
        }
    }

    // calculate the terms in the sum
    double term_1 = proposal_val * std::log(proposal_val);
    double term_2 = m_nu_val * m_log_nu_val;
    double term_3 = R::lgammafn(proposal_val);
    double term_4 = R::lgammafn(m_nu_val);
    double term_5 = (proposal_val - m_nu_val) * sum_log_lik;

    return K * (term_1 + term_2 + term_3 + term_4) + term_5;
}




// calculate
//
//         p(proposal_val | a, b)
//     log ---------------------
//         p(m_nu_val | a, b)
//
//         = (a - 1) * (log(proposal_val) - log(m_nu_val))
//
//             - b * (proposal_val - m_nu_val)

double Nu::calc_log_lik_nu_term(double proposal_val, double log_proposal_val) const {

    double term1 = m_alpha_0_minus_1 * (log_proposal_val - m_log_nu_val);
    double term2 = m_beta_0 * (proposal_val - m_nu_val);

    return term1 - term2;
}
