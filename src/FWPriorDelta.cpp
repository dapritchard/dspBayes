#include <cmath>
#include "Rcpp.h"

#include "CoefGen.h"
#include "GammaFWDay.h"
#include "FWPriors.h"
#include "ProposalFcns.h"


Delta::Delta(int n_samp, bool record_status, double proposal_dispersion) :
    MHCont(n_samp, record_status, proposal_dispersion),
    m_alpha_0_minus_1 {1.0},
    m_beta_0_minus_1  {1.0},
    m_delta_val       {0.7},
    m_log_delta_val   {std::log(0.7)}
{
    *m_vals = m_delta_val;
}




// update the value of `m_delta_val` using a Metropolis step
void Delta::sample(const CoefGen& coefs,
                   const MDay& mday,
                   const Mu& mu,
                   const Nu& nu) {

    // sample the proposal value for Metropolis step.  TODO: generalize proposal function
    const double proposal_val     = ProposalFcns::abs_unif(*m_vals, m_prp_disp);
    const double log_proposal_val = std::log(proposal_val);

    // calculate `log(r)` where `r` is the acceptance ratio for the Metropolis
    // step
    double log_r = calc_log_r(coefs, mday, mu, nu, proposal_val, log_proposal_val);

    // sample the updated value by either accepting the proposal value or by
    // keeping the current value
    double new_val = update(log_r, proposal_val);
    if (new_val != m_delta_val) {
        m_delta_val = new_val;
        m_log_delta_val = std::log(m_delta_val);
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
double Delta::calc_log_r(const CoefGen& coefs,
                         const MDay& mday,
                         const Mu& mu,
                         const Nu& nu,
                         double proposal_val,
                         double log_proposal_val) const {

    return calc_log_lik_gamma_term(coefs, mday, mu, nu, proposal_val, log_proposal_val)
        + calc_log_lik_nu_term(proposal_val, log_proposal_val);
}




// calculate
//
//         p(gamma_1, ..., gamma_K | m, mu, nu, delta_proposal)
//     log ----------------------------------------------------
//         p(gamma_1, ..., gamma_K | m, mu, nu, delta_current)
//
//                           p(gamma_k | m, mu, nu, delta_proposal)
//         = sum_{k=1}^K log --------------------------------------.
//                           p(gamma_k | m, mu, nu, delta_current)

// Furthermore, the k-th term of the sum simplifies to
//
//         = -nu * (log(delta_proposal^{|k-m|}) - log(delta_current^{|k-m|}))
//
//               nu * gamma_k
//             - ------------ ((1 / delta_proposal^{|k-m|}) - (1 / delta_current^{|k-m|})).
//                    mu
//
// Additionally, we note that terms not indexed by `k` can be factored out of
// each term in the sum for efficiency.

double Delta::calc_log_lik_gamma_term(const CoefGen& coefs,
                                      const MDay& mday,
                                      const Mu& mu,
                                      const Nu& nu,
                                      double proposal_val,
                                      double log_proposal_val) const {

    // extract priors for convenience
    const double mday_val = mday.val();
    const double mu_val   = mu.val();
    const double nu_val   = nu.val();

    // used to collect the sum log-likelihood ratio of the individual gamma
    // terms (not including the terms that aren't indexed by the sum)
    double sum_term1 = 0;
    double sum_term2 = 0;

    // each iteration conditionally calculates the log-likelihood ratio of the
    // j-th gamma term and adds it to `sum_log-lik`
    for (int k = 0; k < coefs.m_n_gamma; ++k) {

        // only consider the gamma terms that are part of the fertile window
        // TODO: can we make an interface for gamma params?
        // TODO: when `day_dist` equals 0 then the k-th term in the sum is 0
        if (coefs.m_gamma[k]->is_fw_day()) {

            // current FW day index and coefficient value
            int curr_day_idx    = dynamic_cast<GammaFWDay*>(coefs.m_gamma[k])->m_day_idx;
            double curr_gam_val = coefs.m_gamma[k]->m_gam_val;

            // calculate `|k-m|` and `-|k-m|`
            double day_dist           = abs(curr_day_idx - mday_val);
            double neg_day_dist       = -day_dist;

            // calculate `1 / delta_proposal^{|k-m|}` and `1 / delta_current^{|k-m|}`
            double inv_delta_pow_prp  = pow(proposal_val, neg_day_dist);
            double inv_delta_pow_curr = pow(m_delta_val, neg_day_dist);

            // add in the terms for the k-th FW day
            sum_term1 += day_dist * (log_proposal_val - m_log_delta_val);
            sum_term2 += curr_gam_val * (inv_delta_pow_prp - inv_delta_pow_curr);
        }
    }

    // multiply by the constant terms and return
    sum_term2 /= mu_val;
    return -nu_val * (sum_term1 + sum_term2);
}




// calculate
//
//         p(proposal_val | a, b)
//     log ---------------------
//         p(m_nu_val | a, b)
//
//         = (a - 1) * (log(delta_proposal) - log(delta_current))
//
//             + (b - 1) * (log(1 - delta_proposal) - log(1 - delta_current))

double Delta::calc_log_lik_nu_term(double proposal_val, double log_proposal_val) const {

    double term1 = m_alpha_0_minus_1 * (log_proposal_val - m_log_delta_val);
    double term2 = m_beta_0_minus_1 * (std::log(1 - proposal_val) - std::log(1 - m_delta_val));

    return term1 - term2;
}
