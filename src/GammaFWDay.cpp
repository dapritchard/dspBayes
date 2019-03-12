#include "Rcpp.h"

#include "GammaGen.h"
#include "GammaFWDay.h"
#include "global_vars.h"
#include "ProposalFcns.h"


GammaFWDay::GammaFWDay(const Rcpp::NumericMatrix& U,
                       const Rcpp::NumericVector& gamma_specs) :
    GammaContMH {U, gamma_specs},
    m_day_idx   {static_cast<int>(gamma_specs["h"]) + 1}  // TODO: make this indep of loc in the design mat?
{}




// TODO: the only difference in this function and the GammaContMH version is the
// signature for `get_log_r`.  Any way to consolidate?
double GammaFWDay::sample(const WGen& W,
                          const XiGen& xi,
                          UProdBeta& ubeta,
                          const int* X,
                          const FWPriors& fw_priors) {

    // std::cout << "beginning of GammaFWDay::sample" << std::endl;
    // std::cout << "start m_proposal_fcn" << std::endl;
    // std::cout << "m_beta_val:  " << m_beta_val << std::endl;
    // std::cout << "m_mh_delta:  " << m_mh_delta << std::endl;
    // const double proposal_beta = m_proposal_fcn(m_beta_val, m_mh_delta);
    const double proposal_beta = (m_beta_val - m_mh_delta) + (2.0 * m_mh_delta * R::unif_rand());
    // std::cout << "start exp(proposal_beta)" << std::endl;
    const double proposal_gam  = exp(proposal_beta);

    // calculate the log acceptance ratio
    const double log_r = get_log_r(W, xi, ubeta, X, proposal_beta, proposal_gam, fw_priors);

    // accept proposal value `min(r, 1)-th` proportion of the time
    if ((log_r >= 0) || (log(R::unif_rand()) < log_r)) {

        // update `U * beta` and `exp(U * beta)` based upon accepting the
        // proposal value
        ubeta.update(m_Uh, proposal_beta, m_beta_val);

        // update member variables to based upon accepting the proposal value
        m_beta_val = proposal_beta;
        m_gam_val = proposal_gam;
        ++m_mh_accept_ctr;
    }

    return m_gam_val;
}




// Calculate log gamma acceptance ratio r.  The acceptance ratio for gamma_h is
// given by:
//
//            p(W | gamma*, xi, data) * p(gamma_h*)
//         -------------------------------------------
//         p(W | gamma^(s), xi, data) * p(gamma_h^(s))
//
// where gamma* denotes the gamma vector with the h-th term replaced by the
// proposal value and similarly for gamma^(s).

inline double GammaFWDay::get_log_r(const WGen& W,
                                    const XiGen& xi,
                                    const UProdBeta& ubeta,
                                    const int* X,
                                    double proposal_beta,
                                    double proposal_gam,
                                    const FWPriors& fw_priors) {

    // // std::cout << "start `get_w_log_lik`" << std::endl;
    // get_w_log_lik(W, xi, ubeta, X, proposal_beta);
    // // std::cout << "start `get_w_log_lik`" << std::endl;

    // // std::cout << "start `get_gam_log_lik`" << std::endl;
    // get_gam_log_lik(proposal_beta, proposal_gam, fw_priors);
    // // std::cout << "start `get_gam_log_lik`" << std::endl;

    return get_w_log_lik(W, xi, ubeta, X, proposal_beta)
        + get_gam_log_lik(proposal_beta, proposal_gam, fw_priors);
}




// calculate
//
//         p(proposal_gam | m, mu, nu, delta)
//     log ----------------------------------,
//         p(current_gam | m, mu, nu, delta)
//
// which simplifies to
//
//     (nu - 1)(log(proposal_gam) - log(current_gam))
//
//                   nu
//         - ------------------ (proposal_gam - current_gam).
//           delta^{|k - m|} mu
//
// Also recall that the logarithm of the h-th element of the gamma vector is the
// h-th element of the beta vector.

double GammaFWDay::get_gam_log_lik(double proposal_beta,
                                   double proposal_gam,
                                   FWPriors fw_priors) {

    // extract priors for convenience
    double mday_val  = fw_priors.mday.val();
    double mu_val    = fw_priors.mu.val();
    double nu_val    = fw_priors.nu.val();
    double delta_val = fw_priors.delta.val();

    // calculate `(nu - 1)(log(proposal_gam) - log(current_gam))`
    double term1 = (nu_val - 1) * (proposal_beta - m_beta_val);

    // calculate `(nu / (delta^{|k - m|} mu)) * (proposal_gam - current_gam)`
    double day_dist         = abs(m_day_idx - mday_val);
    double ar_delta_pow     = pow(delta_val, day_dist);
    double term2_multiplier = nu_val / (ar_delta_pow * mu_val);
    double term2_diff       = proposal_gam - m_gam_val;
    double term2            = term2_multiplier * term2_diff;

    return term1 - term2;
}
