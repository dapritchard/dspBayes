#include "Rcpp.h"

#include "GammaGen.h"
#include "global_vars.h"

using R::lgammafn_sign;



// GammaFWDay::GammaFWDay() {

// }




double GammaFWDay::sample(Mday mday, Mu mu, Nu nu, Delta delta) {

    double proposal_val, new_val, log_r;

    // calculate the log acceptance ratio
    const double log_r = get_log_r(W, xi, ubeta, X, proposal_beta, proposal_gam);

    // accept proposal value `min(r, 1)-th` of the time
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
// where gamma* denotes the kgamma vector with the h-th term replaced by the
// proposal value and similarly for gamma^(s).

inline double GammaFWDay::get_log_r(const WGen& W,
				    const XiGen& xi,
				    const UProdBeta& ubeta,
				    const int* X,
				    double proposal_beta,
				    double proposal_gam) {

    return get_w_log_lik(W, xi, ubeta, X, proposal_beta)
	+ get_gam_log_lik(proposal_beta, proposal_gam);
}




// calculate
//
//         p(proposal_gam | m, mu, nu, delta)
//     log ----------------------------------
//         p(current_gam | m, mu, nu, delta)
//
// which simplifies to
//
//     (nu - 1)(log(proposal_gam) - log(current_gam))
//
//                   nu
//         - ------------------ (proposal_gam - current_gam)
//           delta^{|k - m|} mu

double GammaFWDay::get_gam_log_lik(double beta_proposal, Mday mday, Mu mu, Nu nu, Delta delta) {

    double term1, term2_multiple, term2_diff;
    double day_dist, ar_delta_pow;

    // calculate `(nu - 1)(log(proposal_gam) - log(current_gam))`
    term1 = (nu.val - 1) * (beta_proposal - m_beta_val);

    // calculate `(nu / delta^{|k - m|} mu)(proposal_gam - current_gam)`
    day_dist = abs(m_day_idx - mday);
    ar_delta_pow = pow(delta.val, day_dist);
    term2_multiple = nu.val / (ar_delta_pow * mu.val);

    return term1 - term2;
}
