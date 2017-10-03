#include "Rcpp.h"

#include "GammaGen.h"
#include "global_vars.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"



GammaContMH::GammaContMH() :

{


}




double GammaContMH::sample(const WGen& W, const XiGen& xi, const int* X) {

    double proposal_beta = samp_proposal_beta();
    double proposal_gam  = exp(proposal_beta);

    double log_r = get_log_r(proposal_gam);

    if ((log_r >= 0) || (log(R::unif_rand()) < log_r)) {
	m_beta_val = proposal_beta;
	m_gam_val = proposal_gam;
	++m_mh_accept_ctr;
    }

    return m_gam_val;
}




inline double GammaContMH::sample_proposal_beta() const {

    return (R::unif_rand() < m_mh_prob_samp_1) ?
	0.0 :
	m_proposal_fcn(m_beta_val, m_mh_delta);
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

inline double GammaContMH::get_log_r(const WGen& W,
				     const XiGen& xi,
				     const UProdBeta& ubeta,
				     double proposal_beta,
				     double proposal_gam) {

    // // the value of `beta_h* - beta_h^(s)`
    // double beta_diff = log(proposal_gam / m_gam_val);

    return (get_w_log_lik(W, xi, ubeta, beta_diff)
	    + get_gam_log_lik(proposal_beta, proposal_gam)
	    + get_proposal_log_lik(proposal_beta));
}




// calculate `p(W | proposal_gam, xi) / p(W | current_gam, xi)`

double GammaContMH::get_w_log_lik(const WGen& W,
				  const XiGen& xi,
				  const UProdBeta& ubeta,
				  double beta_diff) {

    const int* w_vals            = W.vals();
    const bool is_preg_day       = W.is_preg_day();
    const double* xi_vals        = xi.vals();
    const double* ubeta_exp_vals = ubeta.exp_vals();

    // tracks the running total of the log-likelihood
    double sum_log_lik = 0;

    // each iteration adds the i-th value of the loglikelihood to the running
    // value of `sum_log_lik`
    for (int i = 0; i < n_days; ++i) {

	// TODO: don't need this snippet below?

	// // case:
	// if ((preg_day_idx < w_n_preg_days) && (w_days_idx[preg_day_idx] == i)) {
	//     term1 = w_vals[preg_day_idx++];
	// }
	// // case:
	// else {
	//     term1 = 0.0;
	// }



	// map the current day to the i-th subject to obtain `xi_i`
	double xi_i = xi_vals[d2s[i]];

	// calculate
	//
	//           [xi_i * exp(u_{ijk}^T beta*)]^{w_{ijk}
	//     log -----------------------------------------
	//         [xi_i * exp(u_{ijk}^T beta^(s))]^{w_{ijk}
	//
	//         = w_{ijk} * u{ijkh} * (beta_h* - beta_h^(s))
	//
	// which is one of the terms in `p(W | proposal) / p(W | current)`.

	double term1 = (is_preg_day[i]) ?
	    w_vals++ * m_Uh[i] * beta_diff :
	    0.0;

	// // calculate
	// //
	// //           exp{ -xi_i * exp(u_{ijk}^T beta*) }
	// //     log --------------------------------------
	// //         exp{ -xi_i * exp(u_{ijk}^T beta^(s)) }
	// //
	// //         = -xi_i * [ exp(u_{ijk}^T beta*) - exp(u_{ijk}^T beta^(s)) ]
	// //
	// //         = -xi_i * (exp{ u_ijkh * (beta_h* - beta_h^(s)) } - 1) * exp(u_{ijk}^T beta^(s))
	// //
	// // which is one of the terms in `p(W | proposal) / p(W | current)`.

	// double term2 = -xi_i * (exp(m_Uh[i] * beta_diff) - 1) * ubeta_exp_vals[i];

	// calculate `-xi_i * [exp(U * beta*) - exp(U * beta)]`, which is one of
	// the terms in `p(W | proposal) / p(W | current)`.
	double term2 = -xi_i * (exp(ubeta_vals[i] + (m_Uh[i] * beta_diff)) - ubeta_exp_vals[i]);

	// add the portion of the log-likelihood from the current day to the
	// running total
	sum_log_lik += term1 + term2;
    }
}




// calculate  `p(gamma_h*) / p(gamma_h)` which is given by
//
//       /   1,                                            gamma_h* = 1,  gamma_h^(s) = 1
//      /
//      |                      p_h
//      |    ---------------------------------------- ,    gamma_h* = 1,  gamma_h^(s) != 1
//      |    (1 - p_h) * Gamma(gamma_h^(s); a_h, b_h)
//      /
//     <
//      \    (1 - p_h) * Gamma(gamma_h*; a_h, b_h)
//      |    ------------------------------------- ,       gamma_h* != 1,  gamma_h^(s) = 1
//      |                    p_h
//      |
//      |      Gamma(gamma_h*; a_h, b_h)
//      \    ---------------------------- ,                gamma_h* != 1,  gamma_h^(s) != 1
//       \   Gamma(gamma_h^(s); a_h, b_h)
//
//
// The expressions below for the middle two cases are due to the equality
//
//      log Gamma(x; a, b) = log(norm_const) + ((a - 1) * x) - (b * x)
//
// plus one extra term for the log of the `p_h / (1 - p_h)` term or vice versa.
// The expression below for the last term is given by
//
//         Gamma(x*; a, b)
//     log --------------- = log{ (x* / x)^(a - 1) * exp(-b * (x* - x)) }
//          Gamma(x; a, b)
//
//                         = (a - 1) * log(x* / x) - (b * (x* - x))
//
// and noting that for our case we have
//
//      log(x* / x) = log(gamma_h* / gamma_h)
//
//                  = log(exp(beta_h*) / exp(beta_h))
//
//                  = beta_h* - beta_h

double GammaContMH::get_gam_log_lik(double proposal_beta, double proposal_gam) {

    double gam_log_lik;

    if ((proposal_gam = 1.0) && (m_gam_val == 1.0)) {
	gam_log_lik = 0.0;
    }
    else if ((proposal_gam == 1.0) && (m_gam_val != 1.0)) {
	gam_log_lik = (m_log_ph_over_1_minus_ph
		       - m_log_norm_const
		       - ((m_hyp_a - 1) * m_beta_val)
		       + (m_hyp_b * m_gam_val));
    }
    else if ((proposal_gam != 1.0) && (m_gam_val == 1.0)) {
	gam_log_lik = (m_log_1_minus_ph_over_ph
		       + m_log_norm_const
		       + ((m_hyp_a - 1) * log(proposal_gam))
		       - (m_hyp_b * proposal_gam));
    }
    else {
	const double beta_diff = proposal_beta - m_beta_val;
	const double gam_diff  = proposal_gam - m_gam_val;
	gam_log_lik = ((m_hyp_a - 1) * beta_diff) - (m_hyp_b * (gam_diff));
    }

    return gam_log_lik;
}




inline double GammaContMH::get_proposal_log_lik(double proposal_beta) const {

    // calc `J(beta^(s) | beta*)` where `J` is the density function of the
    // proposal distribution
    double log_numer = (m_beta_val == 1) ?
	m_mh_log_prob_samp_1 :
	m_mh_log_1_minus_prob_samp_1 + m_log_proposal_den(m_beta_val, proposal_beta, m_mh_delta);

    // calc `J(beta* | beta^(s))` where `J` is the density function of the
    // proposal distribution
    double log_denom = (proposal_beta  == 1) ?
	m_mh_log_prob_samp_1 :
	m_mh_log_1_minus_prob_samp_1 + m_log_proposal_den(proposal_beta, m_beta_val, m_mh_delta);

    return log_numer - log_denom;
}
