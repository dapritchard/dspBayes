#include "Rcpp.h"

#include "GammaGen.h"
#include "global_vars.h"
#include "ProposalFcns.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"




GammaContMH::GammaContMH(const Rcpp::NumericMatrix& U,
			 const Rcpp::NumericVector& gamma_specs) :
    GammaGen(U, gamma_specs),
    m_log_norm_const(log_dgamma_trunc_norm_const()),
    m_log_p_over_1_minus_p(log(m_hyp_p / (1 - m_hyp_p))),
    m_log_1_minus_p_over_p(-m_log_p_over_1_minus_p),
    m_mh_p(gamma_specs["mh_p"]),
    m_mh_log_p(log(m_mh_p)),
    m_mh_log_1_minus_p(log(1 - m_mh_p)),
    m_mh_delta(gamma_specs["mh_delta"]),
    m_mh_accept_ctr(0),
    m_proposal_fcn(ProposalFcns::unif),
    m_log_proposal_den(ProposalFcns::log_den_unif) {
}
// TODO: provide a way to choose the proposal functions




double GammaContMH::sample(const WGen& W, const XiGen& xi, UProdBeta& ubeta, const int* X) {

    const double proposal_beta = sample_proposal_beta();
    const double proposal_gam  = exp(proposal_beta);

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




inline double GammaContMH::sample_proposal_beta() const {

    return (R::unif_rand() < m_mh_p) ?
	0.0 :
	m_proposal_fcn(m_beta_val, m_mh_delta);
}




// Calculate log gamma acceptance ratio r.  The acceptance ratio for gamma_h is
// given by:
//
//            p(W | gamma*, xi, data) * p(gamma_h*) / J(beta* | beta^(s))
//         -----------------------------------------------------------------
//         p(W | gamma^(s), xi, data) * p(gamma_h^(s)) / J(beta^(s) | beta*)
//
// where gamma* denotes the gamma vector with the h-th term replaced by the
// proposal value and similarly for gamma^(s).

inline double GammaContMH::get_log_r(const WGen& W,
				     const XiGen& xi,
				     const UProdBeta& ubeta,
				     const int* X,
				     double proposal_beta,
				     double proposal_gam) {

    return (get_w_log_lik(W, xi, ubeta, X, proposal_beta)
	    + get_gam_log_lik(proposal_beta, proposal_gam)
	    + get_proposal_log_lik(proposal_beta));
}




// calculate `p(W | proposal_gam, xi) / p(W | current_gam, xi)`

double GammaContMH::get_w_log_lik(const WGen& W,
				  const XiGen& xi,
				  const UProdBeta& ubeta,
				  const int* X,
				  double proposal_beta) const {

    const int* w_vals            = W.vals();
    const int* w_days_idx        = W.days_idx();
    const double* xi_vals        = xi.vals();
    const double* ubeta_vals     = ubeta.vals();
    const double* ubeta_exp_vals = ubeta.exp_vals();

    // the value of `beta_h* - beta_h^(s)`
    double beta_diff = proposal_beta - m_beta_val;

    // tracks the running total of the log-likelihood
    double sum_log_lik = 0;

    // each iteration adds the i-th value of the loglikelihood to the running
    // value of `sum_log_lik`
    for (int i = 0; i < m_n_days; ++i) {

	// if intercourse did not occur on this day then `W` is non-random and
	// the log ratio is 0
	if (! X[i]) {
	    continue;
	}

	double term1, term2;

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
	//
	// IMPORTANT: note that `w_days_idx` has a sentinal value appended to
	// the end of the data so that we need not worry about reading past the
	// end of the array

	if (*w_days_idx == i) {
	    term1 = *w_vals * m_Uh[i] * beta_diff;
	    ++w_vals;
	    ++w_days_idx;
	}
	else {
	    term1 = 0.0;
	}

	// calculate `-xi_i * [exp(U * beta*) - exp(U * beta)]`, which is one of
	// the terms in `p(W | proposal) / p(W | current)`.
	term2 = -xi_i * (exp(ubeta_vals[i] + (m_Uh[i] * beta_diff)) - ubeta_exp_vals[i]);

	// add the portion of the log-likelihood from the current day to the
	// running total
	sum_log_lik += term1 + term2;
    }

    return sum_log_lik;
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
//      log Gamma(x; a, b) = log(norm_const) + ((a - 1) * log(x)) - (b * x)
//
// plus one extra term for the log of the `p_h / (1 - p_h)` term or its inverse.
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

double GammaContMH::get_gam_log_lik(double proposal_beta, double proposal_gam) const {

    double gam_log_lik;

    if ((proposal_gam == 1.0) && (m_gam_val == 1.0)) {
	gam_log_lik = 0.0;
    }
    else if ((proposal_gam == 1.0) && (m_gam_val != 1.0)) {
	gam_log_lik = (m_log_p_over_1_minus_p
		       - m_log_norm_const
		       - ((m_hyp_a - 1) * m_beta_val)
		       + (m_hyp_b * m_gam_val));
    }
    else if ((proposal_gam != 1.0) && (m_gam_val == 1.0)) {
	gam_log_lik = (m_log_1_minus_p_over_p
		       + m_log_norm_const
		       + ((m_hyp_a - 1) * proposal_beta)
		       - (m_hyp_b * proposal_gam));
    }
    else {
	const double beta_diff = proposal_beta - m_beta_val;
	const double gam_diff  = proposal_gam - m_gam_val;
	gam_log_lik = ((m_hyp_a - 1) * beta_diff) - (m_hyp_b * (gam_diff));
    }

    return gam_log_lik;
}




// calculate J(beta^(s) | beta*) / J(beta* | beta^(s)).

double GammaContMH::get_proposal_log_lik(double proposal_beta) const {

    // case: both proposal and current sampled 0 or both sampled from the
    // continuous part of the distribution.  In either case the ratio is J(beta^(s) | beta*)1 (the
    // latter case is b/c the proposal distributions are symmetric) so that the
    // log is 0.
    if (((proposal_beta != 0.0) && (m_beta_val != 0.0))
	|| ((proposal_beta == 0.0) && (m_beta_val == 0.0))) {

	return 0.0;
    }

    // case: proposal is 0 and current is continuous.  This gives us (for `J`
    // the proposal distribution)
    //
    //         (1 - pi) * J(beta^(s) | beta*)
    //     log ------------------------------
    //                     pi
    //
    else if (proposal_beta == 0.0) {

	double log_dgamma_curr = m_log_proposal_den(m_beta_val, proposal_beta, m_mh_delta);

	return m_mh_log_1_minus_p + log_dgamma_curr - m_mh_log_p;
    }

    // case: proposal is continuous and current is 0.  This gives us (for `J`
    // the proposal distribution)
    //
    //                      pi
    //     log ------------------------------
    //         (1 - pi) * J(beta* | beta^(s))
    //
    else {

	double log_dgamma_proposal = m_log_proposal_den(proposal_beta, m_beta_val, m_mh_delta);

	return m_mh_log_p - m_mh_log_1_minus_p - log_dgamma_proposal;
    }
}




// calculates the log norming constant for a possibly truncated Gamma(a, b)
// distribution, given by
//
//                           1                        b^a
//     log  -------------------------------------  ----------
//          int_{bnd_l}^{bnd_u} Gamma(x; a, b) dx  gammafn(a)

double GammaContMH::log_dgamma_trunc_norm_const() const {

    // if the upper bound is infinity then F(infinity) = 1
    double F_upp = (m_bnd_u == R_PosInf) ?
	1.0 :
	R::pgamma(m_bnd_u, m_hyp_a, 1.0 / m_hyp_b, 1, 0);

    // if the lower bound is 0 then F(0) = 0
    double F_low = (m_bnd_l == 0.0) ?
	0.0 :
	R::pgamma(m_bnd_l, m_hyp_a, 1.0 / m_hyp_b, 1, 0);

    return (m_hyp_a * log(m_hyp_b)) - R::lgammafn_sign(m_hyp_a, NULL) - log(F_upp - F_low);
}




// double (*(*get_proposal_fcn)(double cond, double delta))(int proposal_code) {

//     // initialize the appropriate subclass of `GammaGen`, as specified by
//     // `curr_gamma_specs`
//     switch((int) proposal_code) {
//     case UNIF:
// 	// TODO
// 	break;
//     }
// }
