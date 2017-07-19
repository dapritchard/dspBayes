#include <cmath>
#include "Rcpp.h"
#include "PhiGen.h"
#include "XiGen.h"
#include "ProposalFcns.h"

using std::log;
using R::lgammafn;




PhiGen::PhiGen(Rcpp::NumericVector phi_hyper, int n_samp) :
    // initialization list
    m_hyp_c1(phi_hyper["c1"]),
    m_hyp_c2(phi_hyper["c2"]),
    m_delta(phi_hyper["delta"]),
    m_phi_val(phi_hyper["mean"]),
    m_vals(new double[n_samp]),
    m_output_start(m_vals),
    m_output_end(m_output_start + n_samp),
    m_accept_ctr(0),
    m_is_same_as_prev(false),
    m_log_norm_const(0) {
}




void PhiGen::sample(const XiGen& xi) {

    double proposal_val, log_r;

    // sample the proposal value for Metropolis step
    proposal_val = ProposalFcns::abs_unif(m_phi_val, m_delta);

    // calculate `log(r)` where `r` is the acceptance ratio for the Metropolis
    // step
    log_r = calc_log_r(xi, proposal_val);

    // sample the updated value of phi by either accepting the proposal value or
    // by keeping the current value
    m_phi_val = *m_vals++ = update_phi(log_r, proposal_val);
}




// calculate `log(r)` where `r` is the ratio of the posterior likelihood of the
// proposed value of phi divided by the posterior likelihood of the current
// value of phi.  This value is given by:
//
//         { prod_i pi(xi_i | phi^{*}) } * pi(phi^{*})
//     log -------------------------------------------
//         { prod_i pi(xi_i | phi^(s)) } * pi(phi^(s))
//
// and where pi refers to the r.v.'s respective density functions, phi^{*}
// refers to the proposal value, and phi^(s) refers to the phi value for the
// (current) s-th scan.
//
// The calculation is subdivided into two smaller calculations, namely (i)
// taking the log of the ratio of the xi likelihoods, and then (ii) the log of
// the ratio of the phi likelihoods.

double PhiGen::calc_log_r(const XiGen& xi, double proposal_val) {

    double log_proportion_dgamma_xi, log_proportion_dgamma_phi;

    log_proportion_dgamma_xi = calc_log_proportion_dgamma_xi(xi, proposal_val);
    log_proportion_dgamma_phi = calc_log_proportion_dgamma_phi(xi, proposal_val);

    return log_proportion_dgamma_xi + log_proportion_dgamma_phi;
}




// randomly update the value of phi by either accepting the proposal value for
// phi or by keeping the current value.  The proposal value is accepted with
// probability min(1, r).

double PhiGen::update_phi(double log_r, double proposal_val) {

    if (log(R::unif_rand()) < log_r) {
	m_phi_val = proposal_val;
	m_is_same_as_prev = false;
	++m_accept_ctr;
    }
    else {
	m_is_same_as_prev = true;
    }

    return m_phi_val;
}




// calculate the log likelihood ratio for the xi terms given by:
//
//         { prod_i pi(xi_i | phi^{*}) }
//     log -----------------------------
//         { prod_i pi(xi_i | phi^(s)) }
//
//
//       =  n * phi^{*} * log(phi^{*}) - log(Gamma(phi^{*})) }
//              - n * { phi^(s) * log(phi^(s)) - log(Gamma(phi^(s))) }
//              + (phi^{*} - phi^(s)) * sum_i { log(xi_i) - xi_i }
//
// and where n is the number of terms in the product.

double PhiGen::calc_log_proportion_dgamma_xi(const XiGen& xi, double proposal_val) {

    double curr_xi, numer_log_norm_const, denom_log_norm_const, log_kernel_ratio;

    const double* xi_vals = xi.vals();
    int n_subj = xi.n_subj();

    // calculate the log of the normalizing constant for the proposal value
    numer_log_norm_const = n_subj * log_dgamma_norm_const(proposal_val);

    // calculate the log of the normalizing constant for the current value.  If
    // the current value is the same as the previous (i.e. the proposal value
    // was not accepted), then the term need not be calculated again.
    if (m_is_same_as_prev) {
	denom_log_norm_const = m_log_norm_const;
    } else {
	denom_log_norm_const = n_subj * log_dgamma_norm_const(m_phi_val);
	m_log_norm_const = denom_log_norm_const;
    }

    // each iteration adds the kernel log likelihood ratio of the current
    // subject to `log_kernel_ratio`
    log_kernel_ratio = 0;
    for (int i = 0; i < n_subj; ++i) {
	curr_xi = xi_vals[i];
	log_kernel_ratio += log(curr_xi) - curr_xi;
    }
    // mutiply in the `phi^{*} - phi^(s)` term
    log_kernel_ratio *= proposal_val - m_phi_val;

    return numer_log_norm_const - denom_log_norm_const + log_kernel_ratio;
}




// calculate the log likelihood ratio for the phi terms given by:
//
//         pi(phi^{*})
//     log -----------
//         pi(phi^(s))
//
//       =  (c_1 - 1) * log(phi^{*} / phi^(s)) - c2 * (phi^{*} - phi^(s))
//
// and where c_1 and c_2 are the hyperparameters for the gamma distribution for
// phi.

double PhiGen::calc_log_proportion_dgamma_phi(const XiGen& xi, double proposal_val) {

    double term1, term2;

    // (c_1 - 1) * log(phi^{*} / phi^(s))
    term1 = (m_hyp_c1 - 1) * log(proposal_val / m_phi_val);
    // c2 * (phi^{*} - phi^(s))
    term2 = m_hyp_c2 * (proposal_val - m_phi_val);

    return term1 - term2;
}




// the logarithm of the normalizing term from a gamma density with parameters
// `a` and `a`, and which is given by:
//
//     log(a^a / gamma(a)) = a * log(a) - log( gamma(a) )

double PhiGen::log_dgamma_norm_const(double a) {
    return a * log(a) + lgammafn(a);
}
