#include <cmath>
#include "PhiGen.h"
#include "ProposalFcns.h"


double PhiGen::sample() {

    proposal_val = ProposalFcns::abs_unif(m_val, m_delta);

    log_r = calc_log_r(xi);

    return update_phi(log_r);
}






double PhiGen::calc_log_r(const XiGen& xi, double proposal_val) {

    double log_proportion_dgamma_xi, log_proportion_dgamma_phi;

    log_proportion_dgamma_xi = calc_log_proportion_dgamma_xi(xi, proposal_val);
    log_proportion_dgamma_phi = calc_log_proportion_dgamma_phi(xi, proposal_val);

    return log_proportion_dgamma_xi + log_proportion_dgamma_phi;
}




double PhiGen::update_phi(double log_r, double proposal_val) {

    if (log(R::unif_rand()) < log_r) {
	m_val = proposal_val;
	m_is_same_as_prev = false;
	++m_accept_ctr;
    }

    return m_val;
}




double calc_log_proportion_dgamma_xi(const XiGen& xi, double proposal_val) {

    double curr_xi, numer_log_norm_const, denom_log_norm_const, log_xi_power,
	log_exp_prod_xi, sum_val;
    double* xi_vals;
    int n;

    xi_vals = xi.vals();
    n = xi.size();
    sum_val = 0;

    phi_diff = proposal_val - m_val;


    //
    for (int i = 0; i < n; ++i) {

	curr_xi = xi_vals[i];

	//
	numer_log_norm_const = log_dgamma_norm_const(m_proposal);
	if (m_is_same_as_prev) {
	    denom_log_norm_const = m_prev_log_norm_const;
	} else {
	    m_prev_log_norm_const = denom_log_norm_const = log_dgamma_norm_const(m_val);
	}

	log_xi_power = phi_diff * log(curr_xi);

	log_exp_prod_xi = -phi_diff * curr_xi;

	sum_val += (numer_log_norm_const
		    - denom_log_norm_const
		    + log_xi_power
		    + log_exp_prod_xi);
    }

    return sum_val;
}




double calc_log_proportion_dgamma_phi(const XiGen& xi, double proposal_val) {

    double log_phi_power, log_exp_prod_phi;

    log_phi_power = (m_hyp_c1 - 1) * log(proposal_val - m_val);
    log_exp_prod_phi = -m_hyp_c2 * phi_diff;

    return log_phi_power + log_exp_prod_phi;
}




// the logarithm of the normalizing term from a gamma density with parameters
// `a` and `a`, and which is given by:
//
//     log(a^a / gamma(a)) = a * log(a) - log( gamma(a) )

double PhiGen::log_dgamma_norm_const(double a) {
    return a * log(a) + lgamma(a);
}
