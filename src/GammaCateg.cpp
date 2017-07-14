#include <cmath>  // log, exp, lgamma
#include "GammaGen.h"
#include "global_vars.h"

using std::log;
using std::exp;
using std::lgamma;




GammaCateg::GammaCateg() :
    m_bnd_l_is_zero(m_bnd_l == 0),
    m_bnd_u_is_inf(m_bnd_u == R_Posinf),
    m_is_trunc(!m_bnd_l_is_zero || !m_bnd_u_is_inf) {}




// sample a new value for gamma_h.  As a side-effect, updates `U_prod_beta`,
// `m_beta_val`, and `m_gam_val` to reflect the newly sampled value of
// gamma_h.

double GammaCateg::sample_gammma(double* U_prod_beta, const double* W, const double* xi) {

    double a_tilde, b_tilde, p_tilde;

    // calculate a_tilde
    a_tilde = calc_a_tilde(W);

    // calculate b_tilde and p_tilde.  Note that `calc_b_tilde(U_prod_beta)` has
    // the side-effect of changing the values of the data pointed to by
    // `U_prod_beta` to instead have the values given by `U * beta - U_h *
    // beta_h`.
    b_tilde = calc_b_tilde(U_prod_beta);
    p_tilde = calc_p_tilde(a_tilde, b_tilde);

    // sample a new value for gamma_h and update beta_h
    m_gam_val = sample_gam(a_tilde, b_tilde, p_tilde);
    m_beta_val = log(m_gam_val);

    // change the values of the data pointed to by `U_prod_beta` to take the
    // values of `U * beta` using the newly sampled value of `gamma_h`
    add_uh_prod_beta_h(U_prod_beta);

    return m_gam_val;
}




// calculate the value of a_h tilde for fixed h, which is defined as:
//
//     a_h + sum_ijk u_ijkh * W_ijk
//

// ******* rework.  most W_ijk are 0.  *************

double GammaCateg::calc_a_tilde(const WGen W) {

    double* w_vals, *w_days_idx;
    int curr_w_day;
    double sum_val;

    w_vals = W.vals();
    w_days_idx = W.days_idx();
    sum_val = m_hyp_a;

    // each iteration adds `u_ijkh * W_ijk` to `sum_val`
    for (int r = 0; r < m_n_obs; ++r) {

	// case: the r-th day corresponds to a random (i.e. possibly nonzero)
	// `W_ijk`, so add the value of `u_ijkh * W_ijk` to `sum_val`
	if (r == *w_days_idx) {

	    // case: `U_ijkh` has value of 1, so add `W_ijk` to the running
	    // total.  Otherwise `U_ijkh` has value of 0, so we can just ignore.
	    if (m_Uh[r]) {
		sum_val += *w_vals++;
	    }

	    ++w_days_idx;
	}
    }

    return sum_val;
}




// calculate the value of b_h tilde for fixed h, which is defined as:
//
//     b_h + sum_{i,j,k: X_ijk == 1, U_ijkh == 1} exp(
//         log(xi_i) + sum_{l != h} U_{ijkl} * beta_l
//     )
//
// Furthermore, note that the s-th element of `sum_{l != h} U_{ijkl} * beta_l` is equal to
// `U*beta - U_h*beta_h`, a fact that is used in the calculations below.
//
// ** important **: this function has the side effect of changing the values of
// the data pointed to by `U_prod_beta` to instead have the values given by `U *
// beta - U_h * beta_h`.

double GammaCateg::calc_b_tilde(double* U_prod_beta, const double* xi) {

    // initialize `sum_val` to take the first term in the expression
    double sum_val = m_hyp_b;

    // each iteration checks whether `r` corresponding to index ijk satisfies
    // the consitions of the outer sum, and if so, adds the value of the expression
    // inside of the outer sum to the running total
    for (int r = 0; r < m_nobs; ++r) {

	// case: U_{ijkh} has a value of 1, so update the ijk-th element of
	// `U_prod_beta` to have the value of the ijk-th element of `U*beta -
	// U_h*beta_h`.  If U_{ijkh} has a value of 0 then no update is needed.
	if (m_Uh[r]) {

	    U_prod_beta[r] -= m_beta_val;

	    // case: the index ijk that `r` corresponds to is one of the terms
	    // included in the outer sum, so add the value of the expression to
	    // the running total
	    if (X[r]) {
		sum_val += exp(log(xi[ d2s[r] ]) + U_prod_beta[r]);
	    }
	}
    }

    return sum_val;
}




// Calculate p_h tilde, where p_h tilde is defined as p / (p + d2), and where d2
// is given by:
//
//
//     (1 - p) * C(a, b) * int{ G(x; a_tilde, b_tilde) }dx * exp(b_tilde - b)
//     ----------------------------------------------------------------------
//                   C(a_tilde, b_tilde) * int{ G(x; a, b) }dx
//
//
// Calculations for d2 are performed on the log scale.  Since the only terms
// that change across MCMC samples are expressions containing `a_tilde` or
// `b_tilde`, we collect the remaining terms and store them as
// `m_log_d2_const_terms`.

double GammaCateg::calc_p_tilde(double a_tilde, double b_tilde) {

    // log( int{ G(x; a_tilde, b_tilde) }dx )
    double log_trunc_const_tilde = m_is_trunc ?
	log_dgamma_trunc_const(a_tilde, b_tilde) :
	0;

    // log( C(a_tilde, b_tilde) )
    double log_norm_const_tilde = log_dgamma_norm_const(a_tilde, b_tilde);

    // `log(d2)` for d2 defined as above
    double log_d2 = (m_log_d2_const_terms
		     + log_trunc_const_tilde
		     + b_tilde
		     - log_norm_const_tilde);

    return p / (p + exp(log_d2));
}




// change the values of the data pointed to by `U_prod_beta_no_h` so that each
// element has the value of the corresponding element of `U_h * beta_h* added to
// it

void add_uh_prod_beta_h(double* U_prod_beta_no_h) {

    // each iteration updates the ijk-th element of `U * beta - U_h * beta_h` to
    // instead have the values given by `U_prod_beta`.
    for (int r = 0; r < m_nobs; ++r) {

	// case: U_{ijkh} has a value of 1, so update the ijk-th element of
	// `U*beta - U_h*beta_h` to have the value of the ijk-th element of
	// `U_prod_beta`.  If U_{ijkh} has a value of 0 then no update is
	// needed.
	if (m_Uh[r]) {
	    U_prod_beta_no_h[r] += m_beta_val;
	}
    }
}




// the logarithm of the normalizing term from a gamma density with parameters
// `a` and `b`, and which is given by:
//
//     log(b^a / gamma(a)) = a * log(b) - log( gamma(b) )

double GammaCateg::log_dgamma_norm_const(double a, double b) {
    return a * log(b) + lgamma(b);
}




// calculates the log of the inverse of the norming constant resulting from
// possibly truncating a Gamma(a, b) distribution.  In more detail, if we define
//
//     C = int_{bnd_l}^{bnd_u} Gamma(x; a, b) dx
//
// then the function returns log(C).  It is the inverse because we divide by
// this value to normalize the truncated distribution.  However, rather than
// directly calculating the integral, we instead appeal to the fundamental
// theorem of calculus and calculate `log(F(bnd_u) - F(bnd_l))`, where F is the
// distribution function of Gamma(x; a, b).

double GammaCateg::log_dgamma_trunc_const(double a, double b) {

    double F_low, F_upp;

    // if the upper bound is infinity then F(infinity) = 1
    F_upp = (m_bnd_u_is_inf) ?
	1 :
	pgamma(m_bnd_u, a, 1/b, 0, 0);

    // if the lower bound  is 0 then F(0) = 0
    F_low = (m_bnd_l_is_zero) ?
	0 :
	pgamma(m_bnd_u, a, 1/b, 0, 0);

    return log(F_upp - F_low);
}




// calculates the log of the constant terms in the expression for d2 as defined
// in `GammaCateg::calc_p_tilde`, as given by the expression:
//
//         (1 - p) * C(a, b)  * exp(-b)
//     log ----------------------------
//              int{ G(x; a, b) }dx
//

double calc_log_d2_const_terms() {

    return (log(1 - m_hyp_p)
	    + log_dgamma_norm_const(m_hyp_a, m_hyp_b)
	    - m_hyp_b
	    - log_dgamma_trunc_const(m_hyp_a, m_hyp_b));
}
