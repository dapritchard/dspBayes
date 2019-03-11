#include "Rcpp.h"

// TODO: missing header files
#include "GammaGen.h"
#include "global_vars.h"
#include "FWPriors.h"

using R::lgammafn_sign;




GammaCateg::GammaCateg(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& gamma_specs) :
    GammaGen(U, gamma_specs),
    m_bnd_l_is_zero(m_bnd_l == 0.0),
    m_bnd_u_is_inf(m_bnd_u == R_PosInf),
    m_is_trunc(!m_bnd_l_is_zero || !m_bnd_u_is_inf),
    m_incl_one(m_hyp_p != 0.0),
    m_log_d2_const_terms(calc_log_d2_const_terms()) {
}




// sample a new value for gamma_h.  As a side-effect, updates `ubeta`,
// `m_beta_val`, and `m_gam_val` to reflect the newly sampled value of gamma_h.

double GammaCateg::sample(const WGen& W,
                          const XiGen& xi,
                          UProdBeta& ubeta,
                          const int* X,
                          const FWPriors& fw_priors) {

    double a_tilde, b_tilde, p_tilde;

    // calculate a_tilde
    a_tilde = calc_a_tilde(W);

    // calculate b_tilde and p_tilde.  Note that `calc_b_tilde(ubeta)` has the
    // side-effect of changing the values of the data pointed to by `ubeta` to
    // instead have the values given by `U * beta - U_h * beta_h`.
    b_tilde = calc_b_tilde(ubeta, xi, X);
    p_tilde = m_incl_one ? calc_p_tilde(a_tilde, b_tilde) : 0;

    // sample a new value for gamma_h and update beta_h
    m_gam_val = sample_gamma(a_tilde, b_tilde, p_tilde);
    m_beta_val = log(m_gam_val);

    // change the values of the data pointed to by `ubeta` to take the values of
    // `U * beta` using the newly sampled value of `gamma_h`
    // TODO: check that gamma isn't 1 before calling
    ubeta.add_uh_prod_beta_h(m_Uh, m_beta_val);

    return m_gam_val;
}




double GammaCateg::sample_gamma(double a_tilde, double b_tilde, double p_tilde) {

    double unif_bnd_l, unif_bnd_u, unif_rv;

    // case: with probability `p_tilde`, sample a value of 1
    if (R::unif_rand() < p_tilde) {
        return 1;
    }

    // case: with probability `1 - p_tilde`, sample from a possibly truncated
    // gamma distribution
    else {

        // if the distribution is not truncated then use the usual routine to
        // sample a gamma variate
        if (! m_is_trunc) {
            return R::rgamma(a_tilde, 1 / b_tilde);
        }
        // case: sample from a truncated gamma distribution
        else {

            // calculate `F(m_bnd_l; a_tilde, b_tilde)` and `F(m_bnd_u; a_tilde,
            // b_tilde)`
            unif_bnd_l = m_bnd_l_is_zero ?
                0.0 :
                R::pgamma(m_bnd_l, a_tilde, 1.0 / b_tilde, 1, 0);
            unif_bnd_u = m_bnd_u_is_inf ?
                1.0 :
                R::pgamma(m_bnd_u, a_tilde, 1.0 / b_tilde, 1, 0);

            // sample a uniform r.v. and return the `unif_rv`-th quantile from
            // the gamma distribution
            unif_rv = R::runif(unif_bnd_l, unif_bnd_u);
            return R::qgamma(unif_rv, a_tilde, 1.0 / b_tilde, 1, 0);
        }
    }
}





// calculate the value of a_h tilde for fixed h, which is defined as:
//
//     a_h + sum_ijk u_ijkh * W_ijk
//
// This requires merely stepping through the indices `ijk` that correspond to
// cycles in which a pregnancy occured and checking whether `u_ijkh` has a value
// of 1 for those indices.

double GammaCateg::calc_a_tilde(const WGen& W) {

    // a pointer to the values of W, the indices that the values correspond to,
    // and one past the end of the values of W
    const int* w_vals = W.vals();
    const int* w_days_idx = W.days_idx();
    const int* w_end = w_days_idx + W.n_preg_days();

    // initialize `sum_val` to the first term in the sum
    double sum_val = m_hyp_a;

    // each iteration adds `u_ijkh * W_ijk` to `sum_val` for the indices `ijk`
    // that occur during cycles in which a pregnancy occured (otherwise `W_ijk`
    // is guaranteed to be 0)
    for ( ; w_days_idx < w_end; ++w_days_idx) {

        // case: `U_ijkh` has value of 1, so add `W_ijk` to the running total.
        // Otherwise `U_ijkh` has value of 0, so we can just ignore.
        if (m_Uh[*w_days_idx]) {
            sum_val += *w_vals;
        }

        ++w_vals;
    }

    return sum_val;
}




// calculate the value of b_h tilde for fixed h, which is defined as:
//
//     b_h + sum_{i,j,k: X_ijk == 1, U_ijkh == 1} exp(
//         log(xi_i) + sum_{l != h} U_{ijkl} * beta_l
//     )
//
// Furthermore, note that the s-th element of `sum_{l != h} U_{ijkl} * beta_l`
// is equal to `U * beta - U_h * beta_h`, a fact that is used in the
// calculations below.
//
// Nota bene: this function has the side effect of changing the values of the
// data pointed to by `U * beta` to instead have the values given by `U * beta -
// U_h * beta_h`.

double GammaCateg::calc_b_tilde(UProdBeta& ubeta, const XiGen& xi, const int* X) {

    double* ubeta_vals = ubeta.vals();
    const double* xi_vals = xi.vals();

    // initialize `sum_val` to take the first term in the expression
    double sum_val = m_hyp_b;

    // each iteration checks whether `r` corresponding to index ijk satisfies
    // the consitions of the outer sum, and if so, adds the value of the
    // expression inside of the outer sum to the running total
    for (int r = 0; r < m_n_days; ++r) {

        // case: U_{ijkh} has a value of 1, so update the ijk-th element of
        // `ubeta_vals` to have the value of the ijk-th element of `U*beta -
        // U_h*beta_h`.  If U_{ijkh} has a value of 0 then no update is needed.
        // TODO: check that gamma_h isn't 1
        if (m_Uh[r]) {

            ubeta_vals[r] -= m_beta_val;

            // case: the index ijk that `r` corresponds to is one of the terms
            // included in the outer sum, so add the value of the expression to
            // the running total
            if (X[r]) {
                sum_val += exp(log(xi_vals[ d2s[r] ]) + ubeta_vals[r]);
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

    // log( int{ G(x; a_tilde, b_tilde) }dx ).  Note that this can become -inf,
    // but the behavior of this case propagates through the function as desired,
    // so no special considerations need to be taken.
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

    return m_hyp_p / (m_hyp_p + exp(log_d2));
}




// the logarithm of the normalizing term from a gamma density with parameters
// `a` and `b`, and which is given by:
//
//     log(b^a / gamma(a)) = a * log(b) - log( gamma(b) )

double GammaCateg::log_dgamma_norm_const(double a, double b) {
    return a * log(b) - lgammafn_sign(a, NULL);
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
        1.0 :
        R::pgamma(m_bnd_u, a, 1.0 / b, 1, 0);

    // if the lower bound  is 0 then F(0) = 0
    F_low = (m_bnd_l_is_zero) ?
        0.0 :
        R::pgamma(m_bnd_l, a, 1.0 / b, 1, 0);

    return log(F_upp - F_low);
}




// calculates the log of the constant terms in the expression for d2 as defined
// in `calc_p_tilde`, as given by the expression:
//
//         (1 - p) * C(a, b)  * exp(-b)
//     log ----------------------------
//              int{ G(x; a, b) }dx
//
// and where C(a, b) is as defined by the `log_dgamma_norm_const` function.

double GammaCateg::calc_log_d2_const_terms() {

    return (log(1 - m_hyp_p)
            + log_dgamma_norm_const(m_hyp_a, m_hyp_b)
            - m_hyp_b
            - log_dgamma_trunc_const(m_hyp_a, m_hyp_b));
}
