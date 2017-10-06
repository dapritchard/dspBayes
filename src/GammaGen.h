#ifndef DSP_BAYES_GAMMA_GEN_H_
#define DSP_BAYES_GAMMA_GEN_H_

#include "Rcpp.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"




class GammaGen {

public:

    // current value of beta_h and gamma_h
    double m_beta_val;
    double m_gam_val;

    // gamma_h hyperparameters
    const double m_hyp_a;
    const double m_hyp_b;
    const double m_hyp_p;
    const double m_bnd_l;
    const double m_bnd_u;

    // points to the beginning of the data for U_h
    const double* m_Uh;

    // number of observations in the data
    const int m_n_days;

    GammaGen(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& coef_specs);
    virtual ~GammaGen() {}

    // TODO: change this to XGen& X
    virtual double sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta, const int* X) = 0;

    static GammaGen** create_arr(const Rcpp::NumericMatrix& U, const Rcpp::List& gamma_specs);
};




class GammaCateg : public GammaGen {

public:

    // whether h-th variable is truncated or not
    const bool m_bnd_l_is_zero;
    const bool m_bnd_u_is_inf;
    const bool m_is_trunc;
    const bool m_incl_one;

    // the log of the constant terms in the expression for d2 as defined in
    // `GammaCateg::calc_p_tilde`.
    const double m_log_d2_const_terms;


    GammaCateg(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& gamma_specs);

    double sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta, const int* X);
    double calc_a_tilde(const WGen& W);
    double calc_b_tilde(UProdBeta& u_prod_beta, const XiGen& xi, const int* X);
    double calc_p_tilde(double a_tilde, double b_tilde);
    double sample_gamma(double a_tilde, double b_tilde, double p_tilde);
    double log_dgamma_norm_const(double a, double b);
    double log_dgamma_trunc_const(double a, double b);
    double init_log_d2_const_terms();
    double calc_log_d2_const_terms();
};




/* class GammaContAux : public GammaGen { */

/* public: */

/*     // // whether M == 1 (i.e no effect in the model) */
/*     // bool M_is_one; */


/* }; */




class GammaContMH : public GammaGen {

public:

    // precomputed calculations:
    //
    //   m_log_norm_const:  log truncation const times the gamma dist const
    //
    //   m_log_p_over_1_minus_p:  log(p_h / (1 - p_h))
    //
    //   m_log_1_minus_p_over_ph:  log((1 - p_h) / p_h)
    //
    const double m_log_norm_const;
    const double m_log_p_over_1_minus_p;
    const double m_log_1_minus_p_over_p;

    // Metropolis-Hastings variables:
    //
    //   m_mh_prob_samp_1:  P(select 1 as the proposal gamma)
    //
    //   m_mh_log_prob_samp_1:  log(m_mh_prob_samp_1)
    //
    //   m_mh_log_1_minus_prob_samp_1:  log(1 - m_mh_prob_samp_1)
    //
    //   m_mh_delta:  tuning parameter passed to proposal generating fcn
    //
    //   m_mh_accept_ctr: number of times the proposal value was accepted
    //
    const double m_mh_p;
    const double m_mh_log_p;
    const double m_mh_log_1_minus_p;
    const double m_mh_delta;
    int m_mh_accept_ctr;

    // the continuous part of the proposal distribution and the corresponding
    // density function
    double (*m_proposal_fcn)(double cond, double delta);
    double (*m_log_proposal_den)(double val, double cond, double delta);

    GammaContMH(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& gamma_specs);
    double sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta, const int* X);
    double sample_proposal_beta() const;
    double get_log_r(const WGen& W,
		     const XiGen& xi,
		     const UProdBeta& ubeta,
		     double proposal_beta,
		     double proposal_gam);
    double get_w_log_lik(const WGen& W,
			 const XiGen& xi,
			 const UProdBeta& ubeta,
			 double proposal_beta) const;
    double get_gam_log_lik(double proposal_beta, double proposal_gam) const;
    double get_proposal_log_lik(double proposal_beta) const;
    double log_dgamma_trunc_norm_const() const;
};




/* class GammaMeanShrink : public GammaGen { */

/* public: */



/* }; */


#endif
