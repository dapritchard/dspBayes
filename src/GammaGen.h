#ifndef DSP_BAYES_GAMMA_GEN_H_
#define DSP_BAYES_GAMMA_GEN_H_

#include "Rcpp.h"




class GammaGen {

public:

    // current value of gamma_h
    double m_beta_val;
    double m_gam_val;

    // gamma_h hyperparameters
    const double m_hyp_a;
    const double m_hyp_b;
    const double m_hyp_p;
    const double m_bnd_l
    const double m_bnd_u;

    // number of observations in the data
    const int m_n_obs;

    // points to the beginning and one past the end of the data for U_h
    const double* m_Uh;
    // onst double* m_Uh_end;

    GammaGen(Rcpp::NumericMatrix& U, int h);

    // virtual double samp_gam(const std::vector<double>& U_prod_beta,
    // 			    const std::vector<double>& W);

    // virtual void set_ar_param(double prev_gam_val);


};




class GammaCateg : public GammaGen {

public:

    // whether h-th variable is truncated or not
    const bool m_bnd_l_is_zero;
    const bool m_bnd_u_is_inf;
    const bool m_is_trunc;

    // the log of the constant terms in the expression for d2 as defined in
    // `GammaCateg::calc_p_tilde`.
    const double m_log_d2_const_terms;

    double sample_gammma(double* U_prod_beta, const double* W, const double* xi);
    double calc_a_tilde(const double* W);
    double calc_b_tilde(double* U_prod_beta, const double* xi);
    double calc_p_tilde(double a_tilde, double b_tilde);
    void add_uh_prod_beta_h(double* U_prod_beta_no_h);
    double log_dgamma_norm_const(double a, double b);
    double log_dgamma_trunc_const(double a, double b);
    double init_log_d2_const_terms();
};




class GammaContAux : public GammaGen {

public:

    // // whether M == 1 (i.e no effect in the model)
    // bool M_is_one;


};




class GammaContMetr : public GammaGen {

public:

    // Metropolis-Hastings tuning parameter.  Specifies the probability of
    // selecting 1 as the proposal value.
    double mh_pi;

};




class GammaMeanShrink : public GammaGen {

public:



};


#endif
