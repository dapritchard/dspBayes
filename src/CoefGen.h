#ifndef DSP_BAYES_SRC_COEF_GEN_H
#define DSP_BAYES_SRC_COEF_GEN_H

#include "Rcpp.h"
#include "GammaGen.h"
#include "XiGen.h"
#include "UProdBeta.h"


class CoefGen {

public:

    GammaGen** m_gamma;

    Rcpp::NumericVector m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    const int m_n_psi;
    const int m_n_gamma;

    const double* m_dsp_ar_vals;
    const int m_n_dsp_ar_vals;

    double m_mu_0;
    double m_nu_0;
    int m_mu_accept_ctr;
    int m_nu_accept_ctr;
    const double m_mu_hyp_a;
    const double m_mu_hyp_b;
    const double m_nu_hyp_a;
    const double m_nu_hyp_b;
    const double m_delta;
    double (*const m_proposal_fcn)(double val, double delta);

    CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp);
    ~CoefGen();

    // TODO: change to XGen
    void sample(const WGen& W, const XiGen& xi, UProdBeta& ubeta, const int* X);
    void sample_mu();
    void sample_nu();
    static double calc_nu_term_gamma_1(double proposal_nu, double nu_diff);
    static double calc_log_nu_term_gamma_remain(double proposal_nu, double nu_diff);
    static double calc_log_nu_term_gamma_prior(double proposal_nu, double nu_diff);

    const double* vals() const { return m_vals; }
};


#endif
