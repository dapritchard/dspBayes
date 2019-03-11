#ifndef DSP_BAYES_SRC_COEF_GEN_H
#define DSP_BAYES_SRC_COEF_GEN_H

#include "Rcpp.h"
#include "GammaGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "FWPriors.h"


class CoefGen {

public:

    GammaGen** m_gamma;

    Rcpp::NumericVector m_vals_rcpp;
    Rcpp::NumericVector::iterator m_vals;

    const int m_n_psi;
    const int m_n_gamma;

    CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp);
    ~CoefGen();

    void sample(const WGen& W,
                const XiGen& xi,
                UProdBeta& ubeta,
                const int* X,
                const FWPriors& fw_priors);

    const double* vals() const { return m_vals; }
};


#endif
