#ifndef DSP_BAYES_SRC_COEF_GEN_H
#define DSP_BAYES_SRC_COEF_GEN_H

#include "GammaGen.h"


class CoefGen {

    GammaGen** m_gamma;

    const int m_n_psi;
    const int m_n_gamma;

    CoefGen(Rcpp::List& gamma_list);
    ~CoefGen();

    void sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta);
};


#endif
