#ifndef DSP_BAYES_SRC_COEF_GEN_H
#define DSP_BAYES_SRC_COEF_GEN_H

#include "GammaGen.h"


class CoefGen {

    const GammaGen** gamma;

    const int m_n_psi;
    const int m_n_gamma;

    void sample();
};


#endif
