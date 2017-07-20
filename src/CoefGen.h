#ifndef DSP_BAYES_SRC_COEF_GEN_H
#define DSP_BAYES_SRC_COEF_GEN_H

#include "Rcpp.h"
#include "GammaGen.h"
#include "XiGen.h"
#include "UProdBeta.h"


class CoefGen {

public:

    GammaGen** m_gamma;

    double* m_vals;
    double* const m_output_start;
    double* const m_output_end;

    const int m_n_psi;
    const int m_n_gamma;

    // tracks whether we are past the burn-in phase
    bool m_record_status;

    CoefGen(Rcpp::NumericMatrix& U, Rcpp::List& gamma_specs, int n_samp);
    ~CoefGen();

    void sample(const WGen& W, const XiGen& xi, UProdBeta& u_prod_beta, const int* X);

    double* output_start() const { return m_output_start; }
    double* output_end() const { return m_output_end; }
};


#endif
