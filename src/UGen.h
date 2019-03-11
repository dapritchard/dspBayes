#ifndef DSP_BAYES_SRC_U_GEN_H
#define DSP_BAYES_SRC_U_GEN_H

#include "Rcpp.h"

#include "CoefGen.h"
#include "UGenVar.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"




class UGen {

public:

    UGenVar** m_vars;
    int m_n_vars;

    bool m_record_status;

    UGen(Rcpp::NumericMatrix& u_rcpp,
         Rcpp::List& miss_info,
         Rcpp::IntegerVector& miss_type,
         Rcpp::IntegerVector& preg_map,
         Rcpp::IntegerVector& sex_map,
         bool record_status);
    ~UGen();

    void sample(const WGen& W,
                const XiGen& xi,
                const CoefGen& coefs,
                const XGen& X,
                UProdBeta& ubeta,
                UProdTau& utau);

    Rcpp::List realized_samples();
};


#endif
