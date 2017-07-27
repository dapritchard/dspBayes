#ifndef DSP_BAYES_UTEST_FACTORY_H
#define DSP_BAYES_UTEST_FACTORY_H


#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"


class UTestFactory {

public:

    UTestFactory() {}

    // constructor.  Copy data to class members.
    UTestFactory(Rcpp::NumericMatrix U,
		 Rcpp::IntegerVector X_rcpp,
		 Rcpp::List w_day_blocks,
		 Rcpp::IntegerVector w_to_days_idx,
		 Rcpp::IntegerVector w_cyc_to_cyc_idx,
		 Rcpp::List subj_day_blocks,
		 Rcpp::IntegerVector day_to_subj_idx,
		 Rcpp::List gamma_specs,
		 Rcpp::NumericVector phi_specs,
		 int fw_len,
		 int n_burn,
		 int n_samp,
		 Rcpp::List test_data);

    // factory methods
    XiGen* xi();
    XiGen* xi_no_rec();
    WGen* W();
    UProdBeta* ubeta();
    PhiGen* phi();
    PhiGen* phi_no_rec();

    // class members
    Rcpp::NumericMatrix U;
    Rcpp::IntegerVector X_rcpp;
    Rcpp::List preg_cyc;
    Rcpp::IntegerVector w_to_days_idx;
    Rcpp::IntegerVector w_cyc_to_cyc_idx;
    Rcpp::List subj_day_blocks;
    Rcpp::IntegerVector day_to_subj_idx;
    Rcpp::List gamma_specs;
    Rcpp::NumericVector phi_specs;
    int fw_len;
    int n_burn;
    int n_samp;
    Rcpp::NumericVector input_xi;
    Rcpp::NumericVector input_w;
    Rcpp::NumericVector input_ubeta;
    Rcpp::NumericVector target_samples_xi;
    Rcpp::NumericVector target_data_phi;
    Rcpp::NumericVector target_samples_phi;

    int n_days;
    int n_subj;
    double phi_init;
};


#endif
