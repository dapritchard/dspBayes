#ifndef DSP_BAYES_UTEST_FACTORY_H
#define DSP_BAYES_UTEST_FACTORY_H

#include "Rcpp.h"
#include "GammaGen.h"
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
		 Rcpp::IntegerVector w_cyc_to_subj_idx,
		 Rcpp::List subj_day_blocks,
		 Rcpp::IntegerVector day_to_subj_idx,
		 Rcpp::List gamma_specs,
		 Rcpp::NumericVector phi_specs,
		 int fw_len,
		 int n_burn,
		 int n_samp,
		 Rcpp::List test_data);

    // factory methods
    GammaCateg* gamma_categ_all();
    GammaCateg* gamma_categ_zero_one();
    GammaCateg* gamma_categ_one_inf();
    GammaCateg* gamma_categ_zero_half();
    XiGen* xi();
    XiGen* xi_no_rec();
    WGen* W();
    UProdBeta* ubeta();
    PhiGen* phi();
    PhiGen* phi_no_rec();
    int* X();
    static bool eq_dbl(double a, double b);

    // usual input
    Rcpp::NumericMatrix U;
    Rcpp::IntegerVector X_rcpp;
    Rcpp::List preg_cyc;
    Rcpp::IntegerVector w_to_days_idx;
    Rcpp::IntegerVector w_cyc_to_subj_idx;
    Rcpp::List subj_day_blocks;
    Rcpp::IntegerVector day_to_subj_idx;
    Rcpp::List gamma_specs;
    Rcpp::NumericVector phi_specs;
    int fw_len;
    int n_burn;
    int n_samp;

    // testing data
    Rcpp::List input_gamma_specs;
    Rcpp::NumericVector input_ubeta;
    Rcpp::NumericVector input_w;
    Rcpp::NumericVector input_xi;
    Rcpp::NumericVector target_data_gamma_categ;
    Rcpp::NumericVector target_data_phi;
    Rcpp::List target_samples_gamma_categ;
    Rcpp::NumericVector target_samples_phi;
    Rcpp::NumericVector target_samples_w;
    Rcpp::NumericVector target_samples_xi;

    // global testing objects
    Rcpp::IntegerVector seed_vals;
    static double epsilon;

    // derived data
    int n_days;
    int n_subj;
};


#endif
