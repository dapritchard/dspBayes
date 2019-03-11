#ifndef DSP_BAYES_UTEST_FACTORY_H
#define DSP_BAYES_UTEST_FACTORY_H

#include "Rcpp.h"

#include "CoefGen.h"
#include "GammaGen.h"
#include "UGenVar.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "XGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"


class UTestFactory {

public:

    UTestFactory() {}

    // constructor.  Copy data to class members.
    UTestFactory(Rcpp::NumericMatrix u_rcpp,
                 Rcpp::IntegerVector x_rcpp,
                 Rcpp::List          w_day_blocks,
                 Rcpp::IntegerVector w_to_days_idx,
                 Rcpp::IntegerVector w_cyc_to_subj_idx,
                 Rcpp::List          subj_day_blocks,
                 Rcpp::IntegerVector day_to_subj_idx,
                 Rcpp::List          gamma_specs,
                 Rcpp::NumericVector phi_specs,
                 Rcpp::List          x_miss_cyc,
                 Rcpp::List          x_miss_day,
                 Rcpp::NumericVector utau_rcpp,
                 Rcpp::List          tau_coefs,
                 Rcpp::List          u_miss_info,
                 Rcpp::IntegerVector u_miss_type,
                 Rcpp::IntegerVector u_preg_map,
                 Rcpp::IntegerVector u_sex_map,
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
    UProdTau* utau();
    PhiGen* phi();
    PhiGen* phi_no_rec();
    XGen* X();
    int** X_temp();
    UGenVarCateg* u_categ(Rcpp::NumericMatrix* u_rcpp_copy);
    CoefGen* coefs();
    XGen::XMissDay** XMissDay();
    static bool eq_dbl(double a, double b);

    // usual input
    Rcpp::NumericMatrix u_rcpp;
    Rcpp::IntegerVector x_rcpp;
    Rcpp::List          preg_cyc;
    Rcpp::IntegerVector w_to_days_idx;
    Rcpp::IntegerVector w_cyc_to_subj_idx;
    Rcpp::List          subj_day_blocks;
    Rcpp::IntegerVector day_to_subj_idx;
    Rcpp::List          gamma_specs;
    Rcpp::NumericVector phi_specs;
    Rcpp::List          x_miss_cyc;
    Rcpp::List          x_miss_day;
    Rcpp::NumericVector utau_rcpp;
    Rcpp::List          tau_coefs;
    Rcpp::List          u_miss_info;
    Rcpp::IntegerVector u_miss_type;
    Rcpp::IntegerVector u_preg_map;
    Rcpp::IntegerVector u_sex_map;
    int fw_len;
    int n_burn;
    int n_samp;

    // testing data
    Rcpp::NumericVector input_gam_coefs;
    Rcpp::List          input_gamma_specs;
    Rcpp::List          input_u_categ;
    Rcpp::NumericVector input_ubeta;
    Rcpp::NumericVector input_w;
    Rcpp::NumericVector input_x;
    Rcpp::NumericVector input_xi;
    Rcpp::NumericVector target_data_gamma_categ;
    Rcpp::NumericVector target_data_phi;
    Rcpp::IntegerVector target_data_u_categ;
    Rcpp::List          target_samples_gamma_categ;
    Rcpp::List          target_samples_u_categ;
    Rcpp::NumericVector target_samples_phi;
    Rcpp::NumericVector target_samples_w;
    Rcpp::List          target_samples_x;
    Rcpp::NumericVector target_samples_xi;

    // global testing objects
    Rcpp::IntegerVector seed_vals;
    static double epsilon;

    // derived data
    int n_days;
    int n_subj;
};


#endif
