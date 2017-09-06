#ifndef DSP_BAYES_UTEST_U_GEN_VAR_CATEG_H
#define DSP_BAYES_UTEST_U_GEN_VAR_CATEG_H

#include "Rcpp.h"
#include "cppunit/extensions/HelperMacros.h"

#include "CoefGen.h"
#include "UGenVar.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"


class UGenVarCategTest : public CppUnit::TestFixture {

public:

    UGenVarCategTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_sample();
    void test_sample_covariate();
    void test_calc_posterior_w();
    void test_calc_posterior_x();

    CPPUNIT_TEST_SUITE(UGenVarCategTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST(test_sample_covariate);
    CPPUNIT_TEST(test_calc_posterior_w);
    CPPUNIT_TEST(test_calc_posterior_x);
    CPPUNIT_TEST_SUITE_END();


private:

    UGenVarCateg* u_var;
    WGen* W;
    XiGen* xi;
    CoefGen* coefs;
    UProdBeta* ubeta;
    UProdTau* utau;
    XGen* X;

    Rcpp::NumericMatrix* u_rcpp_copy;
    double* utau_vals_copy;

    // testing data
    int var_idx;
    Rcpp::NumericMatrix u_rcpp;
    int n_days;
    Rcpp::IntegerVector w_idx;
    Rcpp::IntegerVector x_idx;
    int col_start;
    int col_end;
    int ref_col;
    int n_categs;
    int max_n_days_miss;
    int max_n_sex_days_miss;
    Rcpp::NumericVector u_prior_probs;
    Rcpp::NumericVector input_w_probs;
    Rcpp::NumericVector input_x_probs;
    int input_block_idx;
    // TODO: `miss_block` data?

    // targets
    Rcpp::NumericVector target_w_probs;
    Rcpp::NumericVector target_x_probs;
    Rcpp::IntegerVector target_sample_covs;
    Rcpp::NumericVector target_alt_exp_ubeta_vals;
    Rcpp::NumericVector target_alt_utau_vals;
    Rcpp::IntegerVector target_categ_update;
    Rcpp::NumericVector target_ubeta_update;
    Rcpp::NumericVector target_utau_update;
    Rcpp::NumericMatrix target_u_update;

    int seed_val;
    double epsilon;

};


#endif
