#ifndef DSP_BAYES_UTEST_U_GEN_VAR_CATEG_H
#define DSP_BAYES_UTEST_U_GEN_VAR_CATEG_H

/* #include "Rcpp.h" */
/* #include "WGen.h" */
/* #include "XGen.h" */
/* #include "XiGen.h" */
/* #include "UProdBeta.h" */
/* #include "UProdTau.h" */
/* #include "UTestFactory.h" */
#include "cppunit/extensions/HelperMacros.h"


class UGenVarCategTest : public CppUnit::TestFixture {

public:

    UGenVarCategTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_sample_covariate();
    void test_calc_posterior_w();
    void test_calc_posterior_x();

    CPPUNIT_TEST_SUITE(UGenVarCategTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample_covariate);
    // CPPUNIT_TEST(test_calc_posterior_w);
    CPPUNIT_TEST(test_calc_posterior_x);
    CPPUNIT_TEST_SUITE_END();


private:

    UGenVarCateg* u_var;
    XGen* X;
    // WGen* W;
    // XiGen* xi;
    // UProdBeta* ubeta;
    UProdTau* utau;
    // Rcpp::IntegerVector* x_rcpp_copy;

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
    // Rcpp::IntegerVector target_u_update;
    // Rcpp::IntegerVector target_ubeta_update;
    Rcpp::NumericVector target_alt_utau_vals;

    int seed_val;
    double epsilon;

};


#endif
