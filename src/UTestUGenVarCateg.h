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

    CPPUNIT_TEST_SUITE(UGenVarCategTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST_SUITE_END();


private:

    UGenVarCateg* u_var;
    // XGen* X;
    // WGen* W;
    // XiGen* xi;
    // UProdBeta* ubeta;
    // UProdTau* utau;
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
    // `miss_block` data?

    // targets
    Rcpp::NumericVector target_w_probs;
    Rcpp::NumericVector target_x_probs;
    Rcpp::IntegerVector target_sample_covs;
    Rcpp::IntegerVector target_u_update;
    Rcpp::IntegerVector target_ubeta_update;
    Rcpp::IntegerVector target_utau_update;

    /* // testing data */
    /* Rcpp::IntegerVector X_rcpp; */
    /* Rcpp::List miss_cyc_rcpp; */
    /* Rcpp::List miss_day_rcpp; */
    /* double cohort_sex_prob; */
    /* double sex_coef; */
    /* int miss_day_idx; */
    /* int day_idx; */
    /* double prior_prob_yes; */
    /* double posterior_prob_yes; */
    /* double xi_i; */
    /* int seed_val; */
    /* double epsilon; */

    /* // targets */
    /* Rcpp::IntegerVector target_x_samples; */
    /* Rcpp::IntegerVector target_x_ijk_samples; */
    /* Rcpp::IntegerVector target_day_before_samples; */
    /* double target_prior_prob_no_prev; */
    /* double target_prior_prob_yes_prev; */
    /* double target_posterior_prob; */
};


#endif
