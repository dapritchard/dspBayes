#ifndef DSP_BAYES_UTEST_X_GEN_H
#define DSP_BAYES_UTEST_X_GEN_H

#include "Rcpp.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"


class XGenTest : public CppUnit::TestFixture {

 public:

    XGenTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_sample();
    void test_calc_prior_prob();
    void test_calc_posterior_prob();
    void test_sample_x_ijk();
    void test_sample_day_before_fw_sex();

    CPPUNIT_TEST_SUITE(XGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST(test_calc_prior_prob);
    CPPUNIT_TEST(test_calc_posterior_prob);
    CPPUNIT_TEST(test_sample_x_ijk);
    CPPUNIT_TEST(test_sample_day_before_fw_sex);
    CPPUNIT_TEST_SUITE_END();


 private:

    XGen* X;
    WGen* W;
    XiGen* xi;
    UProdBeta* ubeta;
    UProdTau* utau;
    Rcpp::IntegerVector* x_rcpp_copy;

    // testing data
    Rcpp::IntegerVector X_rcpp;
    // Rcpp::List miss_cyc;
    // XGen::XMissDay** miss_day;
    double cohort_sex_prob;
    double sex_coef;
    int miss_day_idx;
    int day_idx;
    double prior_prob_yes;
    double posterior_prob_yes;
    double xi_i;
    int seed_val;
    double epsilon;

    // targets
    Rcpp::IntegerVector target_x_samples;
    Rcpp::IntegerVector target_x_ijk_samples;
    Rcpp::IntegerVector target_day_before_samples;
    double target_prior_prob_no_prev;
    double target_prior_prob_yes_prev;
    double target_posterior_prob;
};


#endif
