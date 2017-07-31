#ifndef DSP_BAYES_UTEST_GAMMA_CATEG_H
#define DSP_BAYES_UTEST_GAMMA_CATEG_H

#include "Rcpp.h"
#include "GammaGen.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"


class GammaCategTest : public CppUnit::TestFixture {

public:

    GammaCategTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_calculations();
    void test_sample_gamma();
    void test_sample();

    CPPUNIT_TEST_SUITE(GammaCategTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_calculations);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST_SUITE_END();


private:

    GammaCateg* gamma_all;
    GammaCateg* gamma_zero_one;
    GammaCateg* gamma_one_inf;
    GammaCateg* gamma_zero_half;
    WGen* W;
    XiGen* xi;
    UProdBeta* ubeta;
    int* X;

    double seed_val;
    double epsilon;
    // double a;
    // double b;
    double hyp_a;
    double hyp_b;
    double m_log_d2_const_terms_all;
    double m_log_d2_const_terms_zero_one;
    double m_log_d2_const_terms_one_inf;
    double m_log_d2_const_terms_zero_half;
    double log_dgamma_norm_const;
    double log_dgamma_trunc_const_all;
    double log_dgamma_trunc_const_zero_one;
    double log_dgamma_trunc_const_one_inf;
    double log_dgamma_trunc_const_zero_half;
    double a_tilde;
    double b_tilde;
    Rcpp::NumericVector target_ubeta_no_h;
    double p_tilde_all;
    double p_tilde_zero_one;
    double p_tilde_one_inf;
    double p_tilde_zero_half;
    Rcpp::NumericVector target_samples_all;
    Rcpp::NumericVector target_samples_zero_one;
    Rcpp::NumericVector target_samples_one_inf;
    Rcpp::NumericVector target_samples_zero_half;
    // double update_all;
    // double update_zero_one;
    // double update_one_inf;
    // double update_zero_half;
    Rcpp::NumericVector target_ubeta_all;
    Rcpp::NumericVector target_ubeta_zero_one;
    Rcpp::NumericVector target_ubeta_one_inf;
    Rcpp::NumericVector target_ubeta_zero_half;
};


#endif
