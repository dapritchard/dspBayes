#ifndef DSP_BAYES_UTEST_GAMMA_CONT_MH_H
#define DSP_BAYES_UTEST_GAMMA_CONT_MH_H


#include "cppunit/extensions/HelperMacros.h"
#include "Rcpp.h"

// #include "GammaGen.h"
#include "UProdBeta.h"
#include "WGen.h"
#include "XiGen.h"

// #include "UTestFactory.h"




class GammaContMHTest : public CppUnit::TestFixture {

public:

    // inner classes
    class FromScratchW;
    class FromScratchXi;
    class FromScratchUbeta;
    class FromScratchGamma;

    // constructor
    GammaContMHTest();

    // setUp, tearDown
    void setUp();
    void tearDown();

    // test methods
    void test_constructor();
    void test_sample();
    void test_sample_proposal_beta_from_0();
    void test_sample_proposal_beta_from_cont();
    void test_get_log_r();
    void test_get_w_log_lik();
    void test_get_gam_log_lik();
    void test_get_proposal_log_lik_prop0_curr0();
    void test_get_proposal_log_lik_prop0_curr2();
    void test_get_proposal_log_lik_prop2_curr0();
    void test_get_proposal_log_lik_prop2_curr2();
    void test_log_dgamma_trunc_norm_const_bnd_0_inf();
    void test_log_dgamma_trunc_norm_const_bnd_05_1();

    // utility functions
    void set_seed(int seed_val);
    GammaContMH gen_gamma_proposal2_propval04();
    GammaContMH gen_gamma_bndl_bndu(double bnd_l,  double bnd_u);
    GammaContMH gen_gamma_curr(double curr);

    CPPUNIT_TEST_SUITE(GammaContMHTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST(test_sample_proposal_beta_from_0);
    CPPUNIT_TEST(test_sample_proposal_beta_from_cont);
    CPPUNIT_TEST(test_get_log_r);
    CPPUNIT_TEST(test_get_w_log_lik);
    CPPUNIT_TEST(test_get_gam_log_lik);
    CPPUNIT_TEST(test_get_proposal_log_lik_prop0_curr0);
    CPPUNIT_TEST(test_get_proposal_log_lik_prop0_curr2);
    CPPUNIT_TEST(test_get_proposal_log_lik_prop2_curr0);
    CPPUNIT_TEST(test_get_proposal_log_lik_prop2_curr2);
    CPPUNIT_TEST(test_log_dgamma_trunc_norm_const_bnd_0_inf);
    CPPUNIT_TEST(test_log_dgamma_trunc_norm_const_bnd_05_1);
    CPPUNIT_TEST_SUITE_END();

private:

    FromScratchGamma* gamma_obj;
    FromScratchW* w_obj;
    FromScratchXi* xi_obj;
    FromScratchUbeta* ubeta_obj;
    Rcpp::NumericMatrix m_Uh;
    int old_d2s[9];
    int new_d2s[9];
};




class GammaContMHTest::FromScratchGamma {

public:

    Rcpp::NumericMatrix U;
    Rcpp::NumericVector gamma_specs;
    GammaContMH gamma_obj;

    FromScratchGamma();
    GammaContMH& gamma() { return gamma_obj; }
};




class GammaContMHTest::FromScratchW {

public:

    Rcpp::IntegerVector subj0_preg_cyc;
    Rcpp::IntegerVector subj1_preg_cyc;
    Rcpp::List preg_cyc;
    Rcpp::IntegerVector w_to_days_idx;
    Rcpp::IntegerVector w_cyc_to_subj_idx;
    WGen w_obj;

    FromScratchW();
    WGen& W() { return w_obj; }
};




class GammaContMHTest::FromScratchXi {

public:

    Rcpp::IntegerVector subj0_day_block;
    Rcpp::IntegerVector subj1_day_block;
    Rcpp::List subj_day_blocks;
    XiGen xi_obj;

    FromScratchXi();
    XiGen& xi() { return xi_obj; }
};




class GammaContMHTest::FromScratchUbeta {

public:

    UProdBeta ubeta_obj;

    FromScratchUbeta();
    UProdBeta& ubeta() { return ubeta_obj; }
};


#endif
