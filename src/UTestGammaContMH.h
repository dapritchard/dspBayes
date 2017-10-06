#ifndef DSP_BAYES_UTEST_GAMMA_CONT_MN_H
#define DSP_BAYES_UTEST_GAMMA_CONT_MN_H

// #include "Rcpp.h"
// #include "GammaGen.h"
// #include "WGen.h"
// #include "XiGen.h"
// #include "UProdBeta.h"
// #include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"


class GammaContMHTest : public CppUnit::TestFixture {

public:

    // void setUp();
    // void tearDown();

    void test_constructor();
    void test_sample();
    void test_sample_proposal_beta_from_0();
    void test_sample_proposal_beta_from_cont();
    void test_get_log_r();
    void test_get_w_log_lik();
    void test_get_gam_log_lik();
    void test_get_proposal_log_lik();
    void test_log_dgamma_trunc_norm_const_bnd_0_inf();
    void test_log_dgamma_trunc_norm_const_bnd_05_1();

    // utility functions
    void set_seed(int seed_val);
    GammaContMH gen_gamma_proposal2_propval04();
    GammaContMH gen_gamma_bndl_bndu(double bnd_l,  double bnd_u);

    CPPUNIT_TEST_SUITE(GammaContMHTest);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST(test_sample_proposal_beta_from_0);
    CPPUNIT_TEST(test_sample_proposal_beta_from_cont);
    CPPUNIT_TEST(test_get_log_r);
    CPPUNIT_TEST(test_get_w_log_lik);
    CPPUNIT_TEST(test_get_gam_log_lik);
    CPPUNIT_TEST(test_get_proposal_log_lik);
    CPPUNIT_TEST(test_log_dgamma_trunc_norm_const);
    CPPUNIT_TEST_SUITE_END();
};


#endif
