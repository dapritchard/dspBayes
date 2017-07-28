#ifndef DSP_BAYES_UTEST_PHI_GEN_H
#define DSP_BAYES_UTEST_PHI_GEN_H

#include "Rcpp.h"
#include "PhiGen.h"
#include "XiGen.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"

extern UTestFactory g_ut_factory;


class PhiGenTest : public CppUnit::TestFixture {

public:

    PhiGenTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_calculations();
    void test_update();
    void test_sample_yes_record();
    void test_sample_no_record();

    CPPUNIT_TEST_SUITE(PhiGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_calculations);
    CPPUNIT_TEST(test_update);
    CPPUNIT_TEST(test_sample_yes_record);
    CPPUNIT_TEST(test_sample_no_record);
    CPPUNIT_TEST_SUITE_END();


private:

    PhiGen* phi;
    PhiGen* phi_no_rec;
    XiGen* xi;

    int seed_val;
    double epsilon;
    int n_samp;

    double c1;
    double c2;
    double mean;
    double delta;
    double phi_init;
    double proposal_val;
    double log_dgamma_norm_const;
    double log_proportion_dgamma_phi;
    double log_proportion_dgamma_xi;
    double m_log_norm_const;
    double calc_log_r;
    int accept_ctr;

    Rcpp::NumericVector target_samples;
};


#endif
