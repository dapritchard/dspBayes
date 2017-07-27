#ifndef DSP_BAYES_UTEST_PHI_GEN_H
#define DSP_BAYES_UTEST_PHI_GEN_H

#include "Rcpp.h"
#include "PhiGen.h"
#include "XiGen.h"
#include "cppunit/extensions/HelperMacros.h"

extern int g_n_samp;
extern double g_eps;

extern Rcpp::NumericVector g_phi_specs;
extern Rcpp::NumericVector g_xi_vals;
extern Rcpp::NumericVector g_test_data_phi;
extern Rcpp::NumericVector g_test_data_phi_samples;
extern Rcpp::List g_subj_day_blocks;


class PhiGenTest : public CppUnit::TestFixture {

public:

    void setUp();
    void tearDown();

    // void init_members(int n_samp,
    // 		      Rcpp::NumericVector phi_specs,
    // 		      Rcpp::NumericVector xi_vals,
    // 		      Rcpp::NumericVector test_data_phi,
    // 		      Rcpp::NumericVector test_data_phi_samples);
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
    XiGen* xi;
};


#endif
