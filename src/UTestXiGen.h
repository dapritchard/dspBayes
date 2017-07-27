#ifndef DSP_BAYES_UTEST_XI_GEN_H
#define DSP_BAYES_UTEST_XI_GEN_H

#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"

// extern int g_n_samp;
// extern double g_eps;

// extern Rcpp::NumericVector g_phi_specs;
// extern Rcpp::NumericVector g_xi_vals;
// extern Rcpp::NumericVector g_test_data_phi;
// extern Rcpp::NumericVector g_test_data_phi_samples;
// extern Rcpp::List g_subj_day_blocks;


class XiGenTest : public CppUnit::TestFixture {

public:

    XiGenTest();

    void setUp();
    void tearDown();

//     // void init_members(int n_samp,
//     // 		      Rcpp::NumericVector phi_specs,
//     // 		      Rcpp::NumericVector xi_vals,
//     // 		      Rcpp::NumericVector test_data_phi,
//     // 		      Rcpp::NumericVector test_data_phi_samples);
    void test_constructor();
    void test_sample_yes_record();
    void test_sample_no_record();

    CPPUNIT_TEST_SUITE(XiGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample_yes_record);
    CPPUNIT_TEST(test_sample_no_record);
    CPPUNIT_TEST_SUITE_END();


private:

    XiGen* xi, * xi_no_rec;
    WGen* W;
    PhiGen* phi;
    UProdBeta* ubeta;

    Rcpp::List subj_day_blocks;
    int n_samp;
//     int m_n_samp;
//     Rcpp::NumericVector m_phi_specs;
//     Rcpp::NumericVector m_xi_vals;
//     Rcpp::NumericVector m_test_data_phi;
//     Rcpp::NumericVector m_test_data_phi_samples;
};


#endif
