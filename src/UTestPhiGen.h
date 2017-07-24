#ifndef DSP_BAYES_UTEST_PHI_GEN_H
#define DSP_BAYES_UTEST_PHI_GEN_H

#include "Rcpp.h"
#include "PhiGen.h"
#include "cppunit/extensions/HelperMacros.h"

extern Rcpp::NumericVector g_phi_specs;
extern int g_n_samp;
extern Rcpp::NumericVector g_xi_vals;


class PhiGenTest : public CppUnit::TestFixture {

public:

    PhiGen* phi;

    void setUp();
    void tearDown();

    void test_constructor();

    CPPUNIT_TEST_SUITE(PhiGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST_SUITE_END();
};


#endif
