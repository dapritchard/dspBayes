#ifndef DSP_BAYES_UTEST_W_GEN_H
#define DSP_BAYES_UTEST_W_GEN_H

#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"


class WGenTest : public CppUnit::TestFixture {

 public:

    WGenTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_sample();

    CPPUNIT_TEST_SUITE(WGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample);
    CPPUNIT_TEST_SUITE_END();


 private:

    WGen* W;
    XiGen* xi;
    UProdBeta* ubeta;

    // testing data
    Rcpp::NumericVector target_samples;
    int n_preg_days;
    int n_preg_cyc;
    int seed_val;
    double epsilon;
};


#endif
