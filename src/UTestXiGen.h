#ifndef DSP_BAYES_UTEST_XI_GEN_H
#define DSP_BAYES_UTEST_XI_GEN_H

#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "cppunit/extensions/HelperMacros.h"


class XiGenTest : public CppUnit::TestFixture {

public:

    XiGenTest();

    void setUp();
    void tearDown();

    void test_constructor();
    void test_sample_yes_record();
    void test_sample_no_record();

    CPPUNIT_TEST_SUITE(XiGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_sample_yes_record);
    CPPUNIT_TEST(test_sample_no_record);
    CPPUNIT_TEST_SUITE_END();


private:

    XiGen* xi;
    XiGen* xi_no_rec;
    WGen* W;
    PhiGen* phi;
    UProdBeta* ubeta;

    // testing data
    Rcpp::NumericVector target_samples;
    int n_subj;
    int n_samp;
    int seed_val;
};


#endif
