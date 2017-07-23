#ifndef DSP_BAYES_UTEST_PHI_GEN_H
#define DSP_BAYES_UTEST_PHI_GEN_H

#include "Rcpp.h"
#include "PhiGen.h"

// #include "cppunit/ui/text/TestRunner.h"
#include "cppunit/extensions/HelperMacros.h"

// #include <vector>
// #include "cppunit/extensions/HelperMacros.h"


// class CopyTest : public CppUnit::TestFixture {

//  private:

//     std::vector<char> *v1;
//     std::vector<char> *v2;

//  public:

//     void setUp();
//     void tearDown();

//     void test_equal_iters();
//     void test_copy();

//     CPPUNIT_TEST_SUITE(CopyTest);
//     CPPUNIT_TEST(test_equal_iters);
//     CPPUNIT_TEST(test_copy);
//     CPPUNIT_TEST_SUITE_END();
// };

extern Rcpp::NumericVector g_phi_specs;
extern int g_n_samp;


class PhiGenTest : public CppUnit::TestFixture {

public:

    PhiGen* phi;

    void setUp();
    void tearDown();
    void test_constructor();

    CPPUNIT_TEST_SUITE(PhiGenTest);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST_SUITE_END();

// private:

//     Rcpp::NumericVector phi_specs;
//     int n_samp;
};


#endif
