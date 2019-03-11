#include <algorithm>
#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestXiGen.h"

using Rcpp::as;

extern UTestFactory g_ut_factory;




XiGenTest::XiGenTest() :
    target_samples(g_ut_factory.target_samples_xi),
    n_subj(g_ut_factory.n_subj),
    n_samp(g_ut_factory.n_samp),
    seed_val(as<double>(g_ut_factory.seed_vals["xi"])) {
}




void XiGenTest::setUp() {

    // construct xi
    xi = g_ut_factory.xi();
    xi_no_rec = g_ut_factory.xi_no_rec();

    // construct W
    W = g_ut_factory.W();

    // construct phi
    phi = g_ut_factory.phi();

    // construct X
    X = g_ut_factory.X();

    // constuct ubeta
    ubeta = g_ut_factory.ubeta();
}




void XiGenTest::tearDown() {
    delete xi;
    delete xi_no_rec;
    delete W;
    delete phi;
    delete X;
    delete ubeta;
}




void XiGenTest::test_constructor() {

    // test record samples variant
    CPPUNIT_ASSERT_EQUAL(n_subj * n_samp, (int) xi->m_vals_rcpp.size());
    CPPUNIT_ASSERT_EQUAL(xi->m_vals_rcpp.begin(), xi->m_vals);
    CPPUNIT_ASSERT_EQUAL(n_subj, xi->m_n_subj);
    CPPUNIT_ASSERT(xi->m_record_status);

    // test non-record samples variant
    CPPUNIT_ASSERT_EQUAL(n_subj, (int) xi_no_rec->m_vals_rcpp.size());
    CPPUNIT_ASSERT_EQUAL(xi_no_rec->m_vals_rcpp.begin(), xi_no_rec->m_vals);
    CPPUNIT_ASSERT_EQUAL(n_subj, xi_no_rec->m_n_subj);
    CPPUNIT_ASSERT(! xi_no_rec->m_record_status);
}




void XiGenTest::test_sample_yes_record() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // two samples using the same seed
    set_seed(seed_val);
    xi->sample(*W, *phi, *ubeta, *X);
    set_seed(seed_val);
    xi->sample(*W, *phi, *ubeta, *X);

    // check that placement of iterator points to beginning of second sample
    CPPUNIT_ASSERT_EQUAL(xi->m_vals_rcpp.begin() + 2 * n_subj, xi->m_vals);
    // check values of samples
    CPPUNIT_ASSERT(std::equal(target_samples.begin(),
                              target_samples.end(),
                              xi->m_vals - n_subj,
                              UTestFactory::eq_dbl));
    CPPUNIT_ASSERT(std::equal(target_samples.begin(),
                              target_samples.end(),
                              xi->m_vals,
                              UTestFactory::eq_dbl));
}




void XiGenTest::test_sample_no_record() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // two samples using the same seed
    set_seed(seed_val);
    xi_no_rec->sample(*W, *phi, *ubeta, *X);
    set_seed(seed_val);
    xi_no_rec->sample(*W, *phi, *ubeta, *X);

    // check that placement of iterator points to beginning of data
    CPPUNIT_ASSERT_EQUAL(xi_no_rec->m_vals_rcpp.begin(), xi_no_rec->m_vals);
    // check values of samples
    CPPUNIT_ASSERT(std::equal(target_samples.begin(),
                              target_samples.end(),
                              xi_no_rec->m_vals,
                              UTestFactory::eq_dbl));
}
