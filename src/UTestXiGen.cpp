#include "Rcpp.h"
#include "XiGen.h"
#include "WGen.h"
#include "PhiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestXiGen.h"

extern UTestFactory g_ut_factory;

// extern Rcpp::List g_subj_day_blocks;
// extern Rcpp::List g_preg_cyc;
// extern Rcpp::NumericVector g_phi_hyper;
// extern Rcpp::IntegerVector g_w_days_idx;
// extern Rcpp::IntegerVector g_w_cyc_idx;
// extern Rcpp::IntegerVector g_w_vals;
// extern Rcpp::IntegerVector g_ubeta;
// extern Rcpp::NumericVector g_phi_specs;
// extern int g_n_days;
// extern int g_n_samp;
// extern int g_fw_len;

XiGenTest::XiGenTest() :
    subj_day_blocks(g_ut_factory.subj_day_blocks),
    n_samp(g_ut_factory.n_samp) {
}


void XiGenTest::setUp() {

    // construct xi
    xi = g_ut_factory.xi();
    xi_no_rec = g_ut_factory.xi_no_rec();

    // construct W
    W = g_ut_factory.W();

    // construct phi
    phi = g_ut_factory.phi();

    // constuct ubeta
    ubeta = g_ut_factory.ubeta();
}




void XiGenTest::tearDown() {
    delete xi;
    delete xi_no_rec;
    delete W;
    delete phi;
    delete ubeta;
}




void XiGenTest::test_constructor() {

    // test record samples variant
    CPPUNIT_ASSERT_EQUAL(subj_day_blocks.size() * n_samp, xi->m_vals_rcpp.size());
    CPPUNIT_ASSERT_EQUAL(xi->m_vals_rcpp.begin(), xi->m_vals);
    CPPUNIT_ASSERT(xi->m_subj != NULL);
    CPPUNIT_ASSERT_EQUAL((int) subj_day_blocks.size(), xi->m_n_subj);
    CPPUNIT_ASSERT(xi->m_record_status);

    // test non-record samples variant
    CPPUNIT_ASSERT_EQUAL(subj_day_blocks.size(), xi->m_vals_rcpp.size());
    CPPUNIT_ASSERT_EQUAL(xi->m_vals_rcpp.begin(), xi->m_vals);
    CPPUNIT_ASSERT(xi->m_subj != NULL);
    CPPUNIT_ASSERT_EQUAL((int) subj_day_blocks.size(), xi->m_n_subj);
    CPPUNIT_ASSERT(! xi_no_rec->m_record_status);
}




void XiGenTest::test_sample_yes_record() {

}




void XiGenTest::test_sample_no_record() {

}
