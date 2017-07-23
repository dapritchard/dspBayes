#include "Rcpp.h"
#include "PhiGen.h"
#include "UTestPhiGen.h"

using Rcpp::as;




void PhiGenTest::setUp() {
    phi = new PhiGen(g_phi_specs, g_n_samp);
}




void PhiGenTest::tearDown() {
    delete phi;
}




void PhiGenTest::test_constructor() {

    // initialize members
    CPPUNIT_ASSERT_EQUAL(phi->m_hyp_c1, as<double>(g_phi_specs["c1"]));
    CPPUNIT_ASSERT_EQUAL(phi->m_hyp_c2, as<double>(g_phi_specs["c2"]));
    CPPUNIT_ASSERT_EQUAL(phi->m_delta, as<double>(g_phi_specs["delta"]));
    CPPUNIT_ASSERT(phi->m_vals != NULL);
    CPPUNIT_ASSERT_EQUAL(phi->m_vals, phi->m_output_start);
    CPPUNIT_ASSERT_EQUAL(phi->m_output_start + g_n_samp, phi->m_output_end);
    CPPUNIT_ASSERT_EQUAL(phi->m_accept_ctr, 0);
    CPPUNIT_ASSERT(! phi->m_record_status);
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_EQUAL(phi->m_log_norm_const, 0.0);

    // set current entry of phi samples as the mean of the prior distribution
    CPPUNIT_ASSERT_EQUAL(*(phi->m_vals), as<double>(g_phi_specs["mean"]));
}
