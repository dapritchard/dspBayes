#include <algorithm>
#include "Rcpp.h"
#include "PhiGen.h"
#include "UTestPhiGen.h"

#define DSP_BAYES_UTEST_PHI_GEN_SEED_VAL 99

using Rcpp::as;




void PhiGenTest::setUp() {

    // construct phi
    PhiGen phi = new PhiGen(g_phi_specs, g_n_samp);

    // construct xi.  We pass it an empty list since we don't need the subject
    // information.
    XiGen xi = new XiGen(Rcpp::List, g_n_samp);
    std::copy(g_xi_vals.begin(), g_xi_vals.end(), xi.output_start());

    // set seed for reproducibility
    set.seed(DSP_BAYES_SEED_VAL);
}




void PhiGenTest::tearDown() {
    delete phi;
    delete xi;
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





void PhiGenTest::test_calculations() {



}
