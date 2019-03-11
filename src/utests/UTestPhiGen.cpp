#include <algorithm>
#include "Rcpp.h"
#include "PhiGen.h"
#include "XiGen.h"
#include "UTestPhiGen.h"
#include "UTestFactory.h"

extern UTestFactory g_ut_factory;

using Rcpp::as;




PhiGenTest::PhiGenTest() :
    seed_val(as<int>(g_ut_factory.seed_vals["phi"])),
    epsilon(UTestFactory::epsilon),
    n_samp(g_ut_factory.n_samp),
    c1(as<double>(g_ut_factory.phi_specs["c1"])),
    c2(as<double>(g_ut_factory.phi_specs["c2"])),
    mean(as<double>(g_ut_factory.phi_specs["mean"])),
    delta(as<double>(g_ut_factory.phi_specs["delta"])),
    phi_init(as<double>(g_ut_factory.target_data_phi["phi_val"])),
    proposal_val(as<double>(g_ut_factory.target_data_phi["proposal_val"])),
    log_dgamma_norm_const(as<double>(g_ut_factory.target_data_phi["log_dgamma_norm_const"])),
    log_proportion_dgamma_phi(as<double>(g_ut_factory.target_data_phi["log_proportion_dgamma_phi"])),
    log_proportion_dgamma_xi(as<double>(g_ut_factory.target_data_phi["log_proportion_dgamma_xi"])),
    m_log_norm_const(as<double>(g_ut_factory.target_data_phi["m_log_norm_const"])),
    calc_log_r(as<double>(g_ut_factory.target_data_phi["calc_log_r"])),
    accept_ctr(as<int>(g_ut_factory.target_data_phi["accept_ctr"])),
    target_samples(g_ut_factory.target_samples_phi) {
}




void PhiGenTest::setUp() {

    // construct phi
    phi = g_ut_factory.phi();
    phi_no_rec = g_ut_factory.phi_no_rec();

    // construct xi
    xi = g_ut_factory.xi();
}




void PhiGenTest::tearDown() {
    delete phi;
    delete phi_no_rec;
    delete xi;
}




void PhiGenTest::test_constructor() {

    // member initialization
    CPPUNIT_ASSERT_EQUAL(c1, phi->m_hyp_c1);
    CPPUNIT_ASSERT_EQUAL(c2, phi->m_hyp_c2);
    CPPUNIT_ASSERT_EQUAL(delta, phi->m_delta);
    CPPUNIT_ASSERT_EQUAL(n_samp, (int) phi->m_vals_rcpp.size());
    CPPUNIT_ASSERT_EQUAL(1, (int) phi_no_rec->m_vals_rcpp.size());
    CPPUNIT_ASSERT(phi->m_vals_rcpp.begin() == phi->m_vals);
    CPPUNIT_ASSERT_EQUAL(0, phi->m_accept_ctr);
    CPPUNIT_ASSERT(phi->m_record_status);
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_EQUAL(0.0, phi->m_log_norm_const);

    // initial value for phi
    CPPUNIT_ASSERT_EQUAL(mean, *(phi->m_vals));
}




void PhiGenTest::test_calculations() {

    // log_dgamma_norm_const()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_norm_const,
                                 phi->log_dgamma_norm_const(proposal_val),
                                 epsilon);

    // calc_log_proportion_dgamma_phi()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_proportion_dgamma_phi,
                                 phi->calc_log_proportion_dgamma_phi(proposal_val),
                                 epsilon);

    // calc_log_proportion_dgamma_xi()
    // run once when `m_is_same_as_prev` is false
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_proportion_dgamma_xi,
                                 phi->calc_log_proportion_dgamma_xi(*xi, proposal_val),
                                 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_norm_const,
                                 phi->m_log_norm_const,
                                 epsilon);
    // now run when `m_is_same_as_prev` is true.  The function should reuse some
    // of the internal calculations from last time
    phi->m_is_same_as_prev = true;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_proportion_dgamma_xi,
                                 phi->calc_log_proportion_dgamma_xi(*xi, proposal_val),
                                 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_norm_const,
                                 phi->m_log_norm_const,
                                 epsilon);
}




void PhiGenTest::test_update() {

    // calc_log_r
    CPPUNIT_ASSERT_DOUBLES_EQUAL(calc_log_r,
                                 phi->calc_log_r(*xi, proposal_val),
                                 epsilon);

    // update_phi
    // reject proposal value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(*(phi->m_vals),
                                 phi->update_phi(R_NegInf, proposal_val),
                                 epsilon);
    CPPUNIT_ASSERT(phi->m_is_same_as_prev);
    // accept proposal value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(proposal_val,
                                 phi->update_phi(0.0, proposal_val),
                                 epsilon);
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_EQUAL(1, phi->m_accept_ctr);
}




// record samples when storing the results
void PhiGenTest::test_sample_yes_record() {

    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // sample
    for (int i = 0; i < target_samples.size(); ++i) {
        phi->sample(*xi);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples[i], phi->val(), epsilon);
    }
    CPPUNIT_ASSERT_EQUAL(target_samples.size(),
                         phi->m_vals - phi->m_vals_rcpp.begin());
    CPPUNIT_ASSERT_EQUAL(accept_ctr, phi->m_accept_ctr);
}




// record samples when not storing the results
void PhiGenTest::test_sample_no_record() {

    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // sample
    for (int i = 0; i < target_samples.size(); ++i) {
        phi_no_rec->sample(*xi);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples[i], phi_no_rec->val(), epsilon);
    }
    CPPUNIT_ASSERT_EQUAL((long int) 0, phi_no_rec->m_vals - phi_no_rec->m_vals_rcpp.begin());
    CPPUNIT_ASSERT_EQUAL(accept_ctr, phi_no_rec->m_accept_ctr);
}
