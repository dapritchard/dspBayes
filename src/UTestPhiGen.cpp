#include <algorithm>
#include "Rcpp.h"
#include "PhiGen.h"
#include "UTestPhiGen.h"

using Rcpp::as;




// void init_members(int n_samp,
// 		  Rcpp::NumericVector phi_specs,
// 		  Rcpp::NumericVector xi_vals,
// 		  Rcpp::NumericVector test_data_phi,
// 		  Rcpp::NumericVector test_data_phi_samples) {

//     m_n_samp = n_samp;
//     m_phi_specs = phi_specs;
//     m_xi_vals = xi_vals;
//     m_test_data_phi = test_data_phi;
//     m_test_data_phi_samples = test_data_phi_samples;
// }




void PhiGenTest::setUp() {

    // construct phi
    phi = new PhiGen(g_phi_specs, g_n_samp, true);

    // construct xi.  We pass it an empty list since we don't need the subject
    // information.
    Rcpp::List empty;
    xi = new XiGen(g_subj_day_blocks, g_n_samp);
    std::copy(g_xi_vals.begin(), g_xi_vals.end(), xi->output_start());
}




void PhiGenTest::tearDown() {
    delete phi;
    delete xi;
}




void PhiGenTest::test_constructor() {

    // initialize members
    CPPUNIT_ASSERT_EQUAL(as<double>(g_phi_specs["c1"]), phi->m_hyp_c1);
    CPPUNIT_ASSERT_EQUAL(as<double>(g_phi_specs["c2"]), phi->m_hyp_c2);
    CPPUNIT_ASSERT_EQUAL(as<double>(g_phi_specs["delta"]), phi->m_delta);
    CPPUNIT_ASSERT_EQUAL(g_n_samp, (int) phi->m_vals_rcpp.size());
    CPPUNIT_ASSERT(phi->m_vals_rcpp.begin() == phi->m_vals);
    CPPUNIT_ASSERT_EQUAL(0, phi->m_accept_ctr);
    CPPUNIT_ASSERT(phi->m_record_status);
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_EQUAL(0.0, phi->m_log_norm_const);

    // set current entry of phi samples as the mean of the prior distribution
    CPPUNIT_ASSERT_EQUAL(as<double>(g_phi_specs["mean"]), *(phi->m_vals));
}




void PhiGenTest::test_calculations() {

    *(phi->m_vals) = as<double>(g_test_data_phi["phi_val"]);
    double proposal_val = as<double>(g_test_data_phi["proposal_val"]);

    // log_dgamma_norm_const
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["log_dgamma_norm_const"]),
				 phi->log_dgamma_norm_const(proposal_val),
				 g_eps);

    // calc_log_proportion_dgamma_phi
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["log_proportion_dgamma_phi"]),
				 phi->calc_log_proportion_dgamma_phi(proposal_val),
				 g_eps);

    // calc_log_proportion_dgamma_xi
    // run once when `m_is_same_as_prev` is false
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["log_proportion_dgamma_xi"]),
				 phi->calc_log_proportion_dgamma_xi(*xi, proposal_val),
				 g_eps);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["m_log_norm_const"]),
				 phi->m_log_norm_const,
				 g_eps);
    // now run when `m_is_same_as_prev` is true.  The function should reuse some
    // of the internal calculations from last time
    phi->m_is_same_as_prev = true;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["calc_log_proportion_dgamma_xi"]),
				 phi->calc_log_proportion_dgamma_xi(*xi, proposal_val),
				 g_eps);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["m_log_norm_const"]),
				 phi->m_log_norm_const,
				 g_eps);
}




void PhiGenTest::test_update() {

    *(phi->m_vals) = as<double>(g_test_data_phi["phi_val"]);
    double proposal_val = as<double>(g_test_data_phi["proposal_val"]);

    // calc_log_r
    CPPUNIT_ASSERT_DOUBLES_EQUAL(as<double>(g_test_data_phi["calc_log_r"]),
				 phi->calc_log_r(*xi, proposal_val),
				 g_eps);

    // update_phi
    // reject proposal value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(*(phi->m_vals),
				 phi->update_phi(R_NegInf, proposal_val),
				 g_eps);
    CPPUNIT_ASSERT(phi->m_is_same_as_prev);
    // accept proposal value
    CPPUNIT_ASSERT_DOUBLES_EQUAL(proposal_val,
				 phi->update_phi(0.0, proposal_val),
				 g_eps);
    CPPUNIT_ASSERT(! phi->m_is_same_as_prev);
    CPPUNIT_ASSERT_EQUAL(1, phi->m_accept_ctr);

}




void PhiGenTest::test_sample() {

    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(as<unsigned int>(g_test_data_phi["seed_val"]));

    // sample
    for (int i = 0; i < g_test_data_phi_samples.size(); ++i) {
	phi->sample(*xi);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(g_test_data_phi_samples[i], phi->val(), g_eps);
    }
    CPPUNIT_ASSERT_EQUAL(g_test_data_phi_samples.size(),
    			 phi->m_vals - phi->m_vals_rcpp.begin());
    CPPUNIT_ASSERT_EQUAL(as<int>(g_test_data_phi["accept_ctr"]),
			 phi->m_accept_ctr);
}
