#include "Rcpp.h"
#include "cppunit/extensions/HelperMacros.h"

#include "GammaGen.h"
#include "UTestGammaContMH.h"

#define EPSILON 0.000000000001
#define SEED_YIELDS_0_29  24
#define SEED_YIELDS_0_91  72




void GammaContMHTest::test_constructor() {

    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_norm_const,       gamma->m_log_norm_const,       epsilon);
    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_p_over_1_minus_p, gamma->m_log_p_over_1_minus_p, epsilon);
    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_1_minus_p_over_p, gamma->m_log_1_minus_p_over_p, epsilon);

    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_mh_p,             gamma->m_mh_p, epsilon);
    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_mh_log_p,         gamma->m_mh_log_p, epsilon);
    // CPPUNIT_ASSERT_DOUBLES_EQUAL(m_mh_log_1_minus_p, gamma->m_mh_log_1_minus_p, epsilon);

    // CPPUNIT_ASSERT_EQUAL(m_mh_accept_ctr, gamma->m_mh_accept_ctr);
    // // TODO: test function pointers?

    CPPUNIT_FAIL("implement test_constructor");
}




void GammaContMHTest::test_sample() {

    CPPUNIT_FAIL("implement test_sample");
}




void GammaContMHTest::test_sample_proposal_beta_from_0() {

    // fixture setup.  Creates an object with a proposal function that returns a
    // value of 2, and samples from the discrete part of the distribution when
    // the random number generator provides a value less than 0.4.
    GammaContMH gamma_proposal2_propval04 = gen_gamma_proposal2_propval04();

    // exercise SUT
    set_seed(SEED_YIELDS_0_29);
    double out = gamma_proposal2_propval04.sample_proposal_beta();

    // verify outcome
    CPPUNIT_ASSERT_EQUAL(0.0, out);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_sample_proposal_beta_from_cont() {

    // fixture setup.  Creates an object with a proposal function that returns a
    // value of 2, and samples from the discrete part of the distribution when
    // the random number generator provides a value less than 0.4.
    GammaContMH gamma_proposal2_propval04 = gen_gamma_proposal2_propval04();

    // exercise SUT
    set_seed(SEED_YIELDS_0_91);
    double out = gamma_proposal2_propval04.sample_proposal_beta();

    // verify outcome
    CPPUNIT_ASSERT_EQUAL(2.0, out);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_get_log_r() {

    CPPUNIT_FAIL("implement test_get_log_r");
}




void GammaContMHTest::test_get_w_log_lik() {

    CPPUNIT_FAIL("implement test_get_w_log_lik");
}




void GammaContMHTest::test_get_gam_log_lik() {

    CPPUNIT_FAIL("implement test_get_gam_log_lik");
}




void GammaContMHTest::test_get_proposal_log_lik() {

    CPPUNIT_FAIL("implement test_get_proposal_log_lik");
}




void GammaContMHTest::test_log_dgamma_trunc_norm_const_bnd_0_inf() {

    // fixture setup.  Instantiates data for various choices of lower and upper
    // bound for the gamma-trucated distribution.  Gamma distribution is set to
    // 1.2 and 0.9 for the shape and rate parameters, respectively.
    GammaContMH gamma_00_inf = gen_gamma_bndl_bndu(0.0, R_PosInf);

    // exercise SUT
    double out_00_inf = gamma_00_inf.log_dgamma_trunc_norm_const();

    // verify outcome
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.0410585287860757, out_00_inf, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_log_dgamma_trunc_norm_const_bnd_05_1() {

    // fixture setup.  Instantiates data for various choices of lower and upper
    // bound for the gamma-trucated distribution.  Gamma distribution is set to
    // 1.2 and 0.9 for the shape and rate parameters, respectively.
    GammaContMH gamma_05_1   = gen_gamma_bndl_bndu(0.5, 1.0);

    // exercise SUT
    double out_05_1   = gamma_05_1.log_dgamma_trunc_norm_const();

    // verify outcome
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.4254311104708031, out_05_1,   EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




// utility functions -----------------------------------------------------------

void GammaContMHTest::set_seed(int seed_val) {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // `gamma_all`
    set_seed(seed_val);
}




// creates an object with a proposal function that returns a value of 2, and
// samples from the discrete part of the distribution when the random number
// generator provides a value less than 0.4.

GammaContMH GammaContMHTest::gen_gamma_proposal2_propval04() {

    // the only values of importance are those for `mh_p` and `mh_delta`
    Rcpp::NumericMatrix U(10, 1);
    Rcpp::NumericVector gamma_specs = Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
								  Rcpp::_["hyp_a"]    = 1.0,
								  Rcpp::_["hyp_b"]    = 1.0,
								  Rcpp::_["hyp_p"]    = 0.5,
								  Rcpp::_["bnd_l"]    = 0.0,
								  Rcpp::_["bnd_u"]    = R_PosInf,
								  Rcpp::_["mh_p"]     = 0.4,
								  Rcpp::_["mh_delta"] = 0.0);

    // specify the current value of `beta_h` and corresponding `gamma_h`
    GammaContMH gamma(U, gamma_specs);
    gamma.m_beta_val = 2.0;
    gamma.m_gam_val  = exp(2.0);

    return gamma;
}




GammaContMH GammaContMHTest::gen_gamma_bndl_bndu(double bnd_l,  double bnd_u) {

    // the only values of importance are those for `bnd_l` and `bnd_u`
    Rcpp::NumericMatrix U(10, 1);
    Rcpp::NumericVector gamma_specs = Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
								  Rcpp::_["hyp_a"]    = 1.2,
								  Rcpp::_["hyp_b"]    = 0.9,
								  Rcpp::_["hyp_p"]    = 0.5,
								  Rcpp::_["bnd_l"]    = bnd_l,
								  Rcpp::_["bnd_u"]    = bnd_u,
								  Rcpp::_["mh_p"]     = 0.4,
								  Rcpp::_["mh_delta"] = 0.1);

    return GammaContMH(U, gamma_specs);
}




// # calculates the log norming constant for a possibly truncated Gamma(a, b)
// # distribution, given by
//
// log_dgamma_trunc_norm_const <- function(bnd_l, bnd_u, a, b) {
//
//     F_upp = pgamma(bnd_u, a, b, log.p = FALSE)
//     F_low = pgamma(bnd_l, a, b, log.p = FALSE)
//
//     return(-log(F_upp - F_low) + (a * log(b)) - lgamma(a))
// }
//
//
// out_00_inf <- log_dgamma_trunc_norm_const(0.0, Inf, 1.2, 0.9)
// out_05_1   <- log_dgamma_trunc_norm_const(0.5, 1.0, 1.2, 0.9)
//
// cat(sprintf("%.16f", out_00_inf), "\n",
//     sprintf("%.16f", out_05_1), "\n")
