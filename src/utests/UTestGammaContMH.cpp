#include <algorithm>
#include "Rcpp.h"
#include "cppunit/extensions/HelperMacros.h"

#include "GammaGen.h"
#include "UTestGammaContMH.h"
#include "WGen.h"
#include "XiGen.h"

#define EPSILON 0.000000000001
#define SEED_YIELDS_0_29  24
#define SEED_YIELDS_0_91  72

extern int* d2s;




GammaContMHTest::GammaContMHTest() :
    m_Uh(Rcpp::NumericMatrix(10, 1))
{
    std::copy(d2s, d2s + 9, old_d2s);
    new_d2s[0] = 0;    new_d2s[1] = 0;    new_d2s[2] = 0;    new_d2s[3] = 0;    new_d2s[4] = 0;
    new_d2s[5] = 1;    new_d2s[6] = 1;    new_d2s[7] = 1;    new_d2s[8] = 1;
}




void GammaContMHTest::setUp() {
    gamma_obj = new GammaContMHTest::FromScratchGamma();
    w_obj     = new GammaContMHTest::FromScratchW();
    xi_obj    = new GammaContMHTest::FromScratchXi();
    ubeta_obj = new GammaContMHTest::FromScratchUbeta();
    X         = new int[9];
    std::fill(X, X + 9, 1);
    std::copy(new_d2s, new_d2s + 9, d2s);
}




void GammaContMHTest::tearDown() {
    delete gamma_obj;
    delete w_obj;
    delete xi_obj;
    delete ubeta_obj;
    delete[] X;
    std::copy(old_d2s, old_d2s + 9, d2s);
}




// TODO: can do a better job here.  Should test for differnt lower bound / upper bound
void GammaContMHTest::test_constructor() {

    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.0410585287860757, gamma_obj->gamma().m_log_norm_const,       EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0,                gamma_obj->gamma().m_log_p_over_1_minus_p, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0,                gamma_obj->gamma().m_log_1_minus_p_over_p, EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1,      gamma_obj->gamma().m_mh_p,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log(0.1), gamma_obj->gamma().m_mh_log_p,         EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log(0.9), gamma_obj->gamma().m_mh_log_1_minus_p, EPSILON);

    CPPUNIT_ASSERT_EQUAL(0, gamma_obj->gamma().m_mh_accept_ctr);
    // TODO: test function pointers?
}




void GammaContMHTest::test_sample() {

    // TODO: test this
    // CPPUNIT_FAIL("implement test_sample");
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

    // implicit setup

    // exercise SUT
    GammaContMH& gamma = gamma_obj->gamma();
    double out = gamma.get_log_r(w_obj->W(), xi_obj->xi(), ubeta_obj->ubeta(), X, 0.6, exp(0.6));

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0891383454025099, out, EPSILON);

    // implicit teardown
}




void GammaContMHTest::test_get_w_log_lik() {

    // implicit setup

    // exercise SUT
    GammaContMH& gamma = gamma_obj->gamma();
    double out = gamma.get_w_log_lik(w_obj->W(), xi_obj->xi(), ubeta_obj->ubeta(), X, 0.6);

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.9530805686811668, out, EPSILON);

    // implicit teardown
}



// 1.0 proposal, 1.0 current
void GammaContMHTest::test_get_gam_log_lik_10_10() {

    // fixture setup.  Instantiates data with current beta value 1.
    GammaContMH gamma_10 = gen_gamma_curr(1.0);

    // exercise SUT
    double out_10_10 = gamma_10.get_gam_log_lik(1.0, exp(1.0));

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,  out_10_10, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}



// 1.0 proposal, 1.2 current
void GammaContMHTest::test_get_gam_log_lik_10_12() {

    // fixture setup.  Instantiates data with current beta value 1.
    GammaContMH gamma_12 = gen_gamma_curr(1.2);

    // exercise SUT
    double out_10_12 = gamma_12.get_gam_log_lik(1.0, exp(1.0));

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5016515848497520, out_10_12, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}



// 1.2 proposal, 1.0 current
void GammaContMHTest::test_get_gam_log_lik_12_10() {

    // fixture setup.  Instantiates data with current beta value 1.
    GammaContMH gamma_10 = gen_gamma_curr(1.0);

    // exercise SUT
    double out_12_10 = gamma_10.get_gam_log_lik(1.2, exp(1.2));

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5016515848497520,  out_12_10, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}



// 1.2 proposal, 1.2 current
void GammaContMHTest::test_get_gam_log_lik_12_12() {

    // fixture setup.  Instantiates data with current beta value 1.
    GammaContMH gamma_12 = gen_gamma_curr(1.2);

    // exercise SUT
    double out_12_12 = gamma_12.get_gam_log_lik(1.2, exp(1.2));

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0,  out_12_12, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_get_proposal_log_lik_prop0_curr0() {

    // fixture setup.  Instantiates data with current beta value 0.
    GammaContMH gamma_0 = gen_gamma_curr(0.0);

    // exercise SUT
    double out_0_0 = gamma_0.get_proposal_log_lik(0.0);

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, out_0_0, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_get_proposal_log_lik_prop0_curr2() {

    // fixture setup.  Instantiates data with current beta value 0.
    GammaContMH gamma_2 = gen_gamma_curr(2.0);

    // exercise SUT
    double out_0_2 = gamma_2.get_proposal_log_lik(0.0);

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4054651081081642, out_0_2, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_get_proposal_log_lik_prop2_curr0() {

    // fixture setup.  Instantiates data with current beta value 0.
    GammaContMH gamma_0 = gen_gamma_curr(0.0);

    // exercise SUT
    double out_2_0 = gamma_0.get_proposal_log_lik(2.0);

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.4054651081081642, out_2_0, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_get_proposal_log_lik_prop2_curr2() {

    // fixture setup.  Instantiates data with current beta value 0.
    GammaContMH gamma_2 = gen_gamma_curr(2.0);

    // exercise SUT
    double out_2_2 = gamma_2.get_proposal_log_lik(2.0);

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, out_2_2, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_log_dgamma_trunc_norm_const_bnd_0_inf() {

    // fixture setup.  Instantiates data for lower bound 0 and upper bound
    // infinity for the gamma-trucated distribution (with hard-coded parameters
    // 1.2 shape and 0.9 rate).
    GammaContMH gamma_00_inf = gen_gamma_bndl_bndu(0.0, R_PosInf);

    // exercise SUT
    double out_00_inf = gamma_00_inf.log_dgamma_trunc_norm_const();

    // verify outcome.  Result is hand-calculated.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.0410585287860757, out_00_inf, EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




void GammaContMHTest::test_log_dgamma_trunc_norm_const_bnd_05_1() {

    // fixture setup.  Instantiates data for lower bound 0.5 and upper bound 1.0
    // for the gamma-trucated distribution (with hard-coded parameters 1.2 shape
    // and 0.9 rate).
    GammaContMH gamma_05_1   = gen_gamma_bndl_bndu(0.5, 1.0);

    // exercise SUT
    double out_05_1   = gamma_05_1.log_dgamma_trunc_norm_const();

    // verify outcome
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.4254311104708031, out_05_1,   EPSILON);

    // teardown automatically occurs due to fixture object going out of scope
}




// construct dataset for testing -----------------------------------------------

GammaContMHTest::FromScratchGamma::FromScratchGamma() :
    U(Rcpp::NumericMatrix(9, 1)),
    gamma_specs(Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
                                            Rcpp::_["hyp_a"]    = 1.2,
                                            Rcpp::_["hyp_b"]    = 0.9,
                                            Rcpp::_["hyp_p"]    = 0.5,
                                            Rcpp::_["bnd_l"]    = 0.0,
                                            Rcpp::_["bnd_u"]    = R_PosInf,
                                            Rcpp::_["mh_p"]     = 0.1,
                                            Rcpp::_["mh_delta"] = 0.2)),
    gamma_obj(U, gamma_specs)
{
    // specify values of `U_h`
    U[0] = 1.0;    U[1] = 0.3;    U[2] = 0.4;    U[3] = 1.2;    U[4] = 1.1;
    U[5] = 1.5;    U[6] = 0.2;    U[7] = 0.6;    U[8] = 1.3;

    // specify current value of beta_h and gamma_h
    gamma_obj.m_beta_val = 0.5;
    gamma_obj.m_gam_val = exp(0.5);
}




GammaContMHTest::FromScratchW::FromScratchW() :
    // create `preg_cyc`, a list with each element providing the location and
    // number of days in a pregnancy cycle
    subj0_preg_cyc(Rcpp::IntegerVector::create(Rcpp::_("beg_idx")  = 3,
                                               Rcpp::_("n_days")   = 3,
                                               Rcpp::_("subj_idx") = 0,
                                               Rcpp::_("cyc_idx")  = 1)),
    subj1_preg_cyc(Rcpp::IntegerVector::create(Rcpp::_("beg_idx")  = 3,
                                               Rcpp::_("n_days")   = 3,
                                               Rcpp::_("subj_idx") = 0,
                                               Rcpp::_("cyc_idx")  = 1)),
    preg_cyc(Rcpp::List::create(subj0_preg_cyc, subj1_preg_cyc)),
    // create a vector giving the indices in the daily data that the `W` days
    // correspond to (plus a sentinal entry to signify the end of the vector)
    w_to_days_idx(Rcpp::IntegerVector::create(2, 3, 4, 7, 8, -1)),
    // create unused vector to pass to constructor
    w_cyc_to_subj_idx(Rcpp::IntegerVector::create(0, 1)),
    // create `WGen` object
    w_obj(WGen(preg_cyc, w_to_days_idx, w_cyc_to_subj_idx, 5))
{
    // set `W` values
    w_obj.m_vals[0] = 1;
    w_obj.m_vals[1] = 0;
    w_obj.m_vals[2] = 2;
    w_obj.m_vals[3] = 1;
    w_obj.m_vals[4] = 0;

    // set `W` cycle sums
    w_obj.m_sums[0] = 3;
    w_obj.m_sums[1] = 1;
}




GammaContMHTest::FromScratchXi::FromScratchXi() :
    // create `subj_day_blocks`, a list with each element providing the location
    // and number of days for a subject
    subj0_day_block(Rcpp::IntegerVector::create(Rcpp::_("beg_idx")  = 0,
                                                Rcpp::_("n_days")   = 5)),
    subj1_day_block(Rcpp::IntegerVector::create(Rcpp::_("beg_idx")  = 5,
                                                Rcpp::_("n_days")   = 4)),
    subj_day_blocks(Rcpp::List::create(subj0_day_block, subj1_day_block)),
    // create `XiGen` object
    xi_obj(subj_day_blocks, 1, false)
{
    // set `xi` values
    xi_obj.m_vals[0] = 1.3;
    xi_obj.m_vals[1] = 0.9;
}




GammaContMHTest::FromScratchUbeta::FromScratchUbeta() :
    ubeta_obj(9)
{
    // set `ubeta_obj` values
    ubeta_obj.m_vals[0] =  0.0;    ubeta_obj.m_exp_vals[0] = exp( 0.0);
    ubeta_obj.m_vals[1] = -0.1;    ubeta_obj.m_exp_vals[1] = exp(-0.1);
    ubeta_obj.m_vals[2] =  0.6;    ubeta_obj.m_exp_vals[2] = exp( 0.6);
    ubeta_obj.m_vals[3] =  0.2;    ubeta_obj.m_exp_vals[3] = exp( 0.2);
    ubeta_obj.m_vals[4] =  1.0;    ubeta_obj.m_exp_vals[4] = exp( 1.0);
    ubeta_obj.m_vals[5] = -0.3;    ubeta_obj.m_exp_vals[5] = exp(-0.3);
    ubeta_obj.m_vals[6] =  0.2;    ubeta_obj.m_exp_vals[6] = exp( 0.2);
    ubeta_obj.m_vals[7] =  0.5;    ubeta_obj.m_exp_vals[7] = exp( 0.5);
    ubeta_obj.m_vals[8] =  0.3;    ubeta_obj.m_exp_vals[8] = exp( 0.3);
}




// utility functions -----------------------------------------------------------

void GammaContMHTest::set_seed(int seed_val) {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // set R's internal seed
    set_seed(seed_val);
}




// creates an object with a proposal function that returns a value of 2, and
// samples from the discrete part of the distribution when the random number
// generator provides a value less than 0.4.

GammaContMH GammaContMHTest::gen_gamma_proposal2_propval04() {

    // the only values of importance are those for `mh_p` and `mh_delta`
    Rcpp::NumericVector gamma_specs = Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
                                                                  Rcpp::_["hyp_a"]    = 1.0,
                                                                  Rcpp::_["hyp_b"]    = 1.0,
                                                                  Rcpp::_["hyp_p"]    = 0.5,
                                                                  Rcpp::_["bnd_l"]    = 0.0,
                                                                  Rcpp::_["bnd_u"]    = R_PosInf,
                                                                  Rcpp::_["mh_p"]     = 0.4,
                                                                  Rcpp::_["mh_delta"] = 0.0);

    // specify the current value of `beta_h` and corresponding `gamma_h`
    GammaContMH gamma(m_Uh, gamma_specs);
    gamma.m_beta_val = 2.0;
    gamma.m_gam_val  = exp(2.0);

    return gamma;
}




GammaContMH GammaContMHTest::gen_gamma_bndl_bndu(double bnd_l,  double bnd_u) {

    // the only values of importance are those for `bnd_l` and `bnd_u`
    Rcpp::NumericVector gamma_specs = Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
                                                                  Rcpp::_["hyp_a"]    = 1.2,
                                                                  Rcpp::_["hyp_b"]    = 0.9,
                                                                  Rcpp::_["hyp_p"]    = 0.5,
                                                                  Rcpp::_["bnd_l"]    = bnd_l,
                                                                  Rcpp::_["bnd_u"]    = bnd_u,
                                                                  Rcpp::_["mh_p"]     = 0.4,
                                                                  Rcpp::_["mh_delta"] = 0.1);

    return GammaContMH(m_Uh, gamma_specs);
}




GammaContMH GammaContMHTest::gen_gamma_curr(double curr) {

    // the only values of importance are those for `bnd_l` and `bnd_u`
    Rcpp::NumericVector gamma_specs = Rcpp::NumericVector::create(Rcpp::_["h"]        = 0.0,
                                                                  Rcpp::_["hyp_a"]    = 1.2,
                                                                  Rcpp::_["hyp_b"]    = 0.9,
                                                                  Rcpp::_["hyp_p"]    = 0.5,
                                                                  Rcpp::_["bnd_l"]    = 0.0,
                                                                  Rcpp::_["bnd_u"]    = R_PosInf,
                                                                  Rcpp::_["mh_p"]     = 0.1,
                                                                  Rcpp::_["mh_delta"] = 3.0);

    // specify the current value of `beta_h` and corresponding `gamma_h`
    GammaContMH gamma(m_Uh, gamma_specs);
    gamma.m_beta_val = curr;
    gamma.m_gam_val  = exp(curr);

    return gamma;
}
