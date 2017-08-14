#include <algorithm>
#include "Rcpp.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestGammaCateg.h"

extern UTestFactory g_ut_factory;

using Rcpp::NumericVector;
using Rcpp::as;




GammaCategTest::GammaCategTest() :
    seed_val(g_ut_factory.seed_vals["gamma_categ"]),
    epsilon(UTestFactory::epsilon),
    hyp_a(g_ut_factory.target_data_gamma_categ["hyp_a"]),
    hyp_b(g_ut_factory.target_data_gamma_categ["hyp_b"]),
    m_log_d2_const_terms_all(g_ut_factory.target_data_gamma_categ["m_log_d2_const_terms_all"]),
    m_log_d2_const_terms_zero_one(g_ut_factory.target_data_gamma_categ["m_log_d2_const_terms_zero_one"]),
    m_log_d2_const_terms_one_inf(g_ut_factory.target_data_gamma_categ["m_log_d2_const_terms_one_inf"]),
    m_log_d2_const_terms_zero_half(g_ut_factory.target_data_gamma_categ["m_log_d2_const_terms_zero_half"]),
    log_dgamma_norm_const(g_ut_factory.target_data_gamma_categ["log_dgamma_norm_const_val"]),
    log_dgamma_trunc_const_all(g_ut_factory.target_data_gamma_categ["log_dgamma_trunc_const_all"]),
    log_dgamma_trunc_const_zero_one(g_ut_factory.target_data_gamma_categ["log_dgamma_trunc_const_zero_one"]),
    log_dgamma_trunc_const_one_inf(g_ut_factory.target_data_gamma_categ["log_dgamma_trunc_const_one_inf"]),
    log_dgamma_trunc_const_zero_half(g_ut_factory.target_data_gamma_categ["log_dgamma_trunc_const_zero_half"]),
    a_tilde(g_ut_factory.target_data_gamma_categ["a_tilde"]),
    b_tilde(g_ut_factory.target_data_gamma_categ["b_tilde"]),
    target_ubeta_no_h(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_ubeta_no_h"])),
    p_tilde_all(g_ut_factory.target_data_gamma_categ["p_tilde_all"]),
    p_tilde_zero_one(g_ut_factory.target_data_gamma_categ["p_tilde_zero_one"]),
    p_tilde_one_inf(g_ut_factory.target_data_gamma_categ["p_tilde_one_inf"]),
    p_tilde_zero_half(g_ut_factory.target_data_gamma_categ["p_tilde_zero_half"]),
    target_samples_all(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_samples_all"])),
    target_samples_zero_one(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_samples_zero_one"])),
    target_samples_one_inf(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_samples_one_inf"])),
    target_samples_zero_half(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_samples_zero_half"])),
    target_ubeta_all(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_ubeta_all"])),
    target_ubeta_zero_one(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_ubeta_zero_one"])),
    target_ubeta_one_inf(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_ubeta_one_inf"])),
    target_ubeta_zero_half(as<NumericVector>(g_ut_factory.target_samples_gamma_categ["target_ubeta_zero_half"])) {
}




void GammaCategTest::setUp() {

    // construct gamma_h.  The objects are constructed with the bounds (0, inf),
    // (0, 1), (1, inf), and (0, 1/2), respectively.
    gamma_all = g_ut_factory.gamma_categ_all();
    gamma_zero_one = g_ut_factory.gamma_categ_zero_one();
    gamma_one_inf = g_ut_factory.gamma_categ_one_inf();
    gamma_zero_half = g_ut_factory.gamma_categ_zero_half();

    // construct W
    W = g_ut_factory.W();

    // construct xi
    xi = g_ut_factory.xi();

    // construct ubeta
    ubeta = g_ut_factory.ubeta();

    // construct X
    X = g_ut_factory.X_temp();
}




void GammaCategTest::tearDown() {
    delete gamma_all;
    delete gamma_zero_one;
    delete gamma_one_inf;
    delete gamma_zero_half;
    delete W;
    delete xi;
    delete ubeta;
    delete X;
}




void GammaCategTest::test_constructor() {

    // `gamma_all`
    CPPUNIT_ASSERT(gamma_all->m_bnd_l_is_zero);
    CPPUNIT_ASSERT(gamma_all->m_bnd_u_is_inf);
    CPPUNIT_ASSERT(! gamma_all->m_is_trunc);
    CPPUNIT_ASSERT(gamma_all->m_incl_one);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_d2_const_terms_all,
    				 gamma_all->m_log_d2_const_terms,
    				 epsilon);

    // `gamma_zero_one`
    CPPUNIT_ASSERT(gamma_zero_one->m_bnd_l_is_zero);
    CPPUNIT_ASSERT(! gamma_zero_one->m_bnd_u_is_inf);
    CPPUNIT_ASSERT(gamma_zero_one->m_is_trunc);
    CPPUNIT_ASSERT(gamma_zero_one->m_incl_one);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_d2_const_terms_zero_one,
    				 gamma_zero_one->m_log_d2_const_terms,
    				 epsilon);

    // `gamma_one_inf`
    CPPUNIT_ASSERT(! gamma_one_inf->m_bnd_l_is_zero);
    CPPUNIT_ASSERT(gamma_one_inf->m_bnd_u_is_inf);
    CPPUNIT_ASSERT(gamma_one_inf->m_is_trunc);
    CPPUNIT_ASSERT(gamma_one_inf->m_incl_one);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_d2_const_terms_one_inf,
    				 gamma_one_inf->m_log_d2_const_terms,
    				 epsilon);

    // `gamma_zero_half`
    CPPUNIT_ASSERT(gamma_zero_half->m_bnd_l_is_zero);
    CPPUNIT_ASSERT(! gamma_zero_half->m_bnd_u_is_inf);
    CPPUNIT_ASSERT(gamma_zero_half->m_is_trunc);
    CPPUNIT_ASSERT(! gamma_zero_half->m_incl_one);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(m_log_d2_const_terms_zero_half,
    				 gamma_zero_half->m_log_d2_const_terms,
    				 epsilon);
}




void GammaCategTest::test_calculations() {

    // calc_log_d2_const_terms()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_norm_const,
    				 gamma_all->log_dgamma_norm_const(hyp_a, hyp_b),
    				 epsilon);

    // log_dgamma_trunc_const()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_trunc_const_all,
    				 gamma_all->log_dgamma_trunc_const(hyp_a, hyp_b),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_trunc_const_zero_one,
    				 gamma_zero_one->log_dgamma_trunc_const(hyp_a, hyp_b),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_trunc_const_one_inf,
    				 gamma_one_inf->log_dgamma_trunc_const(hyp_a, hyp_b),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(log_dgamma_trunc_const_zero_half,
    				 gamma_zero_half->log_dgamma_trunc_const(hyp_a, hyp_b),
    				 epsilon);

    // calc_a_tilde()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(a_tilde,
    				 gamma_all->calc_a_tilde(*W),
    				 epsilon);

    // calc_b_tilde(). Have to check that `U * beta` was updated to assume the
    // values of `U * beta - U_h * beta_h`
    CPPUNIT_ASSERT_DOUBLES_EQUAL(b_tilde,
    				 gamma_all->calc_b_tilde(*ubeta, *xi, *X),
    				 epsilon);
    CPPUNIT_ASSERT(std::equal(target_ubeta_no_h.begin(),
    			      target_ubeta_no_h.end(),
    			      ubeta->vals(),
    			      UTestFactory::eq_dbl));

    // calc_p_tilde()
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p_tilde_all,
    				 gamma_all->calc_p_tilde(a_tilde, b_tilde),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p_tilde_zero_one,
    				 gamma_zero_one->calc_p_tilde(a_tilde, b_tilde),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p_tilde_one_inf,
    				 gamma_one_inf->calc_p_tilde(a_tilde, b_tilde),
    				 epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(p_tilde_zero_half,
    				 gamma_zero_half->calc_p_tilde(a_tilde, b_tilde),
    				 epsilon);
}




void GammaCategTest::test_sample_gamma() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // `gamma_all`
    set_seed(seed_val);
    for (int i = 0; i < target_samples_all.size(); ++i) {
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_all[i],
    				     gamma_all->sample_gamma(a_tilde, b_tilde, p_tilde_all),
    				     epsilon);
    }

    // `gamma_zero_one`
    set_seed(seed_val);
    for (int i = 0; i < target_samples_zero_one.size(); ++i) {
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_zero_one[i],
    				     gamma_zero_one->sample_gamma(a_tilde, b_tilde, p_tilde_zero_one),
    				     epsilon);
    }

    // `gamma_one_inf`
    set_seed(seed_val);
    for (int i = 0; i < target_samples_one_inf.size(); ++i) {
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_one_inf[i],
    				     gamma_one_inf->sample_gamma(a_tilde, b_tilde, p_tilde_one_inf),
    				     epsilon);
    }

    // `gamma_zero_half`
    set_seed(seed_val);
    for (int i = 0; i < target_samples_zero_half.size(); ++i) {
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_zero_half[i],
    				     gamma_zero_half->sample_gamma(a_tilde, b_tilde, p_tilde_zero_half),
    				     epsilon);
    }
}




void GammaCategTest::test_sample() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];

    // obtain independent copies of `U * beta`
    UProdBeta* ubeta_all       = g_ut_factory.ubeta();
    UProdBeta* ubeta_zero_one  = g_ut_factory.ubeta();
    UProdBeta* ubeta_one_inf   = g_ut_factory.ubeta();
    UProdBeta* ubeta_zero_half = g_ut_factory.ubeta();

    // `gamma_all`
    set_seed(seed_val);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_all[0],
    				 gamma_all->sample(*W, *xi, *ubeta_all, *X),
    				 epsilon);
    CPPUNIT_ASSERT(std::equal(target_ubeta_all.begin(),
    			      target_ubeta_all.end(),
    			      ubeta_all->vals(),
    			      UTestFactory::eq_dbl));

    // `gamma_zero_one`
    set_seed(seed_val);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_zero_one[0],
    				 gamma_zero_one->sample(*W, *xi, *ubeta_zero_one, *X),
    				 epsilon);
    CPPUNIT_ASSERT(std::equal(target_ubeta_zero_one.begin(),
    			      target_ubeta_zero_one.end(),
    			      ubeta_zero_one->vals(),
    			      UTestFactory::eq_dbl));

    // `gamma_one_inf`
    set_seed(seed_val);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_one_inf[0],
    				 gamma_one_inf->sample(*W, *xi, *ubeta_one_inf, *X),
    				 epsilon);
    CPPUNIT_ASSERT(std::equal(target_ubeta_one_inf.begin(),
    			      target_ubeta_one_inf.end(),
    			      ubeta_one_inf->vals(),
    			      UTestFactory::eq_dbl));

    // `gamma_zero_half`
    set_seed(seed_val);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_samples_zero_half[0],
    				 gamma_zero_half->sample(*W, *xi, *ubeta_zero_half, *X),
    				 epsilon);
    CPPUNIT_ASSERT(std::equal(target_ubeta_zero_half.begin(),
    			      target_ubeta_zero_half.end(),
    			      ubeta_zero_half->vals(),
    			      UTestFactory::eq_dbl));

    delete ubeta_all;
    delete ubeta_zero_one;
    delete ubeta_one_inf;
    delete ubeta_zero_half;
}
