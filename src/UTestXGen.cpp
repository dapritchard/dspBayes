#include <algorithm>
#include "Rcpp.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestXGen.h"

using Rcpp::as;

extern UTestFactory g_ut_factory;




XGenTest::XGenTest()
{}




void XGenTest::setUp() {
    X     = g_ut_factory.X();
    W     = g_ut_factory.W();
    xi    = g_ut_factory.xi_no_rec();
    ubeta = g_ut_factory.ubeta();
    utau  = g_ut_factory.utau();
}




void XGenTest::tearDown() {
    delete X;
    delete W;
    delete xi;
    delete ubeta;
    delete utau;
}




void XGenTest::test_constructor() {
    CPPUNIT_ASSERT(Rcpp::is_true(Rcpp::all(X_rcpp == X->m_x_rcpp)));
    CPPUNIT_ASSERT_EQUAL((int) miss_cyc.size(), X->m_n_miss_cyc);
    CPPUNIT_ASSERT_EQUAL(cohort_sex_prob, X->m_cohort_sex_prob);
    CPPUNIT_ASSERT_EQUAL(sex_coef, X->m_sex_coef);
}




void XGenTest::test_sample() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    X->sample(*W, *xi, *ubeta, *utau);
    CPPUNIT_ASSERT(std::equal(target_x_samples.begin(),
			      target_x_samples.end(),
			      X->vals()));
}




// void XGenTest::test_sample_cycle() {

//     // register seed function
//     Rcpp::Environment base("package:base");
//     Rcpp::Function set_seed = base["set.seed"];
//     set_seed(seed_val);

//     X->sample_cycle(miss_cyc, W, xi, ubeta, utau);
//     CPPUNIT_ASSERT(std::equal(target_sample_cycle.begin(),
// 			      target_sample_cycle.end(),
// 			      X.vals() + miss_day[miss_cyc->beg_idx].idx));
// }




void XGenTest::test_calc_prior_prob() {

    double out_no = X->calc_prior_prob(*miss_day, *utau, 0);
    double out_yes = X->calc_prior_prob(*miss_day, *utau, 1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_prior_prob_no_prev, out_no, epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_prior_prob_yes_prev, out_yes, epsilon);
}




void XGenTest::test_calc_posterior_prob() {

    double out = XGen::calc_posterior_prob(*miss_day, *ubeta, xi_i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_posterior_prob, out, epsilon);
}




void XGenTest::test_sample_x_ijk() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // test "sex yesterday" samples
    for (Rcpp::IntegerVector::const_iterator curr = target_x_ijk_samples.begin();
	 Rcpp::IntegerVector::const_iterator end = target_x_ijk_samples.end();
	 ++curr) {

	CPPUNIT_ASSERT_EQUAL(*curr, XGen::sample_x_ijk(prior_prob_yes, posterior_prob_yes));
    }
}




void XGenTest::test_sample_day_before_fw_sex() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // test "sex yesterday" samples
    for (Rcpp::IntegerVector::const_iterator curr = target_day_before_samples.begin();
	 Rcpp::IntegerVector::const_iterator end = target_day_before_samples.end();
	 ++curr) {

	CPPUNIT_ASSERT_EQUAL(*curr, X->sample_day_before_fw_sex());
    }
}
