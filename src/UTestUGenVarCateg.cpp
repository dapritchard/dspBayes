// #include <algorithm>
#include "Rcpp.h"

#include "CoefGen.h"
#include "UGenVar.h"
#include "UProdBeta.h"
#include "UProdTau.h"
#include "UTestFactory.h"
#include "UTestUGenVarCateg.h"
#include "WGen.h"
#include "XGen.h"
#include "XiGen.h"


// using Rcpp::IntegerVector;
// using Rcpp::NumericVector;
// using Rcpp::as;

extern UTestFactory g_ut_factory;




// UGenVarCategTest::UGenVarCategTest() :
//     X_rcpp(g_ut_factory.X_rcpp),
//     miss_cyc_rcpp(g_ut_factory.x_miss_cyc),
//     miss_day_rcpp(g_ut_factory.x_miss_day),
//     cohort_sex_prob(g_ut_factory.tau_coefs["cohort_sex_prob"]),
//     sex_coef(g_ut_factory.tau_coefs["sex_coef"]),
//     miss_day_idx((int) g_ut_factory.input_x["miss_day_idx"]),
//     day_idx((int) g_ut_factory.input_x["day_idx"]),
//     prior_prob_yes(g_ut_factory.input_x["prior_prob"]),
//     posterior_prob_yes(g_ut_factory.input_x["posterior_prob"]),
//     xi_i(g_ut_factory.input_x["xi_i"]),
//     seed_val(g_ut_factory.seed_vals["X"]),
//     epsilon(UTestFactory::epsilon),
//     target_x_samples(as<IntegerVector>(g_ut_factory.target_samples_x["x_samples"])),
//     target_x_ijk_samples(as<IntegerVector>(g_ut_factory.target_samples_x["x_ijk_samples"])),
//     target_day_before_samples(as<IntegerVector>(g_ut_factory.target_samples_x["day_before_samples"])),
//     target_prior_prob_no_prev(g_ut_factory.target_samples_x["prior_prob_no_prev"]),
//     target_prior_prob_yes_prev(g_ut_factory.target_samples_x["prior_prob_yes_prev"]),
//     target_posterior_prob(g_ut_factory.target_samples_x["posterior_prob"]) {
// }




void UGenVarCategTest::setUp() {

    // W     = g_ut_factory.W();
    // xi    = g_ut_factory.xi_no_rec();
    // ubeta = g_ut_factory.ubeta();
    // utau  = g_ut_factory.utau();

    // // copy X data so that tests don't cause persistent changes
    // x_rcpp_copy = new Rcpp::IntegerVector(X_rcpp.begin(), X_rcpp.end());
    // X = new XGen(*x_rcpp_copy, miss_cyc_rcpp, miss_day_rcpp, cohort_sex_prob, sex_coef);

}




void UGenVarCategTest::tearDown() {
    // delete X;
    // delete W;
    // delete xi;
    // delete ubeta;
    // delete utau;
    // delete x_rcpp_copy;
}




void UGenVarCategTest::test_constructor() {

    // PregCyc* miss_cyc = PregCyc::list_to_arr(miss_cyc_rcpp);
    // XGen::XMissDay* miss_day = XGen::XMissDay::list_to_arr(miss_day_rcpp);

    // CPPUNIT_ASSERT(Rcpp::is_true(Rcpp::all(X_rcpp == X->m_x_rcpp)));
    // CPPUNIT_ASSERT_EQUAL(X->m_vals, X->m_x_rcpp.begin());
    // // TODO: test m_miss_cyc
    // CPPUNIT_ASSERT_EQUAL((int) miss_cyc_rcpp.size(), X->m_n_miss_cyc);
    // // TODO: test m_miss_day
    // CPPUNIT_ASSERT_EQUAL(cohort_sex_prob, X->m_cohort_sex_prob);
    // CPPUNIT_ASSERT_EQUAL(sex_coef, X->m_sex_coef);

    // delete[] miss_cyc;
    // delete[] miss_day;

}
