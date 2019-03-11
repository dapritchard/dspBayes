#include <algorithm>
#include "Rcpp.h"
#include "WGen.h"
#include "XiGen.h"
#include "UProdBeta.h"
#include "UTestFactory.h"
#include "UTestXGen.h"


// #include <iostream>


using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::as;

extern UTestFactory g_ut_factory;




XGenTest::XGenTest() :
    x_rcpp(g_ut_factory.x_rcpp),
    miss_cyc_rcpp(g_ut_factory.x_miss_cyc),
    miss_day_rcpp(g_ut_factory.x_miss_day),
    cohort_sex_prob(g_ut_factory.tau_coefs["cohort_sex_prob"]),
    sex_coef(g_ut_factory.tau_coefs["sex_coef"]),
    miss_day_idx((int) g_ut_factory.input_x["miss_day_idx"]),
    day_idx((int) g_ut_factory.input_x["day_idx"]),
    prior_prob_yes(g_ut_factory.input_x["prior_prob"]),
    posterior_prob_yes(g_ut_factory.input_x["posterior_prob"]),
    xi_i(g_ut_factory.input_x["xi_i"]),
    seed_val(g_ut_factory.seed_vals["X"]),
    epsilon(UTestFactory::epsilon),
    target_x_samples(as<IntegerVector>(g_ut_factory.target_samples_x["x_samples"])),
    target_x_ijk_samples(as<IntegerVector>(g_ut_factory.target_samples_x["x_ijk_samples"])),
    target_day_before_samples(as<IntegerVector>(g_ut_factory.target_samples_x["day_before_samples"])),
    target_prior_prob_no_prev(g_ut_factory.target_samples_x["prior_prob_no_prev"]),
    target_prior_prob_yes_prev(g_ut_factory.target_samples_x["prior_prob_yes_prev"]),
    target_posterior_prob(g_ut_factory.target_samples_x["posterior_prob"]) {
}




void XGenTest::setUp() {

    W     = g_ut_factory.W();
    xi    = g_ut_factory.xi_no_rec();
    ubeta = g_ut_factory.ubeta();
    utau  = g_ut_factory.utau();

    // copy X data so that tests don't cause persistent changes
    x_rcpp_copy = new Rcpp::IntegerVector(x_rcpp.begin(), x_rcpp.end());
    X = new XGen(*x_rcpp_copy, miss_cyc_rcpp, miss_day_rcpp, cohort_sex_prob, sex_coef);

}




void XGenTest::tearDown() {
    delete X;
    delete W;
    delete xi;
    delete ubeta;
    delete utau;
    delete x_rcpp_copy;
}




void XGenTest::test_constructor() {

    PregCyc* miss_cyc = PregCyc::list_to_arr(miss_cyc_rcpp);
    XGen::XMissDay* miss_day = XGen::XMissDay::list_to_arr(miss_day_rcpp);

    CPPUNIT_ASSERT(Rcpp::is_true(Rcpp::all(x_rcpp == X->m_x_rcpp)));
    CPPUNIT_ASSERT_EQUAL(X->m_vals, X->m_x_rcpp.begin());
    // TODO: test m_miss_cyc
    CPPUNIT_ASSERT_EQUAL((int) miss_cyc_rcpp.size(), X->m_n_miss_cyc);
    // TODO: test m_miss_day
    CPPUNIT_ASSERT_EQUAL(cohort_sex_prob, X->m_cohort_sex_prob);
    CPPUNIT_ASSERT_EQUAL(sex_coef, X->m_sex_coef);

    delete[] miss_cyc;
    delete[] miss_day;
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




// TODO: test imputed previous day sex recording



// // TODO: test a single cycle?
// void XGenTest::test_sample_cycle() {

//     PregCyc* miss_cyc = PregCyc::list_to_arr(miss_cyc_rcpp);

//     // register seed function
//     Rcpp::Environment base("package:base");
//     Rcpp::Function set_seed = base["set.seed"];
//     set_seed(seed_val);

//     X->sample_cycle(miss_cyc + , W, xi, ubeta, utau);
//     CPPUNIT_ASSERT(std::equal(target_sample_cycle.begin(),
//                            target_sample_cycle.end(),
//                            X.vals() + miss_day[miss_cyc->beg_idx].idx));

//     delete[] miss_cyc;
// }




void XGenTest::test_calc_prior_prob() {

    double out_no = X->calc_prior_prob(*utau, miss_day_idx, 0);
    double out_yes = X->calc_prior_prob(*utau, miss_day_idx, 1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_prior_prob_no_prev, out_no, epsilon);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_prior_prob_yes_prev, out_yes, epsilon);
}




void XGenTest::test_calc_posterior_prob() {

    double out = XGen::calc_posterior_prob(*ubeta, xi_i, day_idx);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(target_posterior_prob, out, epsilon);
}




void XGenTest::test_sample_x_ijk() {

    // register seed function
    Rcpp::Environment base("package:base");
    Rcpp::Function set_seed = base["set.seed"];
    set_seed(seed_val);

    // test "sex yesterday" samples
    for (Rcpp::IntegerVector::const_iterator curr = target_x_ijk_samples.begin();
         curr < target_x_ijk_samples.end();
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
         curr < target_day_before_samples.end();
         ++curr) {

        CPPUNIT_ASSERT_EQUAL(*curr, X->sample_day_before_fw_sex());
    }
}
